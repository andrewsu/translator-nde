# translator-nde

Exploring the interface between the [NCATS Biomedical Data
Translator](https://doi.org/10.1111/cts.70284) and the [NIAID Data
Ecosystem](https://doi.org/10.1128/msystems.01270-25).

## The idea

Ask Translator *"what drugs might treat disease X"* and the answers come back as **paths** —
often `drug → gene/protein → disease`. Each such path is a mechanistic hypothesis, but
Translator's evidence is assertion-level: it tells you an edge *was asserted*, and by whom, not
whether data supports it. The Translator paper names this limitation directly — provenance
suffers an *"inability to accurately trace complex and incomplete trails of ownership and
algorithmic transformations of primary data."*

NDE indexes 14.2M dataset records from 48 repositories. So: for each `drug → gene` hop, can we
find data that **tests** it?

```
  MONDO disease
       │  ARS creative-mode (knowledge_type: "inferred", biolink:treats)
       ▼
  drug ──► gene ──► disease
   │        │
   ▼        ▼
 ROUTE A  GXA Inference          precomputed DE contrasts, direction-aware
 ROUTE B  GEO reanalysis         compute DE ourselves where GXA has nothing
```

## Why this is not just a text search

The disease axis is already joined: NDE's `healthCondition` is MONDO-mapped, done *via Translator
KP APIs*. But **NDE has no chemical and no gene field** anywhere in its 1030-field mapping, so the
`drug → gene` hop — the one we most want — has nothing structured to join on.

Naive co-mention does not work. Measured against the live API:

| Query | Datasets |
|---|--:|
| `metformin AND PRKAA1` | **1** |
| `hydroxychloroquine AND TLR9` | 12 |
| `imatinib AND ABL1` | 52 |

Metformin→AMPK is textbook pharmacology and returns one dataset. Descriptions name the drug; they
almost never name the gene.

**Route A** works around this. NDE staging carries 10.4M `@type:Inference` records — one per
differential-expression contrast, mostly from the EBI Gene Expression Atlas. Each has the gene
structured (`observationAbout.identifier`, Ensembl), an effect size (`value`, log2FC), a
significance (`marginOfError.value`), and **Biolink-typed direction and aspect qualifiers**. Those
line up directly with Translator's edge qualifiers:

| Translator drug→gene edge | GXA Inference record |
|---|---|
| `object_aspect_qualifier: abundance` | `tripleSubjectQualifier.aspect = "abundance"` |
| `object_direction_qualifier: increased` | `tripleSubjectQualifier.direction = "increase"` |

So the question upgrades from *"does a dataset exist?"* to ***"does the measured direction agree
with the asserted one?"*** — with an effect size attached.

## Four routes

Routes A and B ask *"is the asserted edge supported by data?"* Three worked examples showed that
question is largely unanswerable as posed, because **Translator's drug→gene edges assert changes in
*activity* while the expression atlases measure *abundance***. Kinase inhibitors, receptor
antagonists and neutralising antibodies do not change their target's transcript — baricitinib →
JAK1/2/3, infliximab → TNF, JQ1 → BRD4 are all textbook-correct and all invisible (or
contradictory) to expression data.

Routes C and D ask questions the data can actually answer.

| route | question | data | status |
|---|---|---|---|
| **A** | does the drug change the gene's *expression*? | GXA `@type:Inference` (staging) | built — **0.9%** coverage, mostly the wrong question |
| **B** | can we *compute* that from raw data? | GEO count matrices via NDE | built + validated |
| **C** | what *other* compounds move this gene? | GXA, queried gene-first | designed, query verified |
| **D** | does *activity* data confirm the inhibition? | PubChem BioAssay, ChEMBL | built — 52% coverage, 20 true negatives |

**Route C** inverts the question. Given `Drug1 –inhibits→ Gene2 –associated_with→ Disease3`, take
the therapeutic hypothesis ("less Gene2 is good for Disease3") and ask GXA which *other* compounds
reduce Gene2. Abundance reduction becomes an alternate modality to activity inhibition, and each
hit is a repurposing candidate Translator did not propose. This plays to GXA's one structured
axis — `observationAbout` carries a gene symbol and Ensembl id, so no text matching is needed.

**Route D** uses assays that measure inhibition directly, and is the highest-coverage route in
the project: **367 of 710 drug→gene edges (52%) have a directly measured compound–target result**,
against 1–9% for Route A on the same edges. The join is exact end to end — Translator emits
`NCBIGene:7124`, PubChem's `Target GeneID` column is `7124`, and Node Normalizer supplies the
compound's CID and ChEMBL id — so no text matching is involved anywhere. PubChem's *Inactive* rows
give 20 genuine negatives, which GXA structurally cannot (it stores only significant results).

**Neither source is in NDE** — its catalog has only LINCS (424) and ReframeDB (408) as
activity-adjacent, no ChEMBL/PubChem/BindingDB. Since activity is the modality Translator's edges
are actually about, that is a real coverage gap.

The headline finding is a structural one about Translator. Of 710 edges, 330 carry a direction
qualifier and 63 have a curated ChEMBL mechanism, but **only 1 has both** — because the qualifiers
come from screening databases (BindingDB, Pharos) whose compounds are research chemicals, while
curated mechanism-of-action exists almost only for approved drugs, which arrive via
DrugBank/DrugCentral/DGIdb emitting a bare `biolink:affects` with no qualifier. Measured: 4% of
directional-edge drugs have any ChEMBL mechanism, versus 79% of the rest. So Route D's value is
not adjudicating directions Translator asserts — it is *supplying* the action type Translator
omits, correctly and for approved drugs, on the very edges Route A scored 1/174. See
[`results/REPORT.md`](results/REPORT.md#worked-example-4--route-d-activity-data).

## Layout

| Path | What |
|---|---|
| `src/translator_nde/nde.py` | NDE Discovery API client (query, facet, scroll) |
| `src/translator_nde/ids.py` | CURIE ↔ NDE bridge: Node Normalizer + Name Resolver |
| `src/translator_nde/translator.py` | ARS creative-mode client + path extraction |
| `src/translator_nde/gxa.py` | Route A: drug→gene edges vs. GXA DE contrasts |
| `src/translator_nde/reanalysis.py` | Route B: GEO reanalysis from NDE sample arm labels |
| `src/translator_nde/activity.py` | Route D: PubChem BioAssay + ChEMBL activity evidence |
| `scripts/run_ars.py` | Submit a disease query, extract paths |
| `scripts/run_route_{a,b,d}.py` | Score a `paths.json` edge set against each route |
| `tests/` | Regression fixtures (see below) |

## Setup

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
PYTHONPATH=src .venv/bin/python -m pytest tests/ -q
```

## Three things that will bite you

**1. `Inference` is staging-only.** It does not exist in production, and an `@type:Inference`
query against prod returns **0 hits with no error**. `NDEClient` raises rather than let that pass
silently. Records were added by nde-crawlers PR #368 (merged 2026-07-10) and had not been promoted
to prod as of 2026-09-01.

**2. Matching a drug naively counts the wrong contrasts.** The compound is free text, and it may
appear in *both* arms of a contrast — meaning it is an experimental tool, not the variable. The
record `gxa_e_geod_54049_g1_g5_at1g03790` contrasts
`'pBeaconRFP_GR::bZIP1; …10 uM dexamethasone…'` vs `'empty vector; …10 uM dexamethasone…'` in
*Arabidopsis*, where dexamethasone is a GR inducer and the real contrast is genotype. So:

```
q=@type:Inference AND species.identifier:9606
  AND variableMeasured.value:<drug>
  AND NOT measurementDenominator.value:<drug>     # <- mandatory
```

The exclusion removes 7,590 dexamethasone contrasts. Both records are pinned in
`tests/test_gxa_fixtures.py`.

**3. …and the query above is still not enough.** Elasticsearch matches a synonym *anywhere* in the
test-arm text, which counts contrasts where the drug name is incidental. Measured, all previously
scored as evidence: `MITF-RFP-HA overexpression` matched **rifampicin** via the synonym "RFP" (red
fluorescent protein); `no response to infliximab treatment` matched **nitric oxide** via "NO";
`A/CA/04/2009 Influenza virus` matched **calcium** via "CA". There is no structured factor type to
filter on — `constraintProperty` is identical for genuine and spurious contrasts — but
`variableMeasured.value` is a comma-separated factor list, so `gxa.factor_supports_drug()` requires
a synonym to **be** a factor or to **lead** one followed only by a dose. This cut Route A's measured
coverage from 19/460 edges to **4/460**; see
[`results/REPORT.md`](results/REPORT.md#the-fix-require-the-drug-to-occupy-the-variable-position).

## Status

Early. Route A works and reproduces known pharmacology — dexamethasone → FKBP5 gives 47/47
direction agreement (median log2FC 3.4, min adj-p 1e-115, 10 experiments), and → TSC22D3 (GILZ)
49/49; both are canonical glucocorticoid-induced genes. Route B is not built yet.

See [`results/REPORT.md`](results/REPORT.md) for the first worked example (asthma).

## License

MIT — see [LICENSE](LICENSE).

`src/translator_nde/_de.py` is vendored from
[andrewsu/DN-meta-analysis](https://github.com/andrewsu/DN-meta-analysis) (same author) and is
covered by the same license here; see the file header for provenance.

Data retrieved through this code belongs to its original sources and carries their terms —
[NIAID Data Ecosystem](https://data.niaid.nih.gov), [NCATS Translator](https://ui.transltr.io),
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), and the
[EBI Gene Expression Atlas](https://www.ebi.ac.uk/gxa/).
