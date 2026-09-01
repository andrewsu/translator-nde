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

## Layout

| Path | What |
|---|---|
| `src/translator_nde/nde.py` | NDE Discovery API client (query, facet, scroll) |
| `src/translator_nde/ids.py` | CURIE ↔ NDE bridge: Node Normalizer + Name Resolver |
| `src/translator_nde/translator.py` | ARS creative-mode client + path extraction |
| `src/translator_nde/gxa.py` | Route A: drug→gene edges vs. GXA DE contrasts |
| `scripts/run_ars.py` | Submit a disease query, extract paths |
| `tests/` | Regression fixtures (see below) |

## Setup

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
PYTHONPATH=src .venv/bin/python -m pytest tests/ -q
```

## Two things that will bite you

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

## Status

Early. Route A works and reproduces known pharmacology — dexamethasone → FKBP5 gives 47/47
direction agreement (median log2FC 3.4, min adj-p 1e-115, 10 experiments), and → TSC22D3 (GILZ)
49/49; both are canonical glucocorticoid-induced genes. Route B is not built yet.

See [`results/REPORT.md`](results/REPORT.md) for the first worked example (asthma).
