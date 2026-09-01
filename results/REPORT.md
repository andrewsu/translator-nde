# Worked example 1 — asthma

**Date:** 2026-09-01 · **Disease:** asthma (`MONDO:0004979`)
**ARS pk:** `5b656c0f-b7da-4db4-ba1f-d3a794b422d4`
**NDE:** staging (`api-staging.data.niaid.nih.gov/v1`)
**Artifacts:** `data/ars/5b656c0f-.../paths.json`, `.../route_a.json`, `results/route_a_asthma.log`

## TL;DR

The pipeline works end to end, and Route A is **precise but narrow**: of 123 evaluable
`drug → gene` edges from a real Translator creative-mode answer, GXA had a matching
differential-expression contrast for **4 (3%)**.

The bottleneck is **GXA's drug coverage, not the matching logic** — Expression Atlas contains no
contrasts at all for the drugs asthma answers are actually built from. This is the finding that
matters for the project: **Route B (GEO reanalysis) is required, not optional.**

## What Translator returned

Creative-mode `biolink:treats` on asthma, 13 ARAs, ~40 s to completion.

| ARA | payload | results | aux graphs | 2-hop drug→gene→disease paths |
|---|--:|--:|--:|--:|
| ara-arax | 10.3 MB | 260 | 296 | **311** |
| ara-unsecret | 0.8 MB | 227 | 248 | **117** |
| ara-bte | 6.9 MB | 500 | 826 | 0 |
| kp-molecular | 1.6 MB | 525 | 0 | 0 |

**428 paths → 220 distinct drug→gene edges → 76 distinct genes, 174 distinct drugs.**

The genes are exactly right for asthma, which is a good sign the extraction is sound:
`P2RX3` (82), `NR3C1` (54, the glucocorticoid receptor), `ADRB2` (32, the β2-agonist target),
`HRH1` (30), `PTGS2` (22), `TNF` (19), `IL4R` (14).

### BTE returns no mechanistic intermediates

Its 0 is not an extractor bug. BTE's knowledge graph for this query contains **0 drug→gene edges
out of 1,666**; its answers are drug→disease directly (`treats`,
`treats_or_applied_or_studied_to_treat`, `in_clinical_trials_for`) plus disease subclass
reasoning. Mechanistic paths come from ARAX and unsecret-agent only.

### Only 19 of 123 edges carry a direction qualifier

**All 311 ARAX paths have `direction: None`.** ARAX does not emit
`object_direction_qualifier` / `object_aspect_qualifier` on these edges; only unsecret-agent does.
So the direction-agreement test — the strongest thing Route A can do — is available on
**15% of edges**. Everything else can only be scored for coverage.

## What Route A found

97 of 220 edges were skipped as not drug-like (IUPAC strings, synthetic peptides such as
`H-D-Phe-His-Leu-Leu-Arg-…`), leaving 123 evaluated in 96 s.

| Verdict | Edges |
|---|--:|
| not covered by any GXA contrast | 119 |
| agrees | 2 |
| ambiguous (covered, no asserted direction) | 2 |
| disagrees | 0 |
| **covered by ≥1 contrast** | **4 / 123 (3%)** |

| Drug → gene | contrasts | agree | median log2FC | verdict |
|---|--:|--:|--:|---|
| Rifampicin → PPARGC1A | 1 | 1 | 3.2 | agrees |
| Cyclic AMP → PPARGC1A | 2 | 2 | 1.85 | agrees |
| Infliximab → TNF | 3 | – | 1.1 | ambiguous |
| Antibodies → TNF | 1 | – | 1.6 | ambiguous |

## Why coverage is so low

Not the matcher. GXA simply does not contain the drugs. Human test-arm contrast counts, measured
directly:

| Drug | any species | human |
|---|--:|--:|
| dexamethasone | 57,679 | 23,314 |
| metformin | 2,818 | 928 |
| rifampicin | 2,312 | 108 |
| prednisolone | 481 | 232 |
| albuterol | 6 | 6 |
| **budesonide** | **0** | 0 |
| **formoterol** | **0** | 0 |
| **salbutamol** | **0** | 0 |
| **fluticasone** | **0** | 0 |
| **montelukast** | **0** | 0 |
| theophylline | 222 | **0** |
| imatinib | 80 | **0** |
| aspirin | 213 | **0** |

The first-line asthma therapeutics have **zero** contrasts. Dexamethasone is an outlier, not the
norm — GXA is a curated re-analysis of selected experiments, not a drug-perturbation atlas.

## Positive control

Route A does work where GXA has data, and reproduces textbook pharmacology:

| Drug → gene | contrasts | agree | median log2FC | min adj-p | experiments |
|---|--:|--:|--:|--:|--:|
| dexamethasone → FKBP5 | 47 | 47 | 3.4 | 1.0e-115 | 10 |
| dexamethasone → TSC22D3 (GILZ) | 49 | 49 | 2.7 | 0.0 | 12 |
| dexamethasone → NR3C1 | 5 | 2 (3 disagree) | – | 9.5e-22 | 3 |

FKBP5 and TSC22D3 are the canonical glucocorticoid-induced genes. NR3C1 coming back **ambiguous**
is correct, not a failure — glucocorticoid-receptor autoregulation is genuinely bidirectional and
context-dependent.

## Conclusions

1. **The bridge is real but the precomputed route is thin.** 3% coverage on a real disease query
   is too low to be useful alone. Route B has to carry the load.
2. **Prioritise Route B on the drugs that matter.** Asthma's actual drugs are absent from GXA but
   present in GEO — production NDE has 5,253 dexamethasone datasets (1,302 GEO), and
   drug-anchored search plus the 8.7M `Sample` records should reach budesonide/formoterol.
3. **Direction-agreement is only testable on 15% of edges** given ARAX omits qualifiers. Either
   restrict the strong claim to unsecret-agent edges, or recover directions from a KP that
   supplies them (DGIdb-derived edges do).
4. **Filter drugs before querying.** 44% of edges were research compounds and peptides. The
   length/IUPAC heuristic is crude; normalising to a real drug vocabulary would be better.

## Reproduce

```bash
PYTHONPATH=src .venv/bin/python scripts/run_ars.py MONDO:0004979 asthma
PYTHONPATH=src .venv/bin/python scripts/run_route_a.py data/ars/<pk>/paths.json
```

⚠️ ARS results change over time, so a fresh run will not match this exactly. The archived
`paths.json` is the citable artifact.

---

# Working backward from GXA

Example 1 asked "here are Translator's answers, does GXA cover them?" and got 3%. This asks the
inverse: **which drugs does GXA actually contain**, and what diseases are they used for? If the
bridge works anywhere, it works there.

## What is in GXA

GXA has no chemical field and no facetable field carrying the compound
(`variableMeasured.value` and `measurementQualifier` are both `text`), so the drug vocabulary has
to be probed rather than listed. Two approaches:

**1. Sample the contrast strings.** 40,000 human `Inference` records yielded only 299 distinct
contrasts — scroll order is grouped by experiment, so this is a biased sample and the counts mean
nothing. It is still informative about *shape*: the top arms are things like
`blood alcohol content 0.04%rising`, `early mesodermal progenitors`, `clear cell renal carcinoma`,
`cultured skin substitute; 28 day`, `SARS coronavirus Urbani`. **Most GXA contrasts are not drug
perturbations at all** — they are disease-vs-normal, cell type, developmental stage, timepoint and
genotype comparisons.

**2. Probe a drug vocabulary.** DrugMechDB supplies 1,652 curated drugs, each already paired with
the diseases it treats — which is exactly what is needed to close the loop back to Translator.

> **66 of 1,643 drugs (4.3%) have ≥1 human GXA test-arm contrast.**

| drug | contrasts | ref-arm | note |
|---|--:|--:|---|
| doxycycline | 29,768 | 3,291 | **Tet-on inducer** — largely a tool compound |
| dexamethasone | 23,314 | 0 | clean |
| valproic acid | 21,056 | 0 | clean |
| cisplatin | 18,766 | 3,209 | cell-line cytotoxicity |
| doxorubicin | 13,244 | 3,370 | ref-arm 25% — mixed |
| infliximab | 12,759 | 375 | clean |
| estradiol | 7,097 | 1,872 | ref-arm 26% — mixed |
| panobinostat | 4,759 | 0 | clean |
| acetaminophen | 4,543 | 0 | clean |
| anakinra | 4,070 | 0 | clean |
| gefitinib | 3,817 | 1,486 | ref-arm 39% — mixed |
| sorafenib | 2,962 | 0 | clean |
| digoxin | 2,817 | 0 | clean |
| methotrexate | 2,289 | 3,966 | **more in ref arm than test arm** |
| prednisone | 2,289 | 0 | clean |
| metformin | 928 | 0 | clean |

The reference-arm column is doing real work here. `doxycycline` topping the list is an artefact —
it is the Tet-on induction agent in a large number of experiments, not a tested therapeutic. And
`methotrexate` appears in *more* reference arms than test arms. Without the
`NOT measurementDenominator.value:<drug>` clause both would be badly overcounted, which is the
same failure mode as the Arabidopsis dexamethasone record in `tests/test_gxa_fixtures.py`.

The covered set skews to three families: **corticosteroids** (dexamethasone, prednisone,
prednisolone, methylprednisolone), **oncology cytotoxics** used as cell-line stressors (cisplatin,
doxorubicin, paclitaxel, docetaxel, methotrexate), and a **small number of biologics**
(infliximab, anakinra, etanercept). Everyday small-molecule therapeutics are largely absent.

## Which diseases does that point to

Ranking DrugMechDB diseases by how many GXA-covered drugs treat them:

| disease | covered drugs | net contrasts | drugs |
|---|--:|--:|---|
| **Rheumatoid arthritis** | **9** | 41,298 | dexamethasone, infliximab, anakinra, methotrexate, etanercept, prednisone… |
| Ankylosing spondylitis | 6 | 39,376 | dexamethasone, infliximab, prednisone, etanercept… |
| Pemphigus | 5 | 28,769 | dexamethasone, prednisone, benzoic acid, methylprednisolone… |
| Psoriasis | 5 | 27,220 | dexamethasone, methotrexate, prednisone, methylprednisolone… |
| Non-small cell lung cancer | 5 | 5,058 | gefitinib, docetaxel, paclitaxel, erlotinib… |
| Asthma | 4 | 27,220 | dexamethasone, prednisone, methylprednisolone, prednisolone |

A caveat on this table: the long tail of 4-drug diseases all sit at exactly 27,220 contrasts
because they share the *same four corticosteroids*. Steroid breadth inflates any
steroid-responsive condition. **Rheumatoid arthritis and ankylosing spondylitis are the genuine
leaders** because they add disease-specific biologics (infliximab, anakinra, etanercept) on top.

Both are immune-mediated and so squarely in NIAID's remit. They are the diseases used for
example 2.

---

# Worked example 2 — rheumatoid arthritis and ankylosing spondylitis

RA and AS were chosen by the backward analysis above. **The selection did not
work**, and why it failed is the most useful result in this report.

## Route A across all three diseases

| disease | edges | GXA-covered | any GXA info | agrees | disagrees | tested-not-sig |
|---|--:|--:|--:|--:|--:|--:|
| asthma (`MONDO:0004979`) | 123 | 4 (3%) | 11 (9%) | 2 | 0 | 7 |
| rheumatoid arthritis (`MONDO:0008383`) | 174 | 1 (**1%**) | 11 (6%) | 0 | 0 | 10 |
| ankylosing spondylitis (`MONDO:0005306`) | 163 | 14 (9%) | 46 (**28%**) | 1 | 0 | 32 |

RA scored **worse than asthma** despite being picked for maximum GXA drug coverage.

## Why the backward selection failed

RA was selected because 9 GXA-covered drugs treat it *according to DrugMechDB*. But
Translator's creative-mode answer proposes a largely disjoint drug set:

- Translator returns **159 distinct drugs** for RA; **3** are GXA-covered
  (Baricitinib, Bevacizumab, Infliximab).
- Of the 9 drugs RA was selected for, **8 never appear in Translator's answer at all** —
  methotrexate, anakinra, etanercept, prednisone, prednisolone, methylprednisolone,
  dexamethasone, doxorubicin.

The bottleneck is a **three-way intersection**, and selecting on only one term does not move
the others:

1. Translator proposes drug D for disease X *with a gene intermediate*
2. GXA tests D in human
3. gene G is *significantly* DE in those contrasts

DrugMechDB indications describe established therapy; Translator's inferred `treats` answers
skew toward mechanism-adjacent and investigational compounds. They are different distributions.

## Splitting "no data" from "tested and null"

GXA stores only *significant* DE results, so an absent record was ambiguous. Separating the two
roughly triples the interpretable fraction (asthma 3%→9%, AS 9%→28%):

- `tested_not_significant` — GXA tests the drug, but the gene is never significantly DE.
  **Evidence against** the edge.
- `no_drug_data` — the drug is absent from GXA. No information either way.

Two cases show why this matters, and both are the same underlying limitation:

| edge | GXA | reading |
|---|---|---|
| Baricitinib → JAK1/2/3, TYK2 | drug covered; JAKs never DE | correct mechanism, invisible to expression |
| Infliximab → TNF | 12,759 drug contrasts, 3 mention TNF, direction *opposes* | correct mechanism, **disagrees** on expression |

Baricitinib inhibits JAK *kinase activity*; infliximab neutralises TNF *protein*. Neither
changes the transcript in the asserted direction. **Expression is not activity** — this was
listed as a caveat in the plan and is now the dominant empirical finding. A `disagrees` verdict
against an `activity`-aspect edge should not be read as contradicting the biology.

Across all three diseases there were **zero genuine `disagrees`** and only 3 `agrees`.

## Route B discovery is far healthier

RA + AS human GEO sequencing series, discovered through NDE, checked against GEO for a deposited
matrix:

| | series |
|---|--:|
| candidates (union of MONDO + free text) | 453 |
| **re-analyzable** | **182 (40%)** |
| — raw counts | 117 |
| — normalized matrix | 65 |
| with NDE sample-level arm labels | **182 (100%)** |

Raw-counts series: median 24 samples, 65 with ≥20, 22 with ≥50.

⚠️ MONDO-only discovery would have found half of this. Only **41.9%** of NDE's GEO datasets carry
any `healthCondition`, so the union with free text is necessary (RA 213→422, AS 15→34).

⚠️ Some hits are single-cell (GSE109449, GSE235508) and need different handling than bulk DE.

## Two bugs worth recording

- **Biologics were being used as gene intermediates.** Infliximab is `biolink:Protein` as well as
  `biolink:Drug`, producing `Etoricoxib → Infliximab → rheumatoid arthritis`. Drug-ness now wins;
  regression test in `tests/test_path_extraction.py`.
- **ARS `result_count` is sometimes a string**, which crashed the poll loop with a TypeError and
  killed the first RA run.

## Where this leaves the project

Route A is a **cheap, precise, very low-recall** filter, and much of what it returns is negative
evidence about transcription rather than support for a mechanism. It is worth keeping as a
first pass, but it cannot carry the project.

Route B has 182 re-analyzable series for two diseases alone, all with usable arm labels. That is
where the remaining effort belongs.

## Route B validated end-to-end on GSE89408

RA synovial biopsies, 218 NDE `Sample` records, `GSE89408_GEO_count_matrix_rename.txt.gz`.
Arms split 152 RA / 28 healthy; 16,284 genes tested in ~1 s.

| gene | logFC | adj. p | role |
|---|--:|--:|---|
| MMP1 | +4.98 | 2.1e-13 | classic RA synovial protease |
| MMP3 | +4.97 | 3.2e-13 | " |
| POSTN | +2.88 | 3.3e-13 | fibroblast activation |
| MS4A1 | +2.67 | 3.4e-09 | B-cell infiltration |
| CXCL13 | +2.32 | 2.0e-08 | ectopic lymphoid follicles |
| IL6 | +1.66 | 1.2e-07 | drug target |
| CD3E | +1.41 | 1.4e-08 | T-cell infiltration |
| TNF | +1.08 | 7.4e-11 | drug target |

All eight known markers up in RA. The pipeline recovers the expected biology.

⚠️ **Arm provenance matters and is reported.** Depositors name matrix columns arbitrarily: this
series uses `normal_tissue_1` / `RA_tissue_148`, while the GEO titles are `healthy tissue 2` —
so neither the GSM accession nor the title matches a column. Matching is layered
(accession → title → column-name regex) and `DEResult.arm_source` records which one was used.
Here it fell through to `matrix_columns`, which loses the link back to NDE sample records; that
is a weaker provenance chain and should be flagged in any result built on it.

---

# Worked example 3 — Route B on three drug-treatment contrasts

Selected on **contrast structure**, not sample count. GSE89408 has 218 samples but contrasts RA
against healthy, which tests no drug→gene edge. Of 87 profiled raw-counts series, only a handful
have a genuine drug arm *and* a matched control arm.

| series | contrast | n | design | genes | sig FDR<0.05 |
|---|---|--:|---|--:|--:|
| GSE97165 | RA synovium, post vs pre triple DMARD | 19+19 | **paired** | 16,452 | 464 (3%) |
| GSE148395 | RA fibroblasts, JQ1 vs DMSO | 12+12 | unpaired | 13,840 | 6,546 (**47%**) |
| GSE141646 | AS whole blood, post vs pre TNF inhibitor | 22+22 | **paired** | 14,089 | 282 (2%) |

## Three pipeline bugs, all found by running only three datasets

1. **Paired designs were analysed unpaired.** GSE97165 and GSE141646 sample the *same patient*
   before and after treatment; between-patient variation in synovium and whole blood dwarfs the
   drug effect. Unpaired gave GSE97165 **0/7** markers and GSE141646 **0** significant genes.
   Paired gives 2/7 and 282. `paired_moderated_ttest` applies the same empirical-Bayes shrinkage
   to within-subject differences.
2. **Arm regexes missed separator variants.** GSE148395 columns use both `ST1359_JQ` and
   `ST1387-JQ`; `r"_JQ"` matched 4 of 12, so the run compared **4 JQ vs 12 DMSO** — unbalanced
   *and* confounded by the IL-1β sub-arm. `run_de` now warns when arms are ≥2:1 unbalanced.
3. **Ensembl-indexed matrices silently defeated gene lookup.** Every GSE141646 marker read
   "not in matrix". Now detected and mapped via MyGene.info (13,256 ids).

## The finding that changes the design

**A generic "expect the drug to decrease the gene" rule is wrong.** For JQ1 the correct
pharmacodynamic readout is *up*:

| gene | logFC | adj p | |
|---|--:|--:|---|
| **HEXIM1** | **+1.77** | 1.8e-14 | canonical BET-inhibitor biomarker — P-TEFb release |
| **BRD2** | **+1.62** | 1.8e-12 | known compensatory BET upregulation |
| IL21R | −3.14 | 1.8e-14 | |
| MMP1 | −1.70 | 4.9e-02 | |
| MYC | +1.01 | 1.6e-06 | JQ1's MYC suppression is cell-type specific, not universal |

The drug demonstrably worked. The panel scored 1/7 only because the "expect decreased"
hypothesis was naive. **Direction must come from the Translator edge's own
`object_direction_qualifier`, per gene** — never from a blanket expectation.

## GSE141646 looks like a technical artifact

Its top hits are *all* small non-coding RNA — SNORD17, SNORA74B, RNU4-2, RNY4, SCARNA1,
SNORA74D, RNU4-1, SNORA73B. That is a library-prep / RNA-degradation signature, not TNF biology,
and none of the TNF-pathway markers move. 282 "significant" genes dominated by housekeeping
non-coding RNA should be treated as a failed run, not a negative result. Worth an automatic
biotype check before trusting any series.

## Power differs by an order of magnitude

In-vitro drug perturbation: **47%** of genes significant. Clinical pre/post: **2–3%**.
For testing drug→gene edges, in-vitro perturbation series are the productive substrate;
clinical pre/post designs are underpowered for anything but the largest effects.

---

# Worked example 4 — Route D: activity data

Routes A and B ask expression data whether a drug changes its target's abundance.
Examples 1–3 established that this is the wrong question for most Translator
edges, which assert changes in **activity**. Route D asks the question the
assertions actually make, of assays that measure it: PubChem BioAssay and ChEMBL.

Run over the same three ARS creative-mode answer sets, deduped to distinct
drug→gene edges. Everything below is measured, 2026-09-01.

## Coverage — 52% of edges have a directly measured compound–target result

```
          edges  measured  mech  bind  inact  nottest  nodata  noid
asthma      220    131(60%)   27   129      2       34       8    20
RA          253    137(54%)   30   131      5       72       8     7
AS          237     99(42%)    6    86     13      116       3    13
TOTAL       710    367(52%)   63   346     20      222      19    40
```

`measured` = a curated mechanism, a measured potency, or a recorded Inactive
outcome for that exact compound against that exact target.

Set against Route A on the same edges — asthma 3%, RA 1%, AS 9% — this is a
**one-to-two order of magnitude** difference in how often the data can say
anything at all. 244 edges carry a potency value; 232 carry a pChEMBL ≥ 6.

The join is exact end to end. Translator emits `NCBIGene:7124`; PubChem's
`Target GeneID` column is `7124`. Node Normalizer supplies the compound's CID
and ChEMBL id from its clique (461/524 and 470/524 of the distinct drugs
respectively; the 32 with neither are mostly biologics). **No text matching
anywhere** — the failure mode that limits Routes A and B does not arise.

## 20 true negatives — the thing Route A structurally could not produce

GXA stores only significant results, so a missing record is uninformative.
PubChem records **Inactive** outcomes, so "this compound was tested against this
target and did nothing" is a real observation. Twenty edges are contradicted
this way, e.g. `Baicalein → CFTR` (4 inactive) against an asserted *increased*
activity, and five separate compounds recorded inactive against TNF in the AS
answer set.

## The finding: Translator's directions and curated mechanism barely intersect

Of 710 edges, **330 carry a direction qualifier and 63 have a curated ChEMBL
mechanism — but only 1 has both** (crofelemer → CFTR, INHIBITOR vs. *decreased*:
agrees). By drug the disjointness is just as sharp: 330 drugs vs. 50 drugs,
overlapping on that same one.

This is not an artifact of the matcher. It follows from provenance:

| | edges with a direction | edges without |
|---|---|---|
| supplying KPs | bindingdb 334, pharos 258, ctd 53, gtopdb 18 | dgidb 737, drugbank 444, semmeddb 296, drugcentral 235 |
| drugs with *any* curated ChEMBL mechanism | **12 / 302 (4%)** | **134 / 169 (79%)** |

Translator's directional qualifiers come from **screening databases**, whose
compounds are overwhelmingly research chemicals; curated mechanism-of-action
exists almost exclusively for **approved drugs**, which reach Translator through
DrugBank/DrugCentral/DGIdb — and those emit a bare `biolink:affects` with no
qualifier at all. The two evidence types are sourced from opposite ends of the
pharmacology pipeline.

Two consequences:

1. **The direction-agreement test is nearly untestable on this data** — n=1. That
   is a negative result about qualifier coverage in Translator, not about
   activity data. Reporting an agreement rate over one edge would be dishonest.
2. **Route D's real value is the inverse**: for the 63 mechanism edges it
   *supplies* the action type Translator omits. Every one is textbook-correct
   and 59/63 are max_phase 4 (approved) — baricitinib→JAK1/JAK2 INHIBITOR,
   infliximab and certolizumab→TNF INHIBITOR, celecoxib/etoricoxib/rofecoxib→
   PTGS2 INHIBITOR, triamcinolone/betamethasone/flunisolide→NR3C1 AGONIST.
   Route A scored 1/174 on the RA set; Route D recovers RA's actual pharmacology.

Action types across all three: INHIBITOR 37, AGONIST 16, ANTAGONIST 8,
BLOCKER 1, OPENER 1.

## Not in NDE

NDE's `includedInDataCatalog` facet has LINCS (424) and ReframeDB (408) as the
only activity-adjacent sources — no ChEMBL, PubChem BioAssay or BindingDB. Route
D therefore reaches outside NDE entirely. Given that activity is the modality
Translator's edges are actually about, that is a real coverage gap for
mechanism-of-action work, and worth reporting to the NDE team as such.

## Four defects found by running it

Each broke a whole class of edges, and each is pinned in `tests/test_activity.py`.

1. **Salt-form mechanisms.** Imatinib (`CHEMBL941`) has *no* mechanism rows of
   its own; all four are filed under the mesylate `CHEMBL1642`. Querying only
   `molecule_chembl_id` downgraded the canonical ABL1 inhibitor to a mere binding
   observation. Fixed by unioning with `parent_molecule_chembl_id`.
2. **Protein-family targets.** Aspirin's mechanism target is `CHEMBL2094253`
   "Cyclooxygenase", a PROTEIN FAMILY, not the PTGS2 single protein. Target-id
   equality missed it; matching on the target's component UniProt accessions
   catches it.
3. **The gene-keyed PubChem endpoint does not scale.** `/gene/geneid/{id}/concise`
   returns 337 MB for CFTR, and for DRD2 a 436 MB body that arrives **truncated
   and unparseable** — so the most heavily screened targets, the interesting ones,
   silently returned nothing and were scored `not_tested`. The compound-keyed
   `/compound/cid/{cid}/assaysummary` view carries a `Target GeneID` column, so
   the join stays an exact integer match at ~300 KB per compound. Spot-checked
   identical: imatinib/ABL1 gives 154 Active + 10 Inactive either way. Cache for
   the whole run fell from 521 MB (9 genes) to 6.2 MB.
4. **A fetch failure was being laundered into evidence.** With (3) in place a
   PubChem error yielded zero rows and a `not_tested` verdict — an error
   presented as a measured negative. Failures now get their own `fetch_failed`
   verdict and a `pubchem_error` flag, and are excluded from denominators.

One hypothesis checked and **rejected**: BindingDB and Pharos directions are
~50/50 increased/decreased, which looked like non-committal hedging. It is not —
of 294 pairs, only 2 carry both directions. The claims are genuine.

## Reproduce

```bash
.venv/bin/python scripts/run_route_d.py data/ars/<pk>/paths.json
PYTHONPATH=src .venv/bin/python -m pytest tests/test_activity.py -q
```
