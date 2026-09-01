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
