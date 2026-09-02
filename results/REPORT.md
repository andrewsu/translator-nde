# Grounding Translator knowledge-graph paths in primary data

**Results report** · all figures measured 2026-09-01 · code and artifacts in this repository.

A Translator `drug → gene → disease` path is a mechanistic hypothesis. This report asks, for each
`drug → gene` hop, whether primary data exists that tests it — and answers that question four
different ways, because the obvious way turns out to be the wrong question.

| route | question it asks | data source |
|---|---|---|
| **A** | does the drug change the gene's **expression**? | NDE `@type:Inference` — Gene Expression Atlas DE contrasts (staging only) |
| **B** | can we **compute** that ourselves from raw data? | GEO count matrices, with treated/control arms from NDE `@type:Sample` |
| **C** | what **other** compounds move this gene? | the same GXA contrasts, queried gene-first |
| **D** | does **activity** data confirm the asserted inhibition? | PubChem BioAssay and ChEMBL (neither indexed by NDE) |

Routes A, B and D are built and run; **Route C is designed but not yet run**, so it has no results
in this report.

Routes A and B interrogate expression. Route D interrogates activity. That distinction is the
central result: **Translator's drug→gene edges overwhelmingly assert changes in *activity*, while
the expression atlases measure *abundance*** — and a kinase inhibitor or a neutralising antibody
does not move its target's transcript.

## Results at a glance

Across three NIAID-relevant diseases — asthma, rheumatoid arthritis, ankylosing spondylitis —
drawn from live ARS creative-mode `treats` queries:

| | Route A (expression) | Route D (activity) |
|---|--:|--:|
| drug→gene edges evaluated | 460 | 710 |
| edges with a usable measurement | **4 (0.9%)** | **367 (52%)** |
| agreements with the asserted direction | 2 | 1 |
| contradictions | 0 | 20 (measured-inactive) |

The denominators differ because Route A skips 252 edges whose "drug" is an IUPAC string or a
synthetic peptide — free-text search cannot match those — whereas Route D joins on identifiers and
can attempt every edge carrying an `NCBIGene` id. Restricting Route D to the same 460 edges would
not change its conclusion; the gap is more than an order of magnitude either way.

Route B is not scored per edge; it is a *capability* result — 182 of 453 RA/AS GEO series (40%)
are re-analyzable from NDE metadata alone, and the pipeline reproduces the expected biology
wherever the underlying series is sound (one of the four series run turned out to be a technical
artifact, which the report treats as a failed run rather than a negative result).

## How to read this

Examples 1–3 are chronological and each one changed the design, so they are worth reading in
order. Example 4 is where the project's most useful evidence comes from.

- **Method** — how a GXA contrast is matched to a Translator edge. Needed to read example 1.
- **Example 1 — asthma.** Route A end to end on a real answer set.
- **Working backward from GXA.** Which drugs GXA actually contains, and which diseases that implies.
- **Example 2 — RA and AS.** The backward selection tested, and why it fails.
- **Example 3 — Route B.** Reanalysis of three drug-treatment series.
- **Example 4 — Route D.** Activity data, and the highest-coverage result in the project.
- **Conclusions.**

---

# Method — matching a drug→gene edge to a GXA contrast

This governs every Route A number in the report.

NDE `@type:Inference` records carry no chemical identifier — the compound appears only as free text
in `variableMeasured.value` (test arm) and `measurementDenominator.value` (reference arm). Three
filters are needed before a contrast counts as evidence, and each is derived from a case that
would otherwise be scored wrongly.

**1. Species.** `species.identifier:9606`. Without it, plant and mouse tool-compound designs
dominate.

**2. The drug must be what *differs* between the arms.** Record
`gxa_e_geod_54049_g1_g5_at1g03790` contrasts `'pBeaconRFP_GR::bZIP1, …, 10 uM dexamethasone in
ethanol'` vs `'empty vector, …, 10 uM dexamethasone in ethanol'` — dexamethasone is a GR *inducer*
given to both arms and the real variable is genotype. Eight human records do the same thing with
disease: `'prednisolone 20 milligram per day, polymyalgia rheumatica'` vs
`'prednisolone 20 milligram per day, normal'`.

The obvious implementation, `AND NOT measurementDenominator.value:<drug>`, is wrong. TG-GATEs
studies (`E-CURD-*`) carry a `cohort` factor holding the compound name on **every** sample, so the
aspirin study's vehicle control reads `'aspirin, liver, 15 day'` while its test arm reads
`'aspirin, aspirin 150 milligram per kilogram, liver, 15 day'` — no aspirin was administered to the
control, but the name is there. Excluding on the name deletes the entire study: **213 aspirin, 222
theophylline, 2,204 rifampicin and 12,528 allopurinol contrasts**, all of them genuine
dose-vs-vehicle designs.

So compare *factors*, not names. `gxa.drug_is_the_variable()` requires the drug to occupy the
variable position in the test arm (filter 3 below) **and** that exact factor not to be a reference
factor too. Where a test arm carries both a bare label and a dosed one, the dosed factor is
preferred, since it is the evidence the compound was actually given. This keeps the TG-GATEs
designs and still excludes both the *Arabidopsis* and the polymyalgia cases, where the drug factor
is identical on each side.

**3. The drug must occupy the variable position.** Elasticsearch matches a synonym *anywhere* in the
test-arm text, which admits contrasts where the drug name is incidental:

| test arm | matched | why it is not evidence |
|---|---|---|
| `MITF-RFP-HA overexpression` | rifampicin, via "RFP" | **RFP is red fluorescent protein** |
| `no response to infliximab treatment, Crohn's disease` | nitric oxide, via "NO" | **"no" is the English word** |
| `A/CA/04/2009 Influenza virus, 30 hour` | calcium, via "CA" | **CA is California** |
| `before first infliximab treatment, …, Crohn's disease` vs `control` | infliximab | the variable is **disease**; two of three contrasts precede treatment |
| `Tr1 cell clone, 6 hour, stimulated with monoclonal antibodies to CD3 and CD28` | "Antibodies" | anti-CD3/CD28 stimulation, not a therapeutic |

There is no structured factor type to filter on — `variableMeasured.constraintProperty` is
`['schema:healthCondition', 'nde:sample']` on genuine and spurious contrasts alike, checked across
eight records. But `variableMeasured.value` is a **comma-separated list of factor values**, and a
real compound factor is the compound name alone or the name followed by a dose:

```
GENUINE   'SK-BR-3, metformin 4 millimolar'          'doxorubicin 0.6 microgram per milliliter'
          'differentiated brown adiopcyte, cyclic AMP'
SPURIOUS  'MITF-RFP-HA overexpression'               'A/CA/04/2009 Influenza virus, 30 hour'
          'before first infliximab treatment, no response to infliximab treatment, Crohn's disease'
```

`gxa.factor_supports_drug()` splits on commas and keeps a contrast only if some synonym **is** a
factor, or **leads** one with nothing but a dose after it. Synonyms under 3 characters and a
stop-word list (`no`, `none`, `control`, `vehicle`, `dmso`, …) are excluded outright. The dose
requirement is load-bearing: `no response to infliximab treatment` *begins* with "no", so a
leading-token rule alone would still admit nitric oxide. The same verification is applied to
`drug_contrast_count`, which text-verifies a 200-record sample rather than trusting the raw Lucene
count. All five collisions above are pinned in `tests/test_gxa_fixtures.py`.

Name Resolver synonym lists are themselves a hazard: azacitidine returns `AZC`, `ac 5`, `5-AC`,
`5 aza`, `AZA-CR`, which match **5-aza-2-deoxycytidine** contrasts — that is decitabine, a
different drug.

---

# Worked example 1 — asthma

**Date:** 2026-09-01 · **Disease:** asthma (`MONDO:0004979`)
**ARS pk:** `5b656c0f-b7da-4db4-ba1f-d3a794b422d4`
**NDE:** staging (`api-staging.data.niaid.nih.gov/v1`)
**Artifacts:** `data/ars/5b656c0f-.../paths.json`, `.../route_a.json`, `results/route_a_asthma.log`

## TL;DR

The pipeline works end to end, and Route A is **precise but very narrow**: of 123 evaluable
`drug → gene` edges from a real Translator creative-mode answer, GXA had a matching
differential-expression contrast — one in which the drug is genuinely the experimental variable —
for **1 (1%)**.

The bottleneck is **GXA's drug coverage, not the matching logic** — Expression Atlas contains no
contrasts at all for the drugs asthma answers are actually built from. Precomputed expression
therefore cannot carry the project on its own, which is what motivates Route B (computing the
contrasts ourselves) in example 3.

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
`H-D-Phe-His-Leu-Leu-Arg-…`), leaving 123 evaluated in 109 s.

| Verdict | Edges |
|---|--:|
| no GXA data for the drug | 116 |
| drug tested, gene never significantly DE | 6 |
| agrees | 1 |
| disagrees | 0 |
| **covered by ≥1 contrast** | **1 / 123 (1%)** |

| Drug → gene | contrasts | agree | median log2FC | verdict | GXA experiment | GEO |
|---|--:|--:|--:|---|---|---|
| Cyclic AMP → PPARGC1A | 2 | 2 | 1.85 | agrees | [E-MTAB-2602](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2602/Results) | — (ENA ERR522174–86) |

The contrast is `'cyclic AMP' vs 'none'` in differentiated brown and white adipocyte. GEO series
for GXA experiments are resolved by looking up each contrast's GSM accessions in NDE production
`@type:Sample` records, not inferred from the `E-GEOD-` naming; E-MTAB-2602 is ArrayExpress-native
and has ENA run accessions instead.

## Why coverage is so low

Covering an edge takes two independent things to go right:

1. GXA has to have **tested that drug in human at all**, and
2. the **gene** has to be significantly differentially expressed in those contrasts.

Almost all the loss happens at step 1. Of the **77 distinct drugs** behind asthma's 123 evaluable
edges, only **6 have even one verified human GXA contrast** (8%): prednisolone, prednisone,
sorafenib, rifampicin, dupilumab and cyclic AMP. Whatever the matcher does at step 2, it can only
ever work on those six drugs. The same measurement for the other two diseases in this report gives
2 of 79 for RA (3%) and 18 of 140 for AS (13%).

So the question becomes: *why* does GXA lack these drugs? The table below spot-checks twelve
compounds to characterise the shape of the gap. It is a hand-assembled list, not a sample — five
are drugs Translator actually proposed for asthma, four are first-line asthma therapeutics
Translator did *not* propose (included to show the gap is in GXA rather than in Translator's
answer), and three are well-studied compounds included as reference points.

Columns: **any species** and **human** are counts of GXA contrasts naming the drug in the test
arm, with and without `species.identifier:9606`. **verified** is how many of up to 200 sampled
human contrasts survive the arm checks above — i.e. how many actually make the drug the variable
rather than merely mentioning it. The 200-record cap means a drug whose genuine contrasts are rare
among many confounded ones can sample to zero, so `verified` is a floor.

| Drug | in asthma answer? | any species | human | verified / sampled |
|---|---|--:|--:|--:|
| dexamethasone | no — reference point | 57,679 | 23,314 | 200/200 |
| metformin | no — reference point | 2,818 | 928 | 200/200 |
| rifampicin | **yes** | 2,312 | 108 | 108/108 |
| prednisolone | **yes** | 481 | 232 | 200/200 |
| theophylline | no — first-line asthma | 222 | **0** | – |
| aspirin | **yes** | 213 | **0** | – |
| imatinib | no — reference point | 80 | **0** | – |
| albuterol / salbutamol | no — first-line asthma | 6 | 6 | **0/6** |
| **budesonide** | **yes** | **0** | 0 | – |
| **formoterol** | **yes** | **0** | 0 | – |
| **fluticasone** | no — first-line asthma | **0** | 0 | – |
| **montelukast** | no — first-line asthma | **0** | 0 | – |

Three things follow.

**The inhaled therapeutics are simply absent.** Budesonide, formoterol, fluticasone and
montelukast have **zero** contrasts in any species. Two of those are drugs Translator proposed for
asthma, so this is a gap in GXA, not a quirk of the answer set. Albuterol's six human contrasts
look like coverage but none survive the factor-position check, so it belongs with the zeros.

**Human coverage is the binding constraint, not coverage in general.** Aspirin (213), theophylline
(222) and imatinib (80) all have contrasts — every one of them in **rat**. Aspirin's come from the
TG-GATEs liver and kidney studies `E-CURD-59` and `E-CURD-50`, which are clean 45/150/450 mg/kg
dose series against untreated controls; they are simply the wrong species for a human drug→gene
edge. GXA is a curated re-analysis of whatever was deposited in ArrayExpress and GEO, not a
drug-perturbation atlas, and its compound coverage tracks use as a *toxicology or laboratory model*
far more than clinical importance in humans.

**A worked case: GSE162120.** The DISARM trial (`PRJNA680616`) is 118 human bronchial-brush
RNA-seq samples from COPD patients randomised to formoterol, formoterol/**budesonide**, or
salmeterol/**fluticasone**, bronchoscoped before and after 12 weeks. It is exactly the experiment
the four zero rows above call for, and it is **not in GXA**: of GXA's 4,562 experiments, none
mentions GSE162120, DISARM, budesonide, formoterol, fluticasone or salmeterol. GXA's entire
human drug-perturbation universe is **363 experiments** (860 of 4,562 carry a `compound` factor;
363 of those are human), so a 2021 respiratory trial not being among them is unremarkable — which
is the point. Route A's blind spot is a selection gap in GXA, not a defect in the matching.

The same series is fully present in **NDE production**: the `Dataset` record, all 118 `Sample`
records, and a deposited counts matrix (`GSE162120_gene_NumReads.txt.gz`). Paired pre/post,
three arms, subject IDs for pairing — a Route B dataset of exactly the shape example 3 wants.

**Dexamethasone is the outlier that makes the atlas look better than it is.** Its 23,314 human
contrasts are two orders of magnitude above the next drug Translator proposed, and they arise
because dexamethasone is a workhorse cell-culture reagent. Reasoning from dexamethasone to "GXA
covers drugs" is the mistake this table is here to prevent.

## Positive control

Route A does work where GXA has data, and reproduces textbook pharmacology:

| Drug → gene | contrasts | agree | median log2FC | min adj-p | experiments |
|---|--:|--:|--:|--:|--:|
| dexamethasone → FKBP5 | 24 | 24 | 3.25 | 1.0e-115 | 9 |
| dexamethasone → TSC22D3 (GILZ) | 26 | 26 | 2.55 | 0.0 | 11 |
| dexamethasone → NR3C1 | 4 | 2 (2 disagree) | 0.25 | 9.5e-22 | 3 |

FKBP5 and TSC22D3 are the canonical glucocorticoid-induced genes. NR3C1 coming back **ambiguous**
is correct, not a failure — glucocorticoid-receptor autoregulation is genuinely bidirectional and
context-dependent.

## Conclusions

1. **The bridge is real but the precomputed route is thin.** 1% coverage on a real disease query
   is too low to be useful alone.
2. **The data exists, just not precomputed.** Asthma's actual drugs are absent from GXA but
   present in GEO — production NDE has 5,253 dexamethasone datasets (1,302 GEO), and
   drug-anchored search plus the 8.7M `Sample` records should reach budesonide and formoterol.
   That is the Route B hypothesis, tested in example 3.
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

Example 1 asked "here are Translator's answers, does GXA cover them?" and got 1%. This asks the
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

These are raw text matches and so an **upper bound**: the factor-position check reduces them
wherever a drug's name or a short synonym also occurs incidentally (albuterol's six human
contrasts, for instance, all fail it). What the table is for — which drugs GXA lacks entirely — is
unaffected.

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

The reference-arm column is a useful flag: `doxycycline` topping the list is an artefact — it is
the Tet-on induction agent in a large number of experiments, not a tested therapeutic — and
`methotrexate` appears in *more* reference arms than test arms. But it is only a flag. A high
reference-arm count can mean the compound is held constant across both arms (the *Arabidopsis*
case) **or** merely that a cohort factor names it (the TG-GATEs case), and those need opposite
treatment. Which one applies is decided per contrast by comparing factor lists, not by this column.

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

RA and AS were chosen by the backward analysis above, on the expectation that picking diseases
for maximum GXA drug coverage would raise Route A's hit rate. **It did not**, and why it failed
says more about the bridge than the coverage numbers do.

## Route A across all three diseases

| disease | edges | GXA-covered | any GXA info | agrees | disagrees | tested-not-sig |
|---|--:|--:|--:|--:|--:|--:|
| asthma (`MONDO:0004979`) | 123 | 1 (1%) | 7 (6%) | 1 | 0 | 6 |
| rheumatoid arthritis (`MONDO:0008383`) | 174 | **0 (0%)** | 6 (3%) | 0 | 0 | 6 |
| ankylosing spondylitis (`MONDO:0005306`) | 163 | 3 (2%) | 24 (**15%**) | 1 | 0 | 21 |

RA scored **worse than asthma**: zero covered edges, from the disease selected precisely because
it had the most GXA-covered drugs.

That leaves **4 covered edges out of 460 (0.9%)** across all three diseases. Every one is a
genuine compound-vs-vehicle design:

| edge | contrast | verdict |
|---|---|---|
| Cyclic AMP → PPARGC1A | `cyclic AMP` vs `none`, brown/white adipocyte (E-MTAB-2602) | agrees |
| Doxorubicin → C3 | `doxorubicin 0.6 µg/mL` vs `none` (E-GEOD-46493, E-MTAB-6045, E-MTAB-9362) | ambiguous |
| Cisplatin → C3 | `Cisplatin` vs `None` (E-MTAB-3645) | agrees |
| Metformin → TNF | `metformin 4 millimolar` vs `none` (E-MTAB-7737) | ambiguous |

Two agreements, zero disagreements. One borderline exclusion worth naming: dinoprostone → TNF is
dropped because its factor is `PGE2-maturation`, a maturation protocol rather than a dosed PGE2
arm. That is defensible but shows the shape of false negative the rule produces.

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
multiplies the interpretable fraction several-fold (asthma 1%→6%, AS 2%→15%):

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

Across all three diseases there were **zero `disagrees`** and only 2 `agrees`.

## Route B discovery is far healthier

RA + AS human GEO sequencing series, discovered through NDE, checked against GEO for a deposited
matrix:

| | series |
|---|--:|
| candidates (union of MONDO + free text) | 453 |
| **re-analyzable** | **182 (40%)** |
| — raw counts | 117 |
| — normalized matrix | 65 |
| with NDE `Sample` records | **182 (100%)** |

⚠️ That last row counts series for which NDE holds per-sample records — **not** series whose
records carry a usable treatment label. Those are different claims, and GSE162120 separates them:
NDE structures its timepoint (`temporalCoverage.duration`, 62 pre / 56 post) and its subject id is
in the sample `name` (`DIS320001-V3`), but the **treatment arm is not a structured field**. There
is no `additionalProperty` for it, and `variableMeasured` lists only the characteristic *names*
(`subject`, `smoking_status`) without values. The arm survives only inside the free-text
`sampleProcess` blob, wedged between the file description and the extraction protocol:

```
… Supplementary_files_format_and_content: Matrix table with number of reads
for every gene and every sample FOR/BUD total RNA RNA was extracted from …
```

A regex over that field does recover all three arms exactly — FOR 41, FOR/BUD 36, SAL/FLU 40, plus
one unlabelled — matching the GEO series matrix. So the labels are retrievable, but by text
scraping rather than field lookup, and Route B's arm assignment should not be assumed to work from
NDE metadata alone until it is checked per series.

Raw-counts series: median 24 samples, 65 with ≥20, 22 with ≥50.

⚠️ MONDO-only discovery would have found half of this. Only **41.9%** of NDE's GEO datasets carry
any `healthCondition`, so the union with free text is necessary (RA 213→422, AS 15→34).

⚠️ Some hits are single-cell (GSE109449, GSE235508) and need different handling than bulk DE.

## Two requirements on the extractor

- **A biologic must not be used as the gene intermediate.** Infliximab carries `biolink:Protein`
  as well as `biolink:Drug`, so a naive category test yields paths like
  `Etoricoxib → Infliximab → rheumatoid arthritis`. Drug-ness wins; regression test in
  `tests/test_path_extraction.py`.
- **ARS `result_count` is sometimes a string**, so the poll loop must coerce it rather than
  compare it numerically.

## Where this leaves Route A

Route A is a **cheap, precise, very low-recall** filter, and much of what it returns is negative
evidence about transcription rather than support for a mechanism. It is worth keeping as a first
pass, but on this evidence it cannot carry the project.

Route B has 182 re-analyzable series for these two diseases alone, all with usable arm labels, so
that is where example 3 goes next. The deeper problem — that expression is the wrong assay for an
activity claim — is not fixed by computing the contrasts ourselves, and is what eventually
motivates Route D.

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

### Relating it back to Translator

Route B as run in example 3 tests `drug → gene` edges, and GSE89408 cannot: it contrasts RA
synovium against healthy, so no drug varies. It can test the *other* hop of the same path. If
Translator proposes `drug → gene → disease` for RA, its gene intermediates ought to be
differentially expressed in RA tissue.

RA's answer set yields **95 distinct gene intermediates**, 68 of them in this matrix. The join is
`scripts/join_route_b_translator.py`, scored against size-matched random gene sets drawn from the
same matrix, because "most of them are DE" means nothing without knowing how many genes are DE
overall — and here that number is enormous.

| threshold | Translator | background | permutation p |
|---|--:|--:|--:|
| adj p < 0.05 | 84% | 80% | 0.27 |
| adj p < 0.05 & \|logFC\| > 0.5 | 76% | 64% | 0.019 |
| adj p < 0.05 & \|logFC\| > 1.0 | 41% | 34% | 0.15 |
| adj p < 0.05 & \|logFC\| > 1.5 | 19% | 12% | 0.056 |
| adj p < 0.05 & \|logFC\| > 2.0 | 6% | 3% | 0.12 |

Median \|logFC\| is 0.92 for Translator's genes against 0.73 for the transcriptome; of the 57
significant ones, 43 are up and 14 down.

**Read this as a negative result.** The shift is in the expected direction and it is small. The
p-values are non-monotonic across five thresholds and two fall nominally below 0.05, which is
about what five tries buys you by chance — picking the 0.019 and reporting that alone would be
cherry-picking. On this contrast, Translator's RA gene intermediates are only marginally
distinguishable from randomly chosen expressed genes.

Three reasons not to read much into it either way:

1. **The contrast is over-powered for this question.** 152 vs 28 samples puts **80% of the
   transcriptome** below FDR 0.05, so DE membership carries almost no information and the test has
   little room to discriminate.
2. **The background is the wrong null.** Translator's intermediates are drug targets — biased
   toward well-studied, well-expressed, immune-related genes. A fair comparison needs an
   expression-matched or druggable-genome background, not all 16,284 genes.
3. **This tests `gene → disease`, not the hop the project is about.** A gene being DE in RA tissue
   says nothing about whether a given drug acts on it, which is the claim Routes A and D address.

The useful conclusion is methodological: **a disease-vs-healthy contrast is the wrong instrument
for validating mechanistic paths**, and the 40% of RA/AS series that are re-analyzable are worth
far more when they carry a drug arm. That is what example 3 goes after.

### Do Route B's perturbagens appear in Translator at all?

Before scoring a drug arm against Translator edges, the prior question is whether Translator ever
proposed those drugs. Checked against the three answer sets:

| Route B perturbagen | series | in Translator's answer? |
|---|---|---|
| **JQ1** | GSE148395 | **absent from all three** |
| methotrexate | GSE97165 | absent from RA; in AS → TNF (SemMedDB) |
| sulfasalazine | GSE97165 | absent from all three |
| hydroxychloroquine | GSE97165 | absent from RA; in AS → TNF |
| infliximab | GSE141646 | present in all three; **6 gene edges** in AS |
| **budesonide** | GSE162120 | **present in asthma, 10 gene edges** |
| **formoterol** | GSE162120 | **present in asthma, 2 gene edges** |
| fluticasone, salmeterol | GSE162120 | absent from all three |

Two things follow, and they cut in opposite directions.

**The best-powered Route B dataset tests a drug Translator never proposed.** GSE148395 (JQ1, 47% of
genes significant) has no counterpart edge to score against. Nor does GSE97165's triple DMARD:
methotrexate, sulfasalazine and hydroxychloroquine are the standard RA regimen and **none of them
appears in Translator's RA answer** — two show up only under AS. This is the three-way-intersection
problem from example 2 again, now biting from the data side: the series worth re-analysing and the
drugs Translator proposes are largely disjoint sets.

**GSE162120 is the exception, and it is a good one.** Both of its ICS/LABA arms — budesonide and
formoterol — are drugs Translator proposes for asthma, with twelve gene edges between them, in a
paired pre/post human design with counts deposited. That makes it the first genuinely joinable
Route B dataset in this project.

But look at what those edges assert:

```
Budesonide → CYP3A4, CYP1A2, CYP2C19, CYP2C9, NR1I2, NR3C2, PGR, CRHR1, EDN1, ANXA1
Formoterol → ADRB1, CYP2C19
```

**Five of budesonide's ten are drug metabolism** — four CYPs plus NR1I2 (PXR), the xenobiotic
sensor that regulates them. ANXA1 is a genuine glucocorticoid effector. **NR3C1, the
glucocorticoid receptor and budesonide's actual target, is not among them** — although NR3C1 *is*
reached in the same answer set by seven other corticosteroids (betamethasone, desonide,
flunisolide, prednisolone, prednisone, triamcinolone, beclomethasone). Formoterol shows the same
shape: ADRB1 rather than ADRB2, though ADRB2 is reached by four other β-agonists in the same
answer. The canonical target is present in the graph for the drug *class* and missing for these
particular members.

Route D corroborates from the activity side: budesonide→NR3C2 and →PGR are binding-confirmed
(3 and 7 active measurements), while budesonide→CYP1A2 is **measured inactive** — four inactive
outcomes, no active ones.

And every one of the twelve edges carries `direction: None, aspect: None`, sourced from DGIdb,
DrugCentral, DrugBank or SemMedDB — precisely the non-directional KP population example 4
identifies. So even on the one joinable dataset, Translator supports the question *"is this gene
differentially expressed under the drug?"* but not *"does it move the way the edge says?"*, because
no edge says.

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

## Three requirements these three datasets impose

1. **Paired designs must be analysed paired.** GSE97165 and GSE141646 sample the *same patient*
   before and after treatment, and between-patient variation in synovium and whole blood dwarfs
   the drug effect — analysed unpaired, neither series yields anything.
   `paired_moderated_ttest` applies the same empirical-Bayes shrinkage to within-subject
   differences, and recovers 464 and 282 significant genes respectively.
2. **Arm assignment must tolerate separator variants.** GSE148395 columns use both `ST1359_JQ`
   and `ST1387-JQ`, so a `r"_JQ"` pattern captures only 4 of 12 — an unbalanced comparison
   confounded by the IL-1β sub-arm. `run_de` warns when arms are ≥2:1 unbalanced.
3. **Ensembl-indexed matrices must be mapped before gene lookup.** GSE141646 is indexed by
   Ensembl gene id, so symbol lookup finds nothing at all; `looks_ensembl` detects this and maps
   via MyGene.info (13,256 ids).

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

Set against Route A on the same edges — asthma 1%, RA 0%, AS 2%, and 4/460 (0.9%) overall — this
is a **fifty-fold** difference in how often the data can say anything at all. 244 edges carry a
potency value; 232 carry a pChEMBL ≥ 6.

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
   Route A covers **zero** of the 174 RA edges; Route D recovers RA's actual pharmacology.

Action types across all three: INHIBITOR 37, AGONIST 16, ANTAGONIST 8,
BLOCKER 1, OPENER 1.

## Not in NDE

NDE's `includedInDataCatalog` facet has LINCS (424) and ReframeDB (408) as the
only activity-adjacent sources — no ChEMBL, PubChem BioAssay or BindingDB. Route
D therefore reaches outside NDE entirely. Given that activity is the modality
Translator's edges are actually about, that is a real coverage gap for
mechanism-of-action work, and worth reporting to the NDE team as such.

## Four things the join has to get right

Each governs a whole class of edges, and each is pinned in `tests/test_activity.py`.

1. **Mechanisms are often filed on the salt.** Imatinib (`CHEMBL941`) has *no*
   mechanism rows of its own; all four are under the mesylate `CHEMBL1642`.
   Querying `molecule_chembl_id` alone reduces the canonical ABL1 inhibitor to a
   mere binding observation, so the lookup unions in `parent_molecule_chembl_id`.
2. **Mechanism targets are often families or complexes.** Aspirin's target is
   `CHEMBL2094253` "Cyclooxygenase", a PROTEIN FAMILY, not the PTGS2 single
   protein. Target-id equality misses these; matching also on the target's
   component UniProt accessions catches them.
3. **PubChem must be queried per compound, not per gene.** Both routes join
   exactly, but `/gene/geneid/{id}/concise` returns 337 MB for CFTR and, for
   DRD2, a 436 MB body that arrives **truncated and unparseable** — so the most
   heavily screened targets, the interesting ones, return nothing. The
   compound-keyed `/compound/cid/{cid}/assaysummary` view carries a
   `Target GeneID` column, keeping the join an exact integer match at ~300 KB per
   compound and 6.2 MB for the whole run. Verified equivalent: imatinib/ABL1
   gives 154 Active + 10 Inactive either way.
4. **A fetch failure must not become a measured negative.** Zero rows from a
   PubChem error is not the same observation as zero rows from a compound never
   tested. Failures get their own `fetch_failed` verdict and a `pubchem_error`
   flag, and are excluded from denominators.

One hypothesis checked and **rejected**: BindingDB and Pharos directions are
~50/50 increased/decreased, which looked like non-committal hedging. It is not —
of 294 pairs, only 2 carry both directions. The claims are genuine.

## Reproduce

```bash
.venv/bin/python scripts/run_route_d.py data/ars/<pk>/paths.json
PYTHONPATH=src .venv/bin/python -m pytest tests/test_activity.py -q
```

---

# Conclusions

## 1. The bridge exists, but not where the sketch put it

The original design assumed the `drug → gene` hop could be grounded in expression data: either
precomputed (Route A) or recomputed (Route B). Measured across 460 edges from three real
Translator answers, **Route A grounds 4 of them (0.9%)**. Two agree with the asserted direction,
none contradict it. That is not a matching failure — the filters were validated against a
dexamethasone positive control that recovers FKBP5 and TSC22D3 exactly as textbook pharmacology
predicts. Expression Atlas simply does not contain the drugs Translator proposes.

## 2. The reason is a modality mismatch, not a coverage gap

Translator's drug→gene edges assert changes in **activity**; the expression atlases measure
**abundance**. Baricitinib inhibits JAK kinase activity and the JAKs are never differentially
expressed. Infliximab neutralises TNF protein and the transcript moves the *other* way. Both are
textbook-correct mechanisms that expression data cannot confirm, and computing the contrasts
ourselves does not change that — Route B reaches the same wall from the other side.

## 3. Activity data answers the question that was actually asked

Route D covers **367 of 710 edges (52%)** — fifty times Route A — because it queries assays that
measure what the edges assert. It also supplies the one evidence class the expression routes
structurally cannot: PubChem records **Inactive** outcomes, so 20 edges have a genuine measured
negative rather than mere absence.

Neither PubChem BioAssay nor ChEMBL is indexed by NDE. Since activity is the modality Translator's
edges are about, that is a substantive coverage gap in NDE for mechanism-of-action work, and the
most actionable thing in this report for the NDE team.

## 4. Translator's qualifiers and its curated mechanisms come from different worlds

Of 710 edges, 330 carry a direction qualifier and 63 have a curated ChEMBL mechanism — **1 has
both**. Directions come from screening databases whose compounds are research chemicals (4% have
any ChEMBL mechanism); curated mechanism exists for approved drugs, which arrive via
DrugBank/DrugCentral/DGIdb emitting a bare `biolink:affects` with no qualifier (79%).

So the direction-agreement test that motivated Route D is not answerable at n=1. Route D's actual
contribution is the inverse — **supplying** the action type Translator omits, correct in all 63
cases and approved-drug-backed in 59.

## 5. Route B is worth keeping, for a different job

40% of RA/AS GEO series are re-analyzable straight from NDE metadata, and the pipeline reproduces
known biology. But **none of its results has yet been joined back to a Translator drug→gene edge.**
The one join the data supports — Translator's RA gene intermediates against RA-vs-healthy synovium
in GSE89408 — tests the `gene → disease` hop instead, and returns a weak, unconvincing enrichment
over random genes. Route B's value is therefore not adjudicating individual edges: clinical
pre/post designs yield 2–3% of genes significant, too underpowered for that, while
disease-vs-healthy designs are so over-powered that 80% of the transcriptome is significant and
nothing discriminates. It is worth keeping to generate fresh contrasts from data nobody has
re-analysed, with in-vitro perturbation series (47% significant) as the productive substrate.

## What to do next

1. **Route C** — invert the GXA query: given `Drug1 –inhibits→ Gene2`, ask which *other* compounds
   move Gene2. This is the one use of GXA that plays to its structured axis (`observationAbout`
   carries symbol and Ensembl id, so no text matching), and it generates repurposing candidates
   Translator did not propose.
2. **Score Route A per `aspect`.** Only `abundance`/`expression` edges should ever be judged by
   expression data; `activity` edges should be out-of-scope, not `disagrees`.
3. **Add a biotype check to Route B**, so a run dominated by snoRNA/snRNA is reported as a failed
   run rather than a result.
4. **Join Route B to Translator on `GSE162120`.** The GSE89408 join tests the wrong hop, and the
   three example-3 series cannot be joined at all — JQ1 and sulfasalazine are absent from every
   answer set, and methotrexate and hydroxychloroquine are absent from RA. GSE162120 is the one
   dataset that overlaps: budesonide and formoterol are both proposed for asthma, in a paired
   pre/post human design. Scoring its twelve edges is a coverage test, not a direction test, since
   all twelve are unqualified.
5. **Report the ChEMBL/PubChem/BindingDB absence to the NDE team** as a coverage gap,
   alongside the GEO `characteristics_ch1` flattening seen in GSE162120.
