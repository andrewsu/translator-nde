"""Test Translator's gene->disease hop against a re-analysed GEO contrast.

Route B as run in example 3 tests drug->gene edges. GSE89408 cannot: it
contrasts RA synovium against healthy, so no drug varies. What it *can* test is
the other hop of the same path -- Translator proposes `drug -> gene -> disease`,
and if a gene is a genuine mechanistic intermediate for RA it ought to be
differentially expressed in RA tissue.

Compares the Translator gene set against a size-matched random background drawn
from the same expression matrix, because "most of them are DE" means nothing
without knowing how many genes are DE overall.
"""
import json, random, sys
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.reanalysis import download_matrix, fetch_samples, run_de

GSE = "GSE89408"
MATRIX = "GSE89408_GEO_count_matrix_rename.txt.gz"
PATHS = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "data/ars/d1f9daae-9330-4c2c-81f5-5c21ee2f44f4/paths.json")
FDR, N_PERM = 0.05, 2000

doc = json.loads(PATHS.read_text())
genes = {}
for p in doc["paths"]:
    if p.get("gene_name"):
        genes.setdefault(p["gene_name"], set()).add(p["drug_name"])
print(f"{PATHS.name}: {len(doc['paths'])} paths -> {len(genes)} distinct gene intermediates")

res = run_de(GSE, download_matrix(GSE, MATRIX), [], [], is_counts=True,
             samples=fetch_samples(GSE), patterns=(r"^RA", r"^normal"))
if res.error:
    sys.exit(f"DE failed: {res.error}")
print(f"{GSE}: {res.matched_treated} RA vs {res.matched_control} healthy, "
      f"{len(res.table):,} genes (arms via {res.arm_source})")

tab = res.table
sig_all = (tab["adj.P.Val"] < FDR).sum()
print(f"  background: {sig_all:,}/{len(tab):,} genes DE at FDR<{FDR} "
      f"({100*sig_all/len(tab):.0f}%)\n")

hits, missing = [], []
for g, drugs in sorted(genes.items()):
    r = res.gene(g)
    (missing if r is None else hits).append(g if r is None else {**r, "n_drugs": len(drugs)})

sig = [h for h in hits if h["adj_p"] < FDR]
print(f"Translator RA genes: {len(genes)} | in matrix {len(hits)} | "
      f"DE at FDR<{FDR} {len(sig)} ({100*len(sig)/max(len(hits),1):.0f}%)")

# Size-matched permutation. A single FDR cut is nearly uninformative here --
# 152 vs 28 samples makes ~80% of the transcriptome "significant" -- so sweep an
# effect-size threshold too, and report the whole sweep rather than its best
# point, since five thresholds will throw a nominal p<0.05 by chance.
idx = {str(i).upper(): i for i in tab.index}
sub = tab.loc[[idx[h["gene"].upper()] for h in hits]]
rng = random.Random(0)
pool = list(tab.index)
sweep = []
for lfc in (0.0, 0.5, 1.0, 1.5, 2.0):
    obs = int(((sub["adj.P.Val"] < FDR) & (sub["logFC"].abs() > lfc)).sum())
    null = []
    for _ in range(N_PERM):
        d = tab.loc[rng.sample(pool, len(sub))]
        null.append(int(((d["adj.P.Val"] < FDR) & (d["logFC"].abs() > lfc)).sum()))
    sweep.append({
        "min_abs_logFC": lfc, "translator": obs,
        "translator_pct": round(100 * obs / len(sub)),
        "background_pct": round(100 * float(
            ((tab["adj.P.Val"] < FDR) & (tab["logFC"].abs() > lfc)).mean())),
        "permutation_median": sorted(null)[N_PERM // 2],
        "p": (sum(1 for n in null if n >= obs) + 1) / (N_PERM + 1),
    })

print(f"\n{'threshold':<26}{'Translator':>12}{'background':>12}{'perm p':>9}")
for r in sweep:
    print(f"  adj_p<{FDR} & |logFC|>{r['min_abs_logFC']:<5}"
          f"{r['translator_pct']:>10}%{r['background_pct']:>11}%{r['p']:>9.3f}")
print(f"\n  median |logFC|: Translator {sub['logFC'].abs().median():.3f} "
      f"vs background {tab['logFC'].abs().median():.3f}")
up = int((sub[sub["adj.P.Val"] < FDR]["logFC"] > 0).sum())
down = int((sub[sub["adj.P.Val"] < FDR]["logFC"] < 0).sum())
print(f"  significant Translator genes: {up} up / {down} down")

out = Path("results/route_b_translator_join.json")
out.write_text(json.dumps({
    "gse": GSE, "paths_file": str(PATHS), "fdr": FDR,
    "n_translator_genes": len(genes), "n_in_matrix": len(hits), "n_significant": len(sig),
    "background_sig": int(sig_all), "background_total": int(len(tab)),
    "sweep": sweep,
    "median_abs_logfc_translator": float(sub["logFC"].abs().median()),
    "median_abs_logfc_background": float(tab["logFC"].abs().median()),
    "n_up": up, "n_down": down,
    "genes": sorted(hits, key=lambda x: x["adj_p"]), "not_in_matrix": sorted(missing),
}, indent=2))
print(f"\nwrote {out}")
