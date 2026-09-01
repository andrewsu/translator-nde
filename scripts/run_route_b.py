"""Route B: differential expression for selected drug-treatment contrasts.

Each entry is a GEO series with a genuine drug arm and a matched control arm --
the thing Route B exists to test. Disease-vs-healthy series are excluded: they
validate the machinery but test no drug->gene edge.
"""
import json, sys, time
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.reanalysis import download_matrix, fetch_samples, run_de

RUNS = [
    {
        "gse": "GSE97165", "file": "GSE97165_GEO_count_matrix.txt.gz",
        "label": "RA synovium, post vs pre triple DMARD (MTX+SSZ+HCQ)",
        "drug": "methotrexate (triple DMARD)", "disease": "rheumatoid arthritis",
        "patterns": (r"_post", r"_pre"),
        "subject_re": r"RA_(?:pre|post)_(\d+)$",
        # Expect inflammation DOWN after effective therapy.
        "genes": ["MMP1", "MMP3", "IL6", "TNF", "CXCL13", "MS4A1", "CD3E", "POSTN"],
        "expect": "decreased",
    },
    {
        "gse": "GSE148395", "file": "GSE148395_final_FLS_JQ1_genecounts.txt.gz",
        "label": "RA synovial fibroblasts, JQ1 (BET inhibitor) vs DMSO",
        "drug": "JQ1", "disease": "rheumatoid arthritis",
        "patterns": (r"[-_]JQ", r"[-_]DMSO"),
        # JQ1 displaces BRD4 from super-enhancers; MYC and IL6 are the canonical
        # readouts and should fall.
        "genes": ["MYC", "IL6", "BRD2", "BRD3", "BRD4", "CXCL8", "MMP1", "MMP3"],
        "expect": "decreased",
    },
    {
        "gse": "GSE141646", "file": "GSE141646_AS-TNF_counts.txt.gz",
        "label": "AS whole blood, post vs pre TNF inhibitor",
        "drug": "TNF inhibitor", "disease": "ankylosing spondylitis",
        "patterns": (r"_POST", r"_PRE"),
        "subject_re": r"^(A\d+)_",
        "genes": ["TNF", "IL6", "IL1B", "CXCL8", "TLR4", "NFKB1", "SOCS3", "MMP3"],
        "expect": "decreased",
    },
]

out = []
for r in RUNS:
    t0 = time.time()
    print(f"\n{'='*78}\n{r['gse']}  {r['label']}\n{'='*78}", flush=True)
    samples = fetch_samples(r["gse"])
    path = download_matrix(r["gse"], r["file"])
    res = run_de(r["gse"], path, [], [], is_counts=True, samples=samples,
                 patterns=r["patterns"], subject_re=r.get("subject_re"))
    if res.error:
        print(f"  ERROR: {res.error}")
        out.append({**{k: r[k] for k in ("gse", "label", "drug", "disease")},
                    "error": res.error})
        continue

    print(f"  arms: treated={res.matched_treated} control={res.matched_control} "
          f"(via {res.arm_source})  design={res.design}"
          f"{f' n_pairs={res.n_pairs}' if res.n_pairs else ''}  "
          f"genes tested={len(res.table)}  [{time.time()-t0:.0f}s]")
    for w in res.warnings:
        print(f"  ! {w}")
    sig = res.table[res.table["adj.P.Val"] < 0.05]
    print(f"  significant at FDR<0.05: {len(sig)} "
          f"({100*len(sig)/len(res.table):.0f}% of genes)")
    print(f"\n  top 8 by p-value:")
    for g, row in res.table.head(8).iterrows():
        print(f"    {str(g)[:14]:16s} logFC={row['logFC']:+6.2f}  "
              f"adj_p={row['adj.P.Val']:.2e}")

    print(f"\n  marker panel (expect {r['expect']} in treated arm):")
    hits = {}
    for g in r["genes"]:
        v = res.gene(g)
        hits[g] = v
        if v is None:
            print(f"    {g:8s} not in matrix")
        else:
            ok = "OK " if v["direction"] == r["expect"] and v["adj_p"] < 0.05 else \
                 ("ns " if v["adj_p"] >= 0.05 else "OPP")
            print(f"    {g:8s} logFC={v['logFC']:+6.2f}  adj_p={v['adj_p']:.2e}  "
                  f"{v['direction']:9s} [{ok}]")
    agree = sum(1 for v in hits.values()
                if v and v["adj_p"] < 0.05 and v["direction"] == r["expect"])
    tested = sum(1 for v in hits.values() if v)
    print(f"\n  panel agreement: {agree}/{tested} significant in expected direction")

    out.append({**{k: r[k] for k in ("gse", "label", "drug", "disease", "expect")},
                "n_treated": res.matched_treated, "n_control": res.matched_control,
                "arm_source": res.arm_source, "design": res.design,
                "n_pairs": res.n_pairs, "warnings": res.warnings,
                "genes_tested": len(res.table),
                "n_sig_fdr05": int(len(sig)), "panel": hits,
                "panel_agree": agree, "panel_tested": tested,
                "top20": [{"gene": str(g), "logFC": float(row["logFC"]),
                           "adj_p": float(row["adj.P.Val"])}
                          for g, row in res.table.head(20).iterrows()]})

Path("results/route_b_runs.json").write_text(json.dumps(out, indent=2))
print(f"\n\nwrote results/route_b_runs.json")
