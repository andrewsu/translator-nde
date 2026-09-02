"""Reanalyse GSE292059: mepolizumab vs placebo in nasal lavage.

Closes the loop the rest of the project only sets up. Translator proposes
mepolizumab for asthma through IL5 and IL5RA; this is a placebo-controlled
trial that treated with it and profiled expression. Two questions:

1. Is the Translator-nominated gene differentially expressed after treatment?
2. What pathways do *all* the differentially expressed genes fall in?

Design. 965 nasal lavage samples over 270 donors, longitudinal (visit01 ...
visit14, plus unscheduled "cold" visits). The primary contrast is the end of
treatment, visit14, mepolizumab against placebo -- a between-arm comparison at
a single timepoint, which avoids mixing the drug effect with season and
intercurrent colds. Baseline (visit01) is analysed identically as a control.

Method. **DESeq2** (via pydeseq2) is primary. Library sizes here span 58-fold
(0.27M to 15.7M reads), and log-CPM with a limma-style moderated t-test is only
recommended when that ratio stays under about 3 -- beyond it the observations
differ too much in precision to treat alike, and either voom weights or a
count-level model is needed. The moderated t-test is still computed, as a
concordance check rather than as the answer: methods that disagree about the
size of a differentially expressed set can still agree about a specific gene,
and it is worth showing which of those holds here.

Confounding. Nasal lavage is a cell mixture, and this trial reports its
composition per sample. Mepolizumab cuts the eosinophil fraction sevenfold, so
a bulk contrast between arms is partly a comparison of *different cell
mixtures* rather than of transcription within a cell type. The contrast is
therefore run twice, `~treatment` and `~treatment + eosinophil fraction`, which
turn out to answer different questions -- see the printed output.
"""
import gzip, io, json, sys, time, urllib.request, warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, "src")
from translator_nde._de import bh_fdr, counts_to_log2cpm, moderated_ttest

GSE = "GSE292059"
CACHE = Path("data/geo") / GSE
BASE = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE[:-3]}nnn/{GSE}"
FILES = {
    "counts": f"{BASE}/suppl/{GSE}_counts_raw.tsv.gz",
    "idmap": f"{BASE}/suppl/{GSE}_hgnc_ensembl_id_mapping.tsv.gz",
    "matrix": f"{BASE}/matrix/{GSE}_series_matrix.txt.gz",
}
TRANSLATOR_GENES = ["IL5", "IL5RA"]
FDR = 0.05
MIN_READS = 1_000_000   # 13 of 152 visit14 libraries fall below this


def fetch(key: str) -> bytes:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / Path(FILES[key]).name
    if not f.exists():
        with urllib.request.urlopen(FILES[key], timeout=600) as r:
            f.write_bytes(r.read())
    return f.read_bytes()


def sample_table() -> pd.DataFrame:
    """Per-sample annotation from the series matrix.

    The counts matrix is keyed by library id (`lib28551`), not GSM, and
    `!Sample_title` is what carries that id -- so the join runs through the
    title rather than the accession.
    """
    text = gzip.decompress(fetch("matrix")).decode("utf-8", "replace")
    head = text.split("!series_matrix_table_begin", 1)[0]

    def row(key):
        for l in head.splitlines():
            if l.startswith(key):
                return [v.strip('"') for v in l.split("\t")[1:]]
        return []

    cols = {"gsm": row("!Sample_geo_accession"), "lib": row("!Sample_title")}
    for l in head.splitlines():
        if not l.startswith("!Sample_characteristics_ch1"):
            continue
        vals = [v.strip('"') for v in l.split("\t")[1:]]
        if ":" not in vals[0]:
            continue
        field = vals[0].split(":", 1)[0].strip()
        cols[field] = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]
    return pd.DataFrame(cols).set_index("lib")


def load_counts() -> pd.DataFrame:
    # Named .tsv.gz but comma-delimited; trusting the extension yields one column.
    raw = gzip.decompress(fetch("counts")).decode("utf-8", "replace")
    df = pd.read_csv(io.StringIO(raw), index_col=0)
    idmap = pd.read_csv(io.BytesIO(gzip.decompress(fetch("idmap"))), sep="\t")
    sym = dict(zip(idmap["ensembl_id"], idmap["hgnc_symbol"]))
    df.index = [sym.get(i, i) for i in df.index]
    df = df[[bool(i) and not str(i).startswith("ENSG") for i in df.index]]
    return df.groupby(level=0).sum()          # collapse duplicate symbols


def arm_subset(counts, samples, visit, min_reads=MIN_READS):
    keep = [k for k in samples.index[(samples["visit"] == visit)
            & samples["treatment"].isin(["Mepolizumab", "Placebo"])]
            if k in counts.columns]
    sub = counts[keep]
    if min_reads:
        ok = sub.sum(axis=0) >= min_reads
        keep = [k for k in keep if ok[k]]
        sub = counts[keep]
    treated = np.array([samples.loc[k, "treatment"] == "Mepolizumab" for k in keep])
    return sub, treated, keep


CELL_TYPES = ["eos total", "neut total", "mac total", "lym total",
              "epi total", "squ total"]


def cell_composition(samples, keep, treated):
    """Differential-cell percentages by arm, with a rank test.

    The most direct evidence that the drug worked, and it needs no expression
    data at all -- which also makes it the explanation for most of what the
    expression analysis finds.
    """
    from scipy import stats
    df = samples.loc[keep, CELL_TYPES].apply(pd.to_numeric, errors="coerce")
    out = {}
    for c in CELL_TYPES:
        a = df.loc[treated, c].dropna()
        b = df.loc[~np.array(treated), c].dropna()
        if len(a) > 5 and len(b) > 5:
            out[c] = {"mepolizumab_median": float(a.median()),
                      "placebo_median": float(b.median()),
                      "p": float(stats.mannwhitneyu(a, b).pvalue)}
    return out


def deseq2(sub, treated, keep, samples=None, adjust_eos=False):
    """Negative-binomial test with median-of-ratios size factors."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    meta = pd.DataFrame(
        {"treatment": np.where(treated, "Mepolizumab", "Placebo")}, index=keep)
    design = "~treatment"
    if adjust_eos:
        eos = pd.to_numeric(samples.loc[keep, "eos total"], errors="coerce")
        meta["eos"] = eos.astype(float).values
        meta = meta.dropna()
        design = "~treatment + eos"
    cts = sub.T.astype(int).loc[meta.index]
    cts = cts.loc[:, cts.sum(axis=0) > 0]
    # pydeseq2's IRLS overflows on a handful of genes and recovers. The warnings
    # come from worker processes, so neither np.seterr nor an errstate context
    # here suppresses them -- redirect stderr when running if they are in the way.
    dds = DeseqDataSet(counts=cts, metadata=meta, design=design, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["treatment", "Mepolizumab", "Placebo"], quiet=True)
    st.summary()
    d = st.results_df.dropna(subset=["padj"])
    return d.rename(columns={"log2FoldChange": "logFC", "pvalue": "P.Value",
                             "padj": "adj.P.Val"})


def limma(sub, treated):
    """limma-style moderated t-test on log2-CPM, for concordance only."""
    expr = counts_to_log2cpm(sub)
    res = moderated_ttest(expr, treated)
    res["adj.P.Val"] = bh_fdr(res["P.Value"].values)
    return res


def enrich(genes, library="GO_Biological_Process_2023", top=12):
    """Enrichr over-representation for one library.

    `addList` takes multipart/form-data; urlencoding it returns HTTP 400.
    """
    if not genes:
        return []
    boundary = "----translatornde"
    body = "".join(
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n"
        for k, v in (("list", "\n".join(genes)), ("description", GSE))
    ) + f"--{boundary}--\r\n"
    req = urllib.request.Request(
        "https://maayanlab.cloud/Enrichr/addList", data=body.encode(),
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    uid = json.loads(urllib.request.urlopen(req, timeout=180).read())["userListId"]
    time.sleep(2)   # Enrichr returns 429 if the calls are back to back
    url = (f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}"
           f"&backgroundType={library}")
    res = json.loads(urllib.request.urlopen(url, timeout=300).read())[library]
    # Enrichr rows: [rank, term, p, odds ratio, combined score, genes, adj p, ...]
    return [{"term": r[1], "p": r[2], "adj_p": r[6], "odds": r[3],
             "n_genes": len(r[5]), "genes": r[5][:8]} for r in res[:top]]


def report(res, name, sig=None):
    sig = res[res["adj.P.Val"] < FDR] if sig is None else sig
    up = sig.index[sig["logFC"] > 0]
    down = sig.index[sig["logFC"] < 0]
    print(f"  {name:12s} {len(res):,} genes tested, {len(sig):,} at FDR<{FDR} "
          f"({len(up)} up, {len(down)} down)")
    for g in TRANSLATOR_GENES:
        if g in res.index:
            r = res.loc[g]
            flag = "DE" if r["adj.P.Val"] < FDR else "not DE"
            print(f"      {g:6s} logFC={r['logFC']:+.3f} p={r['P.Value']:.2g} "
                  f"adj_p={r['adj.P.Val']:.2g}  [{flag}]")
        else:
            print(f"      {g:6s} filtered out")
    return sig, list(up), list(down)


def main():
    samples = sample_table()
    counts = load_counts()
    print(f"{GSE}: {counts.shape[0]:,} genes x {counts.shape[1]:,} libraries")
    print(f"  {len(samples):,} samples, {samples['donorid'].nunique()} donors, "
          f"tissue={samples['tissue'].iloc[0]}")

    out = {"gse": GSE, "fdr": FDR, "min_reads": MIN_READS,
           "primary_method": "DESeq2 (pydeseq2)", "contrasts": {}}

    for visit, label in [("visit14", "END OF TREATMENT (visit14)"),
                         ("visit01", "BASELINE (visit01) — control")]:
        sub, treated, keep = arm_subset(counts, samples, visit)
        ls = sub.sum(axis=0)
        print(f"\n{label}: {int(treated.sum())} mepolizumab vs "
              f"{int((~treated).sum())} placebo")
        print(f"  library sizes {ls.min()/1e6:.2f}–{ls.max()/1e6:.2f}M "
              f"(ratio {ls.max()/ls.min():.1f}x) — why DESeq2 rather than log-CPM")

        comp = cell_composition(samples, keep, treated)
        print("  differential cell % (median, mepolizumab vs placebo):")
        for c, v in comp.items():
            flag = "  <-- differs" if v["p"] < 0.05 else ""
            print(f"      {c:11s} {v['mepolizumab_median']:6.2f} vs "
                  f"{v['placebo_median']:6.2f}   p={v['p']:.2g}{flag}")

        d = deseq2(sub, treated, keep)
        d_sig, d_up, d_down = report(d, "DESeq2")
        l = limma(sub, treated)
        l_sig, _, _ = report(l, "limma (check)")

        shared = set(d_sig.index) & set(l_sig.index)
        print(f"  concordance: {len(shared)} genes significant under both "
              f"({100*len(shared)/max(len(l_sig),1):.0f}% of the limma set)")

        # Does anything survive once the cell mixture is accounted for?
        adj = deseq2(sub, treated, keep, samples=samples, adjust_eos=True)
        adj_sig, _, _ = report(adj, "+eos covar")

        rec = {
            "n_treated": int(treated.sum()), "n_control": int((~treated).sum()),
            "lib_min": float(ls.min()), "lib_max": float(ls.max()),
            "n_tested": int(len(d)), "n_significant": int(len(d_sig)),
            "n_up": len(d_up), "n_down": len(d_down),
            "limma_n_significant": int(len(l_sig)), "n_shared": len(shared),
            "cell_composition": comp,
            "eos_adjusted": {
                "n_significant": int(len(adj_sig)),
                "translator_genes": {
                    g: ({"logFC": float(adj.loc[g, "logFC"]),
                         "adj_p": float(adj.loc[g, "adj.P.Val"])}
                        if g in adj.index else None) for g in TRANSLATOR_GENES}},
            "translator_genes": {
                g: ({"logFC": float(d.loc[g, "logFC"]),
                     "p": float(d.loc[g, "P.Value"]),
                     "adj_p": float(d.loc[g, "adj.P.Val"]),
                     "limma_logFC": float(l.loc[g, "logFC"]) if g in l.index else None,
                     "limma_adj_p": float(l.loc[g, "adj.P.Val"]) if g in l.index else None}
                    if g in d.index else None)
                for g in TRANSLATOR_GENES},
            "top_de": [{"gene": g, "logFC": float(d.loc[g, "logFC"]),
                        "adj_p": float(d.loc[g, "adj.P.Val"])}
                       for g in d_sig.sort_values("adj.P.Val").index[:40]],
        }
        if len(d_sig):
            print(f"  top DE: " + ", ".join(
                f"{g}({d.loc[g,'logFC']:+.2f})"
                for g in d_sig.sort_values("adj.P.Val").index[:10]))
        if visit == "visit14" and len(d_sig):
            for name, glist in (("down", d_down), ("up", d_up)):
                for lib in ("Reactome_2022", "GO_Biological_Process_2023"):
                    try:
                        terms = enrich(glist, library=lib)
                    except Exception as exc:
                        print(f"  {name}/{lib}: enrichment failed ({exc})")
                        time.sleep(5)
                        continue
                    time.sleep(2)
                    rec.setdefault("enrichment", {})[f"{name}:{lib}"] = terms
                    n_hit = sum(1 for t in terms if t["adj_p"] < FDR)
                    print(f"\n  {name} ({len(glist)} genes) — {lib} "
                          f"[{n_hit} terms at adj p<{FDR}]")
                    for t in terms[:6]:
                        print(f"   {'*' if t['adj_p'] < FDR else ' '} "
                              f"adj_p={t['adj_p']:.2e} n={t['n_genes']:<3d} {t['term'][:58]}")
        out["contrasts"][visit] = rec

    Path("results/gse292059_reanalysis.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/gse292059_reanalysis.json")


if __name__ == "__main__":
    main()
