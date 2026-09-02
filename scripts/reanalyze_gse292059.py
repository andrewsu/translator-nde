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
intercurrent colds. Baseline (visit01) is analysed the same way as a negative
control: before treatment the arms should not differ.
"""
import gzip, io, json, sys, time, urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

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


def contrast(counts, samples, visit, label):
    keep = samples.index[(samples["visit"] == visit)
                         & samples["treatment"].isin(["Mepolizumab", "Placebo"])]
    keep = [k for k in keep if k in counts.columns]
    sub = counts[keep]
    is_treated = np.array([samples.loc[k, "treatment"] == "Mepolizumab" for k in keep])
    expr = counts_to_log2cpm(sub)
    res = moderated_ttest(expr, is_treated)   # positive logFC = up in mepolizumab
    res["adj.P.Val"] = bh_fdr(res["P.Value"].values)
    n_sig = int((res["adj.P.Val"] < FDR).sum())
    print(f"\n{label}: {is_treated.sum()} mepolizumab vs {(~is_treated).sum()} placebo, "
          f"{len(expr):,} genes tested")
    print(f"  significant at FDR<{FDR}: {n_sig:,} ({100*n_sig/len(res):.1f}%)")
    return res, int(is_treated.sum()), int((~is_treated).sum())


def enrich(genes, library="GO_Biological_Process_2023", top=12):
    """Enrichr over-representation for one library.

    `addList` takes multipart/form-data; urlencoding it returns HTTP 400.
    """
    if not genes:
        return []
    boundary = "----translatornde"
    body = "".join(
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n"
        for k, v in (("list", "\n".join(genes)), ("description", "GSE292059"))
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


def main():
    samples = sample_table()
    counts = load_counts()
    print(f"{GSE}: {counts.shape[0]:,} genes x {counts.shape[1]:,} libraries")
    print(f"  annotation: {len(samples):,} samples, "
          f"{samples['donorid'].nunique()} donors, tissue={samples['tissue'].iloc[0]}")

    out = {"gse": GSE, "fdr": FDR, "contrasts": {}}
    for visit, label in [("visit14", "END OF TREATMENT (visit14)"),
                         ("visit01", "BASELINE (visit01) — negative control")]:
        res, n_t, n_c = contrast(counts, samples, visit, label)
        sig = res[res["adj.P.Val"] < FDR].sort_values("adj.P.Val")

        print(f"  Translator-nominated genes:")
        for g in TRANSLATOR_GENES:
            if g in res.index:
                r = res.loc[g]
                mark = "DE" if r["adj.P.Val"] < FDR else "not DE"
                print(f"    {g:7s} logFC={r['logFC']:+.3f}  p={r['P.Value']:.3g}  "
                      f"adj_p={r['adj.P.Val']:.3g}  [{mark}]")
            else:
                print(f"    {g:7s} not detected above the expression filter")

        up = sig.index[sig["logFC"] > 0]
        down = sig.index[sig["logFC"] < 0]
        if len(sig):
            print(f"  direction: {len(up)} up, {len(down)} down in mepolizumab")
            print(f"  top DE genes: " + ", ".join(
                f"{g}({res.loc[g,'logFC']:+.2f})" for g in sig.index[:10]))
        out["contrasts"][visit] = {
            "n_treated": n_t, "n_control": n_c, "n_tested": int(len(res)),
            "n_significant": int(len(sig)),
            "translator_genes": {g: ({"logFC": float(res.loc[g, "logFC"]),
                                      "p": float(res.loc[g, "P.Value"]),
                                      "adj_p": float(res.loc[g, "adj.P.Val"])}
                                     if g in res.index else None)
                                 for g in TRANSLATOR_GENES},
            "top_de": [{"gene": g, "logFC": float(res.loc[g, "logFC"]),
                        "adj_p": float(res.loc[g, "adj.P.Val"])}
                       for g in sig.index[:40]],
        }
        out["contrasts"][visit]["n_up"] = int(len(up))
        out["contrasts"][visit]["n_down"] = int(len(down))
        if visit == "visit14" and len(sig):
            # Enrich up and down separately: a mixed list dilutes both, and here
            # the biology is almost entirely one-directional.
            for name, glist in (("all", list(sig.index)),
                                ("down", list(down)), ("up", list(up))):
                if not glist:
                    continue
                for lib in ("Reactome_2022", "GO_Biological_Process_2023"):
                    try:
                        terms = enrich(glist, library=lib)
                    except Exception as exc:
                        print(f"\n  {name}/{lib}: enrichment failed ({exc})")
                        time.sleep(5)
                        continue
                    time.sleep(2)
                    out["contrasts"][visit].setdefault("enrichment", {})[
                        f"{name}:{lib}"] = terms
                    n_hit = sum(1 for t in terms if t["adj_p"] < 0.05)
                    print(f"\n  {name} ({len(glist)} genes) — {lib} "
                          f"[{n_hit} terms at adj p<0.05]")
                    for t in terms[:6]:
                        flag = "*" if t["adj_p"] < 0.05 else " "
                        print(f"   {flag} adj_p={t['adj_p']:.2e} n={t['n_genes']:<3d} "
                              f"{t['term'][:60]}")

    Path("results/gse292059_reanalysis.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/gse292059_reanalysis.json")


if __name__ == "__main__":
    main()
