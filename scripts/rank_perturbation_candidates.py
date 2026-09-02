"""Rank the confirmed drug x series pairs into a shortlist worth re-analysing.

Reads results/perturbation_series.json and keeps pairs that have BOTH a level
naming the drug and a control-like level, then checks GEO for a deposited
expression matrix. Ranked by the Translator answer's own score for that drug,
so the top of the list is what Translator most confidently proposed *and* what
somebody has actually run an experiment on.
"""
import gzip, json, re, sys, urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.geo import inspect


def series_matrix_rows(gse: str) -> int:
    """Expression rows in the GEO series matrix itself.

    `geo.inspect` only classifies files under `suppl/`, so microarray series --
    whose expression values live in the series matrix and nowhere else -- come
    back as `none` despite being perfectly re-analysable. Checked: GSE59671 has
    22,278 rows, GSE45468 16,963, GSE73578 54,676.
    """
    pre = re.sub(r"\d{1,3}$", "nnn", gse)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            text = gzip.decompress(r.read()).decode("utf-8", "replace")
    except Exception:
        return 0
    body = text.split("!series_matrix_table_begin", 1)
    if len(body) < 2:
        return 0
    return max(0, body[1].split("!series_matrix_table_end")[0].count("\n") - 2)

DISEASES = {
    "asthma": "5b656c0f-b7da-4db4-ba1f-d3a794b422d4",
    "rheumatoid arthritis": "d1f9daae-9330-4c2c-81f5-5c21ee2f44f4",
    "ankylosing spondylitis": "c3efea15-1a8f-4c09-a521-7462603bec08",
}
# \b before the 0 matters: without it "0 uM" matches inside "50 uM" and a
# dose-only arm is mistaken for its own vehicle control.
CONTROL = re.compile(
    r"vehicle|placebo|untreated|\bcontrol\b|dmso|\bnone\b|no treat|baseline|"
    r"saline|\bpbs\b|\b0\s?(u|n|m|µ)M\b|\bmock\b", re.I)
TOP_PER_DISEASE = 6

prom = defaultdict(lambda: [0, set(), 0.0])
for dis, pk in DISEASES.items():
    for p in json.loads(Path(f"data/ars/{pk}/paths.json").read_text())["paths"]:
        e = prom[(dis, p.get("drug_name"))]
        e[0] += 1
        e[1].add(p["gene_name"])
        e[2] = max(e[2], p.get("score") or 0)

confirmed = json.loads(Path("results/perturbation_series.json").read_text())["confirmed"]
rows = []
for c in confirmed:
    vals = list(c["values"])
    arms = [v for v in vals if c["drug"].lower() in v.lower()]
    ctrls = [v for v in vals if CONTROL.search(v) and c["drug"].lower() not in v.lower()]
    if not (arms and ctrls):
        continue
    for dis in c["diseases"]:
        n_paths, genes, score = prom.get((dis, c["drug"]), [0, set(), 0.0])
        if not n_paths:
            continue
        rows.append({"disease": dis, "drug": c["drug"], "gse": c["gse"],
                     "field": c["field"], "arm": arms[0], "control": ctrls[0],
                     "n_samples_mentioning": c["n_samples_mentioning"],
                     "translator_paths": n_paths, "translator_genes": sorted(genes),
                     "translator_score": round(score, 3)})

rows.sort(key=lambda r: (-r["translator_score"], -r["translator_paths"]))
print(f"{len(rows)} clean drug x series pairs "
      f"({len({r['drug'] for r in rows})} drugs, {len({r['gse'] for r in rows})} series)\n")

shortlist, seen = [], defaultdict(int)
for r in rows:
    if seen[r["disease"]] >= TOP_PER_DISEASE or (r["disease"], r["drug"]) in shortlist:
        continue
    g = inspect(r["gse"])
    r["matrix_kind"] = g.matrix_kind
    r["matrix_files"] = g.matrix_files[:2]
    r["series_matrix_rows"] = series_matrix_rows(r["gse"]) if g.matrix_kind == "none" else 0
    r["usable"] = g.matrix_kind != "none" or r["series_matrix_rows"] > 1000
    shortlist.append(r)
    seen[r["disease"]] += 1

for dis in DISEASES:
    sub = [r for r in shortlist if r["disease"] == dis]
    print(f"\n=== {dis.upper()} ===")
    for r in sub:
        src = (r["matrix_kind"] if r["matrix_kind"] != "none"
               else (f"series_matrix({r['series_matrix_rows']:,})"
                     if r["usable"] else "NONE"))
        print(f"  {'OK ' if r['usable'] else '-- '}{r['drug'][:20]:20s} {r['gse']:11s} "
              f"{src:22s} n~{r['n_samples_mentioning']:<4d} "
              f"score={r['translator_score']:.2f} "
              f"genes={','.join(r['translator_genes'][:4])}")
        print(f"      {r['arm'][:64]}  VS  {r['control'][:40]}")

out = Path("results/perturbation_shortlist.json")
out.write_text(json.dumps({"all_clean_pairs": rows, "shortlist": shortlist}, indent=2))
print(f"\nwrote {out}")
