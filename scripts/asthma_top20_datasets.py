"""For asthma's top-ranked Translator answers, find GEO series that perturb with them.

Answers the discovery question across the whole head of the answer list rather
than a hand-picked drug: of the 20 highest-scoring drugs, which have a public
expression experiment that actually treated with them?

Reads the two-stage sweep in results/perturbation_series.json (NDE sample-level
search, then GEO `!Sample_characteristics_ch1` confirmation) and joins it back
to the archived Translator answer.
"""
import collections, gzip, json, re, urllib.request
from collections import defaultdict
from pathlib import Path

HEADER_CACHE = Path("data/geo/series_headers")


def series_header(gse: str) -> str:
    """The metadata block of a GEO series matrix, cached."""
    HEADER_CACHE.mkdir(parents=True, exist_ok=True)
    f = HEADER_CACHE / f"{gse}.txt"
    if f.exists():
        return f.read_text()
    pre = re.sub(r"\d{1,3}$", "nnn", gse)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            text = gzip.decompress(r.read()).decode("utf-8", "replace")
    except Exception:
        return ""
    head = text.split("!series_matrix_table_begin", 1)[0]
    f.write_text(head)
    return head


def _rows(header: str, key: str) -> list[list[str]]:
    return [[v.strip('"') for v in l.split("\t")[1:]]
            for l in header.splitlines() if l.startswith(key)]


def describe_system(gse: str) -> str:
    """A short phrase for the experimental system: cell line, tissue, cohort.

    The drug being right does not make the model right -- several series below
    treat with an asthma drug in an unrelated system -- so this belongs in the
    table rather than in a footnote.
    """
    header = series_header(gse)
    if not header:
        return "—"
    chars: dict[str, str] = {}
    for vals in _rows(header, "!Sample_characteristics_ch1"):
        named = [v for v in vals if ":" in v]
        if not named:
            continue
        field = named[0].split(":", 1)[0].strip().lower()
        common = collections.Counter(
            v.split(":", 1)[1].strip() for v in named).most_common(1)[0][0]
        chars.setdefault(field, common)
    src = _rows(header, "!Sample_source_name_ch1")
    src_common = (collections.Counter(v for v in src[0] if v).most_common(1)[0][0]
                  if src and any(src[0]) else "")

    # Keep both when a series records a line and what kind of cell it is
    # ("cell line: WTC-11" + "cell type: human iPSCs").
    model = ", ".join(dict.fromkeys(
        chars[k] for k in ("cell line", "cell type") if chars.get(k)))
    tissue = next((chars[k] for k in ("tissue", "organism part") if k in chars), "")
    cohort = next((chars[k] for k in ("disease status", "diagnosis", "disease state",
                                      "disease", "condition") if k in chars), "")
    parts = [x for x in (model, tissue, cohort) if x]
    if not parts:
        parts = [src_common]
    seen, out = set(), []
    for x in parts:                      # "tissue: Nasal lavage" often repeats source
        k = x.lower()
        if k not in seen:
            seen.add(k)
            out.append(x)
    desc = ", ".join(out)
    # Some series abbreviate the characteristic ("tissue: LS" for lesional skin);
    # the free-text source name is more informative when that happens.
    if len(desc) < 8 and src_common:
        desc = src_common
    return desc[:70] or "—"

PK = "5b656c0f-b7da-4db4-ba1f-d3a794b422d4"
TOP_N = 20
# \b before the 0 matters: "0 uM" otherwise matches inside "50 uM".
CONTROL = re.compile(
    r"vehicle|placebo|untreated|\bcontrol\b|dmso|\bnone\b|no treat|baseline|"
    r"saline|\bpbs\b|\b0\s?(u|n|m|µ)M\b|\bmock\b", re.I)
# A dose-series control reads "<drug>_dosis: 0" -- it names the drug, so it looks
# like a treatment arm unless zero and absent levels are excluded explicitly.
ZERO_LEVEL = re.compile(r":\s*(0|0\.0|NA|-|FALSE|no|none|neg\w*)\s*$", re.I)

paths = json.loads(Path(f"data/ars/{PK}/paths.json").read_text())["paths"]
score, genes = {}, defaultdict(set)
for p in paths:
    score[p["drug_name"]] = max(score.get(p["drug_name"], 0), p.get("score") or 0)
    genes[p["drug_name"]].add(p["gene_name"])
top = sorted(score, key=lambda d: -score[d])[:TOP_N]

sweep = json.loads(Path("results/perturbation_series.json").read_text())
stage1 = {r["drug"].lower(): r for r in sweep["stage1"]}
confirmed = defaultdict(list)
for c in sweep["confirmed"]:
    confirmed[c["drug"].lower()].append(c)

rows = []
for rank, drug in enumerate(top, 1):
    s1 = stage1.get(drug.lower())
    hits = confirmed.get(drug.lower(), [])
    # Prefer a series with an explicit control arm; those are directly analysable.
    def has_control(c):
        return any(CONTROL.search(v) and drug.lower() not in v.lower() for v in c["values"])
    hits = sorted(hits, key=lambda c: (not has_control(c), -c["n_samples_mentioning"]))
    best = hits[0] if hits else None
    arm = ctrl = None
    if best:
        named = [v for v in best["values"] if drug.lower() in v.lower()]
        # Single agent beats a combination, and a real dose beats a zero level.
        dosed = sorted((v for v in named if not ZERO_LEVEL.search(v)),
                       key=lambda v: ("+" in v, len(v)))
        arm = dosed[0] if dosed else (named[0] if named else None)
        ctrl = next(
            (v for v in best["values"]
             if v != arm and (ZERO_LEVEL.search(v)
                              or (CONTROL.search(v) and drug.lower() not in v.lower()))),
            None)
    rows.append({
        "rank": rank, "drug": drug, "score": round(score[drug], 3),
        "system": describe_system(best["gse"]) if best else None,
        "genes": sorted(genes[drug]),
        "nde_samples": s1["n_samples"] if s1 else 0,
        "candidate_series": s1["n_series"] if s1 else 0,
        "confirmed_series": len(hits),
        "best_series": best["gse"] if best else None,
        "arm": arm, "control": ctrl,
        "all_confirmed": [c["gse"] for c in hits],
    })

out = Path("results/asthma_top20_datasets.json")
out.write_text(json.dumps({"pk": PK, "top_n": TOP_N, "rows": rows}, indent=2))

with_data = [r for r in rows if r["confirmed_series"]]
with_ctrl = [r for r in with_data if r["control"]]
print(f"top {TOP_N} asthma answers: {len(with_data)} have >=1 confirmed perturbation series, "
      f"{len(with_ctrl)} with an explicit control arm")
print(f"{'#':<3}{'drug':24s}{'score':>6}{'samples':>9}{'series':>7}  best / arm vs control")
for r in rows:
    print(f"{r['rank']:<3}{r['drug'][:22]:24s}{r['score']:>6.2f}"
          f"{r['nde_samples']:>9,}{r['confirmed_series']:>7}"
          f"  {r['best_series'] or '—'}   {r['system'] or ''}")
    if r["arm"]:
        print(f"      {r['arm'][:58]}  VS  {(r['control'] or '(no control level)')[:34]}")
print(f"\nwrote {out}")
