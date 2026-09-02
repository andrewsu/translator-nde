"""For asthma's top-ranked Translator answers, find GEO series that perturb with them.

Answers the discovery question across the whole head of the answer list rather
than a hand-picked drug: of the 20 highest-scoring drugs, which have a public
expression experiment that actually treated with them?

Reads the two-stage sweep in results/perturbation_series.json (NDE sample-level
search, then GEO `!Sample_characteristics_ch1` confirmation) and joins it back
to the archived Translator answer.
"""
import collections, gzip, json, re, sys, urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.nde import NDEClient, PROD

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


# GEO's own declaration of what the experiment is. Worth carrying into the table:
# GSE136034 is a triamcinolone-vs-DMSO perturbation with clean arms, so it passes
# arm confirmation, but its only data file is a 4C-seq chromatin-capture table --
# it cannot answer a differential-expression question at all.
_TYPE_SHORT = {
    "expression profiling by high throughput sequencing": "RNA-seq",
    "expression profiling by array": "array",
    "genome binding/occupancy profiling by high throughput sequencing": "ChIP-seq",
    "methylation profiling by high throughput sequencing": "methylation",
    "non-coding rna profiling by high throughput sequencing": "ncRNA-seq",
}


def experiment_type(gse: str) -> tuple[str, bool]:
    """(short label, is it expression data) as declared by GEO."""
    header = series_header(gse)
    types = [v.strip('"') for l in header.splitlines()
             if l.startswith("!Series_type") for v in l.split("\t")[1:]]
    strat = {v.strip('"').upper() for l in header.splitlines()
             if l.startswith("!Sample_library_strategy") for v in l.split("\t")[1:]}
    labels, is_expr = [], False
    for t in dict.fromkeys(types):
        short = _TYPE_SHORT.get(t.lower())
        if short:
            is_expr = True
            labels.append(short)
        else:
            # "Other" is uninformative on its own; the library strategy says more.
            extra = ", ".join(sorted(s for s in strat if s and s != "OTHER"))
            labels.append(f"other ({extra})" if extra else "other")
    return (", ".join(dict.fromkeys(labels)) or "—"), is_expr


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
    org = _rows(header, "!Sample_organism_ch1")
    organism = (collections.Counter(v for v in org[0] if v).most_common(1)[0][0]
                if org and any(org[0]) else "")
    # Species belongs in the description: a mouse organoid screen and a human
    # bronchial brushing are not interchangeable evidence for a human disease.
    prefix = "" if organism == "Homo sapiens" else (f"{organism}: " if organism else "")
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
    return (prefix + desc)[:74] or "—"

PK = "5b656c0f-b7da-4db4-ba1f-d3a794b422d4"
TOP_N = 20
# \b before the 0 matters: "0 uM" otherwise matches inside "50 uM".
CONTROL = re.compile(
    r"vehicle|placebo|untreated|\bcontrol\b|dmso|\bnone\b|no treat|baseline|"
    r"saline|\bpbs\b|\b0\s?(u|n|m|µ)M\b|\bmock\b|unstim\w*|"
    r"no stimulation|\bDMSO\b", re.I)
# A dose-series control reads "<drug>_dosis: 0" -- it names the drug, so it looks
# like a treatment arm unless zero and absent levels are excluded explicitly.
ZERO_LEVEL = re.compile(r":\s*(0|0\.0|NA|-|FALSE|no|none|neg\w*)\s*$", re.I)

paths = json.loads(Path(f"data/ars/{PK}/paths.json").read_text())["paths"]
score, genes = {}, defaultdict(set)
for p in paths:
    score[p["drug_name"]] = max(score.get(p["drug_name"], 0), p.get("score") or 0)
    genes[p["drug_name"]].add(p["gene_name"])
top = sorted(score, key=lambda d: -score[d])[:TOP_N]

# Fields a depositor may use to encode the treatment arm. `characteristics_ch1`
# is the conventional place, but GSE20297 -- a real terbutaline experiment --
# keeps its arms only in the sample title ("HaCaT_TNF+IFNg+terbutaline") and
# leaves characteristics constant, so confirming on characteristics alone
# discards it.
ARM_FIELDS = ("!Sample_characteristics_ch1", "!Sample_title",
              "!Sample_source_name_ch1")


def confirm_arm(gse: str, drug: str) -> dict | None:
    """A per-sample field that names the drug and takes more than one value."""
    header = series_header(gse)
    for key in ARM_FIELDS:
        for vals in _rows(header, key):
            counts = collections.Counter(v for v in vals if v)
            if len(counts) < 2 or not any(drug.lower() in v.lower() for v in counts):
                continue
            first = next(iter(counts))
            field = (first.split(":", 1)[0].strip() if ":" in first
                     else key.replace("!Sample_", ""))
            org = _rows(header, "!Sample_organism_ch1")
            organism = (collections.Counter(v for v in org[0] if v).most_common(1)[0][0]
                        if org and any(org[0]) else "")
            # Always retain the levels that name the drug. GSE157167 is a
            # ~300-compound organoid screen, so pemirolast's 4-sample level
            # never survives a plain most_common(8) and the arm looks absent.
            keep = dict(counts.most_common(8))
            keep.update({v: n for v, n in counts.items()
                         if drug.lower() in v.lower()})
            return {"gse": gse, "field": field, "values": keep,
                    "organism": organism, "arm_field": key.replace("!Sample_", ""),
                    "n_samples_mentioning": sum(
                        n for v, n in counts.items() if drug.lower() in v.lower())}
    return None


MIN_PER_SERIES = 3   # below this a mention is more likely incidental than an arm
SCROLL_CAP = 1500
nde = NDEClient(base_url=PROD)


def nde_candidates(drug: str) -> tuple[int, int, list[tuple[str, int]]]:
    """(human samples, samples any species, candidate series) naming the drug.

    Both counts are reported because the human filter is not free: pemirolast's
    only perturbation data is GSE157167, a mouse intestinal-organoid screen, and
    filtering to Homo sapiens made the drug look like it had no data at all. The
    series search is left unfiltered and the organism is surfaced in the
    experimental-system column instead, so the reader can judge.
    """
    q_any = f'@type:Sample AND "{drug}"'
    n_any = nde.count(q_any)
    n_human = nde.count(f'{q_any} AND species.name:"Homo sapiens"')
    by = collections.Counter()
    if n_any:
        for h in nde.scroll(q_any, fields="isBasisFor.identifier", max_records=SCROLL_CAP):
            b = h.get("isBasisFor") or {}
            if isinstance(b, list):
                b = b[0] if b else {}
            if str(b.get("identifier", "")).startswith("GSE"):
                by[b["identifier"]] += 1
    return n_human, n_any, [(g, k) for g, k in by.most_common() if k >= MIN_PER_SERIES]


rows = []
for rank, drug in enumerate(top, 1):
    n_human, n_any, cands = nde_candidates(drug)
    hits = [h for h in (confirm_arm(g, drug) for g, _ in cands[:8]) if h]
    # Prefer a series with an explicit control arm; those are directly analysable.
    def has_control(c):
        return any(CONTROL.search(v) and drug.lower() not in v.lower() for v in c["values"])

    # An arm declared in `characteristics_ch1` is a deliberate annotation; one
    # recovered from a sample title is inferred from a naming convention, so
    # prefer the former when a drug has both.
    # Human first. A mouse or zebrafish series is worth showing when it is all
    # that exists -- pemirolast has only GSE157167 -- but it should never
    # displace a human one, and unranked it does: roflumilast's best series
    # moved from BEAS-2B airway epithelium to mouse pro-B cells.
    def rank_key(c):
        return (c.get("organism") != "Homo sapiens",
                c["arm_field"] != "characteristics_ch1", not has_control(c),
                -c["n_samples_mentioning"])

    hits = sorted(hits, key=rank_key)
    best = hits[0] if hits else None
    arm = ctrl = None
    n_arm = n_ctrl = 0
    etype, is_expr = experiment_type(best["gse"]) if best else ("—", False)
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
        n_arm = best["values"].get(arm, 0)
        n_ctrl = best["values"].get(ctrl, 0) if ctrl else 0
    rows.append({
        "rank": rank, "drug": drug, "score": round(score[drug], 3),
        "system": describe_system(best["gse"]) if best else None,
        "genes": sorted(genes[drug]),
        "nde_samples_human": n_human, "nde_samples_any": n_any,
        "candidate_series": len(cands),
        "confirmed_series": len(hits),
        "best_series": best["gse"] if best else None,
        "arm_field": best["arm_field"] if best else None,
        "experiment_type": etype, "is_expression": is_expr,
        "n_arm": n_arm, "n_control": n_ctrl,
        "arm": arm, "control": ctrl,
        "all_confirmed": [c["gse"] for c in hits],
    })

out = Path("results/asthma_top20_datasets.json")
out.write_text(json.dumps({"pk": PK, "top_n": TOP_N, "rows": rows}, indent=2))

with_data = [r for r in rows if r["confirmed_series"]]
with_ctrl = [r for r in with_data if r["control"]]
print(f"top {TOP_N} asthma answers: {len(with_data)} have >=1 confirmed perturbation series, "
      f"{len(with_ctrl)} with an explicit control arm")
print(f"{'#':<3}{'drug':22s}{'score':>5}{'human':>7}{'any':>7}  {'series':11s}"
      f"{'type':>16s}{'n':>8s}  system")
for r in rows:
    print(f"{r['rank']:<3}{r['drug'][:20]:22s}{r['score']:>5.2f}"
          f"{r['nde_samples_human']:>7,}{r['nde_samples_any']:>7,}"
          f"  {r['best_series'] or '—':11s}{r['experiment_type']:>16s}"
          f"{('  ' + str(r['n_arm']) + 'v' + str(r['n_control'])) if r['n_arm'] else '':>8s}"
          f"  {(r['system'] or '')[:40]}")
    if r["arm"]:
        print(f"      {r['arm'][:58]}  VS  {(r['control'] or '(no control level)')[:34]}")
print(f"\nwrote {out}")
