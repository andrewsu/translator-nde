"""Route B: differential expression from a deposited GEO count matrix.

Three steps, each of which can fail independently and is reported separately:

1. **Arms** -- get per-GSM labels from NDE ``@type:Sample`` records. This is the
   step that used to be hand-written per series (see ``cls_30528()`` etc. in
   DN-meta-analysis); NDE now carries the labels as queryable metadata.
2. **Matrix** -- download the supplementary counts/expression matrix that
   ``geo.inspect()`` located.
3. **DE** -- moderated t-test via the vendored DN-meta-analysis library.

The statistics are deliberately not reimplemented; ``_de`` is that library
verbatim so results stay comparable to the published DN analysis.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests

from scipy import stats

from . import _de
from .geo import suppl_url
from .nde import PROD, NDEClient


@dataclass
class Sample:
    gsm: str
    name: str | None
    description: str | None
    properties: dict[str, str] = field(default_factory=dict)

    def text(self) -> str:
        """Everything we could match an arm label against, lowercased."""
        bits = [self.name or "", self.description or ""]
        bits += [f"{k}: {v}" for k, v in self.properties.items()]
        return " | ".join(bits).lower()


def fetch_samples(gse: str, client: NDEClient | None = None) -> list[Sample]:
    """Per-GSM metadata for a series, from NDE rather than GEO SOFT."""
    c = client or NDEClient(base_url=PROD)
    q = f'@type:Sample AND isBasisFor.identifier:"{gse}"'
    out: list[Sample] = []
    for hit in c.scroll(q, fields="identifier,name,description,additionalProperty"):
        ident = hit.get("identifier")
        ident = ident[0] if isinstance(ident, list) else ident
        props = {
            p.get("propertyID"): p.get("value")
            for p in (hit.get("additionalProperty") or [])
            if isinstance(p, dict) and p.get("propertyID")
        }
        out.append(Sample(gsm=ident, name=hit.get("name"),
                          description=hit.get("description"), properties=props))
    return out


def assign_arms(
    samples: list[Sample], treated: str, control: str
) -> tuple[list[str], list[str], list[str]]:
    """Split samples by regex against their NDE metadata text.

    Returns (treated GSMs, control GSMs, ambiguous GSMs). A sample matching both
    patterns is ambiguous, not treated -- the commonest real failure is a
    control arm whose label also names the drug (vehicle-for-X, or the tool
    compound present in both arms).
    """
    t_re, c_re = re.compile(treated, re.I), re.compile(control, re.I)
    t, c, amb = [], [], []
    for s in samples:
        txt = s.text()
        in_t, in_c = bool(t_re.search(txt)), bool(c_re.search(txt))
        if in_t and in_c:
            amb.append(s.gsm)
        elif in_t:
            t.append(s.gsm)
        elif in_c:
            c.append(s.gsm)
    return t, c, amb


def download_matrix(gse: str, filename: str, cache: Path | str = "data/geo") -> Path:
    cache = Path(cache) / gse
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / filename
    if dest.exists() and dest.stat().st_size:
        return dest
    r = requests.get(suppl_url(gse) + filename, timeout=300, stream=True)
    r.raise_for_status()
    with dest.open("wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    return dest


def load_matrix(path: Path) -> pd.DataFrame:
    """Read a gene x sample matrix, sniffing the delimiter."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        head = fh.readline()
    sep = "\t" if head.count("\t") >= head.count(",") else ","
    df = pd.read_csv(path, sep=sep, index_col=0, low_memory=False)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce").dropna(how="all")


def _norm(s: str) -> str:
    """Collapse to comparable form: lowercase alphanumerics only."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def match_columns(
    df: pd.DataFrame, gsms: list[str], samples: list[Sample] | None = None
) -> list[str]:
    """Map GSMs onto matrix columns by accession, then by sample title.

    Depositors name supplementary-matrix columns however they like, so the
    accession is frequently absent. GSE89408 for instance uses
    ``normal_tissue_1`` / ``RA_tissue_148`` while the GEO titles are
    ``healthy tissue 2`` -- neither accession nor title matches.
    """
    cols = [str(c) for c in df.columns]
    norm_cols = {_norm(c): c for c in cols}
    titles = {s.gsm: s.name for s in (samples or []) if s.name}

    out = []
    for g in gsms:
        hit = next((c for c in cols if g.lower() in c.lower()), None)
        if hit is None and g in titles:
            hit = norm_cols.get(_norm(titles[g]))
        if hit:
            out.append(hit)
    return out


def arms_from_columns(
    df: pd.DataFrame, treated: str, control: str
) -> tuple[list[str], list[str]]:
    """Last-resort arm assignment from the matrix column names themselves.

    Used when GSM/title matching fails. This loses the link back to NDE sample
    records, so callers must report that the weaker strategy was used.
    """
    t_re, c_re = re.compile(treated, re.I), re.compile(control, re.I)
    t, c = [], []
    for col in (str(x) for x in df.columns):
        in_t, in_c = bool(t_re.search(col)), bool(c_re.search(col))
        if in_t and not in_c:
            t.append(col)
        elif in_c and not in_t:
            c.append(col)
    return t, c


@dataclass
class DEResult:
    gse: str
    n_treated: int
    n_control: int
    matched_treated: int
    matched_control: int
    arm_source: str = "gsm"   # "gsm" | "title" | "matrix_columns"
    design: str = "unpaired"  # "unpaired" | "paired"
    n_pairs: int = 0
    warnings: list[str] = field(default_factory=list)
    table: pd.DataFrame | None = None
    error: str | None = None

    def gene(self, symbol: str) -> dict[str, Any] | None:
        """Look up one gene's result, tolerating case and Ensembl-versioned ids."""
        if self.table is None:
            return None
        idx = {str(i).upper().split(".")[0]: i for i in self.table.index}
        key = idx.get(symbol.upper())
        if key is None:
            return None
        row = self.table.loc[key]
        return {"gene": symbol, "logFC": float(row["logFC"]),
                "p": float(row["P.Value"]), "adj_p": float(row["adj.P.Val"]),
                "direction": "increased" if row["logFC"] > 0 else "decreased"}


def run_de(
    gse: str, matrix_path: Path, treated: list[str], control: list[str],
    *, is_counts: bool = True, min_per_group: int = 2,
    samples: list[Sample] | None = None,
    patterns: tuple[str, str] | None = None,
    subject_re: str | None = None,
    map_symbols: bool = True,
) -> DEResult:
    """Moderated t-test of treated vs control, positive logFC = up in treated."""
    try:
        df = load_matrix(matrix_path)
    except Exception as exc:
        return DEResult(gse, len(treated), len(control), 0, 0,
                        error=f"load: {exc}"[:150])

    t_cols = match_columns(df, treated, samples)
    c_cols = match_columns(df, control, samples)
    arm_source = "gsm"
    if (len(t_cols) < min_per_group or len(c_cols) < min_per_group) and patterns:
        t_cols, c_cols = arms_from_columns(df, *patterns)
        arm_source = "matrix_columns"
    if len(t_cols) < min_per_group or len(c_cols) < min_per_group:
        return DEResult(gse, len(treated), len(control), len(t_cols), len(c_cols),
                        error="too few samples matched to matrix columns")

    warnings: list[str] = []
    ratio = max(len(t_cols), len(c_cols)) / max(min(len(t_cols), len(c_cols)), 1)
    if ratio >= 2:
        warnings.append(
            f"arms badly unbalanced ({len(t_cols)} vs {len(c_cols)}); "
            f"check the arm patterns match every column naming variant")

    sub = df[t_cols + c_cols]
    # Counts need CPM filtering + log2; already-normalized matrices only need log2
    # if they still look linear.
    genes = _de.counts_to_log2cpm(sub) if is_counts else _de.maybe_log2(sub)

    if map_symbols and looks_ensembl(genes.index):
        mapping = map_ensembl_to_symbol(list(genes.index))
        if mapping:
            genes = genes.rename(index=lambda i: mapping.get(str(i).split(".")[0], str(i)))
            warnings.append(f"mapped {len(mapping)} Ensembl ids to symbols")
        else:
            warnings.append("index looks Ensembl but symbol mapping returned nothing")
    genes = _de.collapse_by_symbol(genes)

    pairs = pair_subjects(t_cols, c_cols, subject_re) if subject_re else []
    try:
        if len(pairs) >= 3:
            table = paired_moderated_ttest(genes, pairs)
            design, n_pairs = "paired", len(pairs)
        else:
            if subject_re:
                warnings.append(
                    f"paired analysis requested but only {len(pairs)} pairs matched; "
                    f"fell back to unpaired")
            is_treated = np.array([True] * len(t_cols) + [False] * len(c_cols))
            table = _de.moderated_ttest(genes, is_treated)
            design, n_pairs = "unpaired", 0
    except Exception as exc:
        return DEResult(gse, len(treated), len(control), len(t_cols), len(c_cols),
                        arm_source, error=f"de: {exc}"[:150], warnings=warnings)
    return DEResult(gse, len(treated), len(control), len(t_cols), len(c_cols),
                    arm_source, design, n_pairs, warnings, table=table)

# ---------------------------------------------------------------- paired designs


def pair_subjects(
    treated_cols: list[str], control_cols: list[str], subject_re: str
) -> list[tuple[str, str]]:
    """Pair treated/control columns by a subject id captured from the name.

    Pre/post designs sample the *same patient* twice. Analyzing them unpaired
    throws away the within-subject blocking and is badly underpowered -- between
    -patient variation in synovium or whole blood dwarfs the treatment effect.
    """
    rx = re.compile(subject_re, re.I)
    def key(c: str) -> str | None:
        m = rx.search(c)
        return (m.group(1) if m.groups() else m.group(0)).lower() if m else None

    ctrl = {key(c): c for c in control_cols if key(c)}
    return [(t, ctrl[key(t)]) for t in treated_cols if key(t) and key(t) in ctrl]


def paired_moderated_ttest(genes: pd.DataFrame, pairs: list[tuple[str, str]]):
    """Empirical-Bayes moderated paired t-test on within-subject differences.

    Same shrinkage as the two-group test in ``_de`` (reusing its ``fit_f_dist``),
    applied to the one-sample problem mean(treated - control) != 0.
    Positive logFC = up in the treated member of each pair.
    """
    diffs = np.column_stack([
        genes[t].to_numpy(float) - genes[c].to_numpy(float) for t, c in pairs
    ])
    n = diffs.shape[1]
    if n < 3:
        raise ValueError(f"need >=3 pairs, got {n}")
    mean_d = diffs.mean(1)
    s2 = diffs.var(1, ddof=1)
    dg = n - 1
    d0, s0_2 = _de.fit_f_dist(s2, dg)
    if np.isinf(d0):
        s2_post, df_tot = s2, dg
    else:
        s2_post = (d0 * s0_2 + dg * s2) / (d0 + dg)
        df_tot = dg + d0
    se = np.sqrt(s2_post / n)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = mean_d / se
    p = 2 * stats.t.sf(np.abs(t), df_tot)
    out = pd.DataFrame({"gene": genes.index, "logFC": mean_d, "t": t,
                        "P.Value": p, "n_pairs": n}).set_index("gene")
    out["adj.P.Val"] = _de.bh_fdr(out["P.Value"].to_numpy())
    return out.sort_values("P.Value")


# ------------------------------------------------------------- gene identifiers

_ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d{6,}", re.I)


def looks_ensembl(index: pd.Index) -> bool:
    sample = [str(i) for i in list(index)[:50]]
    return sum(bool(_ENSEMBL_RE.match(i)) for i in sample) > len(sample) / 2


def map_ensembl_to_symbol(ids: list[str], *, batch: int = 900) -> dict[str, str]:
    """Ensembl gene id -> HGNC symbol via MyGene.info.

    Needed because several GEO count matrices are indexed on Ensembl ids, which
    silently defeats symbol-based gene lookup (every marker reads "not in
    matrix").
    """
    stripped = sorted({str(i).split(".")[0] for i in ids})
    out: dict[str, str] = {}
    for i in range(0, len(stripped), batch):
        chunk = stripped[i : i + batch]
        r = requests.post("https://mygene.info/v3/query",
                          data={"q": ",".join(chunk), "scopes": "ensembl.gene",
                                "fields": "symbol", "species": "human"},
                          timeout=120)
        if not r.ok:
            continue
        for rec in r.json():
            if rec.get("symbol") and rec.get("query"):
                out.setdefault(rec["query"], rec["symbol"])
    return out
