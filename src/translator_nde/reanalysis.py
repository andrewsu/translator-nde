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


def match_columns(df: pd.DataFrame, gsms: list[str]) -> list[str]:
    """Map GSM accessions onto matrix columns, which rarely use them verbatim."""
    cols = {str(c): str(c) for c in df.columns}
    out = []
    for g in gsms:
        if g in cols:
            out.append(g)
            continue
        hit = next((c for c in cols if g.lower() in c.lower()), None)
        if hit:
            out.append(hit)
    return out


@dataclass
class DEResult:
    gse: str
    n_treated: int
    n_control: int
    matched_treated: int
    matched_control: int
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
) -> DEResult:
    """Moderated t-test of treated vs control, positive logFC = up in treated."""
    try:
        df = load_matrix(matrix_path)
    except Exception as exc:
        return DEResult(gse, len(treated), len(control), 0, 0, error=f"load: {exc}"[:150])

    t_cols, c_cols = match_columns(df, treated), match_columns(df, control)
    if len(t_cols) < min_per_group or len(c_cols) < min_per_group:
        return DEResult(gse, len(treated), len(control), len(t_cols), len(c_cols),
                        error="too few samples matched to matrix columns")

    sub = df[t_cols + c_cols]
    # Counts need CPM filtering + log2; already-normalized matrices only need log2
    # if they still look linear.
    genes = _de.counts_to_log2cpm(sub) if is_counts else _de.maybe_log2(sub)
    genes = _de.collapse_by_symbol(genes)
    is_treated = np.array([True] * len(t_cols) + [False] * len(c_cols))
    try:
        table = _de.moderated_ttest(genes, is_treated)
    except Exception as exc:
        return DEResult(gse, len(treated), len(control), len(t_cols), len(c_cols),
                        error=f"de: {exc}"[:150])
    return DEResult(gse, len(treated), len(control), len(t_cols), len(c_cols),
                    table=table)
