"""GEO supplementary-file inspection for Route B (reanalysis).

NDE tells us a GEO series exists and what it is about, but its ``distribution``
field only points at the accession page -- it does not list supplementary files.
So whether a series is *re-analyzable without realigning raw reads* has to be
answered against GEO's FTP directory listing.

We want series that already ship a gene-level matrix. Raw counts are preferred
(proper count-based statistics); normalized matrices (TPM/FPKM/CPM) are usable
on a log scale but weaker; anything else means realignment, which is out of
scope here.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Literal

import requests

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

MatrixKind = Literal["raw_counts", "normalized", "processed_other", "none"]

# Ordered: first match wins, so raw counts beat normalized when both appear.
_RAW_COUNTS = re.compile(
    r"(raw[_.-]?count|gene[_.-]?count|count[_.-]?matrix|counts?[_.-]?table"
    r"|featurecounts|htseq|[_.-]counts?\.|rsem.*expected|expected[_.-]count)",
    re.IGNORECASE,
)
_NORMALIZED = re.compile(
    # (?<!non-) so "non-normalized" (raw array intensities) is not mislabelled.
    r"(tpm|fpkm|rpkm|[_.-]cpm[_.-]|(?<!non-)(?<!non_)normali[sz]ed|vst|rlog|deseq"
    r"|edger|expression[_.-]?matrix|gene[_.-]?expression|abundance)",
    re.IGNORECASE,
)
_PROCESSED = re.compile(r"\.(txt|tsv|csv|xlsx?)(\.gz)?$", re.IGNORECASE)

# Files that are never a usable gene-level matrix.
_NOT_MATRIX = re.compile(
    r"(\.bed|\.bw|\.bigwig|\.wig|\.bam|\.sam|\.fastq|\.fq|\.vcf|\.mtx"
    r"|barcodes|features\.tsv|readme|filelist|\.pdf|\.png|_RAW\.tar)",
    re.IGNORECASE,
)

_FILE_RE = re.compile(r'href="([^"?/][^"]*)"')


def suppl_url(gse: str) -> str:
    """GEO stores series under a 1000-series bucket: GSE89408 -> GSE89nnn/."""
    return f"{FTP}/{gse[:-3]}nnn/{gse}/suppl/"


@dataclass
class GEOSupplementary:
    gse: str
    files: list[str] = field(default_factory=list)
    matrix_kind: MatrixKind = "none"
    matrix_files: list[str] = field(default_factory=list)
    error: str | None = None

    @property
    def reanalyzable(self) -> bool:
        return self.matrix_kind in ("raw_counts", "normalized")


def classify(files: list[str]) -> tuple[MatrixKind, list[str]]:
    """Best gene-level matrix available, preferring raw counts."""
    usable = [f for f in files if not _NOT_MATRIX.search(f)]
    for kind, pattern in (("raw_counts", _RAW_COUNTS), ("normalized", _NORMALIZED)):
        hits = [f for f in usable if pattern.search(f)]
        if hits:
            return kind, sorted(hits)
    other = [f for f in usable if _PROCESSED.search(f)]
    return ("processed_other", sorted(other)) if other else ("none", [])


def inspect(
    gse: str, *, session: requests.Session | None = None, timeout: int = 30,
    pause: float = 0.34,
) -> GEOSupplementary:
    """List a series' supplementary files and classify what is there.

    ``pause`` defaults to ~3 req/s, NCBI's unauthenticated rate limit.
    """
    s = session or requests
    try:
        r = s.get(suppl_url(gse), timeout=timeout)
        time.sleep(pause)
        if r.status_code == 404:
            return GEOSupplementary(gse, error="no suppl directory")
        r.raise_for_status()
    except requests.RequestException as exc:
        return GEOSupplementary(gse, error=str(exc)[:120])

    files = sorted({
        f for f in _FILE_RE.findall(r.text)
        if not f.startswith(("/", "http")) and "." in f
    })
    kind, matrices = classify(files)
    return GEOSupplementary(gse, files=files, matrix_kind=kind, matrix_files=matrices)
