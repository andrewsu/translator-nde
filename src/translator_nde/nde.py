"""Thin client for the NIAID Data Ecosystem Discovery API.

No existing Python client exists for NDE (everything upstream is curl/requests
one-offs), so this is deliberately small: a query wrapper, faceting, and scroll
pagination, with the two response gotchas handled.

Two deployments matter here:

* ``PROD``    -- ``Dataset`` (5.4M), ``Sample`` (8.7M), ``ComputationalTool``.
* ``STAGING`` -- additionally ``Inference`` (10.4M GXA/ImmuneSpace DE contrasts)
  and ``DataCollection``. These do **not** exist in prod and an ``@type:Inference``
  query there returns 0 hits *with no error*, so the base URL is always explicit.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Iterator

import requests

PROD = "https://api.data.niaid.nih.gov/v1"
STAGING = "https://api-staging.data.niaid.nih.gov/v1"

# @type values that only exist on staging (verified 2026-09-01, prod build 20260830).
STAGING_ONLY_TYPES = frozenset({"Inference", "DataCollection"})

# The API rejects from > 10000 with HTTP 400; deep paging must use scroll.
MAX_FROM = 10_000


class NDEError(RuntimeError):
    pass


@dataclass
class NDEClient:
    base_url: str = PROD
    timeout: int = 120
    retries: int = 3
    pause: float = 0.1
    session: requests.Session = field(default_factory=requests.Session)

    def _get(self, path: str, params: dict[str, Any]) -> dict:
        # Drop Nones so callers can pass optional params unconditionally.
        params = {k: v for k, v in params.items() if v is not None}
        last: Exception | None = None
        for attempt in range(self.retries):
            try:
                r = self.session.get(
                    f"{self.base_url}{path}", params=params, timeout=self.timeout
                )
                if r.status_code == 400:
                    raise NDEError(f"400 from NDE: {r.text[:300]}")
                r.raise_for_status()
                time.sleep(self.pause)
                return r.json()
            except (requests.RequestException, ValueError) as exc:  # retry transient
                last = exc
                time.sleep(2**attempt)
        raise NDEError(f"NDE request failed after {self.retries} tries: {last}")

    def query(
        self,
        q: str,
        *,
        fields: str | None = None,
        size: int = 10,
        frm: int = 0,
        sort: str | None = None,
        facets: str | None = None,
        facet_size: int | None = None,
        extra_filter: str | None = None,
    ) -> dict:
        """Single /query call. Returns the raw response dict."""
        if frm > MAX_FROM:
            raise NDEError(
                f"from={frm} exceeds the API cap of {MAX_FROM}; use scroll() instead"
            )
        self._warn_if_staging_only(q)
        return self._get(
            "/query",
            {
                "q": q,
                "fields": fields,
                "size": size,
                "from": frm or None,
                "sort": sort,
                "facets": facets,
                "facet_size": facet_size,
                "extra_filter": extra_filter,
            },
        )

    def count(self, q: str, *, extra_filter: str | None = None) -> int:
        """Total hits for a query, without fetching any."""
        return int(self.query(q, size=0, extra_filter=extra_filter)["total"])

    def facet(self, q: str, field_name: str, *, facet_size: int = 100) -> dict[str, int]:
        """Facet counts as an ordered {term: count} dict.

        Note NDE facet buckets use ``term``/``count``, not Elasticsearch's
        ``key``/``doc_count``. Faceting requires a keyword-mapped field; text
        fields silently return an empty bucket list rather than erroring.
        """
        resp = self.query(q, size=0, facets=field_name, facet_size=facet_size)
        buckets = resp.get("facets", {}).get(field_name, {}).get("terms", [])
        return {b["term"]: b["count"] for b in buckets}

    def scroll(
        self, q: str, *, fields: str | None = None, max_records: int | None = None
    ) -> Iterator[dict]:
        """Iterate all hits past the 10k `from` cap, 500 per page."""
        self._warn_if_staging_only(q)
        resp = self._get(
            "/query", {"q": q, "fields": fields, "fetch_all": "true"}
        )
        seen = 0
        while True:
            hits = resp.get("hits", [])
            if not hits:
                return
            for hit in hits:
                yield hit
                seen += 1
                if max_records is not None and seen >= max_records:
                    return
            scroll_id = resp.get("_scroll_id")
            if not scroll_id:
                return
            resp = self._get("/query", {"scroll_id": scroll_id})

    def build_info(self) -> dict:
        """Build date/version — pin this in result files for reproducibility.

        Best-effort: staging's /metadata carries a large per-source ``src`` map
        and has been observed to truncate mid-response, so a failure here must
        not block a run.
        """
        try:
            meta = self._get("/metadata", {})
        except NDEError as exc:
            return {"base_url": self.base_url, "error": str(exc)[:120]}
        info = {k: meta.get(k) for k in ("build_date", "build_version", "biothing_type")}
        info["base_url"] = self.base_url
        return info

    def _warn_if_staging_only(self, q: str) -> None:
        if self.base_url != STAGING:
            for t in STAGING_ONLY_TYPES:
                if f"@type:{t}" in q:
                    raise NDEError(
                        f"@type:{t} exists only on staging ({STAGING}); querying "
                        f"{self.base_url} would silently return 0 hits. "
                        f"Use NDEClient(base_url=STAGING)."
                    )


def lucene_or(field_name: str, values: list[str]) -> str:
    """OR-clause over one field, phrase-quoting multi-word values.

    Escapes the Lucene specials that actually appear in drug/gene synonyms
    (parentheses, brackets, colons, slashes and the like).
    """
    parts = []
    for v in values:
        v = v.strip()
        if not v:
            continue
        if any(c in v for c in ' \t()[]{}:/\\+-!^"~*?|&'):
            parts.append(f'{field_name}:"{v}"')
        else:
            parts.append(f"{field_name}:{v}")
    if not parts:
        raise ValueError("no usable values")
    return "(" + " OR ".join(parts) + ")"
