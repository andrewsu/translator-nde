"""Identifier bridge between Translator CURIEs and NDE's text/ontology fields.

Translator speaks CURIEs; NDE speaks MONDO local IDs (for disease) and free text
(for everything else). Three jobs:

1. disease CURIE  -> NDE ``healthCondition.identifier`` (bare local ID, no prefix)
2. drug CURIE     -> synonym list for free-text matching
3. gene CURIE     -> Ensembl ID + symbol, which is how GXA keys genes
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Any, Iterable

import requests

NODENORM = "https://nodenormalization-sri.renci.org/1.5/get_normalized_nodes"
NAMERES = "https://name-resolution-sri.renci.org/lookup"

# Node Normalizer labels are noisy -- PubChem CID stubs and bare registry numbers
# are not usable as free-text search terms.
_JUNK_LABEL = re.compile(
    r"^(CID\s*\d+|SID\s*\d+|\d+[-\d]*|[A-Z0-9]{8,}|UNII[-:\s].*)$", re.IGNORECASE
)


def mondo_local_id(curie: str) -> str:
    """``MONDO:0018076`` -> ``0018076``.

    NDE stores ``healthCondition.identifier`` as the bare local ID and keeps the
    prefix separately in ``inDefinedTermSet``, so the prefix must be stripped
    before joining.
    """
    return curie.split(":", 1)[1] if ":" in curie else curie


@dataclass
class IdResolver:
    timeout: int = 120
    pause: float = 0.2
    session: requests.Session = field(default_factory=requests.Session)
    _norm_cache: dict[str, Any] = field(default_factory=dict)

    # ---------------------------------------------------------------- nodenorm

    def normalize(
        self,
        curies: Iterable[str],
        *,
        conflate: bool = True,
        drug_chemical_conflate: bool = True,
        batch: int = 500,
    ) -> dict[str, Any]:
        """Batch Node Normalizer lookup. Unresolvable CURIEs map to None.

        ``conflate`` merges gene/protein cliques (so an NCBIGene and its
        UniProtKB counterpart normalize together) -- necessary because Translator
        conflates genes and proteins and GXA keys on Ensembl.
        """
        curies = [c for c in dict.fromkeys(curies) if c]
        todo = [c for c in curies if c not in self._norm_cache]
        for i in range(0, len(todo), batch):
            chunk = todo[i : i + batch]
            r = self.session.post(
                NODENORM,
                json={
                    "curies": chunk,
                    "conflate": conflate,
                    "drug_chemical_conflate": drug_chemical_conflate,
                },
                timeout=self.timeout,
            )
            r.raise_for_status()
            self._norm_cache.update(r.json())
            time.sleep(self.pause)
        return {c: self._norm_cache.get(c) for c in curies}

    def clique_members(self, curie: str, prefix: str) -> list[str]:
        """All equivalent identifiers of a given prefix, e.g. ENSEMBL for a gene."""
        node = self.normalize([curie]).get(curie)
        if not node:
            return []
        return [
            e["identifier"]
            for e in node.get("equivalent_identifiers", [])
            if e["identifier"].upper().startswith(prefix.upper() + ":")
        ]

    def canonical_label(self, curie: str) -> str | None:
        node = self.normalize([curie]).get(curie)
        return (node or {}).get("id", {}).get("label")

    # --------------------------------------------------------------- nameres

    def synonyms(self, curie: str, *, limit: int = 1) -> list[str]:
        """Free-text synonyms for a CURIE, for NDE text matching.

        Name Resolver's ``synonyms`` list is much cleaner than Node Normalizer's
        equivalent-identifier labels (which include ``CID 152743144``-style
        stubs), so it is preferred and the normalizer is only the fallback.
        """
        label = self.canonical_label(curie)
        out: list[str] = []
        if label:
            out.append(label)
            r = self.session.get(
                NAMERES,
                params={"string": label, "autocomplete": "false", "limit": limit},
                timeout=self.timeout,
            )
            if r.ok:
                for rec in r.json():
                    if rec.get("curie") == curie or not out[1:]:
                        out.extend(rec.get("synonyms", []))
            time.sleep(self.pause)
        if not out:  # fall back to normalizer labels
            node = self.normalize([curie]).get(curie) or {}
            out = [
                e["label"] for e in node.get("equivalent_identifiers", []) if e.get("label")
            ]
        return _dedupe_clean(out)

    def lookup(self, text: str, *, biolink_type: str | None = None) -> str | None:
        """Free text -> best-match CURIE (e.g. ``"tuberculosis"`` -> MONDO:0018076)."""
        r = self.session.get(
            NAMERES,
            params={
                "string": text,
                "autocomplete": "false",
                "limit": 1,
                "biolink_type": biolink_type,
            },
            timeout=self.timeout,
        )
        r.raise_for_status()
        hits = r.json()
        time.sleep(self.pause)
        return hits[0]["curie"] if hits else None

    # ------------------------------------------------------------------ genes

    def gene_keys(self, curie: str) -> dict[str, Any]:
        """Gene CURIE -> the keys GXA can be queried on.

        GXA's ``observationAbout`` carries an Ensembl ID and a symbol, so those
        are what we need out of whatever Translator hands us (NCBIGene, HGNC or
        UniProtKB).
        """
        return {
            "curie": curie,
            "symbol": self.canonical_label(curie),
            "ensembl": [c.split(":", 1)[1] for c in self.clique_members(curie, "ENSEMBL")],
        }


def _dedupe_clean(values: Iterable[str]) -> list[str]:
    """Case-insensitive dedupe, dropping junk labels and over-long strings."""
    seen: dict[str, str] = {}
    for v in values:
        v = (v or "").strip()
        if not v or len(v) > 60 or _JUNK_LABEL.match(v):
            continue
        seen.setdefault(v.lower(), v)
    return list(seen.values())
