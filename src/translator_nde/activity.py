"""Route D: test a Translator drug->gene edge against *activity* data.

Routes A and B ask expression data whether a drug changes a gene's abundance.
For the majority of Translator drug->gene edges that is the wrong question --
the edges assert changes in **activity**, and a kinase inhibitor or receptor
antagonist does not move its target's transcript. This module asks the question
the assertion actually makes, of assays that actually measure it:

* **PubChem BioAssay** -- screening and dose-response data, keyed on NCBI Gene
  ID. Translator emits ``NCBIGene:7124``; PubChem takes the bare ``7124``. No
  text matching anywhere. Crucially it reports **Inactive** as well as Active
  outcomes, so a compound tested against a target and found inert is a genuine
  negative -- something the expression routes structurally never had, because
  GXA only stores significant results.

* **ChEMBL** -- manually curated mechanism of action, *typed by action*
  (INHIBITOR / AGONIST / ANTAGONIST / ...). That vocabulary carries the same
  semantics as Translator's ``object_direction_qualifier``, which makes this a
  like-for-like check rather than a proxy.

Neither source is indexed by NDE (its catalog has LINCS and ReframeDB but no
ChEMBL, PubChem or BindingDB), so Route D deliberately reaches outside NDE.
That gap is itself a finding, not an inconvenience.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Literal

import requests

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# ChEMBL action types grouped by the direction of effect on target activity, so
# they can be compared against Translator's object_direction_qualifier.
_DECREASING = frozenset({
    "INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR",
    "NEGATIVE MODULATOR", "DISRUPTING AGENT", "INVERSE AGONIST", "DOWNREGULATOR",
    "SEQUESTERING AGENT", "CHELATING AGENT", "DEGRADER",
})
_INCREASING = frozenset({
    "AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR", "OPENER", "STABILISER", "RELEASING AGENT", "UPREGULATOR",
})
# MODULATOR, BINDING AGENT, SUBSTRATE, OTHER, HYDROLYTIC ENZYME etc. assert an
# interaction without a direction -- they confirm engagement, not sign.

Verdict = Literal[
    "mechanism_agrees",     # curated mechanism, action type matches the asserted direction
    "mechanism_disagrees",  # curated mechanism, action type opposes it
    "mechanism_untyped",    # curated mechanism, but no directional claim on one side
    "binding_confirmed",    # measured potency against the target, no curated mechanism
    "measured_inactive",    # compound WAS tested against the target and was inert
    "not_tested",           # target has activity data; this compound is not in it
    "no_activity_data",     # no PubChem assays and no ChEMBL target for the gene
    "no_compound_id",       # drug has no CID and no ChEMBL id (typically a biologic)
    "fetch_failed",         # PubChem could not be reached -- not an observation
]


@dataclass
class Assay:
    """One PubChem BioAssay activity row, flattened."""

    aid: str
    cid: str | None
    gene_id: str | None          # NCBI Gene ID of the assay target
    outcome: str | None          # Active / Inactive / Inconclusive / Unspecified
    activity_name: str | None    # IC50 / Potency / Kd / ...
    value_um: float | None
    assay_name: str | None
    target_accession: str | None
    pubmed_id: str | None


@dataclass
class Mechanism:
    """One curated ChEMBL mechanism-of-action row."""

    molecule_chembl_id: str
    target_chembl_id: str | None
    action_type: str | None
    mechanism_of_action: str | None
    direct_interaction: bool
    max_phase: float | None
    refs: list[str] = field(default_factory=list)
    # UniProt accessions of the mechanism target's components. ChEMBL often
    # records a mechanism against a PROTEIN FAMILY or COMPLEX rather than the
    # single protein -- aspirin's target is "Cyclooxygenase" (CHEMBL2094253),
    # not PTGS2 -- so target_chembl_id equality alone silently misses them.
    target_accessions: list[str] = field(default_factory=list)

    @property
    def implied_direction(self) -> str | None:
        """ChEMBL action type -> Biolink direction qualifier vocabulary."""
        a = (self.action_type or "").upper()
        if a in _DECREASING:
            return "decreased"
        if a in _INCREASING:
            return "increased"
        return None


@dataclass
class ActivityEvidence:
    """Aggregated activity evidence for one Translator drug->gene edge."""

    drug: str
    drug_name: str
    gene: str
    gene_name: str
    asserted_direction: str | None
    asserted_aspect: str | None
    verdict: Verdict

    drug_cid: str | None = None
    drug_chembl: str | None = None
    target_chembl: str | None = None

    # PubChem
    n_active: int = 0
    n_inactive: int = 0
    n_inconclusive: int = 0
    best_potency_um: float | None = None
    best_potency_type: str | None = None
    assay_ids: list[str] = field(default_factory=list)

    # ChEMBL
    chembl_action_type: str | None = None
    chembl_moa: str | None = None
    max_phase: float | None = None
    pchembl_max: float | None = None
    n_bioactivities: int = 0

    # True when PubChem could not be reached for this compound, so absence of
    # rows is a fetch failure rather than evidence the compound was untested.
    pubchem_error: bool = False

    assays: list[Assay] = field(default_factory=list)

    def to_dict(self, *, include_assays: bool = False) -> dict:
        d = asdict(self)
        if not include_assays:
            d.pop("assays")
        return d


class _Cached:
    """Shared on-disk JSON cache. PubChem gene pulls run to several MB, and both
    APIs are rate-limited, so nothing is fetched twice."""

    def __init__(self, cache_dir: str | Path, pause: float = 0.25, timeout: int = 120):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.pause = pause
        self.timeout = timeout
        self.session = requests.Session()

    # PubChem answers a large gene pull with PUGREST.Timeout (504) or
    # PUGREST.ServerBusy (503) under load, intermittently and for genes that
    # succeed on a retry, so a single failure must not be read as "no data".
    _RETRY_STATUS = frozenset({429, 500, 502, 503, 504})

    def _fetch(self, url: str, params: dict | None = None, *, tries: int = 4) -> Any:
        delay = 2.0
        for attempt in range(tries):
            r = self.session.get(url, params=params, timeout=self.timeout)
            if r.status_code == 404:
                time.sleep(self.pause)
                return None
            if r.status_code in self._RETRY_STATUS and attempt < tries - 1:
                time.sleep(delay)
                delay *= 2
                continue
            time.sleep(self.pause)
            r.raise_for_status()
            try:
                return r.json()
            except ValueError:
                # A 200 with an unparseable body means the response was cut
                # short -- indistinguishable from a 504 in effect, so retry.
                if attempt == tries - 1:
                    raise
                time.sleep(delay)
                delay *= 2
        raise RuntimeError("unreachable")

    def _get(self, key: str, url: str, params: dict | None = None) -> Any:
        path = self.cache_dir / f"{key}.json"
        if path.exists():
            return json.loads(path.read_text())
        payload = self._fetch(url, params)
        path.write_text(json.dumps(payload))
        return payload


class PubChemBioAssay(_Cached):
    """Compound activity rows, keyed on PubChem CID and filtered by NCBI Gene ID.

    Queried **per compound**, not per gene, though PubChem offers both. The
    gene-keyed ``/gene/geneid/{id}/concise`` route joins just as exactly but is
    unusable at scale: CFTR returns 337 MB and DRD2 returns a 436 MB body that
    arrives truncated and fails to parse, so heavily screened targets -- exactly
    the interesting ones -- silently yield nothing. The compound-keyed
    ``assaysummary`` view carries a ``Target GeneID`` column, so the join is
    still an exact integer match with no text anywhere, at ~300 KB per compound.
    Spot-checked identical: imatinib/ABL1 154 active + 10 inactive either way.
    """

    # Columns kept, in the order they are stored in the compact cache.
    _COLS = ("AID", "CID", "Target GeneID", "Activity Outcome", "Activity Name",
             "Activity Value [uM]", "Target Accession", "PubMed ID")

    def compound_assays(self, cid: str) -> list[Assay]:
        path = self.cache_dir / f"pubchem_cid_{cid}.json"
        if path.exists():
            doc = json.loads(path.read_text())
        else:
            doc = self._compact(
                self._fetch(f"{PUBCHEM}/compound/cid/{cid}/assaysummary/JSON")
            )
            path.write_text(json.dumps(doc, separators=(",", ":")))
        names = doc["aid_names"]
        return [
            Assay(
                aid=r[0], cid=r[1], gene_id=r[2], outcome=r[3], activity_name=r[4],
                value_um=r[5], target_accession=r[6], pubmed_id=r[7],
                assay_name=names.get(r[0]),
            )
            for r in doc["rows"]
        ]

    def assays_for(self, cid: str, ncbigene: str) -> list[Assay]:
        """Rows where this compound was assayed against this gene's product."""
        gid = ncbigene.split(":")[-1]
        return [a for a in self.compound_assays(cid) if a.gene_id == gid]

    def potency_by_type(self, cid: str, ncbigene: str) -> dict[str, float]:
        """Most potent Active measurement per assay type, in micromolar.

        Reporting a single minimum across all types is misleading: it can hand
        back a Kd for a compound that also has an IC50, and the meaningful
        readout depends on the mechanism -- EC50/AC50 for an agonist, IC50/Ki
        for an inhibitor.
        """
        out: dict[str, float] = {}
        for a in self.assays_for(cid, ncbigene):
            if a.value_um is None or (a.outcome or "").lower() != "active":
                continue
            t = a.activity_name or "unspecified"
            if t not in out or a.value_um < out[t]:
                out[t] = a.value_um
        return out

    @classmethod
    def _compact(cls, payload: Any) -> dict:
        """Shrink the dump before caching: lift the repeated assay-name string
        into an AID lookup and drop RNAi knockdown rows, which test no compound
        and would otherwise inflate the "tested" denominator."""
        table = (payload or {}).get("Table") or {}
        cols = (table.get("Columns") or {}).get("Column") or []
        idx = {c: i for i, c in enumerate(cols)}
        rnai_i, name_i = idx.get("RNAi"), idx.get("Assay Name")
        take = [idx.get(c) for c in cls._COLS]
        rows, names = [], {}
        for row in table.get("Row") or []:
            c = row.get("Cell") or []

            def cell(i):
                return (c[i] or None) if i is not None and i < len(c) else None

            if rnai_i is not None and cell(rnai_i):
                continue
            out = [cell(i) for i in take]
            if not out[2]:  # no target gene -- nothing to join on
                continue
            out[5] = _num(out[5])  # Activity Value [uM]
            rows.append(out)
            if name_i is not None and out[0] not in names:
                names[out[0]] = cell(name_i)
        return {"aid_names": names, "rows": rows}


# Assay types that answer "how strongly does it engage the target", split by
# the mechanism they belong to. An AGONIST measured by IC50 is usually a
# counter-screen, not the readout of interest.
_AGONIST_ASSAYS = ("EC50", "AC50", "Potency")
_ANTAGONIST_ASSAYS = ("IC50", "Ki", "Kd", "Potency")


def preferred_potency(
    by_type: dict[str, float], action_type: str | None
) -> tuple[str, float] | None:
    """Pick the assay type that matches the mechanism, falling back to any."""
    a = (action_type or "").upper()
    order = _AGONIST_ASSAYS if a in _INCREASING else _ANTAGONIST_ASSAYS
    for t in order:
        if t in by_type:
            return t, by_type[t]
    if not by_type:
        return None
    t = min(by_type, key=by_type.get)
    return t, by_type[t]


def pchembl_to_um(pchembl: float | None) -> float | None:
    """pChEMBL is -log10(molar); convert to micromolar for comparability."""
    return None if pchembl is None else 10 ** (6 - pchembl)


class ChEMBLClient(_Cached):
    """Curated mechanism and potency, joined on UniProt accession or gene symbol."""

    def target_for_gene(
        self, *, uniprot: str | None = None, symbol: str | None = None
    ) -> str | None:
        """UniProt accession is an exact join and is tried first; the gene-symbol
        synonym lookup is the fallback for genes Node Normalizer gives no
        UniProtKB member for."""
        if uniprot:
            hit = self._target_query(
                f"chembl_target_acc_{uniprot}", {"target_components__accession": uniprot}
            )
            if hit:
                return hit
        if symbol:
            return self._target_query(
                f"chembl_target_sym_{symbol}",
                {
                    "target_components__target_component_synonyms__component_synonym":
                        symbol,
                    "organism": "Homo sapiens",
                },
            )
        return None

    def _target_query(self, key: str, params: dict) -> str | None:
        payload = self._get(
            key, f"{CHEMBL}/target.json",
            {**params, "target_type": "SINGLE PROTEIN", "limit": 5},
        )
        targets = (payload or {}).get("targets") or []
        return targets[0]["target_chembl_id"] if targets else None

    def target_components(self, target_chembl_id: str) -> list[str]:
        """UniProt accessions making up a ChEMBL target (1 for a single protein,
        several for a family or complex)."""
        payload = self._get(
            f"chembl_target_{target_chembl_id}",
            f"{CHEMBL}/target/{target_chembl_id}.json",
        )
        return [
            c["accession"]
            for c in (payload or {}).get("target_components") or []
            if c.get("accession")
        ]

    def mechanisms(self, molecule_chembl_id: str) -> list[Mechanism]:
        """All curated mechanisms for a compound, across every target.

        Fetched per *molecule* rather than per target: a molecule has a handful
        of mechanisms, whereas a well-studied target has hundreds of them.

        Queried on both ``molecule_chembl_id`` and ``parent_molecule_chembl_id``
        because ChEMBL attaches mechanisms to whichever form was studied --
        imatinib (CHEMBL941) has *no* mechanism rows of its own; all four are
        filed under the mesylate salt CHEMBL1642.
        """
        rows: dict[tuple, dict] = {}
        for key in ("molecule_chembl_id", "parent_molecule_chembl_id"):
            payload = self._get(
                f"chembl_mech_{key}_{molecule_chembl_id}", f"{CHEMBL}/mechanism.json",
                {key: molecule_chembl_id, "limit": 100},
            )
            for m in (payload or {}).get("mechanisms") or []:
                rows.setdefault((m.get("target_chembl_id"), m.get("action_type")), m)
        out = []
        for m in rows.values():
            target = m.get("target_chembl_id")
            out.append(
                Mechanism(
                    molecule_chembl_id=molecule_chembl_id,
                    target_chembl_id=target,
                    action_type=m.get("action_type"),
                    mechanism_of_action=m.get("mechanism_of_action"),
                    direct_interaction=bool(m.get("direct_interaction")),
                    max_phase=_num(m.get("max_phase")),
                    refs=[
                        r.get("ref_id", "") for r in (m.get("mechanism_refs") or [])
                    ],
                    target_accessions=self.target_components(target) if target else [],
                )
            )
        return out

    def pchembl(self, molecule_chembl_id: str, target_chembl_id: str) -> tuple[float | None, int]:
        """Best pChEMBL value for a compound-target pair, and how many were measured.

        pChEMBL is -log10 of a standardised IC50/Ki/EC50, so it is comparable
        across assay types; >=6 (1 uM) is the usual "meaningfully potent" cut.
        """
        payload = self._get(
            f"chembl_act_{molecule_chembl_id}_{target_chembl_id}",
            f"{CHEMBL}/activity.json",
            {
                "molecule_chembl_id": molecule_chembl_id,
                "target_chembl_id": target_chembl_id,
                "pchembl_value__isnull": "false",
                "limit": 100,
            },
        )
        vals = [
            _num(a.get("pchembl_value"))
            for a in (payload or {}).get("activities") or []
        ]
        vals = [v for v in vals if v is not None]
        return (max(vals) if vals else None), len(vals)


class ActivityMatcher:
    """Scores Translator drug->gene edges against PubChem and ChEMBL."""

    def __init__(
        self,
        cache_dir: str | Path = "data/activity",
        *,
        resolver: Any = None,
        pause: float = 0.25,
    ):
        from .ids import IdResolver

        self.pubchem = PubChemBioAssay(cache_dir, pause=pause)
        self.chembl = ChEMBLClient(cache_dir, pause=pause)
        self.resolver = resolver or IdResolver()
        self._gene_cache: dict[str, tuple] = {}
        self._mech_cache: dict[str, list[Mechanism]] = {}

    # ------------------------------------------------------------- resolution

    def compound_ids(self, drug_curie: str) -> tuple[str | None, str | None]:
        """Drug CURIE -> (PubChem CID, ChEMBL id) via Node Normalizer cliques."""
        cids = self.resolver.clique_members(drug_curie, "PUBCHEM.COMPOUND")
        chembls = self.resolver.clique_members(drug_curie, "CHEMBL.COMPOUND")
        return (
            cids[0].split(":", 1)[1] if cids else None,
            chembls[0].split(":", 1)[1] if chembls else None,
        )

    def gene_context(self, gene_curie: str, gene_name: str | None = None):
        """Gene CURIE -> (ChEMBL target id, UniProt accession)."""
        if gene_curie in self._gene_cache:
            return self._gene_cache[gene_curie]
        uniprot = [
            c.split(":", 1)[1]
            for c in self.resolver.clique_members(gene_curie, "UniProtKB")
        ]
        symbol = gene_name or self.resolver.canonical_label(gene_curie)
        acc = uniprot[0] if uniprot else None
        target = self.chembl.target_for_gene(uniprot=acc, symbol=symbol)
        self._gene_cache[gene_curie] = (target, acc)
        return target, acc

    def _mechanisms(self, chembl_id: str) -> list[Mechanism]:
        if chembl_id not in self._mech_cache:
            self._mech_cache[chembl_id] = self.chembl.mechanisms(chembl_id)
        return self._mech_cache[chembl_id]

    # ---------------------------------------------------------------- scoring

    def evaluate(
        self,
        *,
        drug: str,
        drug_name: str,
        gene: str,
        gene_name: str,
        asserted_direction: str | None = None,
        asserted_aspect: str | None = None,
    ) -> ActivityEvidence:
        cid, chembl_id = self.compound_ids(drug)
        ev = ActivityEvidence(
            drug=drug, drug_name=drug_name, gene=gene, gene_name=gene_name,
            asserted_direction=asserted_direction, asserted_aspect=asserted_aspect,
            verdict="no_compound_id", drug_cid=cid, drug_chembl=chembl_id,
        )
        if not cid and not chembl_id:
            return ev

        target, accession = self.gene_context(gene, gene_name)
        ev.target_chembl = target

        # ---- PubChem: was this exact compound tested against this exact gene?
        rows: list[Assay] = []
        if cid:
            try:
                rows = self.pubchem.assays_for(cid, gene)
            except Exception:
                ev.pubchem_error = True
        ev.assays = rows
        ev.assay_ids = sorted({a.aid for a in rows})
        for a in rows:
            o = (a.outcome or "").lower()
            if o == "active":
                ev.n_active += 1
            elif o == "inactive":
                ev.n_inactive += 1
            elif o == "inconclusive":
                ev.n_inconclusive += 1
        potent = [
            a for a in rows
            if a.value_um is not None and (a.outcome or "").lower() == "active"
        ]
        if potent:
            best = min(potent, key=lambda a: a.value_um)
            ev.best_potency_um = best.value_um
            ev.best_potency_type = best.activity_name

        # ---- ChEMBL: is there a curated, action-typed mechanism for the pair?
        mech = None
        if chembl_id:
            if target:
                ev.pchembl_max, ev.n_bioactivities = self.chembl.pchembl(chembl_id, target)
            mech = _match_mechanism(self._mechanisms(chembl_id), target, accession)
        if mech:
            ev.chembl_action_type = mech.action_type
            ev.chembl_moa = mech.mechanism_of_action
            ev.max_phase = mech.max_phase

        ev.verdict = _verdict(ev, mech, has_target_data=bool(rows or target))
        return ev

    def evaluate_paths(
        self, paths: Iterable[dict], *, dedupe: bool = True
    ) -> list[ActivityEvidence]:
        """Score every distinct drug->gene edge in a ``paths.json`` path list."""
        seen: set[tuple] = set()
        out = []
        for p in paths:
            key = (p["drug"], p["gene"], p.get("direction"), p.get("aspect"))
            if dedupe and key in seen:
                continue
            seen.add(key)
            out.append(
                self.evaluate(
                    drug=p["drug"], drug_name=p.get("drug_name") or p["drug"],
                    gene=p["gene"], gene_name=p.get("gene_name") or "",
                    asserted_direction=p.get("direction"),
                    asserted_aspect=p.get("aspect"),
                )
            )
        return out


def _match_mechanism(
    mechs: list[Mechanism], target: str | None, accession: str | None
) -> Mechanism | None:
    """Find the mechanism aimed at this gene.

    Exact target identity first; then component membership, which is what
    catches mechanisms filed against a protein family or complex. A directional
    action type is preferred over an untyped one when both match.
    """
    hits = [
        m for m in mechs
        if (target and m.target_chembl_id == target)
        or (accession and accession in m.target_accessions)
    ]
    if not hits:
        return None
    return next((m for m in hits if m.implied_direction), hits[0])


def _verdict(
    ev: ActivityEvidence, mech: Mechanism | None, *, has_target_data: bool
) -> Verdict:
    if mech:
        implied = mech.implied_direction
        if implied and ev.asserted_direction:
            return (
                "mechanism_agrees"
                if implied == ev.asserted_direction
                else "mechanism_disagrees"
            )
        return "mechanism_untyped"
    if ev.pubchem_error and not ev.n_active and ev.pchembl_max is None:
        # Absence of rows here is a fetch failure, not a measured negative;
        # calling it "not tested" would launder an error into evidence.
        return "fetch_failed"
    if ev.n_active or (ev.pchembl_max is not None):
        return "binding_confirmed"
    if ev.n_inactive:
        # Tested and inert. Unlike a missing GXA record, this is real evidence
        # against the edge -- the negative control Route A never had.
        return "measured_inactive"
    if has_target_data:
        return "not_tested"
    return "no_activity_data"


def _num(v: Any) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
