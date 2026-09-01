"""Phase 1: ARS creative-mode queries and drug->gene->disease path extraction.

Submits an inferred ``biolink:treats`` query for a disease, then reconstructs the
mechanistic paths that justify each answer.

Supporting paths hang off a TRAPI message in **two** places and both must be
followed:

1. ``results[].analyses[].support_graphs``  -> keys into ``auxiliary_graphs``
2. a knowledge-graph edge attribute ``biolink:support_graphs`` -> ditto.
   This is how an *inferred* ``treats`` edge points at the concrete chain that
   produced it, and it is the one that actually carries drug->gene->disease.

Auxiliary-graph ``edges`` arrays are explicitly unordered by the spec, so paths
are re-chained on subject/object rather than trusted as written.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import requests

ARS_SUBMIT = "https://ars-prod.transltr.io/ars/api/submit"
ARS_MESSAGES = "https://ars-prod.transltr.io/ars/api/messages"

DRUG_CATEGORIES = {"biolink:ChemicalEntity", "biolink:SmallMolecule", "biolink:Drug",
                   "biolink:MolecularEntity", "biolink:ChemicalMixture"}
GENE_CATEGORIES = {"biolink:Gene", "biolink:Protein", "biolink:GeneOrGeneProduct"}
DISEASE_CATEGORIES = {"biolink:Disease", "biolink:DiseaseOrPhenotypicFeature",
                      "biolink:PhenotypicFeature"}

TERMINAL_STATUSES = {"Done", "Error"}


def creative_treats_query(disease_curie: str) -> dict:
    """Creative-mode 'what chemicals treat this disease' TRAPI message."""
    return {
        "message": {
            "query_graph": {
                "nodes": {
                    "on": {"categories": ["biolink:Disease"], "ids": [disease_curie]},
                    "sn": {"categories": ["biolink:ChemicalEntity"]},
                },
                "edges": {
                    "t_edge": {
                        "subject": "sn",
                        "object": "on",
                        "predicates": ["biolink:treats"],
                        "knowledge_type": "inferred",
                    }
                },
            }
        }
    }


@dataclass
class Hop:
    subject: str
    subject_name: str | None
    predicate: str
    object: str
    object_name: str | None
    qualifiers: dict[str, str] = field(default_factory=dict)
    sources: list[str] = field(default_factory=list)

    @property
    def direction(self) -> str | None:
        """Biolink object_direction_qualifier, e.g. 'increased' / 'decreased'."""
        return self.qualifiers.get("biolink:object_direction_qualifier")

    @property
    def aspect(self) -> str | None:
        return self.qualifiers.get("biolink:object_aspect_qualifier")


@dataclass
class MechanisticPath:
    """A drug -> gene -> disease chain pulled out of a creative-mode answer."""

    drug: str
    drug_name: str | None
    gene: str
    gene_name: str | None
    disease: str
    disease_name: str | None
    drug_gene: Hop
    gene_disease: Hop
    score: float | None = None

    def to_dict(self) -> dict:
        return {
            "drug": self.drug, "drug_name": self.drug_name,
            "gene": self.gene, "gene_name": self.gene_name,
            "disease": self.disease, "disease_name": self.disease_name,
            "drug_gene_predicate": self.drug_gene.predicate,
            "direction": self.drug_gene.direction,
            "aspect": self.drug_gene.aspect,
            "drug_gene_sources": self.drug_gene.sources,
            "gene_disease_predicate": self.gene_disease.predicate,
            "score": self.score,
        }


class ARSClient:
    """Submit to and read from the Translator ARS.

    Per-ARA payloads run to tens of megabytes, so responses are written straight
    to disk and parsed from there rather than held in memory.
    """

    def __init__(self, cache_dir: Path | str = "data/ars", timeout: int = 300):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.session = requests.Session()

    def submit(self, disease_curie: str) -> str:
        r = self.session.post(
            ARS_SUBMIT, json=creative_treats_query(disease_curie), timeout=self.timeout
        )
        r.raise_for_status()
        pk = r.json()["pk"]
        (self.cache_dir / pk).mkdir(exist_ok=True)
        (self.cache_dir / pk / "query.json").write_text(
            json.dumps({"disease": disease_curie, "pk": pk}, indent=2)
        )
        return pk

    def trace(self, pk: str) -> dict:
        r = self.session.get(f"{ARS_MESSAGES}/{pk}", params={"trace": "y"},
                             timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def wait(self, pk: str, *, poll: int = 15, max_wait: int = 900,
             verbose: bool = True) -> dict:
        """Poll the parent trace until every child ARA is in a terminal state."""
        deadline = time.time() + max_wait
        while True:
            tr = self.trace(pk)
            children = tr.get("children", [])
            done = [c for c in children if c.get("status") in TERMINAL_STATUSES]
            if verbose:
                ready = sum(1 for c in done if (c.get("result_count") or 0) > 0)
                print(f"  [{int(time.time()-deadline+max_wait):>4}s] "
                      f"{len(done)}/{len(children)} ARAs finished, {ready} with results",
                      flush=True)
            if children and len(done) == len(children):
                return tr
            if time.time() > deadline:
                if verbose:
                    print("  timed out; returning partial trace", flush=True)
                return tr
            time.sleep(poll)

    def fetch_child(self, child_pk: str, *, agent: str = "") -> Path | None:
        """Download one ARA's TRAPI payload to disk. Returns the path."""
        out = self.cache_dir / f"{child_pk}.json"
        if out.exists():
            return out
        r = self.session.get(f"{ARS_MESSAGES}/{child_pk}", timeout=self.timeout,
                             stream=True)
        if not r.ok:
            return None
        with out.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        return out

    def fetch_all_children(self, pk: str, *, min_results: int = 1) -> dict[str, Path]:
        """Download every child ARA response that actually returned results."""
        tr = self.trace(pk)
        out: dict[str, Path] = {}
        for child in tr.get("children", []):
            if (child.get("result_count") or 0) < min_results:
                continue
            agent = (child.get("actor") or {}).get("agent", "?")
            path = self.fetch_child(child["message"], agent=agent)
            if path:
                out[agent] = path
        return out


# --------------------------------------------------------------------- parsing


def load_message(path: Path) -> dict:
    """Extract the TRAPI message from an ARS child record on disk."""
    doc = json.loads(Path(path).read_text())
    fields = doc.get("fields", doc)
    data = fields.get("data") or {}
    return data.get("message", {})


def _node_categories(node: dict) -> set[str]:
    return set(node.get("categories") or [])


def _qualifiers(edge: dict) -> dict[str, str]:
    return {
        q["qualifier_type_id"]: q["qualifier_value"]
        for q in edge.get("qualifiers") or []
        if "qualifier_type_id" in q
    }


def _sources(edge: dict) -> list[str]:
    return [
        s["resource_id"]
        for s in edge.get("sources") or []
        if s.get("resource_role") == "primary_knowledge_source"
    ]


def _support_graph_ids(edge: dict) -> list[str]:
    out: list[str] = []
    for attr in edge.get("attributes") or []:
        if attr.get("attribute_type_id") == "biolink:support_graphs":
            val = attr.get("value")
            out.extend(val if isinstance(val, list) else [val])
    return [v for v in out if v]


def collect_edge_ids(message: dict, result: dict) -> set[str]:
    """All KG edge ids reachable from one result, following both attach points."""
    kg_edges = (message.get("knowledge_graph") or {}).get("edges") or {}
    aux = message.get("auxiliary_graphs") or {}

    edge_ids: set[str] = set()
    aux_todo: list[str] = []

    for analysis in result.get("analyses") or []:
        for binding in (analysis.get("edge_bindings") or {}).values():
            for b in binding:
                if b.get("id"):
                    edge_ids.add(b["id"])
        aux_todo.extend(analysis.get("support_graphs") or [])
        # Pathfinder results bind paths rather than edges.
        for binding in (analysis.get("path_bindings") or {}).values():
            for b in binding:
                if b.get("id"):
                    aux_todo.append(b["id"])

    seen_aux: set[str] = set()
    # Edges can reference aux graphs, whose edges can reference more aux graphs.
    while aux_todo:
        aid = aux_todo.pop()
        if aid in seen_aux or aid not in aux:
            continue
        seen_aux.add(aid)
        for eid in aux[aid].get("edges") or []:
            edge_ids.add(eid)
            aux_todo.extend(_support_graph_ids(kg_edges.get(eid, {})))

    for eid in list(edge_ids):
        aux_todo.extend(_support_graph_ids(kg_edges.get(eid, {})))
    while aux_todo:
        aid = aux_todo.pop()
        if aid in seen_aux or aid not in aux:
            continue
        seen_aux.add(aid)
        for eid in aux[aid].get("edges") or []:
            edge_ids.add(eid)
            aux_todo.extend(_support_graph_ids(kg_edges.get(eid, {})))

    return edge_ids


def extract_paths(message: dict, disease_curie: str) -> Iterator[MechanisticPath]:
    """Yield drug -> gene -> disease chains from a creative-mode message."""
    kg = message.get("knowledge_graph") or {}
    nodes, edges = kg.get("nodes") or {}, kg.get("edges") or {}

    for result in message.get("results") or []:
        edge_ids = collect_edge_ids(message, result)
        # Re-chain rather than trusting order: aux-graph edge lists are unordered.
        by_subject: dict[str, list[tuple[str, dict]]] = {}
        for eid in edge_ids:
            e = edges.get(eid)
            if e:
                by_subject.setdefault(e["subject"], []).append((eid, e))

        score = None
        for analysis in result.get("analyses") or []:
            if analysis.get("score") is not None:
                score = analysis["score"]
                break

        for drug_id, first_hops in by_subject.items():
            if not (_node_categories(nodes.get(drug_id, {})) & DRUG_CATEGORIES):
                continue
            for _, e1 in first_hops:
                gene_id = e1["object"]
                if not (_node_categories(nodes.get(gene_id, {})) & GENE_CATEGORIES):
                    continue
                for _, e2 in by_subject.get(gene_id, []):
                    dis_id = e2["object"]
                    if dis_id != disease_curie and not (
                        _node_categories(nodes.get(dis_id, {})) & DISEASE_CATEGORIES
                    ):
                        continue
                    yield MechanisticPath(
                        drug=drug_id,
                        drug_name=nodes.get(drug_id, {}).get("name"),
                        gene=gene_id,
                        gene_name=nodes.get(gene_id, {}).get("name"),
                        disease=dis_id,
                        disease_name=nodes.get(dis_id, {}).get("name"),
                        drug_gene=Hop(
                            drug_id, nodes.get(drug_id, {}).get("name"),
                            e1.get("predicate", ""), gene_id,
                            nodes.get(gene_id, {}).get("name"),
                            _qualifiers(e1), _sources(e1),
                        ),
                        gene_disease=Hop(
                            gene_id, nodes.get(gene_id, {}).get("name"),
                            e2.get("predicate", ""), dis_id,
                            nodes.get(dis_id, {}).get("name"),
                            _qualifiers(e2), _sources(e2),
                        ),
                        score=score,
                    )


def dedupe_paths(paths: Iterator[MechanisticPath]) -> list[MechanisticPath]:
    """Collapse identical drug/gene/disease/predicate chains, keeping best score."""
    best: dict[tuple, MechanisticPath] = {}
    for p in paths:
        key = (p.drug, p.gene, p.disease, p.drug_gene.predicate,
               p.drug_gene.direction, p.gene_disease.predicate)
        cur = best.get(key)
        if cur is None or (p.score or 0) > (cur.score or 0):
            best[key] = p
    return list(best.values())
