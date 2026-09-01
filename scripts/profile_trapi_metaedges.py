"""Meta-edge profile of a TRAPI knowledge graph from an ARS child payload.

Reduces every KG edge to (subject category, predicate, object category) and
counts the unique combinations. Works on any ARS payload so ARAs can be compared.

Usage:
    python scripts/profile_trapi_metaedges.py <payload.json> [-o out.tsv] [--top N]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.translator import load_message


def primary_category(node: dict) -> str:
    """Most specific declared category, else the CURIE prefix, else Unknown."""
    cats = node.get("categories") or []
    if not cats:
        return "Unknown"
    # TRAPI lists ancestors too; the shortest name is usually the most specific
    # leaf (e.g. SmallMolecule over ChemicalEntity/NamedThing).
    return sorted(cats, key=lambda c: (c == "biolink:NamedThing", len(c)))[0]


def profile(message: dict) -> tuple[Counter, Counter, Counter, int]:
    kg = message.get("knowledge_graph") or {}
    nodes, edges = kg.get("nodes") or {}, kg.get("edges") or {}
    cat = {nid: primary_category(n) for nid, n in nodes.items()}

    meta = Counter()
    for e in edges.values():
        meta[(cat.get(e.get("subject"), "Unknown"),
              e.get("predicate", "?"),
              cat.get(e.get("object"), "Unknown"))] += 1

    return meta, Counter(cat.values()), Counter(
        e.get("predicate", "?") for e in edges.values()
    ), len(edges)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("payload")
    ap.add_argument("-o", "--out")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    msg = load_message(Path(args.payload))
    meta, types, preds, n_edges = profile(msg)

    print(f"{Path(args.payload).name}")
    print(f"  nodes={len((msg.get('knowledge_graph') or {}).get('nodes') or {})} "
          f"edges={n_edges} unique meta-edges={len(meta)} "
          f"node types={len(types)} predicates={len(preds)}")

    short = lambda s: s.replace("biolink:", "")
    print(f"\n  {'count':>6}  meta-edge")
    for (s, p, o), n in meta.most_common(args.top):
        print(f"  {n:>6}  {short(s)} -[{short(p)}]-> {short(o)}")

    print(f"\n  node types: " + ", ".join(
        f"{short(t)} ({n})" for t, n in types.most_common(12)))

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("subject_category\tpredicate\tobject_category\tcount\n")
            for (s, p, o), n in meta.most_common():
                fh.write(f"{s}\t{p}\t{o}\t{n}\n")
        print(f"\n  wrote {args.out}")


if __name__ == "__main__":
    main()
