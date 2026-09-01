"""Test Translator drug->gene edges against activity data (Route D).

Reads a paths.json from run_ars.py and, for each distinct drug->gene edge, asks
PubChem BioAssay and ChEMBL whether that exact compound was measured against
that exact target -- and whether ChEMBL types the interaction the same way
Translator's qualifier does.

Unlike Routes A and B this needs no text matching: the gene join is the bare
NCBI Gene ID Translator already emits, and the compound join is a Node
Normalizer clique member (PubChem CID / ChEMBL id).
"""
import json, sys, time
from collections import Counter
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.activity import ActivityMatcher

paths_file = Path(sys.argv[1])
out_file = paths_file.parent / "route_d.json"
limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0

doc = json.loads(paths_file.read_text())
paths = doc["paths"]

# ars-ars-agent is the ARS *merged* result set and overlaps the individual ARAs;
# edges are deduped by (drug, gene) but the agent list should stay honest.
edges: dict[tuple, dict] = {}
for p in paths:
    if not p.get("gene", "").startswith("NCBIGene:"):
        continue  # PubChem's gene endpoint keys on NCBI Gene ID
    key = (p["drug"], p["gene"])
    e = edges.setdefault(key, {
        "drug": p["drug"], "drug_name": p.get("drug_name") or p["drug"],
        "gene": p["gene"], "gene_name": p.get("gene_name") or "",
        "direction": None, "aspect": None, "agents": set(), "predicates": set(),
    })
    e["agents"].add(p["agent"])
    e["predicates"].add(p["drug_gene_predicate"])
    if p.get("direction") and not e["direction"]:
        e["direction"], e["aspect"] = p["direction"], p["aspect"]

items = list(edges.values())
if limit:
    items = items[:limit]
print(f"{len(paths)} paths -> {len(edges)} distinct drug->gene edges"
      f"{f' (evaluating first {limit})' if limit else ''}\n", flush=True)

matcher = ActivityMatcher()
results, t0 = [], time.time()
for i, e in enumerate(items, 1):
    try:
        ev = matcher.evaluate(
            drug=e["drug"], drug_name=e["drug_name"],
            gene=e["gene"], gene_name=e["gene_name"],
            asserted_direction=e["direction"], asserted_aspect=e["aspect"],
        )
    except Exception as exc:
        print(f"  [{i}/{len(items)}] {e['drug_name'][:28]:28s} -> {e['gene_name']:8s} "
              f"ERROR {type(exc).__name__}: {exc}", flush=True)
        continue
    rec = ev.to_dict()
    rec.update(agents=sorted(e["agents"]), predicates=sorted(e["predicates"]))
    results.append(rec)
    mark = "*" if ev.verdict.startswith("mechanism") else (
        "+" if ev.verdict == "binding_confirmed" else
        "-" if ev.verdict == "measured_inactive" else " ")
    print(f" {mark}[{i}/{len(items)}] {e['drug_name'][:28]:28s} -> {e['gene_name']:8s} "
          f"{ev.verdict:19s} act/inact={ev.n_active:4d}/{ev.n_inactive:4d} "
          f"pChEMBL={ev.pchembl_max} {ev.chembl_action_type or ''}", flush=True)

out_file.write_text(json.dumps(
    {"source": str(paths_file), "disease": doc.get("disease"), "results": results},
    indent=2))

n = max(len(results), 1)
print(f"\n--- {len(results)} edges evaluated in {time.time()-t0:.0f}s ---")
print("verdicts:", dict(Counter(r["verdict"] for r in results).most_common()))
informative = [r for r in results if r["verdict"] not in
               ("no_activity_data", "no_compound_id")]
evidence = [r for r in results if r["verdict"] in
            ("mechanism_agrees", "mechanism_disagrees", "binding_confirmed",
             "measured_inactive")]
print(f"edges with ANY activity information: {len(informative)}/{len(results)}"
      f" ({100*len(informative)/n:.0f}%)")
print(f"edges with a measured compound-target result: {len(evidence)}/{len(results)}"
      f" ({100*len(evidence)/n:.0f}%)")
directional = [r for r in results if r["asserted_direction"]]
scored = [r for r in directional if r["verdict"].startswith("mechanism_")]
print(f"edges carrying an asserted direction: {len(directional)}"
      f"  of which curated-mechanism-scored: {len(scored)}")
print(f"wrote {out_file}")
