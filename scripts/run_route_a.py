"""Join Translator drug->gene edges against GXA DE contrasts (Route A).

Reads a paths.json produced by run_ars.py and, for each distinct drug->gene
edge, asks NDE staging whether any differential-expression contrast measured
that gene under that drug -- and if so, whether the direction agrees.
"""
import json, sys, time
from collections import Counter
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.gxa import GXAMatcher
from translator_nde.ids import IdResolver

paths_file = Path(sys.argv[1])
out_file = paths_file.parent / "route_a.json"
limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0

doc = json.loads(paths_file.read_text())
paths = doc["paths"]

# One evaluation per distinct drug->gene edge; keep an asserted direction if any
# path for that pair carries one (ARAX omits qualifiers entirely).
# ars-ars-agent is the ARS *merged* result set; it overlaps the individual ARAs,
# so counting it alongside them double-counts provenance. Edges are deduped by
# (drug, gene) regardless, but the agent list should stay honest.
MERGED_AGENT = "ars-ars-agent"

edges: dict[tuple, dict] = {}
for p in paths:
    if not p.get("drug_name") or not p.get("gene_name"):
        continue
    key = (p["drug_name"], p["gene_name"])
    e = edges.setdefault(key, {"drug": p["drug_name"], "gene": p["gene_name"],
                               "drug_curie": p["drug"], "gene_curie": p["gene"],
                               "direction": None, "aspect": None, "agents": set(),
                               "predicates": set()})
    e["agents"].add(p["agent"])
    e["predicates"].add(p["drug_gene_predicate"])
    if p.get("direction") and not e["direction"]:
        e["direction"], e["aspect"] = p["direction"], p["aspect"]

items = list(edges.values())
if limit:
    items = items[:limit]
print(f"{len(paths)} paths -> {len(edges)} distinct drug->gene edges"
      f"{f' (evaluating first {limit})' if limit else ''}\n", flush=True)

def drug_like(name: str) -> bool:
    """Cheap filter: IUPAC strings and synthetic peptides never match NDE text."""
    return len(name) <= 30 and name.count("-") <= 3 and "(" not in name


matcher, resolver = GXAMatcher(), IdResolver()
results, skipped, t0 = [], [], time.time()
for i, e in enumerate(items, 1):
    syns = [e["drug"]]
    if not drug_like(e["drug"]):
        skipped.append(e["drug"])
        print(f"  [{i}/{len(items)}] {e['drug'][:28]:28s} -> {e['gene']:8s} "
              f"skip (not drug-like)", flush=True)
        continue
    try:
        # Only pay for synonym expansion if the plain label finds nothing.
        if matcher.client.count(matcher.build_query(syns, [e["gene"]])) == 0:
            extra = resolver.synonyms(e["drug_curie"])
            syns = list(dict.fromkeys(syns + extra))[:12]
    except Exception as exc:
        print(f"  [{i}/{len(items)}] {e['drug'][:28]:28s} synonym lookup failed: {exc}",
              flush=True)
    try:
        ev = matcher.evaluate(
            drug_label=e["drug"], gene_label=e["gene"], drug_synonyms=syns,
            gene_terms=[e["gene"]], asserted_direction=e["direction"], max_records=300,
        )
    except Exception as exc:
        print(f"  [{i}/{len(items)}] {e['drug'][:28]:28s} -> {e['gene']:8s} ERROR {exc}",
              flush=True)
        continue
    rec = ev.to_dict()
    rec.update(agents=sorted(e["agents"]), predicates=sorted(e["predicates"]),
               drug_curie=e["drug_curie"], gene_curie=e["gene_curie"],
               n_synonyms_used=len(syns))
    results.append(rec)
    mark = "*" if ev.n_contrasts else " "
    if True:
        print(f" {mark}[{i}/{len(items)}] {e['drug'][:28]:28s} -> {e['gene']:8s} "
              f"n={ev.n_contrasts:4d} agree={ev.n_agree:4d} dis={ev.n_disagree:4d} "
              f"{ev.verdict:12s} log2FC={ev.median_log2fc}", flush=True)

out_file.write_text(json.dumps(
    {"source": str(paths_file), "disease": doc.get("disease"),
     "nde_build": matcher.client.build_info(), "results": results}, indent=2))

print(f"\n--- {len(results)} edges evaluated in {time.time()-t0:.0f}s ---")
print("verdicts:", dict(Counter(r["verdict"] for r in results)))
covered = [r for r in results if r["n_contrasts"] > 0]
informative = [r for r in results
               if r["n_contrasts"] > 0 or r["verdict"] == "tested_not_significant"]
print(f"covered by >=1 GXA contrast: {len(covered)}/{len(results)}"
      f" ({100*len(covered)/max(len(results),1):.0f}%)")
testable = [r for r in results if r["asserted_direction"]]
print(f"edges with ANY GXA information (incl. tested-not-significant): "
      f"{len(informative)}/{len(results)} "
      f"({100*len(informative)/max(len(results),1):.0f}%)")
print(f"edges carrying an asserted direction: {len(testable)}/{len(results)}")
print(f"skipped as not drug-like: {len(skipped)}")
print(f"wrote {out_file}")
