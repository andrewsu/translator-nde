"""Route C: which *other* compounds move a Translator gene the useful way?

Routes A, B and D all ask whether an asserted edge is supported. Route C asks a
question Translator did not: given `Drug1 -inhibits-> Gene2 -associated_with->
Disease3`, the therapeutic hypothesis is *less Gene2 helps Disease3*. So query
GXA gene-first for any compound that lowers Gene2, and every hit Translator did
not already propose is a repurposing candidate.

This plays to GXA's one structured axis. `observationAbout` carries a gene
symbol and an Ensembl id, so the gene side needs no text matching -- which is
where Routes A and B lose most of their recall.
"""
import json, sys, time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.gxa import (
    GXAMatcher, _GOAL_DIRECTION, _VEHICLE_TERMS, perturbagen_kind,
)
from translator_nde.nde import lucene_or

paths_file = Path(sys.argv[1])
out_file = paths_file.parent / "route_c.json"
doc = json.loads(paths_file.read_text())
paths = doc["paths"]

# A gene's therapeutic goal comes from the Translator edge's own direction
# qualifier. Most edges carry none (they arrive via DGIdb/DrugBank/DrugCentral),
# so those genes are run in both directions and flagged: without a qualifier we
# know the gene matters for the disease but not which way to push it.
goal: dict[str, set[str]] = defaultdict(set)
proposed: dict[str, set[str]] = defaultdict(set)
for p in paths:
    g = p.get("gene_name")
    if not g:
        continue
    proposed[g].add((p.get("drug_name") or "").lower())
    d = _GOAL_DIRECTION.get(p.get("direction") or "")
    if d:
        goal[g].add(d)
genes = sorted(proposed)
for g in genes:
    if not goal[g]:
        goal[g] = {"Downregulated", "Upregulated"}
directed = sum(1 for g in genes if len(goal[g]) == 1)
print(f"{len(paths)} paths -> {len(genes)} genes "
      f"({directed} with a Translator direction, {len(genes)-directed} run both ways)\n",
      flush=True)

# Compounds Translator already proposed for this disease are not "alternates".
all_proposed = {d for s in proposed.values() for d in s if d}

m = GXAMatcher()
rows, t0 = [], time.time()
for i, g in enumerate(sorted(genes), 1):
    for direction in sorted(goal[g]):
        try:
            found = m.alternate_compounds(g, direction, max_records=400)
        except Exception as exc:
            print(f"  [{i}/{len(genes)}] {g:10s} {direction:14s} ERROR {exc}", flush=True)
            continue
        if not found:
            continue
        novel = [a for a in found if a.compound.lower() not in all_proposed]
        print(f"  [{i}/{len(genes)}] {g:10s} {direction:14s} "
              f"{len(found):4d} contrasts, {len({a.compound for a in novel}):3d} novel compounds",
              flush=True)
        for a in novel:
            rows.append({
                "gene": g, "goal_direction": direction,
                "translator_direction_known": len(goal[g]) == 1,
                "compound": a.compound, "dose": a.dose,
                "log2fc": a.log2fc, "adj_p": a.adj_p,
                "experiment": a.experiment, "comparison": a.comparison,
                "test_arm": a.test_arm, "reference_arm": a.reference_arm,
                "n_translator_drugs_for_gene": len([d for d in proposed[g] if d]),
            })

# Ranking by raw gene count selects for promiscuity, not relevance: valproic
# acid and trichostatin A are HDAC inhibitors that move thousands of genes, and
# doxycycline is the Tet-on inducer. Normalise by each compound's own footprint
# -- how many distinct genes it moves in GXA at all -- so the statistic becomes
# "what fraction of this compound's activity lands on the disease's genes".
VEH = lucene_or("measurementDenominator.value", list(_VEHICLE_TERMS))
NOT_VEH = lucene_or("variableMeasured.value", list(_VEHICLE_TERMS))
FACET_CAP = 1000  # NDE's maximum; saturation means "promiscuous", not exactly 1000


def footprint(compound: str) -> tuple[int, int]:
    """(distinct genes moved, contrasts) for a compound across all of GXA."""
    q = (f'@type:Inference AND species.identifier:9606 '
         f'AND variableMeasured.value:"{compound}" AND {VEH} AND NOT {NOT_VEH}')
    n = m.client.count(q)
    if not n:
        return 0, 0
    return len(m.client.facet(q, "observationAbout.name", facet_size=FACET_CAP)), n


by_compound = defaultdict(list)
for r in rows:
    by_compound[r["compound"]].append(r)

n_genes_queried = len(genes)
ranked = []
for c, v in by_compound.items():
    hit_genes = sorted({r["gene"] for r in v})
    n_all_genes, n_contrasts = footprint(c)
    saturated = n_all_genes >= FACET_CAP
    ranked.append({
        "compound": c, "kind": perturbagen_kind(c),
        "n_genes": len(hit_genes), "genes": hit_genes, "n_contrasts": len(v),
        "gxa_genes_moved": n_all_genes, "gxa_genes_saturated": saturated,
        "gxa_contrasts": n_contrasts,
        "specificity": round(len(hit_genes) / n_all_genes, 4) if n_all_genes else None,
        "best_log2fc": max((abs(r["log2fc"] or 0) for r in v), default=0),
        "experiments": sorted({r["experiment"] for r in v if r["experiment"]})[:4],
        "directions_known": any(r["translator_direction_known"] for r in v),
    })
# Specificity, then gene count. A compound that moves 6 of the disease's genes
# out of 40 genes it touches beats one that moves 51 out of 1,000+.
ranked.sort(key=lambda x: (-(x["specificity"] or 0), -x["n_genes"]))

print(f"\n--- {len(rows)} gene x compound contrasts in {time.time()-t0:.0f}s ---")
print(f"distinct novel compounds: {len(ranked)}   "
      f"({m.n_rejected:,} contrasts rejected as non-compound arms)")
drugs_only = [r for r in ranked if r["kind"] == "compound"]
kinds = Counter(r["kind"] for r in ranked)
print(f"perturbagen kinds: {dict(kinds)}")
print(f"\ntop alternates -- small-molecule compounds, ranked by specificity")
print(f"{'compound':28s}{'genes':>6s}{'GXA genes':>11s}{'spec':>8s}{'|log2FC|':>10s}")
for r in drugs_only[:15]:
    sat = "+" if r["gxa_genes_saturated"] else " "
    print(f"  {r['compound'][:26]:26s}{r['n_genes']:>5d}"
          f"{r['gxa_genes_moved']:>10,}{sat}{(r['specificity'] or 0):>8.3f}"
          f"{r['best_log2fc']:>10.1f}  {','.join(r['genes'][:4])}")
print("\n  (+ = hit NDE's 1,000-bucket facet cap, i.e. promiscuous; specificity is a ceiling)")

print(f"\nfor contrast, the raw gene-count ranking that specificity replaces:")
for r in sorted(drugs_only, key=lambda x: -x["n_genes"])[:6]:
    print(f"  {r['compound'][:26]:26s}{r['n_genes']:>5d} genes  "
          f"but moves {r['gxa_genes_moved']:,}{'+' if r['gxa_genes_saturated'] else ''} "
          f"genes overall -> specificity {(r['specificity'] or 0):.3f}")

out_file.write_text(json.dumps(
    {"source": str(paths_file), "disease": doc.get("disease"),
     "n_genes": len(genes), "n_directed": directed,
     "compounds": ranked, "contrasts": rows}, indent=2))
print(f"\nwrote {out_file}")
