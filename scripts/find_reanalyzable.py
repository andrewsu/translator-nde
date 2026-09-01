"""Route B discovery: RA/AS GEO series that already ship a gene-level matrix.

Pulls candidate series from NDE (disease via MONDO, human, sequencing assay),
then asks GEO which of them deposited a counts or expression matrix -- those are
re-analyzable without realigning raw reads.

Also records how many NDE `@type:Sample` records each series has, since those
carry the per-GSM arm labels the DE step needs for grouping.
"""
import json, sys, time
from collections import Counter
from pathlib import Path

import requests

sys.path.insert(0, "src")
from translator_nde.nde import NDEClient, PROD
from translator_nde.geo import inspect

# MONDO id -> (label, free-text synonyms). Only 41.9% of NDE's GEO datasets carry
# any healthCondition, so MONDO-only discovery misses about half the relevant
# series; the union with free text roughly doubles recall (RA 213 -> 422).
DISEASES = {
    "0008383": ("rheumatoid arthritis", ['"rheumatoid arthritis"']),
    "0005306": ("ankylosing spondylitis",
                ['"ankylosing spondylitis"', '"axial spondyloarthritis"']),
}
# Sequencing assays only -- arrays don't produce counts.
SEQ = ('("assay by high throughput sequencer" OR '
       '"high-throughput expression assay" OR "rna-seq of coding rna")')

c = NDEClient(base_url=PROD, pause=0.05)
session = requests.Session()

candidates: dict[str, dict] = {}
for mondo, (label, synonyms) in DISEASES.items():
    disease_clause = " OR ".join([f"healthCondition.identifier:{mondo}"] + synonyms)
    q = (f'@type:Dataset AND ({disease_clause}) '
         f'AND includedInDataCatalog.name:"NCBI GEO" '
         f'AND species.name:"homo sapiens" '
         f'AND measurementTechnique.name:{SEQ}')
    n = c.count(q)
    print(f"{label}: {n} human GEO sequencing datasets", flush=True)
    for hit in c.scroll(q, fields="identifier,name,description,measurementTechnique.name,url",
                        max_records=n):
        ident = hit.get("identifier")
        ident = ident[0] if isinstance(ident, list) else ident
        gse = ident if isinstance(ident, str) and ident.startswith("GSE") else None
        if not gse:
            continue
        e = candidates.setdefault(gse, {"gse": gse, "name": hit.get("name"),
                                        "diseases": [], "url": hit.get("url")})
        if label not in e["diseases"]:
            e["diseases"].append(label)

print(f"\n{len(candidates)} distinct GSE candidates; checking GEO for matrices\n", flush=True)

rows, t0 = [], time.time()
for i, (gse, e) in enumerate(sorted(candidates.items()), 1):
    sup = inspect(gse, session=session)
    n_samples = c.count(f'@type:Sample AND isBasisFor.identifier:"{gse}"')
    rows.append({**e, "matrix_kind": sup.matrix_kind, "matrix_files": sup.matrix_files,
                 "n_suppl_files": len(sup.files), "nde_samples": n_samples,
                 "error": sup.error})
    if sup.reanalyzable:
        print(f"  [{i}/{len(candidates)}] {gse:12s} {sup.matrix_kind:12s} "
              f"samples={n_samples:4d}  {(e['name'] or '')[:52]}", flush=True)
    elif i % 25 == 0:
        print(f"  ... {i}/{len(candidates)} checked, "
              f"{sum(r['matrix_kind'] in ('raw_counts','normalized') for r in rows)} usable, "
              f"{time.time()-t0:.0f}s", flush=True)

out = Path("results/reanalyzable_ra_as.json")
out.write_text(json.dumps({"diseases": {k: v[0] for k, v in DISEASES.items()}, "n_candidates": len(candidates),
                           "datasets": rows}, indent=2))

kinds = Counter(r["matrix_kind"] for r in rows)
usable = [r for r in rows if r["matrix_kind"] in ("raw_counts", "normalized")]
print(f"\n--- {len(rows)} series checked in {time.time()-t0:.0f}s ---")
for k, n in kinds.most_common():
    print(f"  {k:18s} {n}")
print(f"\nre-analyzable: {len(usable)}/{len(rows)} "
      f"({100*len(usable)/max(len(rows),1):.0f}%)")
print(f"  with NDE sample-level arm labels: "
      f"{sum(1 for r in usable if r['nde_samples'] > 0)}")
print(f"wrote {out}")
