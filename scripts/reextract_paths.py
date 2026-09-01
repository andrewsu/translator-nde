"""Re-extract paths from cached ARA payloads, without re-querying the ARS.

The raw per-ARA payloads are on disk under data/ars/, so extractor changes can
be replayed against them -- important because ARS results are not reproducible
over time.
"""
import json, sys
from pathlib import Path
sys.path.insert(0, "src")
from translator_nde.translator import ARSClient, load_message, extract_paths, dedupe_paths

for pk_dir in sorted(Path("data/ars").iterdir()):
    if not pk_dir.is_dir() or not (pk_dir / "query.json").exists():
        continue
    meta = json.loads((pk_dir / "query.json").read_text())
    disease, pk = meta["disease"], meta["pk"]
    c = ARSClient()
    children = c.fetch_all_children(pk)
    all_paths = []
    for agent, path in children.items():
        try:
            paths = dedupe_paths(extract_paths(load_message(path), disease))
        except Exception as exc:
            print(f"  {agent}: ERROR {exc}"); continue
        for p in paths:
            d = p.to_dict(); d["agent"] = agent; all_paths.append(d)
    (pk_dir / "paths.json").write_text(
        json.dumps({"disease": disease, "pk": pk, "paths": all_paths}, indent=2))
    print(f"{disease:22s} pk={pk[:8]}  {len(all_paths):5d} paths "
          f"({len({p['gene_name'] for p in all_paths})} genes, "
          f"{len({p['drug_name'] for p in all_paths})} drugs)")
