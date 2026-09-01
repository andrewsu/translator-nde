import json, sys, time
sys.path.insert(0, "src")
from pathlib import Path
from translator_nde.translator import ARSClient, load_message, extract_paths, dedupe_paths

DISEASE = sys.argv[1] if len(sys.argv) > 1 else "MONDO:0004979"  # asthma
NAME = sys.argv[2] if len(sys.argv) > 2 else "asthma"

c = ARSClient()
pk = c.submit(DISEASE)
print(f"submitted {NAME} ({DISEASE}) -> pk={pk}", flush=True)
print(f"UI: https://ui.transltr.io/main/results?l={NAME}&i={DISEASE}&t=0&r=0&q={pk}", flush=True)
c.wait(pk, poll=20, max_wait=900)

children = c.fetch_all_children(pk)
print(f"\ndownloaded {len(children)} ARA payloads", flush=True)
all_paths = []
for agent, path in children.items():
    mb = path.stat().st_size / 1e6
    try:
        msg = load_message(path)
        paths = dedupe_paths(extract_paths(msg, DISEASE))
    except Exception as e:
        print(f"  {agent:32s} {mb:7.1f} MB  ERROR {e}", flush=True)
        continue
    n_res = len(msg.get("results") or [])
    n_aux = len(msg.get("auxiliary_graphs") or {})
    print(f"  {agent:32s} {mb:7.1f} MB  results={n_res:5d} aux={n_aux:6d} 2hop_paths={len(paths)}", flush=True)
    for p in paths:
        d = p.to_dict(); d["agent"] = agent; all_paths.append(d)

out = Path("data/ars") / pk / "paths.json"
out.write_text(json.dumps({"disease": DISEASE, "pk": pk, "paths": all_paths}, indent=2))
print(f"\ntotal drug->gene->disease paths: {len(all_paths)}  -> {out}", flush=True)
for p in all_paths[:15]:
    print(f"  {p['drug_name']} --{p['drug_gene_predicate'].replace('biolink:','')}"
          f"[{p['direction']}/{p['aspect']}]--> {p['gene_name']} --> {p['disease_name']}", flush=True)
