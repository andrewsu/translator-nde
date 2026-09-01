"""Which drugs actually appear in GXA human DE contrasts?

Works backward from the data instead of from Translator: takes DrugMechDB's
drug vocabulary (1,652 drugs, each already paired with the diseases it treats)
and asks GXA how many human test-arm contrasts exist for each.

Drugs with coverage are the ones where Route A can say anything at all, so their
indications are the diseases worth querying Translator with.
"""
import json, sys, time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.nde import NDEClient, STAGING, lucene_or

DMDB = "/home/asu/Science/DrugMechDB/indication_paths.json"

paths = json.load(open(DMDB))
drug_dis = defaultdict(set)     # drug name -> {(disease name, disease mesh)}
drug_mesh = {}
for p in paths:
    g = p["graph"]
    d, dm = g.get("drug"), g.get("drug_mesh")
    dis, dsm = g.get("disease"), g.get("disease_mesh")
    if d and dis:
        drug_dis[d].add((dis, dsm))
        drug_mesh[d] = dm

drugs = sorted(drug_dis)
print(f"testing {len(drugs)} DrugMechDB drugs against GXA human contrasts\n", flush=True)

c = NDEClient(base_url=STAGING, pause=0.02)
BASE = "@type:Inference AND species.identifier:9606"

rows, t0 = [], time.time()
for i, d in enumerate(drugs, 1):
    try:
        test = lucene_or("variableMeasured.value", [d])
        ref = lucene_or("measurementDenominator.value", [d])
        n = c.count(f"{BASE} AND {test} AND NOT {ref}")
        n_ref = c.count(f"{BASE} AND {ref}") if n else 0
    except Exception as exc:
        print(f"  [{i}/{len(drugs)}] {d[:34]:34s} ERROR {str(exc)[:60]}", flush=True)
        continue
    if n:
        rows.append({"drug": d, "drug_mesh": drug_mesh.get(d), "contrasts": n,
                     "reference_arm_contrasts": n_ref,
                     "diseases": sorted(drug_dis[d])})
        print(f"  [{i}/{len(drugs)}] {d[:34]:34s} contrasts={n:6d} "
              f"(ref-arm {n_ref})  indications={len(drug_dis[d])}", flush=True)
    if i % 200 == 0:
        print(f"  ... {i}/{len(drugs)} tested, {len(rows)} with coverage, "
              f"{time.time()-t0:.0f}s", flush=True)

rows.sort(key=lambda r: -r["contrasts"])
out = Path("results/gxa_drug_coverage.json")
out.write_text(json.dumps({"n_tested": len(drugs), "n_covered": len(rows),
                           "drugs": rows}, indent=2))
print(f"\n--- {len(rows)}/{len(drugs)} DrugMechDB drugs have >=1 human GXA contrast "
      f"({100*len(rows)/len(drugs):.1f}%) in {time.time()-t0:.0f}s ---")
print(f"wrote {out}")
