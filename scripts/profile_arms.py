"""Characterize the arm structure of candidate GEO series.

Sample count alone is a bad selector: GSE89408 has 218 samples but contrasts RA
against healthy, which tests no drug->gene edge. What Route B needs is a series
with a *drug-treatment* arm and a matched control arm.
"""
import json, re, sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.nde import NDEClient, PROD
from translator_nde.reanalysis import fetch_samples

MIN_SAMPLES = int(sys.argv[1]) if len(sys.argv) > 1 else 12

CONTROL = re.compile(
    r"\b(vehicle|untreated|unstim|control|ctrl|dmso|mock|baseline|placebo|"
    r"pre[- ]?treatment|day ?0|0 ?h(r|our)?\b|resting|naive)\b", re.I)
# Suffixes shared by most drug INNs, plus explicit treatment language.
TREATMENT = re.compile(
    r"\b(treated|treatment|stimulat|\+ ?drug|"
    r"\w+(mab|nib|cept|statin|mycin|cycline|prazole|oxib|olol|sartan|azole|"
    r"cillin|micin|pine|zumab|ximab|tinib|ciclib|parib)\b)", re.I)
SINGLE_CELL = re.compile(r"single[- ]?cell|scrna|10x|smart-?seq|drop-?seq", re.I)

cat = json.loads(Path("results/reanalyzable_ra_as.json").read_text())
raw = [d for d in cat["datasets"]
       if d["matrix_kind"] == "raw_counts" and d["nde_samples"] >= MIN_SAMPLES]
raw.sort(key=lambda d: -d["nde_samples"])
print(f"{len(raw)} raw-counts series with >={MIN_SAMPLES} samples\n", flush=True)

client = NDEClient(base_url=PROD, pause=0.02)
rows = []
for i, d in enumerate(raw, 1):
    gse = d["gse"]
    if SINGLE_CELL.search((d["name"] or "")):
        print(f"  [{i}/{len(raw)}] {gse:11s} SKIP single-cell", flush=True)
        continue
    try:
        samples = fetch_samples(gse, client)
    except Exception as exc:
        print(f"  [{i}/{len(raw)}] {gse:11s} ERROR {str(exc)[:50]}", flush=True)
        continue
    texts = [s.text() for s in samples]
    n_ctrl = sum(bool(CONTROL.search(t)) for t in texts)
    n_trt = sum(bool(TREATMENT.search(t)) for t in texts)
    labels = Counter(re.sub(r"\d+", "#", (s.name or "").strip().lower())
                     for s in samples)
    rows.append({"gse": gse, "name": d["name"], "n": len(samples),
                 "n_control": n_ctrl, "n_treatment": n_trt,
                 "drug_contrast": n_ctrl >= 2 and n_trt >= 2,
                 "matrix_files": d["matrix_files"], "diseases": d["diseases"],
                 "top_labels": labels.most_common(6)})
    flag = "DRUG" if rows[-1]["drug_contrast"] else "    "
    print(f"  [{i}/{len(raw)}] {gse:11s} n={len(samples):4d} ctrl={n_ctrl:3d} "
          f"trt={n_trt:3d} {flag}  {(d['name'] or '')[:44]}", flush=True)

rows.sort(key=lambda r: (-r["drug_contrast"], -min(r["n_control"], r["n_treatment"]), -r["n"]))
Path("results/arm_profiles.json").write_text(json.dumps(rows, indent=2))
n_drug = sum(r["drug_contrast"] for r in rows)
print(f"\n--- {len(rows)} profiled, {n_drug} with a drug-treatment contrast ---")
print("wrote results/arm_profiles.json")
