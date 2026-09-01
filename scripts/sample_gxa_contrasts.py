"""Sample human GXA Inference records to enumerate the distinct contrasts.

GXA Inference holds only *significant* DE results, so there is no gene that
appears in every contrast; the contrast vocabulary has to be sampled.
"""
import json, sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, "src")
from translator_nde.nde import NDEClient, STAGING

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
c = NDEClient(base_url=STAGING)
q = "@type:Inference AND species.identifier:9606"

contrasts, experiments, test_arms = Counter(), Counter(), Counter()
seen = 0
for hit in c.scroll(q, fields="measurementQualifier,subjectOf.identifier,variableMeasured.value",
                    max_records=N):
    seen += 1
    mq = hit.get("measurementQualifier")
    if mq:
        contrasts[mq] += 1
        # test arm is everything before the ' vs ' separator
        test_arms[mq.split("' vs '")[0].strip("'")] += 1
    so = hit.get("subjectOf") or {}
    ids = so.get("identifier") or []
    if isinstance(ids, str): ids = [ids]
    for i in ids:
        if not i.startswith("GSM"):
            experiments[i] += 1
    if seen % 5000 == 0:
        print(f"  {seen} records -> {len(contrasts)} distinct contrasts, "
              f"{len(experiments)} experiments", flush=True)

out = Path("results/gxa_contrast_sample.json")
out.write_text(json.dumps({
    "sampled_records": seen,
    "distinct_contrasts": len(contrasts),
    "distinct_experiments": len(experiments),
    "test_arms": test_arms.most_common(400),
    "contrasts": contrasts.most_common(400),
}, indent=2))
print(f"\nsampled {seen} records")
print(f"distinct contrasts: {len(contrasts)}   distinct experiments: {len(experiments)}")
print(f"wrote {out}")
