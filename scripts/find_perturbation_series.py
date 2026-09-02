"""Find GEO series that perturb with a Translator-proposed drug and measure expression.

The question is dataset *discovery*, not edge adjudication: for each drug
Translator proposes for a disease, does a GEO experiment exist that treats with
that drug and profiles expression against a control arm?

Two stages, because the cheap signal is not trustworthy on its own:

1. **NDE sample-level search.** A drug named on individual `@type:Sample`
   records is a better hint than one named on the `Dataset`, which may only be
   an abstract mention. Group the hits by parent series.
2. **GEO characteristics confirmation.** Sample-level mentions are still
   ambiguous -- a study-wide descriptor ("asthma patients on ICS") mentions the
   drug on every sample without perturbing anything. The authoritative signal is
   a `!Sample_characteristics_ch1` row whose values name the drug *and* take at
   least one other value, i.e. a real treatment arm with a control.
"""
import collections, gzip, io, json, re, sys, time, urllib.request
from pathlib import Path

sys.path.insert(0, "src")
from translator_nde.nde import NDEClient, PROD

MIN_SAMPLES = 3          # per series, to call it an arm rather than a stray mention
MAX_SERIES_PER_DRUG = 6  # cap the GEO confirmation fetches
SCROLL_CAP = 1200


def drug_like(name: str) -> bool:
    """Skip IUPAC strings and synthetic peptides; they never match free text."""
    return bool(name) and len(name) <= 30 and name.count("-") <= 3 and "(" not in name


def series_matrix_header(gse: str, timeout: int = 60) -> str | None:
    pre = re.sub(r"\d{1,3}$", "nnn", gse)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            raw = r.read()
    except Exception:
        return None
    try:
        return gzip.decompress(raw).decode("utf-8", "replace")
    except Exception:
        return None


def confirm_arm(header: str, drug: str) -> dict | None:
    """A characteristics row naming the drug and holding >1 distinct value."""
    for line in header.splitlines():
        if not line.startswith("!Sample_characteristics_ch1"):
            continue
        vals = [v.strip('"') for v in line.split("\t")[1:]]
        counts = collections.Counter(v for v in vals if v)
        if len(counts) < 2:
            continue
        if any(drug.lower() in v.lower() for v in counts):
            key = vals[0].split(":")[0] if ":" in vals[0] else "characteristic"
            return {"field": key, "values": dict(counts.most_common(6)),
                    "n_levels": len(counts)}
    return None


def main() -> None:
    disease_sets = {
        "asthma": "5b656c0f-b7da-4db4-ba1f-d3a794b422d4",
        "rheumatoid arthritis": "d1f9daae-9330-4c2c-81f5-5c21ee2f44f4",
        "ankylosing spondylitis": "c3efea15-1a8f-4c09-a521-7462603bec08",
    }
    drugs: dict[str, set[str]] = {}
    for disease, pk in disease_sets.items():
        for p in json.loads(Path(f"data/ars/{pk}/paths.json").read_text())["paths"]:
            n = p.get("drug_name")
            if drug_like(n):
                drugs.setdefault(n, set()).add(disease)
    print(f"{len(drugs)} drug-like Translator drugs across {len(disease_sets)} diseases\n")

    c = NDEClient(base_url=PROD)
    t0 = time.time()
    ranked = []
    for i, drug in enumerate(sorted(drugs), 1):
        q = f'@type:Sample AND "{drug}" AND species.name:"Homo sapiens"'
        try:
            n = c.count(q)
        except Exception as exc:
            print(f"  [{i}/{len(drugs)}] {drug[:28]:28s} count failed: {exc}", flush=True)
            continue
        if n < MIN_SAMPLES:
            continue
        by = collections.Counter()
        for h in c.scroll(q, fields="isBasisFor.identifier", max_records=SCROLL_CAP):
            b = h.get("isBasisFor") or {}
            if isinstance(b, list):
                b = b[0] if b else {}
            if b.get("identifier", "").startswith("GSE"):
                by[b["identifier"]] += 1
        cand = [(g, k) for g, k in by.most_common() if k >= MIN_SAMPLES]
        if cand:
            ranked.append({"drug": drug, "diseases": sorted(drugs[drug]),
                           "n_samples": n, "n_series": len(cand),
                           "series": cand[:MAX_SERIES_PER_DRUG]})
            print(f"  [{i}/{len(drugs)}] {drug[:28]:28s} {n:6,} samples, "
                  f"{len(cand):3d} candidate series", flush=True)

    ranked.sort(key=lambda r: -r["n_series"])
    print(f"\n--- stage 1: {len(ranked)} drugs with >=1 candidate series "
          f"in {time.time()-t0:.0f}s ---")

    print("\n--- stage 2: confirming treatment arms from GEO characteristics ---")
    confirmed = []
    for r in ranked:
        for gse, k in r["series"]:
            hdr = series_matrix_header(gse)
            if not hdr:
                continue
            arm = confirm_arm(hdr, r["drug"])
            if arm:
                rec = {"drug": r["drug"], "diseases": r["diseases"], "gse": gse,
                       "n_samples_mentioning": k, **arm}
                confirmed.append(rec)
                print(f"  ✓ {r['drug'][:22]:22s} {gse:12s} {arm['field'][:24]:24s} "
                      f"{list(arm['values'])[:3]}", flush=True)

    out = Path("results/perturbation_series.json")
    out.write_text(json.dumps({"stage1": ranked, "confirmed": confirmed}, indent=2))
    print(f"\nconfirmed drug x series pairs: {len(confirmed)}")
    print(f"distinct drugs with a confirmed arm: "
          f"{len({r['drug'] for r in confirmed})}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
