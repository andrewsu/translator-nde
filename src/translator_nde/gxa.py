"""Route A: test a Translator drug->gene edge against GXA Inference records.

Each NDE ``@type:Inference`` record is one differential-expression contrast:
a gene, a log2 fold change, an adjusted p-value, and Biolink-typed direction and
aspect qualifiers. That lines up directly with a Translator drug->gene edge:

    Translator  object_aspect_qualifier: abundance   <->  tripleSubjectQualifier aspect
    Translator  object_direction_qualifier: increased <->  tripleSubjectQualifier direction

so the question we can answer is not merely "does data exist?" but "does the
measured direction agree with the asserted one?".

STAGING ONLY -- these records do not exist in production.
"""

from __future__ import annotations

import re
import statistics
from dataclasses import asdict, dataclass, field
from typing import Any, Literal

from .nde import STAGING, NDEClient, lucene_or

HUMAN_TAXON = "9606"

# A dose-like remainder: "4 millimolar", "0.6 microgram per milliliter", "5 uM".
_DOSE = re.compile(r"^[\d.]+\s*(e-?\d+\s*)?[a-z%/ ]*$", re.I)

# Terms that are legitimate chemical names or abbreviations but overwhelmingly
# occur in GXA as ordinary English or as reference-arm labels.
_STOPWORDS = frozenset({
    "no", "none", "control", "untreated", "vehicle", "na", "normal",
    "wild type", "mock", "as", "in", "it", "at", "dmso", "placebo",
})

# Below this length a synonym collides with unrelated text far more often than
# it identifies a compound ("RFP", "NO", "CA", "RIF").
_MIN_SYNONYM_LEN = 3


def factor_supports_drug(variable_measured: str | None, synonyms: list[str]) -> str | None:
    """Is the drug the experimental *variable*, or does its name merely occur?

    Elasticsearch matches a synonym anywhere in the test-arm text, which admits
    contrasts where the term is incidental. Four measured examples, all of which
    Route A originally counted as evidence:

    * ``MITF-RFP-HA overexpression`` -- "RFP" is red fluorescent protein, not
      rifampicin;
    * ``no response to infliximab treatment`` -- "NO" is the English word, not
      nitric oxide;
    * ``A/CA/04/2009 Influenza virus`` -- "CA" is California, not calcium;
    * ``before first infliximab treatment, ..., Crohn's disease`` -- infliximab
      describes the patient group; the variable is disease.

    ``variableMeasured.value`` is a comma-separated list of *factor values*, and
    a compound factor is the compound name alone or the name followed by a dose.
    So require the synonym to be a whole factor, or to lead one with nothing but
    a dose after it. Returns the matching factor, or None.
    """
    for factor in re.split(r"[,;]", variable_measured or ""):
        f = factor.strip()
        low = f.lower()
        for syn in synonyms:
            t = (syn or "").strip().lower()
            if len(t) < _MIN_SYNONYM_LEN or t in _STOPWORDS:
                continue
            if low == t:
                return f
            if low.startswith(t + " ") and _DOSE.match(low[len(t):].strip()):
                return f
    return None

# Biolink direction qualifier values <-> GXA's tripleSubjectQualifier direction.
_DIRECTION_MAP = {"increase": "increased", "decrease": "decreased"}

Verdict = Literal[
    "agrees",                 # gene is DE under the drug, direction matches
    "disagrees",              # gene is DE under the drug, direction opposes
    "ambiguous",              # gene is DE but direction is mixed or unasserted
    "tested_not_significant", # drug HAS contrasts; gene never significantly DE
    "no_drug_data",           # drug absent from GXA -- no information either way
]


@dataclass
class Contrast:
    """One GXA DE contrast, flattened to the fields we reason over."""

    record_id: str
    gene_symbol: str | None
    gene_ensembl: str | None
    log2fc: float | None
    adj_p: float | None
    direction: str | None          # "increase" / "decrease"
    aspect: str | None             # "abundance"
    comparison: str | None         # human-readable contrast string
    test_arm: str | None           # variableMeasured.value, the factor list
    matched_factor: str | None     # the factor the drug name matched, if any
    experiment: str | None         # e.g. E-MTAB-7745
    gsm_ids: list[str] = field(default_factory=list)
    url: str | None = None

    @classmethod
    def from_hit(cls, hit: dict[str, Any]) -> "Contrast":
        about = hit.get("observationAbout") or {}
        if isinstance(about, list):
            about = about[0] if about else {}
        sem = hit.get("semanticMapping") or {}
        quals = {
            q.get("name"): q.get("value")
            for q in _as_list(sem.get("tripleSubjectQualifier"))
            if isinstance(q, dict)
        }
        var = hit.get("variableMeasured") or {}
        if isinstance(var, list):
            var = var[0] if var else {}
        subject_of = hit.get("subjectOf") or {}
        if isinstance(subject_of, list):
            subject_of = subject_of[0] if subject_of else {}
        ids = _as_list(subject_of.get("identifier"))
        return cls(
            record_id=hit.get("_id", ""),
            gene_symbol=about.get("name"),
            gene_ensembl=about.get("identifier"),
            log2fc=_num(hit.get("value")),
            adj_p=_num((hit.get("marginOfError") or {}).get("value")),
            direction=quals.get("direction"),
            aspect=quals.get("aspect"),
            comparison=hit.get("measurementQualifier"),
            test_arm=var.get("value"),
            matched_factor=None,
            experiment=next((i for i in ids if not i.startswith("GSM")), None),
            gsm_ids=[i for i in ids if i.startswith("GSM")],
            url=subject_of.get("url"),
        )


@dataclass
class EdgeEvidence:
    """Aggregated GXA evidence for one Translator drug->gene edge."""

    drug: str
    gene: str
    asserted_direction: str | None
    n_contrasts: int
    n_drug_contrasts: int      # contrasts for the drug across ALL genes
    n_agree: int
    n_disagree: int
    verdict: Verdict
    median_log2fc: float | None
    min_adj_p: float | None
    experiments: list[str]
    contrasts: list[Contrast] = field(default_factory=list)

    def to_dict(self, *, include_contrasts: bool = False) -> dict:
        d = asdict(self)
        if not include_contrasts:
            d.pop("contrasts")
        return d


class GXAMatcher:
    """Queries NDE staging for DE contrasts supporting a drug->gene edge."""

    def __init__(self, client: NDEClient | None = None, *, taxon: str = HUMAN_TAXON):
        self.client = client or NDEClient(base_url=STAGING)
        if self.client.base_url != STAGING:
            raise ValueError("GXA Inference records exist only on NDE staging")
        self.taxon = taxon
        # Contrasts the Lucene query returned but whose test arm does not put
        # the drug in the variable position. Reported, not silently discarded.
        self.n_rejected = 0

    def build_query(self, drug_synonyms: list[str], gene_terms: list[str]) -> str:
        """Compose the query, including the two filters that make it correct.

        The drug must appear in the **test** arm and be absent from the
        **reference** arm. Without the exclusion, contrasts where the compound is
        present in both arms -- i.e. where it is an experimental tool rather than
        the variable -- are wrongly counted as drug->gene evidence.
        """
        test_arm = lucene_or("variableMeasured.value", drug_synonyms)
        ref_arm = lucene_or("measurementDenominator.value", drug_synonyms)
        gene = lucene_or("observationAbout.name", gene_terms)
        return (
            f"@type:Inference AND species.identifier:{self.taxon} "
            f"AND {test_arm} AND NOT {ref_arm} AND {gene}"
        )

    def drug_contrast_count(self, drug_synonyms: list[str], *, sample: int = 200) -> int:
        """How many contrasts test this drug at all, across every gene.

        GXA stores only *significant* DE results, so a gene having no record is
        not the same as the drug being absent. This separates the two.
        """
        test_arm = lucene_or("variableMeasured.value", drug_synonyms)
        ref_arm = lucene_or("measurementDenominator.value", drug_synonyms)
        q = (f"@type:Inference AND species.identifier:{self.taxon} "
             f"AND {test_arm} AND NOT {ref_arm}")
        if self.client.count(q) == 0:
            return 0
        # Text-verify a capped sample rather than trusting the raw count, which
        # is inflated by the same incidental matches factor_supports_drug drops.
        return sum(
            1
            for h in self.client.scroll(
                q, fields="variableMeasured.value", max_records=sample
            )
            if factor_supports_drug(
                (lambda v: (v[0] if isinstance(v, list) else v) or {})(
                    h.get("variableMeasured") or {}
                ).get("value"),
                drug_synonyms,
            )
        )

    def contrasts(
        self, drug_synonyms: list[str], gene_terms: list[str], *, max_records: int = 500
    ) -> list[Contrast]:
        q = self.build_query(drug_synonyms, gene_terms)
        fields = ",".join(
            [
                "observationAbout",
                "value",
                "unitText",
                "marginOfError.value",
                "measurementQualifier",
                "variableMeasured.value",
                "semanticMapping.tripleSubjectQualifier",
                "subjectOf.identifier",
                "subjectOf.url",
            ]
        )
        if self.client.count(q) == 0:
            return []
        kept = []
        for h in self.client.scroll(q, fields=fields, max_records=max_records):
            c = Contrast.from_hit(h)
            factor = factor_supports_drug(c.test_arm, drug_synonyms)
            if factor is None:
                self.n_rejected += 1
                continue
            c.matched_factor = factor
            kept.append(c)
        return kept

    def evaluate(
        self,
        *,
        drug_label: str,
        gene_label: str,
        drug_synonyms: list[str],
        gene_terms: list[str],
        asserted_direction: str | None,
        max_records: int = 500,
    ) -> EdgeEvidence:
        """Compare measured direction against the direction Translator asserts.

        ``asserted_direction`` is the Biolink ``object_direction_qualifier``
        ("increased"/"decreased"), or None when the edge carries no direction --
        in which case we can report coverage but not agreement.
        """
        found = self.contrasts(drug_synonyms, gene_terms, max_records=max_records)
        n_drug = len(found) and -1  # only pay for the extra query when needed
        if not found:
            n_drug = self.drug_contrast_count(drug_synonyms)
        agree = disagree = 0
        for c in found:
            mapped = _DIRECTION_MAP.get((c.direction or "").lower())
            if asserted_direction and mapped:
                if mapped == asserted_direction:
                    agree += 1
                else:
                    disagree += 1

        fcs = [c.log2fc for c in found if c.log2fc is not None]
        ps = [c.adj_p for c in found if c.adj_p is not None]
        return EdgeEvidence(
            drug=drug_label,
            gene=gene_label,
            asserted_direction=asserted_direction,
            n_contrasts=len(found),
            n_drug_contrasts=n_drug if n_drug >= 0 else len(found),
            n_agree=agree,
            n_disagree=disagree,
            verdict=_verdict(found, asserted_direction, agree, disagree, n_drug),
            median_log2fc=statistics.median(fcs) if fcs else None,
            min_adj_p=min(ps) if ps else None,
            experiments=sorted({c.experiment for c in found if c.experiment}),
            contrasts=found,
        )


def _verdict(
    found: list[Contrast], asserted: str | None, agree: int, disagree: int,
    n_drug_contrasts: int = 0,
) -> Verdict:
    if not found:
        # The drug being present in GXA but the gene never turning up means the
        # gene was measured genome-wide and was not significantly DE -- that is
        # evidence against the edge, not absence of evidence.
        return "tested_not_significant" if n_drug_contrasts > 0 else "no_drug_data"
    if not asserted or (agree == 0 and disagree == 0):
        return "ambiguous"  # covered by data, but no directional claim to test
    if agree and not disagree:
        return "agrees"
    if disagree and not agree:
        return "disagrees"
    # Mixed across tissues/doses is the common real case, so report the majority
    # only when it is decisive.
    total = agree + disagree
    if agree / total >= 0.8:
        return "agrees"
    if disagree / total >= 0.8:
        return "disagrees"
    return "ambiguous"


def _as_list(v: Any) -> list:
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


def _num(v: Any) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
