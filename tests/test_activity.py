"""Regression fixtures for Route D's compound-target matching.

Two live-data cases each broke the first implementation, and each represents a
whole class of ChEMBL records, so both are pinned:

* imatinib files *no* mechanism under its own ChEMBL id -- all four are on the
  mesylate salt, reachable only via ``parent_molecule_chembl_id``;
* aspirin's mechanism target is the "Cyclooxygenase" PROTEIN FAMILY, not the
  PTGS2 single protein, so target-id equality alone misses it.

Both silently downgraded a curated, action-typed mechanism to a mere binding
observation -- the exact distinction Route D exists to make.
"""

import pytest

from translator_nde.activity import (
    ActivityMatcher, Mechanism, _match_mechanism, _DECREASING, _INCREASING,
)

IMATINIB = "CHEBI:45783"
ASPIRIN = "CHEBI:15365"
BARICITINIB = "CHEBI:95341"
ABL1, PTGS2, JAK1 = "NCBIGene:25", "NCBIGene:5743", "NCBIGene:3716"


@pytest.fixture(scope="module")
def matcher():
    return ActivityMatcher()


def test_salt_form_mechanism_is_found(matcher):
    """Imatinib -> ABL1 must score as a curated mechanism, not merely binding."""
    ev = matcher.evaluate(
        drug=IMATINIB, drug_name="Imatinib", gene=ABL1, gene_name="ABL1",
        asserted_direction="decreased", asserted_aspect="activity",
    )
    assert ev.chembl_action_type == "INHIBITOR"
    assert ev.verdict == "mechanism_agrees"


def test_protein_family_target_is_matched(matcher):
    """Aspirin -> PTGS2 is filed against the Cyclooxygenase family."""
    ev = matcher.evaluate(
        drug=ASPIRIN, drug_name="Aspirin", gene=PTGS2, gene_name="PTGS2",
        asserted_direction="decreased", asserted_aspect="activity",
    )
    assert ev.chembl_action_type == "INHIBITOR"
    assert ev.verdict == "mechanism_agrees"


def test_edge_route_a_could_not_score(matcher):
    """Baricitinib -> JAK1: GXA has the drug but the JAKs are never DE.

    Route A can say nothing here; Route D confirms it outright. This is the
    worked example that motivated the whole route.
    """
    ev = matcher.evaluate(
        drug=BARICITINIB, drug_name="Baricitinib", gene=JAK1, gene_name="JAK1",
        asserted_direction="decreased", asserted_aspect="activity",
    )
    assert ev.verdict == "mechanism_agrees"
    assert ev.max_phase == 4.0
    assert ev.pchembl_max and ev.pchembl_max > 8


def test_unrelated_pair_is_not_tested(matcher):
    """A drug never assayed against a target must not be scored as evidence."""
    ev = matcher.evaluate(
        drug=ASPIRIN, drug_name="Aspirin", gene=JAK1, gene_name="JAK1",
        asserted_direction="decreased", asserted_aspect="activity",
    )
    assert ev.verdict == "not_tested"
    assert ev.n_active == 0


def test_direction_vocabularies_are_disjoint():
    """An action type must never imply both directions."""
    assert not (_DECREASING & _INCREASING)


def test_directional_mechanism_wins_over_untyped():
    """When a compound has both a typed and an untyped mechanism on the same
    target, the typed one carries the evidence and must be preferred."""
    untyped = Mechanism("M", "T1", "BINDING AGENT", "binds", True, 4.0)
    typed = Mechanism("M", "T1", "INHIBITOR", "inhibits", True, 4.0)
    assert _match_mechanism([untyped, typed], "T1", None) is typed


def test_accession_matching_when_target_id_differs():
    """Family/complex targets are matched through their component accessions."""
    fam = Mechanism("M", "CHEMBL2094253", "INHIBITOR", "COX inhibitor", True, 4.0,
                    target_accessions=["P35354", "P23219"])
    assert _match_mechanism([fam], "CHEMBL230", "P35354") is fam
    assert _match_mechanism([fam], "CHEMBL230", "Q99999") is None


def test_heavily_screened_target_is_reachable(matcher):
    """DRD2 is the case that forced the compound-keyed PubChem endpoint.

    Its gene-keyed ``concise`` dump is a 436 MB body that arrives truncated and
    fails to parse, so the gene-first client returned nothing and the pair was
    scored ``not_tested`` -- a fetch failure laundered into a measured negative.
    The compound-keyed ``assaysummary`` view returns the same rows in ~300 KB.
    """
    ev = matcher.evaluate(
        drug="CHEBI:10454", drug_name="Flupentixol",
        gene="NCBIGene:1813", gene_name="DRD2",
    )
    assert not ev.pubchem_error
    assert ev.n_active > 0
    assert ev.verdict == "binding_confirmed"


def test_fetch_failure_is_not_scored_as_a_negative():
    """A PubChem error must never be reported as evidence against an edge."""
    from translator_nde.activity import ActivityEvidence, _verdict

    ev = ActivityEvidence(
        drug="X", drug_name="X", gene="NCBIGene:1", gene_name="G",
        asserted_direction="decreased", asserted_aspect="activity",
        verdict="not_tested", drug_cid="1", pubchem_error=True,
    )
    assert _verdict(ev, None, has_target_data=True) == "fetch_failed"
