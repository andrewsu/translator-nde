"""Regression fixtures for the drug-arm discrimination rule.

These two records are the difference between a correct and an incorrect
drug->gene evidence count, so they are pinned as tests.
"""

import pytest

from translator_nde.nde import NDEClient, STAGING
from translator_nde.gxa import (
    GXAMatcher, factor_supports_drug, drug_is_the_variable,
)

# Arabidopsis; dexamethasone is a GR *inducer* present in BOTH arms, so the real
# contrast is genotype (bZIP1 vs empty vector), not drug. Must be excluded.
TOOL_COMPOUND_RECORD = "gxa_e_geod_54049_g1_g5_at1g03790"

# Human; 'dexamethasone; 1 micromolar' vs 'vehicle'. A genuine drug contrast.
GENUINE_DRUG_RECORD = "gxa_e_mtab_7745_g2_g1_ensg00000142319"


@pytest.fixture(scope="module")
def client():
    return NDEClient(base_url=STAGING)


def test_tool_compound_record_is_excluded(client):
    """The both-arms record must not survive the matcher's filters."""
    matcher = GXAMatcher(client)
    q = matcher.build_query(["dexamethasone"], ["SOM"])
    assert client.count(f'({q}) AND _id:"{TOOL_COMPOUND_RECORD}"') == 0


def test_tool_compound_record_would_match_naive_query(client):
    """...and it *would* be caught by the naive query, so the rule is load-bearing."""
    naive = f'@type:Inference AND dexamethasone AND _id:"{TOOL_COMPOUND_RECORD}"'
    assert client.count(naive) == 1


def test_genuine_drug_record_is_kept(client):
    matcher = GXAMatcher(client)
    q = matcher.build_query(["dexamethasone"], ["SLC6A3"])
    assert client.count(f'({q}) AND _id:"{GENUINE_DRUG_RECORD}"') == 1


def test_reference_arm_exclusion_removes_records(client):
    """The NOT clause must actually remove a meaningful number of contrasts."""
    base = "@type:Inference AND variableMeasured.value:dexamethasone"
    with_excl = f"{base} AND NOT measurementDenominator.value:dexamethasone"
    assert client.count(base) > client.count(with_excl) > 0


# --- incidental text matches -------------------------------------------------
# Elasticsearch matches a synonym anywhere in the test-arm text. Each of these
# was counted as real drug->gene evidence by Route A before the factor-position
# filter, and each represents a distinct way the match can be incidental.

def test_rejects_synonym_colliding_with_a_protein_tag():
    """'RFP' is red fluorescent protein here, not rifampicin."""
    assert factor_supports_drug(
        "MITF-RFP-HA overexpression", ["Rifampicin", "RIF", "RFP"]
    ) is None


def test_rejects_synonym_colliding_with_english():
    """'NO' is the word, not nitric oxide."""
    assert factor_supports_drug(
        "no response to infliximab treatment, Crohn's disease", ["Nitric Oxide", "NO"]
    ) is None


def test_rejects_synonym_inside_a_strain_name():
    """'CA' is California, not calcium."""
    assert factor_supports_drug(
        "A/CA/04/2009 Influenza virus, 30 hour", ["Calcium", "CA"]
    ) is None


def test_rejects_drug_naming_the_patient_group():
    """The variable is disease; infliximab only describes the cohort, and this
    contrast is taken *before* any treatment was given."""
    assert factor_supports_drug(
        "before first infliximab treatment, no response to infliximab treatment, "
        "Crohn's disease, colon",
        ["Infliximab", "Remicade"],
    ) is None


def test_rejects_drug_named_only_as_a_stimulus_description():
    assert factor_supports_drug(
        "Tr1 cell clone, 6 hour, stimulated with monoclonal antibodies to CD3 and CD28",
        ["Antibodies"],
    ) is None


def test_keeps_bare_compound_factor():
    assert factor_supports_drug("A2780cis, Cisplatin, Normoxia", ["Cisplatin"]) == "Cisplatin"
    assert factor_supports_drug(
        "differentiated brown adiopcyte, cyclic AMP", ["Cyclic AMP"]
    ) == "cyclic AMP"


def test_keeps_compound_with_a_dose():
    assert factor_supports_drug("SK-BR-3, metformin 4 millimolar", ["Metformin"]) \
        == "metformin 4 millimolar"
    assert factor_supports_drug(
        "doxorubicin 0.6 microgram per milliliter", ["Doxorubicin"]
    ) == "doxorubicin 0.6 microgram per milliliter"
    assert factor_supports_drug(
        "5-aza-deoxy-cytidine 5 micromolar", ["Azacitidine", "5-aza-deoxy-cytidine"]
    ) == "5-aza-deoxy-cytidine 5 micromolar"


# --- arm discrimination: the drug must be what *differs* between the arms -----
# `NOT measurementDenominator.value:<drug>` looks like the right rule and is far
# too blunt. These four pin the distinction it gets wrong.

def test_keeps_tggates_cohort_labelled_control():
    """TG-GATEs names the compound on every sample via a `cohort` factor, so the
    aspirin study's vehicle control reads 'aspirin, liver, 15 day' without any
    aspirin having been given. Excluding on the name deletes the whole study."""
    assert drug_is_the_variable(
        "aspirin, aspirin 150 milligram per kilogram, liver, 15 day",
        "aspirin, liver, 15 day",
        ["aspirin"],
    ) == "aspirin 150 milligram per kilogram"


def test_excludes_identical_drug_factor_in_both_arms():
    """The Arabidopsis case: dexamethasone is a GR inducer given to both arms
    and the variable is genotype."""
    assert drug_is_the_variable(
        "pBeaconRFP_GR::bZIP1, Treated with cycloheximide, "
        "10 uM dexamethasone in ethanol",
        "empty vector, Treated with cycloheximide, 10 uM dexamethasone in ethanol",
        ["dexamethasone"],
    ) is None


def test_excludes_drug_held_constant_across_a_disease_contrast():
    """'prednisolone 20 milligram per day, polymyalgia rheumatica' vs the same
    dose in 'normal' -- the variable is disease. Eight real human records."""
    assert drug_is_the_variable(
        "prednisolone 20 milligram per day, polymyalgia rheumatica",
        "prednisolone 20 milligram per day, normal",
        ["prednisolone"],
    ) is None


def test_prefers_the_dosed_factor_over_a_bare_label():
    """When a test arm carries both a bare and a dosed occurrence, the dosed one
    is the evidence that the compound was administered."""
    from translator_nde.gxa import _matching_factors

    assert _matching_factors(
        "aspirin, aspirin 45 milligram per kilogram", ["aspirin"]
    )[0] == "aspirin 45 milligram per kilogram"
