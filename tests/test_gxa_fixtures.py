"""Regression fixtures for the drug-arm discrimination rule.

These two records are the difference between a correct and an incorrect
drug->gene evidence count, so they are pinned as tests.
"""

import pytest

from translator_nde.nde import NDEClient, STAGING
from translator_nde.gxa import GXAMatcher

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
