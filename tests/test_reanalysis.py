"""Route B unit tests -- arm assignment and column matching, no network."""

import numpy as np
import pandas as pd
import pytest

from translator_nde.reanalysis import (
    Sample, arms_from_columns, assign_arms, match_columns, run_de,
)


def s(gsm, name, **props):
    return Sample(gsm=gsm, name=name, description=None, properties=props)


def test_assign_arms_from_nde_metadata():
    samples = [s("GSM1", "healthy tissue 2"), s("GSM2", "RA tissue 7"),
               s("GSM3", "treated with metformin", dose="2 mM")]
    t, c, amb = assign_arms(samples, treated=r"\bRA\b|metformin", control=r"healthy|normal")
    assert set(t) == {"GSM2", "GSM3"} and c == ["GSM1"] and amb == []


def test_sample_matching_both_arms_is_ambiguous():
    """A control arm naming the drug is the commonest real failure."""
    samples = [s("GSM1", "vehicle control for dexamethasone")]
    t, c, amb = assign_arms(samples, treated="dexamethasone", control="vehicle|control")
    assert amb == ["GSM1"] and not t and not c


def test_match_columns_falls_back_to_sample_title():
    df = pd.DataFrame({"healthy_tissue_2": [1], "RA_tissue_7": [2]})
    samples = [s("GSM1", "healthy tissue 2"), s("GSM2", "RA tissue 7")]
    # Accessions appear nowhere in the columns; titles normalise onto them.
    assert match_columns(df, ["GSM1"], samples) == ["healthy_tissue_2"]
    assert match_columns(df, ["GSM2"], samples) == ["RA_tissue_7"]


def test_arms_from_columns_last_resort():
    """GSE89408's real shape: columns encode the arm, titles do not match."""
    df = pd.DataFrame({f"normal_tissue_{i}": [1] for i in range(3)} |
                      {f"RA_tissue_{i}": [2] for i in range(4)})
    t, c = arms_from_columns(df, r"^RA_", r"^normal")
    assert len(t) == 4 and len(c) == 3


def test_run_de_reports_arm_source_and_recovers_signal(tmp_path):
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(200)]
    ctrl = rng.poisson(100, (200, 6))
    trt = rng.poisson(100, (200, 6))
    trt[:5] *= 8                                   # 5 genes strongly up
    df = pd.DataFrame(np.hstack([ctrl, trt]), index=genes,
                      columns=[f"normal_{i}" for i in range(6)] +
                              [f"RA_{i}" for i in range(6)])
    path = tmp_path / "m.csv"
    df.to_csv(path)

    res = run_de("GSEX", path, [], [], is_counts=True, patterns=(r"^RA_", r"^normal"))
    assert res.error is None
    assert res.arm_source == "matrix_columns"
    assert res.matched_treated == 6 and res.matched_control == 6
    top = res.gene("G0")
    assert top is not None and top["direction"] == "increased" and top["adj_p"] < 0.01
