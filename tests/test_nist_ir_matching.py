from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir.matching import (  # noqa: E402
    build_calculated_modes,
    build_experimental_peaks,
    generate_match_candidates,
    infer_assignment_class,
    infer_region,
    match_reference_to_orcaveda_v2,
    solve_peak_matching,
)


def test_infer_region_and_assignment_class():
    assert infer_region(420.0) == "low"
    assert infer_region(1260.0) == "fingerprint"
    assert infer_region(1710.0) == "double_bond"
    assert infer_region(3200.0) == "xh_stretch"
    assert infer_assignment_class("aryl-conjugated C=O stretch") == "carbonyl"
    assert infer_assignment_class("O-H stretch") == "oh_stretch"
    assert infer_assignment_class("torsion") == "torsion"


def test_scored_matching_builds_one_to_one_pairs():
    reference_peaks = pd.DataFrame(
        [
            {"wavenumber_cm-1": 1643.0, "intensity": 1.0},
            {"wavenumber_cm-1": 1215.0, "intensity": 0.7},
        ]
    )
    assignment_audit = pd.DataFrame(
        [
            {
                "mode": 10,
                "scaled_frequency_cm-1": 1641.0,
                "IR_intensity": 300.0,
                "functional_group_assignment": "C=O stretch",
                "warnings": "",
            },
            {
                "mode": 11,
                "scaled_frequency_cm-1": 1214.0,
                "IR_intensity": 120.0,
                "functional_group_assignment": "C-O stretch",
                "warnings": "",
            },
            {
                "mode": 12,
                "scaled_frequency_cm-1": 1645.0,
                "IR_intensity": 10.0,
                "functional_group_assignment": "angle bend",
                "warnings": "",
            },
        ]
    )

    exp_peaks = build_experimental_peaks(reference_peaks)
    calc_modes = build_calculated_modes(assignment_audit)
    candidates = generate_match_candidates(exp_peaks, calc_modes)
    solved = solve_peak_matching(candidates, exp_peaks, calc_modes)
    matched_pairs = solved["matched_pairs"]

    assert len(matched_pairs) == 2
    assert {pair.mode for pair in matched_pairs} == {10, 11}


def test_match_reference_to_orcaveda_v2_dataframe():
    reference_peaks = pd.DataFrame(
        [
            {"wavenumber_cm-1": 3060.0, "intensity": 0.9},
            {"wavenumber_cm-1": 1702.0, "intensity": 1.0},
        ]
    )
    assignment_audit = pd.DataFrame(
        [
            {
                "mode": 5,
                "scaled_frequency_cm-1": 3058.0,
                "IR_intensity": 80.0,
                "functional_group_assignment": "aromatic C-H stretch",
                "warnings": "",
            },
            {
                "mode": 6,
                "scaled_frequency_cm-1": 1700.0,
                "IR_intensity": 200.0,
                "functional_group_assignment": "C=O stretch",
                "warnings": "",
            },
        ]
    )
    matched = match_reference_to_orcaveda_v2(reference_peaks, assignment_audit)
    assert set(matched["orcaveda_mode"]) == {5, 6}
    assert {"match_confidence", "orcaveda_assignment_class", "total_cost"} <= set(matched.columns)
