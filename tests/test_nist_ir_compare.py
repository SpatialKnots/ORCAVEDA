from __future__ import annotations

import sys
from pathlib import Path
import shutil

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir.compare import (  # noqa: E402
    assignment_modes_to_dataframe,
    build_matched_peak_pairs,
    build_scale_engine_payload,
    compare_scale_engines_on_matched_peaks,
    load_orcaveda_assignments,
    match_reference_to_orcaveda,
    pick_reference_peaks,
    reference_points_to_peaks,
)


def test_pick_reference_peaks_simple():
    spectrum = pd.DataFrame(
        {
            "wavenumber_cm-1": [1000, 1004, 1008, 1012, 1016, 1020, 1024],
            "intensity": [0.1, 0.5, 0.2, 0.8, 0.3, 0.6, 0.1],
        }
    )
    peaks = pick_reference_peaks(spectrum, top_n=2, min_separation_cm1=4.0)
    assert len(peaks) == 2
    assert set(round(v, 1) for v in peaks["wavenumber_cm-1"]) == {1012.0, 1020.0}


def test_match_reference_to_orcaveda_simple():
    tmp_path = ROOT / "outputs" / "pytest_nist_ir_compare"
    if tmp_path.exists():
        shutil.rmtree(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    audit_path = tmp_path / "audit.csv"
    pd.DataFrame(
        [
            {
                "mode": 10,
                "frequency_cm-1": 1710.0,
                "IR_intensity": 300.0,
                "functional_group_assignment": "C=O stretch",
            },
            {
                "mode": 11,
                "frequency_cm-1": 1260.0,
                "IR_intensity": 120.0,
                "functional_group_assignment": "C-O stretch",
            },
        ]
    ).to_csv(audit_path, index=False, encoding="utf-8")

    audit = load_orcaveda_assignments(audit_path, scale_factor=0.96)
    reference = pd.DataFrame(
        [
            {"wavenumber_cm-1": 1640.0, "intensity": 0.8},
            {"wavenumber_cm-1": 1210.0, "intensity": 0.4},
        ]
    )
    matched = match_reference_to_orcaveda(reference, audit)

    assert list(matched["orcaveda_mode"]) == [10, 11]
    assert list(matched["orcaveda_assignment"]) == ["C=O stretch", "C-O stretch"]


def test_assignment_modes_to_dataframe_and_reference_points_to_peaks():
    modes = [
        {"mode": 3, "frequency_cm1": 1700.0, "intensity": 200.0, "assignment": "C=O stretch"},
        {"mode": 4, "frequency_cm1": -50.0, "intensity": 10.0, "assignment": "ignore"},
        {"mode": 5, "frequency_cm1": 1260.0, "intensity": 120.0, "assignment": "C-O stretch"},
    ]
    audit = assignment_modes_to_dataframe(modes, scale_factor=0.96)
    assert list(audit["mode"]) == [3, 5]
    assert list(round(v, 1) for v in audit["scaled_frequency_cm-1"]) == [1632.0, 1209.6]

    points = [
        {"x": 1200.0, "y": 0.2},
        {"x": 1210.0, "y": 0.8},
        {"x": 1220.0, "y": 0.3},
        {"x": 1625.0, "y": 0.4},
        {"x": 1632.0, "y": 1.0},
        {"x": 1640.0, "y": 0.5},
    ]
    peaks = reference_points_to_peaks(points, top_n=2, min_separation_cm1=8.0)
    assert len(peaks) == 2
    assert set(round(v, 1) for v in peaks["wavenumber_cm-1"]) == {1210.0, 1632.0}


def test_reference_points_to_peaks_condensed_phase_is_more_selective():
    points = [
        {"x": 980.0, "y": 0.10},
        {"x": 1000.0, "y": 0.80},
        {"x": 1010.0, "y": 0.72},
        {"x": 1020.0, "y": 0.15},
        {"x": 1600.0, "y": 0.95},
        {"x": 1610.0, "y": 0.84},
        {"x": 1620.0, "y": 0.20},
    ]
    gas_peaks = reference_points_to_peaks(points, top_n=4, min_separation_cm1=8.0)
    condensed_peaks = reference_points_to_peaks(
        points,
        top_n=4,
        min_separation_cm1=8.0,
        reference_context={"phase_tag": "liquid_neat", "state": "LIQUID"},
    )

    assert len(gas_peaks) >= len(condensed_peaks)
    assert len(condensed_peaks) == 2
    assert set(round(v, 1) for v in condensed_peaks["wavenumber_cm-1"]) == {1000.0, 1600.0}


def test_reference_points_to_peaks_transmittance_uses_minima():
    points = [
        {"x": 1000.0, "y": 92.0},
        {"x": 1010.0, "y": 70.0},
        {"x": 1020.0, "y": 90.0},
        {"x": 1600.0, "y": 88.0},
        {"x": 1610.0, "y": 60.0},
        {"x": 1620.0, "y": 85.0},
    ]
    peaks = reference_points_to_peaks(
        points,
        top_n=2,
        min_separation_cm1=8.0,
        reference_context={"y_units": "TRANSMITTANCE"},
    )
    assert len(peaks) == 2
    assert set(round(v, 1) for v in peaks["wavenumber_cm-1"]) == {1010.0, 1610.0}


def test_build_matched_peak_pairs_and_engine_table():
    assignment_audit = pd.DataFrame(
        [
            {
                "mode": 10,
                "frequency_cm-1": 1700.0,
                "IR_intensity": 300.0,
                "functional_group_assignment": "C=O stretch",
                "scaled_frequency_cm-1": 1632.0,
            },
            {
                "mode": 11,
                "frequency_cm-1": 1260.0,
                "IR_intensity": 120.0,
                "functional_group_assignment": "C-O stretch",
                "scaled_frequency_cm-1": 1209.6,
            },
            {
                "mode": 12,
                "frequency_cm-1": 3050.0,
                "IR_intensity": 80.0,
                "functional_group_assignment": "C-H stretch",
                "scaled_frequency_cm-1": 2928.0,
            },
        ]
    )
    reference_peaks = pd.DataFrame(
        [
            {"wavenumber_cm-1": 1635.0, "intensity": 1.0},
            {"wavenumber_cm-1": 1212.0, "intensity": 0.7},
            {"wavenumber_cm-1": 2930.0, "intensity": 0.3},
        ]
    )
    matched_pairs = build_matched_peak_pairs(reference_peaks, assignment_audit)
    assert list(matched_pairs["mode"]) == [10, 11, 12]
    assert list(matched_pairs["assignment"]) == ["C=O stretch", "C-O stretch", "C-H stretch"]

    matched_pairs_scored = build_matched_peak_pairs(reference_peaks, assignment_audit, method="scored")
    assert set(matched_pairs_scored["mode"]) == {10, 11, 12}
    matched_pairs_hc = build_matched_peak_pairs(reference_peaks, assignment_audit, method="scored_high_confidence")
    assert set(matched_pairs_hc["mode"]) <= {10, 11, 12}
    assert len(matched_pairs_hc) <= len(matched_pairs_scored)

    engine_table = compare_scale_engines_on_matched_peaks(matched_pairs)
    assert not engine_table.empty
    assert {"engine", "parameters_json", "mean_percent_deviation", "rmse_percent_deviation", "max_percent_deviation", "matched_count"} <= set(engine_table.columns)
    assert {"global_ls", "global_weighted_ls", "global_huber", "piecewise_region", "power_law"} <= set(engine_table["engine"])


def test_build_scale_engine_payload_contains_engine_fits():
    assignment_audit = pd.DataFrame(
        [
            {
                "mode": 1,
                "frequency_cm-1": 1710.0,
                "IR_intensity": 250.0,
                "functional_group_assignment": "C=O stretch",
                "scaled_frequency_cm-1": 1641.6,
            },
            {
                "mode": 2,
                "frequency_cm-1": 1265.0,
                "IR_intensity": 110.0,
                "functional_group_assignment": "C-O stretch",
                "scaled_frequency_cm-1": 1214.4,
            },
        ]
    )
    reference_peaks = pd.DataFrame(
        [
            {"wavenumber_cm-1": 1643.0, "intensity": 0.9},
            {"wavenumber_cm-1": 1215.0, "intensity": 0.6},
        ]
    )
    payload = build_scale_engine_payload(reference_peaks, assignment_audit)
    assert "matched_pairs" in payload
    assert "high_confidence_matched_pairs" in payload
    assert "extended_matched_pairs" in payload
    assert "nearest_matched_pairs" in payload
    assert "engine_table" in payload
    assert "engine_fits" in payload
    assert "matching_layers" in payload
    assert "matching_layer_overview" in payload
    assert "engine_layer_matrix" in payload
    assert payload["default_manual_scale"] > 0.0
    assert "global_ls" in payload["engine_fits"]
    assert "power_law" in payload["engine_fits"]
    assert "high_confidence" in payload["matching_layers"]
    assert "extended" in payload["matching_layers"]
    assert "nearest" in payload["matching_layers"]
    assert "mean_percent_deviation" in payload["engine_fits"]["global_ls"]["metrics"]
    assert "max_percent_deviation" in payload["engine_fits"]["global_ls"]["metrics"]
