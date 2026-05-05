from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks" / "vibrational_assignments"
if str(BENCH) not in sys.path:
    sys.path.insert(0, str(BENCH))

from compare_orcaveda_assignments import compare  # noqa: E402


OUTROOT = ROOT / "outputs" / "pytest_vibrational_assignment_multiscale"


def _write_csv(path: Path, rows: list[dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_acetophenone_like_carbonyl_multiscale_raw_and_scaled():
    outdir = OUTROOT / "acetophenone_like"
    benchmark = _write_csv(
        outdir / "benchmark.csv",
        [
            {
                "molecule": "acetophenone",
                "hess_file": "acetophenone.hess",
                "observed_ir_cm1": 1685.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "aryl ketone C=O stretch",
                "mode_family": "stretch",
                "functional_group": "aryl ketone carbonyl",
                "confidence": "gold",
            }
        ],
    )
    audit = _write_csv(
        outdir / "audit.csv",
        [
            {
                "Filename": "acetophenone.hess",
                "mode": 20,
                "frequency_cm-1": 1690.0,
                "IR_intensity": 80.0,
                "functional_group_assignment": "aromatic ring stretch",
                "top_internal_coordinates": "r(C1-C2)=70.0%",
                "top1_coord": "r(C1-C2)",
                "assignment_confidence": "medium",
            },
            {
                "Filename": "acetophenone.hess",
                "mode": 21,
                "frequency_cm-1": 1755.0,
                "IR_intensity": 100.0,
                "functional_group_assignment": "aryl-conjugated C=O stretch",
                "top_internal_coordinates": "carbonyl_CO_stretch(C1=O2)=90.0%",
                "top1_coord": "carbonyl_CO_stretch(C1=O2)",
                "assignment_confidence": "high",
            },
        ],
    )
    ped = _write_csv(
        outdir / "ped.csv",
        [
            {
                "Filename": "acetophenone.hess",
                "mode": 20,
                "frequency_cm-1": 1690.0,
                "ped_rank": 1,
                "coordinate_family": "aromatic ring stretch",
                "internal_coordinate": "r(C1-C2)",
                "coordinate_class": "stretch",
                "contribution_percent": 70.0,
                "ped_warnings": "",
            },
            {
                "Filename": "acetophenone.hess",
                "mode": 21,
                "frequency_cm-1": 1755.0,
                "ped_rank": 1,
                "coordinate_family": "aryl-conjugated C=O stretch",
                "internal_coordinate": "carbonyl_CO_stretch(C1=O2)",
                "coordinate_class": "stretch",
                "contribution_percent": 88.0,
                "ped_warnings": "",
            },
        ],
    )

    raw = compare(
        benchmark,
        audit,
        outdir / "raw_out.csv",
        max_delta_cm1=80.0,
        windows_cm1=[50.0, 100.0, 200.0],
        scale_factor=0.96,
        primary_frequency="raw",
        ped_audit_csv=ped,
    ).iloc[0]

    assert int(raw["nearest_mode"]) == 20
    assert int(raw["orcaveda_mode"]) == 21
    assert raw["match_strategy"] == "semantic_within_window"
    assert raw["raw_window_50_status"] == "FAIL"
    assert raw["raw_window_100_status"] == "FAIL"
    assert raw["raw_window_200_status"] == "PASS"
    assert raw["scaled_window_50_status"] == "PASS"
    assert raw["scaled_window_50_assignment"] == "aryl-conjugated C=O stretch"
    assert raw["ped_semantic_status"] == "PASS"
    assert raw["stage3d_ped_warning"].startswith("ped_confirms_stage3d")
    assert "nearest_frequency_differs_from_semantic_mode" in raw["stage3d_ped_warning"]
    assert raw["raw_ped_window_200_status"] == "PASS"
    assert raw["scaled_ped_window_50_status"] == "PASS"
    assert "C=O stretch" in raw["ped_top_contributors"]

    scaled = compare(
        benchmark,
        audit,
        outdir / "scaled_out.csv",
        max_delta_cm1=80.0,
        windows_cm1=[50.0, 100.0, 200.0],
        scale_factor=0.96,
        primary_frequency="scaled",
    ).iloc[0]

    assert int(scaled["nearest_mode"]) == 21
    assert int(scaled["orcaveda_mode"]) == 21
    assert abs(float(scaled["delta_cm-1"])) < 1.0
    assert scaled["status"] == "PASS"


def test_broad_ch_stretch_windows_keep_semantic_and_frequency_separate():
    outdir = OUTROOT / "broad_ch"
    benchmark = _write_csv(
        outdir / "benchmark.csv",
        [
            {
                "molecule": "synthetic",
                "hess_file": "synthetic.hess",
                "observed_ir_cm1": 3060.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "aromatic C-H stretch",
                "mode_family": "stretch",
                "functional_group": "aromatic C-H",
                "confidence": "gold",
            }
        ],
    )
    audit = _write_csv(
        outdir / "audit.csv",
        [
            {
                "Filename": "synthetic.hess",
                "mode": 30,
                "frequency_cm-1": 3000.0,
                "IR_intensity": 20.0,
                "functional_group_assignment": "ring bend",
                "top_internal_coordinates": "ang(C1-C2-C3)=60.0%",
                "top1_coord": "ang(C1-C2-C3)",
                "assignment_confidence": "low",
            },
            {
                "Filename": "synthetic.hess",
                "mode": 31,
                "frequency_cm-1": 3188.0,
                "IR_intensity": 90.0,
                "functional_group_assignment": "aromatic C-H stretch",
                "top_internal_coordinates": "aromatic_CH_stretch(C1-H2)=95.0%",
                "top1_coord": "aromatic_CH_stretch(C1-H2)",
                "assignment_confidence": "high",
            },
        ],
    )

    result = compare(
        benchmark,
        audit,
        outdir / "out.csv",
        max_delta_cm1=80.0,
        windows_cm1=[50.0, 100.0, 200.0, 500.0],
        scale_factor=0.96,
        primary_frequency="raw",
    ).iloc[0]

    assert result["raw_window_50_status"] == "FAIL"
    assert result["raw_window_500_status"] == "PASS"
    assert result["scaled_window_50_status"] == "PASS"
    assert result["scaled_window_50_assignment"] == "aromatic C-H stretch"


def test_ped_aware_comparator_accepts_wilson_rank_schema():
    outdir = OUTROOT / "wilson_schema"
    benchmark = _write_csv(
        outdir / "benchmark.csv",
        [
            {
                "molecule": "water",
                "hess_file": "water.hess",
                "observed_ir_cm1": 1600.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "H-O-H bend",
                "mode_family": "bend",
                "functional_group": "H-O-H",
                "confidence": "gold",
            }
        ],
    )
    audit = _write_csv(
        outdir / "audit.csv",
        [
            {
                "Filename": "water.hess",
                "mode": 6,
                "frequency_cm-1": 1590.0,
                "IR_intensity": 100.0,
                "functional_group_assignment": "H-O-H bend",
                "top_internal_coordinates": "ang(H2-O1-H3)=99.0%",
                "top1_coord": "ang(H2-O1-H3)",
                "assignment_confidence": "high",
            }
        ],
    )
    wilson = _write_csv(
        outdir / "wilson.csv",
        [
            {
                "Filename": "water.hess",
                "mode": 6,
                "frequency_cm-1": 1590.0,
                "wilson_rank": 1,
                "coordinate_family": "H-O-H bend",
                "internal_coordinate": "ang(H2-O1-H3)",
                "coordinate_class": "bend",
                "contribution_percent": 99.0,
                "wilson_ped_warnings": "",
            }
        ],
    )

    result = compare(
        benchmark,
        audit,
        outdir / "out.csv",
        max_delta_cm1=80.0,
        windows_cm1=[50.0],
        scale_factor=1.0,
        primary_frequency="raw",
        ped_audit_csv=wilson,
    ).iloc[0]

    assert result["ped_semantic_status"] == "PASS"
    assert "H-O-H bend" in result["ped_top_contributors"]
    assert result["raw_ped_window_50_status"] == "PASS"
