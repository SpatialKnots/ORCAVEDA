from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks" / "vibrational_assignments"
if str(BENCH) not in sys.path:
    sys.path.insert(0, str(BENCH))

from compare_orcaveda_assignments import compare, ped_coverage_audit  # noqa: E402


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


def test_ped_final_assignment_comparator_and_coverage():
    outdir = OUTROOT / "ped_final"
    benchmark = _write_csv(
        outdir / "benchmark.csv",
        [
            {
                "molecule": "synthetic",
                "hess_file": "synthetic.hess",
                "observed_ir_cm1": 1700.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "C=O stretch",
                "mode_family": "stretch",
                "functional_group": "carbonyl",
                "confidence": "gold",
            },
            {
                "molecule": "synthetic",
                "hess_file": "synthetic.hess",
                "observed_ir_cm1": 1000.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "O-H stretch",
                "mode_family": "stretch",
                "functional_group": "O-H",
                "confidence": "gold",
            },
        ],
    )
    legacy_audit = _write_csv(
        outdir / "audit.csv",
        [
            {
                "Filename": "synthetic.hess",
                "mode": 1,
                "frequency_cm-1": 1700.0,
                "functional_group_assignment": "aromatic ring stretch",
            },
            {
                "Filename": "synthetic.hess",
                "mode": 2,
                "frequency_cm-1": 1000.0,
                "functional_group_assignment": "O-H stretch",
            },
        ],
    )
    ped_final = _write_csv(
        outdir / "ped_final.csv",
        [
            {
                "Filename": "synthetic.hess",
                "mode": 1,
                "frequency_cm-1": 1700.0,
                "final_assignment": "C=O stretch",
                "final_assignment_source": "Wilson GF-style PED audit",
                "final_assignment_policy": "ped_confirms_stage3d",
                "final_assignment_warning": "ped_diagnostic_basis",
                "stage3d_assignment": "aromatic ring stretch",
                "ped_assignment": "C=O stretch",
                "ped_source": "Wilson GF-style PED audit",
                "ped_agreement_status": "confirms",
                "ped_policy_warning": "ped_diagnostic_basis",
                "ped_top_family": "C=O stretch",
                "ped_top_percent": 76.0,
                "ped_top_contributors": "C=O stretch [carbonyl_CO_stretch(C1=O2)] 76.0%",
            },
            {
                "Filename": "synthetic.hess",
                "mode": 2,
                "frequency_cm-1": 1000.0,
                "final_assignment": "O-H stretch",
                "final_assignment_source": "Stage 3D assignment audit",
                "final_assignment_policy": "stage3d_fallback_due_to_ped_disagreement",
                "final_assignment_warning": "ped_stage3d_semantic_disagreement; final_label_kept_stage3d",
                "stage3d_assignment": "O-H stretch",
                "ped_assignment": "H-O-H bend",
                "ped_source": "Wilson GF-style PED audit",
                "ped_agreement_status": "disagrees",
                "ped_policy_warning": "ped_stage3d_semantic_disagreement",
                "ped_top_family": "H-O-H bend",
                "ped_top_percent": 90.0,
                "ped_top_contributors": "H-O-H bend [ang(H2-O1-H3)] 90.0%",
            },
        ],
    )

    result = compare(
        benchmark,
        legacy_audit,
        outdir / "out.csv",
        max_delta_cm1=80.0,
        windows_cm1=[50.0],
        scale_factor=1.0,
        primary_frequency="raw",
        ped_final_assignment_csv=ped_final,
    )

    by_mode = {int(row["orcaveda_mode"]): row for _, row in result.iterrows()}
    assert by_mode[1]["status"] == "PASS"
    assert by_mode[1]["orcaveda_assignment"] == "C=O stretch"
    assert by_mode[1]["final_assignment_policy"] == "ped_confirms_stage3d"
    assert by_mode[2]["status"] == "PASS"
    assert by_mode[2]["orcaveda_assignment"] == "O-H stretch"
    assert by_mode[2]["final_assignment_policy"] == "stage3d_fallback_due_to_ped_disagreement"

    coverage, detail = ped_coverage_audit(pd.read_csv(ped_final))
    assert coverage.iloc[0]["modes"] == 2
    assert coverage.iloc[0]["ped_final"] == 1
    assert coverage.iloc[0]["stage3d_fallback"] == 1
    assert coverage.iloc[0]["disagrees"] == 1
    assert len(detail) == 2


def test_carboxylic_acid_mixed_rows_accept_acid_context_without_explicit_oh():
    outdir = OUTROOT / "carboxylic_acid_context"
    benchmark = _write_csv(
        outdir / "benchmark.csv",
        [
            {
                "molecule": "benzoic acid",
                "hess_file": "benzoic_acid.hess",
                "observed_ir_cm1": 1767.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "mixed mode",
                "mode_family": "carbonyl stretch / bend",
                "functional_group": "carboxylic acid C=O; carboxylic acid O-H",
                "confidence": "gold",
            },
            {
                "molecule": "benzoic acid",
                "hess_file": "benzoic_acid.hess",
                "observed_ir_cm1": 1162.0,
                "observed_raman_cm1": "",
                "calculated_cm1": "",
                "assignment_normalized": "mixed mode",
                "mode_family": "mixed in-plane bend",
                "functional_group": "aromatic C-H; carboxylic acid O-H",
                "confidence": "gold",
            },
        ],
    )
    legacy_audit = _write_csv(
        outdir / "audit.csv",
        [
            {
                "Filename": "benzoic_acid.hess",
                "mode": 38,
                "frequency_cm-1": 1767.0,
                "functional_group_assignment": "legacy placeholder",
            },
            {
                "Filename": "benzoic_acid.hess",
                "mode": 28,
                "frequency_cm-1": 1162.0,
                "functional_group_assignment": "legacy placeholder",
            },
        ],
    )
    ped_final = _write_csv(
        outdir / "ped_final.csv",
        [
            {
                "Filename": "benzoic_acid.hess",
                "mode": 38,
                "frequency_cm-1": 1767.0,
                "final_assignment": "carboxylic C=O stretch mixed with C-C-O bend",
                "final_assignment_source": "Wilson GF-style PED audit",
                "final_assignment_policy": "ped_adds_context",
                "stage3d_assignment": "carboxylic C=O stretch mixed with carboxylic O-H bend",
                "ped_assignment": "carboxylic C=O stretch mixed with C-C-O bend",
                "ped_agreement_status": "adds_context",
                "ped_top_family": "carboxylic C=O stretch",
                "ped_top_percent": 52.0,
                "ped_top_contributors": "carboxylic C=O stretch 52.0%; C-C-O bend 12.0%; O-C-O bend 9.0%",
            },
            {
                "Filename": "benzoic_acid.hess",
                "mode": 28,
                "frequency_cm-1": 1162.0,
                "final_assignment": "C-C-H bend",
                "final_assignment_source": "Wilson GF-style PED audit",
                "final_assignment_policy": "ped_confirms_stage3d",
                "stage3d_assignment": "C-C-H bend",
                "ped_assignment": "C-C-H bend",
                "ped_agreement_status": "confirms",
                "ped_top_family": "C-C-H bend",
                "ped_top_percent": 31.0,
                "ped_top_contributors": "C-C-H bend 31.0%; C-C-H bend 23.9%",
            },
        ],
    )

    result = compare(
        benchmark,
        legacy_audit,
        outdir / "out.csv",
        max_delta_cm1=80.0,
        windows_cm1=[50.0],
        scale_factor=1.0,
        primary_frequency="raw",
        ped_final_assignment_csv=ped_final,
    )

    by_freq = {float(row["target_frequency_cm-1"]): row for _, row in result.iterrows()}
    assert by_freq[1767.0]["status"] == "WARN"
    assert by_freq[1767.0]["reason"] == "acid_context_without_explicit_oh"
    assert by_freq[1162.0]["status"] == "FAIL"
    assert by_freq[1162.0]["reason"] == "missing_oh"
