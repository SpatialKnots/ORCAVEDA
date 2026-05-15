from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from validate_veda_like_outputs import validate_veda_like_outputs  # noqa: E402


def _write_minimal_veda_like_set(
    outdir: Path,
    *,
    status: str = "PASS",
    warning: str = "",
    contribution_sum: float = 100.0,
) -> None:
    prefix = outdir / "sample"
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "basis_size": 1,
                "expected_vibrational_rank": 1,
                "validation_status": status,
                "warnings": warning,
            }
        ]
    ).to_csv(prefix.with_name(prefix.name + "__veda_like_basis_diagnostics.csv"), index=False)
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "validation_status": status,
                "warnings": warning,
            }
        ]
    ).to_csv(prefix.with_name(prefix.name + "__veda_like_mode_correspondence.csv"), index=False)
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "coord_index": 0,
                "contribution_percent": contribution_sum,
                "validation_status": status,
                "warnings": warning,
            }
        ]
    ).to_csv(prefix.with_name(prefix.name + "__veda_like_ped_matrix.csv"), index=False)
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "veda_like_rank": 1,
                "contribution_percent": contribution_sum,
                "validation_status": status,
                "warnings": warning,
            }
        ]
    ).to_csv(prefix.with_name(prefix.name + "__veda_like_ped_audit.csv"), index=False)


def test_validate_veda_like_outputs_reports_pass_for_clean_artifacts(tmp_path: Path):
    _write_minimal_veda_like_set(tmp_path)

    summary = validate_veda_like_outputs(tmp_path)

    assert summary["validation_status"] == "PASS"
    assert summary["file_count"] == 1
    assert summary["basis_status_counts"] == {"PASS": 1}
    assert summary["mode_correspondence_status_counts"] == {"PASS": 1}
    assert summary["normalization_failure_count"] == 0
    assert summary["warning_token_counts"] == {}


def test_validate_veda_like_outputs_reports_warn_tokens_without_normalization_failure(tmp_path: Path):
    _write_minimal_veda_like_set(
        tmp_path,
        status="WARN",
        warning="fixed_conversion_failed; linear_or_near_linear_fixed_conversion_review; empirical_ratio_only",
    )

    summary = validate_veda_like_outputs(tmp_path)

    assert summary["validation_status"] == "WARN"
    assert summary["basis_status_counts"] == {"WARN": 1}
    assert summary["mode_correspondence_status_counts"] == {"WARN": 1}
    assert summary["normalization_failure_count"] == 0
    assert summary["warning_token_counts"]["fixed_conversion_failed"] == 4
    assert summary["warning_token_counts"]["linear_or_near_linear_fixed_conversion_review"] == 4


def test_validate_veda_like_outputs_reports_fail_for_bad_normalization(tmp_path: Path):
    _write_minimal_veda_like_set(
        tmp_path,
        warning="high_frequency_hbond_dominates_xh_stretch_secondary",
        contribution_sum=95.0,
    )

    summary = validate_veda_like_outputs(tmp_path)

    assert summary["validation_status"] == "FAIL"
    assert summary["normalization_failure_count"] == 1
    assert summary["normalization_failures"][0]["Filename"] == "sample.hess"
    assert summary["normalization_failures"][0]["mode"] == "6"
    assert summary["normalization_failures"][0]["contribution_percent_sum"] == 95.0
    assert summary["warning_token_counts"]["high_frequency_hbond_dominates_xh_stretch_secondary"] == 4
