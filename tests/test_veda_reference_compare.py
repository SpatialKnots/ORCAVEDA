from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "benchmarks" / "veda_compare"
if str(BENCHMARK) not in sys.path:
    sys.path.insert(0, str(BENCHMARK))

from compare_veda_outputs import compare_veda_outputs  # noqa: E402


def _write_orcaveda_veda_like(outdir: Path, *, percent: float = 100.0, dominant: str = "r(O1-H2)") -> None:
    prefix = outdir / "sample"
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "internal_coordinate": "r(O1-H2)",
                "contribution_percent": percent,
            }
        ]
    ).to_csv(prefix.with_name(prefix.name + "__veda_like_ped_matrix.csv"), index=False)
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "veda_like_rank": 1,
                "internal_coordinate": dominant,
                "coordinate_family": "O-H stretch",
                "contribution_percent": percent,
            }
        ]
    ).to_csv(prefix.with_name(prefix.name + "__veda_like_ped_audit.csv"), index=False)


def _write_reference(reference_dir: Path, *, percent: float = 100.0, dominant: str = "r(O1-H2)") -> None:
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "internal_coordinate": "r(O1-H2)",
                "contribution_percent": percent,
            }
        ]
    ).to_csv(reference_dir / "veda_reference_ped_matrix.csv", index=False)
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "internal_coordinate": dominant,
            }
        ]
    ).to_csv(reference_dir / "veda_reference_dominant_assignments.csv", index=False)


def test_veda_reference_compare_skips_when_reference_directory_missing(tmp_path: Path):
    orcaveda_dir = tmp_path / "orcaveda"
    out_dir = tmp_path / "comparison"
    orcaveda_dir.mkdir()
    _write_orcaveda_veda_like(orcaveda_dir)

    summary = compare_veda_outputs(orcaveda_dir, tmp_path / "missing_reference", out_dir)

    assert summary["comparison_status"] == "SKIP"
    assert summary["acceptance_status"] == "SKIP"
    assert summary["reason"] == "veda_reference_directory_missing"
    assert (out_dir / "veda_reference_comparison_summary.json").is_file()


def test_veda_reference_compare_passes_for_matching_synthetic_rows(tmp_path: Path):
    orcaveda_dir = tmp_path / "orcaveda"
    reference_dir = tmp_path / "reference"
    orcaveda_dir.mkdir()
    reference_dir.mkdir()
    _write_orcaveda_veda_like(orcaveda_dir, percent=99.0)
    _write_reference(reference_dir, percent=100.0)

    summary = compare_veda_outputs(orcaveda_dir, reference_dir, tmp_path / "comparison", tolerance_percent=2.0)

    assert summary["comparison_status"] == "PASS"
    assert summary["acceptance_status"] == "PASS"
    assert summary["matched_matrix_rows"] == 1
    assert summary["max_abs_delta_percent"] == 1.0
    assert summary["dominant_mismatch_count"] == 0


def test_veda_reference_compare_fails_for_out_of_tolerance_percent(tmp_path: Path):
    orcaveda_dir = tmp_path / "orcaveda"
    reference_dir = tmp_path / "reference"
    orcaveda_dir.mkdir()
    reference_dir.mkdir()
    _write_orcaveda_veda_like(orcaveda_dir, percent=90.0)
    _write_reference(reference_dir, percent=100.0)

    summary = compare_veda_outputs(orcaveda_dir, reference_dir, tmp_path / "comparison", tolerance_percent=2.0)

    assert summary["comparison_status"] == "FAIL"
    assert summary["acceptance_status"] == "FAIL"
    assert summary["out_of_tolerance_count"] == 1
    assert "ped_percent_delta_out_of_tolerance" in summary["reason"]


def test_veda_reference_compare_fails_for_dominant_mismatch(tmp_path: Path):
    orcaveda_dir = tmp_path / "orcaveda"
    reference_dir = tmp_path / "reference"
    orcaveda_dir.mkdir()
    reference_dir.mkdir()
    _write_orcaveda_veda_like(orcaveda_dir, dominant="ang(H2-O1-H3)")
    _write_reference(reference_dir, dominant="r(O1-H2)")

    summary = compare_veda_outputs(orcaveda_dir, reference_dir, tmp_path / "comparison")

    assert summary["comparison_status"] == "FAIL"
    assert summary["dominant_mismatch_count"] == 1
    assert "dominant_coordinate_mismatch" in summary["reason"]
