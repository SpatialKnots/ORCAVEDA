from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

COMPARE = ROOT / "benchmarks" / "bmatrix_compare"
if str(COMPARE) not in sys.path:
    sys.path.insert(0, str(COMPARE))

from compare_bmatrix_methods import compare_hess_file, compare_hess_files, main  # noqa: E402


def test_bmatrix_method_compare_h2o_reports_no_delta_or_rank_change():
    summary, rows, selection_differences = compare_hess_file(ROOT / "data" / "hess" / "H2O_freq.hess")

    assert summary["Filename"] == "H2O_freq.hess"
    assert summary["internal_count"] == len(rows)
    assert summary["rows_above_tolerance"] == 0
    assert summary["selected_basis_difference_count"] == 0
    assert summary["redundant_finite_rank"] == summary["redundant_hybrid_rank"]
    assert summary["finite_selected_rank"] == summary["hybrid_selected_rank"]
    assert summary["method_counts"] == {"analytical_distance": 2, "analytical_angle": 1}
    assert selection_differences == []
    assert rows[0]["atoms0"] == "[0,1]"
    assert rows[-1]["angle_degrees"] != ""
    assert rows[-1]["angle_sine"] != ""


def test_bmatrix_method_compare_cli_writes_summary_and_row_artifacts(tmp_path):
    outdir = tmp_path / "bmatrix_compare"
    exit_code = main(
        [
            "--hess",
            "data/hess/H2O_freq.hess",
            "data/hess/NH3.hess",
            "--out",
            str(outdir),
        ]
    )

    assert exit_code == 0
    summary_path = outdir / "bmatrix_method_comparison_summary.json"
    row_path = outdir / "bmatrix_method_comparison_rows.csv"
    selected_path = outdir / "bmatrix_method_comparison_selected_basis_differences.csv"
    assert summary_path.is_file()
    assert row_path.is_file()
    assert selected_path.is_file()

    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["comparison_scope"] == "finite_difference_B_vs_hybrid_analytical_B"
    assert payload["method_boundary"] == "diagnostic only; production pipeline remains finite_difference_B"
    assert payload["file_count"] == 2
    assert payload["files_with_selected_rank_change"] == 0
    assert payload["selected_basis_difference_count"] == 0
    assert payload["rows_above_tolerance_count"] == 0
    assert payload["selected_basis_replacement_rank_loss_count"] == 0


def test_bmatrix_method_compare_reports_selection_replacement_metrics():
    summary, _rows, selection_differences = compare_hess_file(ROOT / "data" / "hess" / "benzene.hess")

    assert summary["selected_basis_difference_count"] == 2
    first = selection_differences[0]
    assert first["replacement_rank_preserved"] is True
    assert first["finite_basis_rank_with_hybrid_row"] == summary["finite_selected_rank"]
    assert first["hybrid_basis_rank_with_finite_row"] == summary["hybrid_selected_rank"]
    assert first["finite_basis_min_singular_with_hybrid_row"] > 0.0
    assert first["hybrid_basis_min_singular_with_finite_row"] > 0.0


def test_bmatrix_method_compare_full_sweep_acceptance_policy_allows_rank_preserving_selection_swaps():
    hess_paths = sorted((ROOT / "data" / "hess").glob("*.hess"))
    summaries, _rows, selection_differences = compare_hess_files(hess_paths)

    assert len(summaries) == 55
    assert sum(int(row["rows_above_tolerance"]) for row in summaries) == 0
    assert all(int(row["redundant_finite_rank"]) == int(row["redundant_hybrid_rank"]) for row in summaries)
    assert all(int(row["finite_selected_rank"]) == int(row["hybrid_selected_rank"]) for row in summaries)
    assert all(bool(row["replacement_rank_preserved"]) for row in selection_differences)

    differing_files = {str(row["Filename"]) for row in selection_differences}
    assert differing_files == {"aniline.hess", "benzene.hess", "benzonitrile.hess", "pyridine.hess"}
    assert len(selection_differences) == 8
