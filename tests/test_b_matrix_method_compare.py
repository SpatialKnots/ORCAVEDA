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

from compare_bmatrix_methods import compare_hess_file, main  # noqa: E402


def test_bmatrix_method_compare_h2o_reports_no_delta_or_rank_change():
    summary, rows = compare_hess_file(ROOT / "data" / "hess" / "H2O_freq.hess")

    assert summary["Filename"] == "H2O_freq.hess"
    assert summary["internal_count"] == len(rows)
    assert summary["rows_above_tolerance"] == 0
    assert summary["redundant_finite_rank"] == summary["redundant_hybrid_rank"]
    assert summary["finite_selected_rank"] == summary["hybrid_selected_rank"]
    assert summary["method_counts"] == {"analytical_distance": 2, "analytical_angle": 1}


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
    assert summary_path.is_file()
    assert row_path.is_file()

    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["comparison_scope"] == "finite_difference_B_vs_hybrid_analytical_B"
    assert payload["method_boundary"] == "diagnostic only; production pipeline remains finite_difference_B"
    assert payload["file_count"] == 2
    assert payload["files_with_selected_rank_change"] == 0
