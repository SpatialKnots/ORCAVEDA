from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "benchmarks" / "veda_compare"
if str(BENCHMARK) not in sys.path:
    sys.path.insert(0, str(BENCHMARK))

from convert_veda_reference import convert_veda_reference  # noqa: E402


def test_veda_reference_ingest_skips_when_raw_directory_missing(tmp_path: Path):
    summary = convert_veda_reference(tmp_path / "missing", tmp_path / "out")

    assert summary["conversion_status"] == "SKIP"
    assert summary["acceptance_status"] == "SKIP"
    assert summary["reason"] == "raw_reference_directory_missing"
    assert (tmp_path / "out" / "veda_reference_ingest_summary.json").is_file()


def test_veda_reference_ingest_copies_normalized_reference_csv(tmp_path: Path):
    raw = tmp_path / "raw"
    out = tmp_path / "out"
    raw.mkdir()
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "internal_coordinate": "r(O1-H2)",
                "contribution_percent": 99.5,
            }
        ]
    ).to_csv(raw / "veda_reference_ped_matrix.csv", index=False)
    pd.DataFrame(
        [
            {
                "Filename": "sample.hess",
                "mode": 6,
                "internal_coordinate": "r(O1-H2)",
            }
        ]
    ).to_csv(raw / "veda_reference_dominant_assignments.csv", index=False)

    summary = convert_veda_reference(raw, out)

    assert summary["conversion_status"] == "PASS"
    assert summary["matrix_rows"] == 1
    assert summary["dominant_rows"] == 1
    matrix = pd.read_csv(out / "veda_reference_ped_matrix.csv")
    assert list(matrix.columns) == ["Filename", "mode", "internal_coordinate", "contribution_percent"]
    assert matrix.loc[0, "contribution_percent"] == 99.5


def test_veda_reference_ingest_uses_explicit_column_mapping(tmp_path: Path):
    raw = tmp_path / "raw"
    out = tmp_path / "out"
    raw.mkdir()
    matrix_csv = raw / "raw_matrix.csv"
    pd.DataFrame(
        [
            {
                "file": "sample.hess",
                "mode_no": "6",
                "veda_coord": "r(O1-H2)",
                "ped_percent": "88.25",
            }
        ]
    ).to_csv(matrix_csv, index=False)

    summary = convert_veda_reference(
        raw,
        out,
        matrix_csv=matrix_csv,
        matrix_mapping={
            "Filename": "file",
            "mode": "mode_no",
            "internal_coordinate": "veda_coord",
            "contribution_percent": "ped_percent",
        },
    )

    assert summary["conversion_status"] == "PASS"
    matrix = pd.read_csv(out / "veda_reference_ped_matrix.csv")
    assert matrix.loc[0, "Filename"] == "sample.hess"
    assert matrix.loc[0, "mode"] == 6
    assert matrix.loc[0, "internal_coordinate"] == "r(O1-H2)"
    assert matrix.loc[0, "contribution_percent"] == 88.25


def test_veda_reference_ingest_fails_for_missing_required_columns(tmp_path: Path):
    raw = tmp_path / "raw"
    out = tmp_path / "out"
    raw.mkdir()
    pd.DataFrame([{"Filename": "sample.hess", "mode": 6}]).to_csv(
        raw / "veda_reference_ped_matrix.csv",
        index=False,
    )

    summary = convert_veda_reference(raw, out)

    assert summary["conversion_status"] == "FAIL"
    assert summary["acceptance_status"] == "FAIL"
    assert "missing columns" in summary["reason"]
