#!/usr/bin/env python3
"""Normalize checked-in original VEDA reference rows for comparison.

This converter is conservative by design. It only accepts an already-normalized
reference matrix or an explicit column mapping for a CSV file. Unknown native
VEDA exports are reported as SKIP/FAIL diagnostics; no reference rows are
fabricated or inferred from ORCAVEDA outputs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
from pandas.errors import EmptyDataError, ParserError


REFERENCE_MATRIX_NAME = "veda_reference_ped_matrix.csv"
REFERENCE_DOMINANT_NAME = "veda_reference_dominant_assignments.csv"
INGEST_SUMMARY_NAME = "veda_reference_ingest_summary.json"

MATRIX_COLUMNS = ["Filename", "mode", "internal_coordinate", "contribution_percent"]
DOMINANT_COLUMNS = ["Filename", "mode", "internal_coordinate"]


def _write_summary(out_dir: Path | None, summary: dict[str, Any]) -> None:
    if out_dir is None:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / INGEST_SUMMARY_NAME).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _normalize_columns(frame: pd.DataFrame, mapping: dict[str, str], required: list[str], label: str) -> pd.DataFrame:
    missing_source = [source for source in mapping.values() if source not in frame.columns]
    if missing_source:
        raise ValueError(f"{label} is missing source columns: {', '.join(sorted(missing_source))}")

    normalized = pd.DataFrame({target: frame[source] for target, source in mapping.items()})
    missing_target = [target for target in required if target not in normalized.columns]
    if missing_target:
        raise ValueError(f"{label} is missing target columns: {', '.join(sorted(missing_target))}")

    normalized = normalized[required].copy()
    normalized["Filename"] = normalized["Filename"].astype(str)
    normalized["mode"] = pd.to_numeric(normalized["mode"], errors="raise").astype(int)
    normalized["internal_coordinate"] = normalized["internal_coordinate"].astype(str)
    if "contribution_percent" in normalized.columns:
        normalized["contribution_percent"] = pd.to_numeric(normalized["contribution_percent"], errors="raise")
    return normalized


def _copy_or_normalize_csv(
    source_path: Path,
    out_path: Path,
    *,
    required: list[str],
    mapping: dict[str, str] | None,
    label: str,
) -> int:
    if mapping is None:
        frame = pd.read_csv(source_path)
        missing = set(required) - set(frame.columns)
        if missing:
            raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")
        normalized = _normalize_columns(frame, {column: column for column in required}, required, label)
    else:
        frame = pd.read_csv(source_path)
        normalized = _normalize_columns(frame, mapping, required, label)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(out_path, index=False)
    return int(len(normalized))


def convert_veda_reference(
    raw_reference_dir: str | Path,
    out_dir: str | Path,
    *,
    matrix_csv: str | Path | None = None,
    dominant_csv: str | Path | None = None,
    matrix_mapping: dict[str, str] | None = None,
    dominant_mapping: dict[str, str] | None = None,
) -> dict[str, Any]:
    raw_root = Path(raw_reference_dir)
    out_root = Path(out_dir)
    summary: dict[str, Any] = {
        "raw_reference_dir": str(raw_root),
        "out_dir": str(out_root),
    }

    if not raw_root.exists():
        summary.update({"conversion_status": "SKIP", "acceptance_status": "SKIP", "reason": "raw_reference_directory_missing"})
        _write_summary(out_root, summary)
        return summary
    if not raw_root.is_dir():
        raise NotADirectoryError(f"not a directory: {raw_root}")

    source_matrix = Path(matrix_csv) if matrix_csv is not None else raw_root / REFERENCE_MATRIX_NAME
    if not source_matrix.exists():
        summary.update({"conversion_status": "SKIP", "acceptance_status": "SKIP", "reason": "missing_reference_matrix_csv"})
        _write_summary(out_root, summary)
        return summary

    try:
        matrix_rows = _copy_or_normalize_csv(
            source_matrix,
            out_root / REFERENCE_MATRIX_NAME,
            required=MATRIX_COLUMNS,
            mapping=matrix_mapping,
            label=source_matrix.name,
        )
        source_dominant = Path(dominant_csv) if dominant_csv is not None else raw_root / REFERENCE_DOMINANT_NAME
        dominant_rows = 0
        dominant_status = "SKIP"
        dominant_reason = f"missing_{REFERENCE_DOMINANT_NAME}"
        if source_dominant.exists():
            dominant_rows = _copy_or_normalize_csv(
                source_dominant,
                out_root / REFERENCE_DOMINANT_NAME,
                required=DOMINANT_COLUMNS,
                mapping=dominant_mapping,
                label=source_dominant.name,
            )
            dominant_status = "PASS"
            dominant_reason = "dominant_reference_normalized"

        summary.update(
            {
                "conversion_status": "PASS",
                "acceptance_status": "PASS",
                "reason": "reference_rows_normalized",
                "matrix_rows": matrix_rows,
                "dominant_rows": dominant_rows,
                "dominant_status": dominant_status,
                "dominant_reason": dominant_reason,
                "matrix_output": str(out_root / REFERENCE_MATRIX_NAME),
            }
        )
        if dominant_status == "PASS":
            summary["dominant_output"] = str(out_root / REFERENCE_DOMINANT_NAME)
    except (EmptyDataError, ParserError, TypeError, ValueError) as exc:
        summary.update({"conversion_status": "FAIL", "acceptance_status": "FAIL", "reason": str(exc)})
        _write_summary(out_root, summary)
        return summary

    _write_summary(out_root, summary)
    return summary


def _mapping_from_args(args: argparse.Namespace, *, dominant: bool = False) -> dict[str, str] | None:
    if dominant:
        values = [args.dominant_filename_column, args.dominant_mode_column, args.dominant_coordinate_column]
        if not any(values):
            return None
        if not all(values):
            raise ValueError("dominant column mapping requires filename, mode, and coordinate columns")
        return {
            "Filename": args.dominant_filename_column,
            "mode": args.dominant_mode_column,
            "internal_coordinate": args.dominant_coordinate_column,
        }

    values = [args.filename_column, args.mode_column, args.coordinate_column, args.percent_column]
    if not any(values):
        return None
    if not all(values):
        raise ValueError("matrix column mapping requires filename, mode, coordinate, and percent columns")
    return {
        "Filename": args.filename_column,
        "mode": args.mode_column,
        "internal_coordinate": args.coordinate_column,
        "contribution_percent": args.percent_column,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-reference", required=True, help="Directory containing raw or normalized original VEDA reference CSV files.")
    parser.add_argument("--out", required=True, help="Output directory for normalized reference CSV files.")
    parser.add_argument("--matrix-csv", help="Explicit source matrix CSV. Defaults to raw-reference/veda_reference_ped_matrix.csv.")
    parser.add_argument("--dominant-csv", help="Explicit source dominant-coordinate CSV. Defaults to raw-reference/veda_reference_dominant_assignments.csv when present.")
    parser.add_argument("--filename-column", help="Source column for Filename in matrix CSV.")
    parser.add_argument("--mode-column", help="Source column for mode in matrix CSV.")
    parser.add_argument("--coordinate-column", help="Source column for internal_coordinate in matrix CSV.")
    parser.add_argument("--percent-column", help="Source column for contribution_percent in matrix CSV.")
    parser.add_argument("--dominant-filename-column", help="Source column for Filename in dominant CSV.")
    parser.add_argument("--dominant-mode-column", help="Source column for mode in dominant CSV.")
    parser.add_argument("--dominant-coordinate-column", help="Source column for internal_coordinate in dominant CSV.")
    parser.add_argument("--fail-on-skip", action="store_true", help="Return nonzero when no convertible reference rows are found.")
    args = parser.parse_args()

    summary = convert_veda_reference(
        args.raw_reference,
        args.out,
        matrix_csv=args.matrix_csv,
        dominant_csv=args.dominant_csv,
        matrix_mapping=_mapping_from_args(args),
        dominant_mapping=_mapping_from_args(args, dominant=True),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["conversion_status"] == "FAIL":
        return 1
    if args.fail_on_skip and summary["conversion_status"] == "SKIP":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
