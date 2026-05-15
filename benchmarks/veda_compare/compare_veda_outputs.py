#!/usr/bin/env python3
"""Compare ORCAVEDA VEDA-like artifacts with checked-in VEDA reference rows.

This harness is deliberately skip-safe. Missing original VEDA reference outputs
produce a SKIP summary, not a PASS. The supported reference schema is documented
in README.md next to this script.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


MATRIX_SUFFIX = "veda_like_ped_matrix.csv"
AUDIT_SUFFIX = "veda_like_ped_audit.csv"
REFERENCE_MATRIX_NAME = "veda_reference_ped_matrix.csv"
REFERENCE_DOMINANT_NAME = "veda_reference_dominant_assignments.csv"


def _find_one(directory: Path, suffix: str) -> Path | None:
    matches = sorted(directory.glob(f"*__{suffix}"))
    if len(matches) == 1:
        return matches[0]
    return None


def _read_required_csv(path: Path, required: set[str], label: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")
    return frame


def _normalize_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy()
    normalized["Filename"] = normalized["Filename"].astype(str)
    normalized["mode"] = pd.to_numeric(normalized["mode"], errors="coerce").astype("Int64")
    normalized["internal_coordinate"] = normalized["internal_coordinate"].astype(str)
    normalized["contribution_percent"] = pd.to_numeric(normalized["contribution_percent"], errors="coerce")
    return normalized


def _matrix_comparison(
    orcaveda_matrix: pd.DataFrame,
    reference_matrix: pd.DataFrame,
    *,
    tolerance_percent: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    left = _normalize_matrix(orcaveda_matrix)
    right = _normalize_matrix(reference_matrix)
    keys = ["Filename", "mode", "internal_coordinate"]
    merged = right.merge(
        left[keys + ["contribution_percent"]],
        on=keys,
        how="outer",
        suffixes=("_reference", "_orcaveda"),
        indicator=True,
    )
    merged["delta_percent"] = merged["contribution_percent_orcaveda"] - merged["contribution_percent_reference"]
    merged["abs_delta_percent"] = merged["delta_percent"].abs()
    comparable = merged[merged["_merge"] == "both"].copy()
    out_of_tolerance = comparable[comparable["abs_delta_percent"] > float(tolerance_percent)]
    summary = {
        "reference_matrix_rows": int(len(right)),
        "orcaveda_matrix_rows": int(len(left)),
        "matched_matrix_rows": int(len(comparable)),
        "missing_reference_row_count": int((merged["_merge"] == "left_only").sum()),
        "extra_orcaveda_row_count": int((merged["_merge"] == "right_only").sum()),
        "out_of_tolerance_count": int(len(out_of_tolerance)),
        "max_abs_delta_percent": 0.0 if comparable.empty else float(comparable["abs_delta_percent"].max()),
    }
    return merged, summary


def _dominant_comparison(
    orcaveda_audit: pd.DataFrame,
    reference_dominant: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {"Filename", "mode", "internal_coordinate"}
    missing = required - set(reference_dominant.columns)
    if missing:
        raise ValueError(f"{REFERENCE_DOMINANT_NAME} is missing columns: {', '.join(sorted(missing))}")

    audit = orcaveda_audit.copy()
    audit["veda_like_rank"] = pd.to_numeric(audit["veda_like_rank"], errors="coerce")
    audit = audit[audit["veda_like_rank"] == 1].copy()
    audit["Filename"] = audit["Filename"].astype(str)
    audit["mode"] = pd.to_numeric(audit["mode"], errors="coerce").astype("Int64")
    audit["internal_coordinate"] = audit["internal_coordinate"].astype(str)

    reference = reference_dominant.copy()
    reference["Filename"] = reference["Filename"].astype(str)
    reference["mode"] = pd.to_numeric(reference["mode"], errors="coerce").astype("Int64")
    reference["internal_coordinate"] = reference["internal_coordinate"].astype(str)

    keys = ["Filename", "mode"]
    merged = reference.merge(
        audit[keys + ["internal_coordinate", "coordinate_family", "contribution_percent"]],
        on=keys,
        how="left",
        suffixes=("_reference", "_orcaveda"),
    )
    merged["dominant_match"] = merged["internal_coordinate_reference"] == merged["internal_coordinate_orcaveda"]
    summary = {
        "reference_dominant_rows": int(len(reference)),
        "dominant_compared_rows": int(len(merged)),
        "dominant_match_count": int(merged["dominant_match"].sum()),
        "dominant_mismatch_count": int((~merged["dominant_match"]).sum()),
    }
    return merged, summary


def compare_veda_outputs(
    orcaveda_dir: str | Path,
    reference_dir: str | Path,
    out_dir: str | Path | None = None,
    *,
    tolerance_percent: float = 5.0,
) -> dict[str, Any]:
    orcaveda_root = Path(orcaveda_dir)
    reference_root = Path(reference_dir)
    out_root = Path(out_dir) if out_dir is not None else None
    summary: dict[str, Any] = {
        "orcaveda_dir": str(orcaveda_root),
        "reference_dir": str(reference_root),
        "tolerance_percent": float(tolerance_percent),
    }

    if not reference_root.exists():
        summary.update(
            {"comparison_status": "SKIP", "acceptance_status": "SKIP", "reason": "veda_reference_directory_missing"}
        )
        _write_outputs(out_root, summary)
        return summary
    if not reference_root.is_dir():
        raise NotADirectoryError(f"not a directory: {reference_root}")

    reference_matrix_path = reference_root / REFERENCE_MATRIX_NAME
    if not reference_matrix_path.exists():
        summary.update({"comparison_status": "SKIP", "acceptance_status": "SKIP", "reason": f"missing_{REFERENCE_MATRIX_NAME}"})
        _write_outputs(out_root, summary)
        return summary

    orcaveda_matrix_path = _find_one(orcaveda_root, MATRIX_SUFFIX)
    orcaveda_audit_path = _find_one(orcaveda_root, AUDIT_SUFFIX)
    if orcaveda_matrix_path is None or orcaveda_audit_path is None:
        summary.update(
            {"comparison_status": "FAIL", "acceptance_status": "FAIL", "reason": "missing_orcaveda_veda_like_artifacts"}
        )
        _write_outputs(out_root, summary)
        return summary

    reference_matrix = _read_required_csv(
        reference_matrix_path,
        {"Filename", "mode", "internal_coordinate", "contribution_percent"},
        REFERENCE_MATRIX_NAME,
    )
    orcaveda_matrix = _read_required_csv(
        orcaveda_matrix_path,
        {"Filename", "mode", "internal_coordinate", "contribution_percent"},
        orcaveda_matrix_path.name,
    )
    matrix_delta, matrix_summary = _matrix_comparison(
        orcaveda_matrix,
        reference_matrix,
        tolerance_percent=tolerance_percent,
    )
    summary.update(matrix_summary)

    dominant_path = reference_root / REFERENCE_DOMINANT_NAME
    dominant_delta: pd.DataFrame | None = None
    if dominant_path.exists():
        reference_dominant = pd.read_csv(dominant_path)
        orcaveda_audit = _read_required_csv(
            orcaveda_audit_path,
            {"Filename", "mode", "veda_like_rank", "internal_coordinate", "coordinate_family", "contribution_percent"},
            orcaveda_audit_path.name,
        )
        dominant_delta, dominant_summary = _dominant_comparison(orcaveda_audit, reference_dominant)
        summary.update(dominant_summary)
    else:
        summary.update(
            {
                "reference_dominant_rows": 0,
                "dominant_compared_rows": 0,
                "dominant_match_count": 0,
                "dominant_mismatch_count": 0,
                "dominant_status": "SKIP",
                "dominant_reason": f"missing_{REFERENCE_DOMINANT_NAME}",
            }
        )

    fail_reasons = []
    if summary["missing_reference_row_count"]:
        fail_reasons.append("reference_rows_missing_in_orcaveda")
    if summary["out_of_tolerance_count"]:
        fail_reasons.append("ped_percent_delta_out_of_tolerance")
    if summary.get("dominant_mismatch_count", 0):
        fail_reasons.append("dominant_coordinate_mismatch")

    if fail_reasons:
        summary["comparison_status"] = "FAIL"
        summary["acceptance_status"] = "FAIL"
        summary["reason"] = ";".join(fail_reasons)
    else:
        summary["comparison_status"] = "PASS"
        summary["acceptance_status"] = "PASS"
        summary["reason"] = "reference_comparison_within_declared_tolerance"

    _write_outputs(out_root, summary, matrix_delta=matrix_delta, dominant_delta=dominant_delta)
    return summary


def _write_outputs(
    out_dir: Path | None,
    summary: dict[str, Any],
    *,
    matrix_delta: pd.DataFrame | None = None,
    dominant_delta: pd.DataFrame | None = None,
) -> None:
    if out_dir is None:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "veda_reference_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if matrix_delta is not None:
        matrix_delta.to_csv(out_dir / "veda_reference_matrix_delta.csv", index=False)
    if dominant_delta is not None:
        dominant_delta.to_csv(out_dir / "veda_reference_dominant_delta.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orcaveda", required=True, help="Directory containing ORCAVEDA *__veda_like_*.csv artifacts.")
    parser.add_argument("--reference", required=True, help="Directory containing original VEDA reference CSV files.")
    parser.add_argument("--out", required=True, help="Directory for comparison summary and delta CSV files.")
    parser.add_argument("--tolerance-percent", type=float, default=5.0, help="Allowed absolute PED percent delta.")
    parser.add_argument("--fail-on-skip", action="store_true", help="Return nonzero when reference outputs are absent.")
    args = parser.parse_args()

    summary = compare_veda_outputs(
        args.orcaveda,
        args.reference,
        args.out,
        tolerance_percent=args.tolerance_percent,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["comparison_status"] == "FAIL":
        return 1
    if args.fail_on_skip and summary["comparison_status"] == "SKIP":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
