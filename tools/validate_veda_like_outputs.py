from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd


WARNING_TOKENS = (
    "fixed_conversion_failed",
    "linear_or_near_linear_fixed_conversion_review",
    "high_frequency_hbond_dominates_xh_stretch_secondary",
    "missing_cartesian_hessian",
)


def _find_one(directory: Path, suffix: str) -> Path:
    matches = sorted(directory.glob(f"*__{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(f"expected exactly one '*__{suffix}' in {directory}, found {len(matches)}")
    return matches[0]


def _read_csv(directory: Path, suffix: str) -> pd.DataFrame:
    return pd.read_csv(_find_one(directory, suffix))


def _status_counts(frame: pd.DataFrame, column: str = "validation_status") -> dict[str, int]:
    if column not in frame.columns:
        return {}
    counts = Counter(str(value) for value in frame[column].fillna(""))
    counts.pop("", None)
    return dict(sorted(counts.items()))


def _count_warning_tokens(*frames: pd.DataFrame) -> dict[str, int]:
    counts = Counter()
    for frame in frames:
        if "warnings" not in frame.columns:
            continue
        for value in frame["warnings"].fillna(""):
            text = str(value)
            for token in WARNING_TOKENS:
                if token in text:
                    counts[token] += 1
    return {token: int(counts[token]) for token in WARNING_TOKENS if counts[token]}


def _normalization_failures(matrix: pd.DataFrame, tolerance: float) -> list[dict[str, Any]]:
    required = {"Filename", "mode", "contribution_percent"}
    missing = required - set(matrix.columns)
    if missing:
        raise ValueError(f"veda_like_ped_matrix is missing columns: {', '.join(sorted(missing))}")

    work = matrix.copy()
    work["contribution_percent"] = pd.to_numeric(work["contribution_percent"], errors="coerce")
    grouped = work.groupby(["Filename", "mode"], dropna=False)["contribution_percent"].sum(min_count=1).reset_index()
    failures: list[dict[str, Any]] = []
    for row in grouped.to_dict("records"):
        total = row["contribution_percent"]
        if pd.isna(total) or abs(float(total) - 100.0) > tolerance:
            failures.append(
                {
                    "Filename": str(row["Filename"]),
                    "mode": str(row["mode"]),
                    "contribution_percent_sum": None if pd.isna(total) else float(total),
                }
            )
    return failures


def validate_veda_like_outputs(directory: str | Path, *, tolerance: float = 1.0e-6) -> dict[str, Any]:
    root = Path(directory)
    if not root.exists():
        raise FileNotFoundError(f"output directory does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"not a directory: {root}")

    basis = _read_csv(root, "veda_like_basis_diagnostics.csv")
    correspondence = _read_csv(root, "veda_like_mode_correspondence.csv")
    matrix = _read_csv(root, "veda_like_ped_matrix.csv")
    audit = _read_csv(root, "veda_like_ped_audit.csv")

    normalization = _normalization_failures(matrix, tolerance)
    files = sorted(str(value) for value in basis.get("Filename", pd.Series(dtype=str)).dropna().unique())
    summary = {
        "directory": str(root),
        "file_count": len(files),
        "files": files,
        "basis_status_counts": _status_counts(basis),
        "mode_correspondence_status_counts": _status_counts(correspondence),
        "ped_audit_status_counts": _status_counts(audit),
        "normalization_tolerance": tolerance,
        "normalization_failure_count": len(normalization),
        "normalization_failures": normalization,
        "warning_token_counts": _count_warning_tokens(basis, correspondence, matrix, audit),
        "artifact_rows": {
            "basis_diagnostics": int(len(basis)),
            "mode_correspondence": int(len(correspondence)),
            "ped_matrix": int(len(matrix)),
            "ped_audit": int(len(audit)),
        },
    }
    status_values = set(summary["basis_status_counts"]) | set(summary["mode_correspondence_status_counts"])
    if normalization or "FAIL" in status_values:
        summary["validation_status"] = "FAIL"
    elif "WARN" in status_values or summary["warning_token_counts"]:
        summary["validation_status"] = "WARN"
    else:
        summary["validation_status"] = "PASS"
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate ORCAVEDA opt-in veda_like_* output artifacts.")
    parser.add_argument("directory", help="Directory containing *__veda_like_*.csv artifacts.")
    parser.add_argument("--tolerance", type=float, default=1.0e-6, help="Per-mode PED percent sum tolerance.")
    parser.add_argument("--json-out", help="Optional path for summary JSON.")
    parser.add_argument("--csv-out", help="Optional path for per-mode normalization failures CSV.")
    args = parser.parse_args()

    summary = validate_veda_like_outputs(args.directory, tolerance=args.tolerance)
    text = json.dumps(summary, indent=2, sort_keys=True)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n", encoding="utf-8")
    if args.csv_out:
        pd.DataFrame(summary["normalization_failures"]).to_csv(args.csv_out, index=False)
    return 0 if summary["validation_status"] != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
