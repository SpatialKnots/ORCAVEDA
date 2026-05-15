from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from b_matrix import analytical_B, finite_difference_B, select_independent_coordinates, svd_rank_condition  # noqa: E402
from chemistry import annotate_chemical_system, expected_vibrational_rank  # noqa: E402
from internal_coordinates import build_internal_coordinates  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402


DEFAULT_HESS_NAMES = (
    "H2O_freq.hess",
    "NH3.hess",
    "formaldehyde.hess",
    "ethene.hess",
)


def _as_int_dict(values: object) -> dict[str, int]:
    if not isinstance(values, dict):
        return {}
    return {str(key): int(value) for key, value in values.items()}


def _build_internals_for_hess(hess_path: Path):
    hess = read_orca_hess(hess_path)
    annotation = annotate_chemical_system(hess.atoms, hess.coords_A)
    internals = build_internal_coordinates(
        hess.atoms,
        hess.coords_A,
        list(annotation.bonds),
        [list(fragment) for fragment in annotation.fragments],
        list(annotation.interfragment_hbonds),
        list(annotation.functional_groups),
    )
    return hess, annotation, internals


def compare_hess_file(
    hess_path: Path,
    *,
    row_tolerance: float = 1.0e-5,
    rank_tolerance: float = 1.0e-6,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    hess, annotation, internals = _build_internals_for_hess(hess_path)
    finite = finite_difference_B(hess.coords_A, internals)
    hybrid, diagnostics = analytical_B(hess.coords_A, internals)

    abs_delta = np.abs(hybrid - finite)
    row_max_delta = np.max(abs_delta, axis=1) if abs_delta.size else np.array([], dtype=float)
    row_methods = list(diagnostics.get("row_methods", []))
    if len(row_methods) != len(internals):
        raise ValueError(
            f"analytical_B row_methods length mismatch for {hess_path.name}: "
            f"{len(row_methods)} != {len(internals)}"
        )

    finite_rank, finite_condition, _ = svd_rank_condition(finite, tol_abs=rank_tolerance)
    hybrid_rank, hybrid_condition, _ = svd_rank_condition(hybrid, tol_abs=rank_tolerance)
    target_rank = expected_vibrational_rank(len(hess.atoms), linear=False)
    finite_selected, finite_selected_rank, finite_selected_condition, _ = select_independent_coordinates(
        finite,
        internals,
        target_rank,
        tol_abs=rank_tolerance,
    )
    hybrid_selected, hybrid_selected_rank, hybrid_selected_condition, _ = select_independent_coordinates(
        hybrid,
        internals,
        target_rank,
        tol_abs=rank_tolerance,
    )

    mismatch_indices = [int(idx) for idx, value in enumerate(row_max_delta) if float(value) > row_tolerance]
    rows: list[dict[str, object]] = []
    for idx, internal in enumerate(internals):
        rows.append(
            {
                "Filename": hess_path.name,
                "row_index": int(idx),
                "internal_coordinate": internal.name,
                "kind": internal.kind,
                "source": internal.source,
                "method": row_methods[idx],
                "max_abs_delta": float(row_max_delta[idx]),
                "above_tolerance": bool(float(row_max_delta[idx]) > row_tolerance),
            }
        )

    summary = {
        "Filename": hess_path.name,
        "natoms": int(len(hess.atoms)),
        "internal_count": int(len(internals)),
        "method_counts": _as_int_dict(diagnostics.get("method_counts")),
        "fallback_reasons": _as_int_dict(diagnostics.get("fallback_reasons")),
        "max_abs_delta": float(np.max(row_max_delta)) if row_max_delta.size else 0.0,
        "mean_row_max_abs_delta": float(np.mean(row_max_delta)) if row_max_delta.size else 0.0,
        "rows_above_tolerance": int(len(mismatch_indices)),
        "row_tolerance": float(row_tolerance),
        "redundant_finite_rank": int(finite_rank),
        "redundant_hybrid_rank": int(hybrid_rank),
        "redundant_finite_condition": float(finite_condition),
        "redundant_hybrid_condition": float(hybrid_condition),
        "target_rank": int(target_rank),
        "finite_selected_rank": int(finite_selected_rank),
        "hybrid_selected_rank": int(hybrid_selected_rank),
        "finite_selected_condition": float(finite_selected_condition),
        "hybrid_selected_condition": float(hybrid_selected_condition),
        "selected_basis_same_indices": bool(tuple(finite_selected) == tuple(hybrid_selected)),
        "finite_selected_indices": [int(idx) for idx in finite_selected],
        "hybrid_selected_indices": [int(idx) for idx in hybrid_selected],
        "formula": annotation.formula,
        "system_type": annotation.system_type,
    }
    return summary, rows


def compare_hess_files(
    hess_paths: Iterable[Path],
    *,
    row_tolerance: float = 1.0e-5,
    rank_tolerance: float = 1.0e-6,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summaries: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    for hess_path in hess_paths:
        summary, row_details = compare_hess_file(
            hess_path,
            row_tolerance=row_tolerance,
            rank_tolerance=rank_tolerance,
        )
        summaries.append(summary)
        rows.extend(row_details)
    return summaries, rows


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _resolve_hess_paths(args: argparse.Namespace) -> list[Path]:
    if args.full_sweep:
        paths = sorted((ROOT / "data" / "hess").glob("*.hess"))
    elif args.hess:
        paths = [Path(item) for item in args.hess]
    else:
        paths = [ROOT / "data" / "hess" / name for name in DEFAULT_HESS_NAMES]

    resolved = []
    for path in paths:
        candidate = path if path.is_absolute() else (ROOT / path)
        if not candidate.is_file():
            raise FileNotFoundError(f"Hess file not found: {candidate}")
        resolved.append(candidate)
    return resolved


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare finite_difference_B with hybrid analytical_B on ORCA .hess inputs. "
            "This is a diagnostic harness only; it does not switch the production pipeline."
        )
    )
    parser.add_argument("--hess", nargs="*", help="Specific .hess files to compare. Defaults to H2O, NH3, formaldehyde, ethene.")
    parser.add_argument("--full-sweep", action="store_true", help="Compare every data/hess/*.hess fixture.")
    parser.add_argument("--out", default="outputs/bmatrix_compare", help="Output directory for JSON and CSV diagnostics.")
    parser.add_argument("--row-tolerance", type=float, default=1.0e-5, help="Row max-absolute-delta tolerance.")
    parser.add_argument("--rank-tolerance", type=float, default=1.0e-6, help="SVD rank tolerance.")
    args = parser.parse_args(argv)

    hess_paths = _resolve_hess_paths(args)
    summaries, rows = compare_hess_files(
        hess_paths,
        row_tolerance=float(args.row_tolerance),
        rank_tolerance=float(args.rank_tolerance),
    )

    outdir = Path(args.out)
    if not outdir.is_absolute():
        outdir = ROOT / outdir
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {
        "comparison_scope": "finite_difference_B_vs_hybrid_analytical_B",
        "method_boundary": "diagnostic only; production pipeline remains finite_difference_B",
        "row_tolerance": float(args.row_tolerance),
        "rank_tolerance": float(args.rank_tolerance),
        "file_count": int(len(summaries)),
        "files_with_rows_above_tolerance": int(sum(1 for row in summaries if int(row["rows_above_tolerance"]) > 0)),
        "files_with_redundant_rank_change": int(
            sum(1 for row in summaries if int(row["redundant_finite_rank"]) != int(row["redundant_hybrid_rank"]))
        ),
        "files_with_selected_rank_change": int(
            sum(1 for row in summaries if int(row["finite_selected_rank"]) != int(row["hybrid_selected_rank"]))
        ),
        "files_with_selected_basis_index_change": int(
            sum(1 for row in summaries if not bool(row["selected_basis_same_indices"]))
        ),
        "summaries": summaries,
    }
    (outdir / "bmatrix_method_comparison_summary.json").write_text(
        json.dumps(_json_ready(payload), indent=2),
        encoding="utf-8",
    )
    _write_csv(outdir / "bmatrix_method_comparison_summary.csv", summaries)
    _write_csv(outdir / "bmatrix_method_comparison_rows.csv", rows)

    print(json.dumps(_json_ready(payload), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
