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


def _angle_geometry(coords_A: np.ndarray, atoms0: Sequence[int]) -> tuple[float | None, float | None]:
    if len(atoms0) != 3:
        return None, None
    i, j, k = (int(idx) for idx in atoms0)
    u = coords_A[i] - coords_A[j]
    v = coords_A[k] - coords_A[j]
    u_norm = float(np.linalg.norm(u))
    v_norm = float(np.linalg.norm(v))
    if u_norm <= 0.0 or v_norm <= 0.0:
        return None, None
    cosine = float(np.dot(u, v) / (u_norm * v_norm))
    cosine = max(-1.0, min(1.0, cosine))
    sine = float(np.sqrt(max(0.0, 1.0 - cosine * cosine)))
    angle_deg = float(np.degrees(np.arccos(cosine)))
    return angle_deg, sine


def _row_descriptor(
    *,
    hess_path: Path,
    coords_A: np.ndarray,
    internals: Sequence[object],
    row_methods: Sequence[str],
    row_max_delta: np.ndarray,
    idx: int | None,
    prefix: str,
) -> dict[str, object]:
    if idx is None:
        return {
            f"{prefix}_row_index": "",
            f"{prefix}_internal_coordinate": "",
            f"{prefix}_kind": "",
            f"{prefix}_source": "",
            f"{prefix}_atoms0": "",
            f"{prefix}_method": "",
            f"{prefix}_max_abs_delta": "",
            f"{prefix}_angle_degrees": "",
            f"{prefix}_angle_sine": "",
        }
    internal = internals[int(idx)]
    atoms0 = tuple(int(atom_idx) for atom_idx in internal.atoms0)
    angle_degrees, angle_sine = _angle_geometry(coords_A, atoms0)
    return {
        f"{prefix}_row_index": int(idx),
        f"{prefix}_internal_coordinate": internal.name,
        f"{prefix}_kind": internal.kind,
        f"{prefix}_source": internal.source,
        f"{prefix}_atoms0": json.dumps(list(atoms0), separators=(",", ":")),
        f"{prefix}_method": row_methods[int(idx)],
        f"{prefix}_max_abs_delta": float(row_max_delta[int(idx)]),
        f"{prefix}_angle_degrees": "" if angle_degrees is None else float(angle_degrees),
        f"{prefix}_angle_sine": "" if angle_sine is None else float(angle_sine),
    }


def _replacement_basis_metrics(
    B: np.ndarray,
    selected: Sequence[int],
    *,
    position: int,
    replacement_idx: int | None,
    rank_tolerance: float,
) -> tuple[int | str, float | str, float | str]:
    if replacement_idx is None or position < 0 or position >= len(selected):
        return "", "", ""
    trial = [int(idx) for idx in selected]
    trial[int(position)] = int(replacement_idx)
    rank, condition, singular_values = svd_rank_condition(B[trial, :], tol_abs=rank_tolerance)
    min_singular = float(singular_values[-1]) if singular_values.size else float("inf")
    return int(rank), float(condition), min_singular


def _selected_basis_differences(
    *,
    hess_path: Path,
    coords_A: np.ndarray,
    internals: Sequence[object],
    finite_B: np.ndarray,
    hybrid_B: np.ndarray,
    row_methods: Sequence[str],
    row_max_delta: np.ndarray,
    finite_selected: Sequence[int],
    hybrid_selected: Sequence[int],
    finite_selected_rank: int,
    hybrid_selected_rank: int,
    rank_tolerance: float,
) -> list[dict[str, object]]:
    max_len = max(len(finite_selected), len(hybrid_selected))
    rows: list[dict[str, object]] = []
    for position in range(max_len):
        finite_idx = int(finite_selected[position]) if position < len(finite_selected) else None
        hybrid_idx = int(hybrid_selected[position]) if position < len(hybrid_selected) else None
        if finite_idx == hybrid_idx:
            continue
        finite_basis_rank_with_hybrid_row, finite_basis_condition_with_hybrid_row, finite_basis_min_singular_with_hybrid_row = (
            _replacement_basis_metrics(
                finite_B,
                finite_selected,
                position=position,
                replacement_idx=hybrid_idx,
                rank_tolerance=rank_tolerance,
            )
        )
        hybrid_basis_rank_with_finite_row, hybrid_basis_condition_with_finite_row, hybrid_basis_min_singular_with_finite_row = (
            _replacement_basis_metrics(
                hybrid_B,
                hybrid_selected,
                position=position,
                replacement_idx=finite_idx,
                rank_tolerance=rank_tolerance,
            )
        )
        row: dict[str, object] = {
            "Filename": hess_path.name,
            "basis_position": int(position),
            "replacement_rank_preserved": bool(
                finite_basis_rank_with_hybrid_row == int(finite_selected_rank)
                and hybrid_basis_rank_with_finite_row == int(hybrid_selected_rank)
            ),
            "finite_basis_rank_with_hybrid_row": finite_basis_rank_with_hybrid_row,
            "finite_basis_condition_with_hybrid_row": finite_basis_condition_with_hybrid_row,
            "finite_basis_min_singular_with_hybrid_row": finite_basis_min_singular_with_hybrid_row,
            "hybrid_basis_rank_with_finite_row": hybrid_basis_rank_with_finite_row,
            "hybrid_basis_condition_with_finite_row": hybrid_basis_condition_with_finite_row,
            "hybrid_basis_min_singular_with_finite_row": hybrid_basis_min_singular_with_finite_row,
        }
        row.update(
            _row_descriptor(
                hess_path=hess_path,
                coords_A=coords_A,
                internals=internals,
                row_methods=row_methods,
                row_max_delta=row_max_delta,
                idx=finite_idx,
                prefix="finite",
            )
        )
        row.update(
            _row_descriptor(
                hess_path=hess_path,
                coords_A=coords_A,
                internals=internals,
                row_methods=row_methods,
                row_max_delta=row_max_delta,
                idx=hybrid_idx,
                prefix="hybrid",
            )
        )
        rows.append(row)
    return rows


def compare_hess_file(
    hess_path: Path,
    *,
    row_tolerance: float = 1.0e-5,
    rank_tolerance: float = 1.0e-6,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
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
        atoms0 = tuple(int(atom_idx) for atom_idx in internal.atoms0)
        angle_degrees, angle_sine = _angle_geometry(hess.coords_A, atoms0)
        rows.append(
            {
                "Filename": hess_path.name,
                "row_index": int(idx),
                "internal_coordinate": internal.name,
                "kind": internal.kind,
                "source": internal.source,
                "atoms0": json.dumps(list(atoms0), separators=(",", ":")),
                "method": row_methods[idx],
                "max_abs_delta": float(row_max_delta[idx]),
                "above_tolerance": bool(float(row_max_delta[idx]) > row_tolerance),
                "angle_degrees": "" if angle_degrees is None else float(angle_degrees),
                "angle_sine": "" if angle_sine is None else float(angle_sine),
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
    selection_differences = _selected_basis_differences(
        hess_path=hess_path,
        coords_A=hess.coords_A,
        internals=internals,
        finite_B=finite,
        hybrid_B=hybrid,
        row_methods=row_methods,
        row_max_delta=row_max_delta,
        finite_selected=finite_selected,
        hybrid_selected=hybrid_selected,
        finite_selected_rank=finite_selected_rank,
        hybrid_selected_rank=hybrid_selected_rank,
        rank_tolerance=rank_tolerance,
    )
    summary["selected_basis_difference_count"] = int(len(selection_differences))
    return summary, rows, selection_differences


def compare_hess_files(
    hess_paths: Iterable[Path],
    *,
    row_tolerance: float = 1.0e-5,
    rank_tolerance: float = 1.0e-6,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    summaries: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    selection_differences: list[dict[str, object]] = []
    for hess_path in hess_paths:
        summary, row_details, selected_details = compare_hess_file(
            hess_path,
            row_tolerance=row_tolerance,
            rank_tolerance=rank_tolerance,
        )
        summaries.append(summary)
        rows.extend(row_details)
        selection_differences.extend(selected_details)
    return summaries, rows, selection_differences


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
    summaries, rows, selection_differences = compare_hess_files(
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
        "rows_above_tolerance_count": int(sum(int(row["rows_above_tolerance"]) for row in summaries)),
        "selected_basis_difference_count": int(len(selection_differences)),
        "selected_basis_replacement_rank_loss_count": int(
            sum(1 for row in selection_differences if not bool(row.get("replacement_rank_preserved")))
        ),
        "summaries": summaries,
    }
    (outdir / "bmatrix_method_comparison_summary.json").write_text(
        json.dumps(_json_ready(payload), indent=2),
        encoding="utf-8",
    )
    _write_csv(outdir / "bmatrix_method_comparison_summary.csv", summaries)
    _write_csv(outdir / "bmatrix_method_comparison_rows.csv", rows)
    _write_csv(outdir / "bmatrix_method_comparison_selected_basis_differences.csv", selection_differences)

    print(json.dumps(_json_ready(payload), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
