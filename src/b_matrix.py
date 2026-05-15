from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

import numpy as np

from orcaveda_models import ComposedCoordinateTerm, InternalCoordinate


EPS_FD_A = 1.0e-4


def finite_difference_B(coords_A: np.ndarray, internals: Sequence[InternalCoordinate], eps: float = EPS_FD_A) -> np.ndarray:
    n = coords_A.size
    B = np.zeros((len(internals), n), dtype=float)
    flat0 = coords_A.reshape(-1)
    for r, internal in enumerate(internals):
        for c in range(n):
            plus = flat0.copy()
            minus = flat0.copy()
            plus[c] += eps
            minus[c] -= eps
            xp = plus.reshape(coords_A.shape)
            xm = minus.reshape(coords_A.shape)
            vp = internal.fn(xp)
            vm = internal.fn(xm)
            if internal.kind == "torsion":
                dv = (vp - vm + math.pi) % (2 * math.pi) - math.pi
                B[r, c] = dv / (2 * eps)
            else:
                B[r, c] = (vp - vm) / (2 * eps)
    B[~np.isfinite(B)] = 0.0
    return B


def analytical_B(
    coords_A: np.ndarray,
    internals: Sequence[InternalCoordinate],
    *,
    eps: float = EPS_FD_A,
    singular_tol: float = 1.0e-12,
    angle_sin_tol: float = 1.0e-3,
) -> tuple[np.ndarray, Dict[str, object]]:
    """
    Build a hybrid analytical B matrix with finite-difference fallback.

    The first GAP 2 implementation is deliberately narrow and additive:
    distance-like two-atom coordinates and regular angle/bend three-atom
    coordinates are analytical; torsions, composed coordinates, linear-bend
    components, and singular or near-linear angle geometries fall back to the
    existing finite difference row. This function does not replace
    `finite_difference_B` in the production pipeline.
    """
    coords = np.asarray(coords_A, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coords_A must have shape (natoms, 3), got {coords.shape}")
    if not np.all(np.isfinite(coords)):
        raise ValueError("coords_A must contain only finite values")

    B = np.zeros((len(internals), coords.size), dtype=float)
    method_counts: Dict[str, int] = {}
    fallback_reasons: Dict[str, int] = {}
    row_methods: List[str] = []
    for row_idx, internal in enumerate(internals):
        row, method, reason = _analytical_internal_coordinate_row(
            coords,
            internal,
            singular_tol=singular_tol,
            angle_sin_tol=angle_sin_tol,
        )
        if row is None:
            row = finite_difference_B(coords, [internal], eps=eps)[0]
            method = "finite_difference_fallback"
            fallback_reasons[reason] = fallback_reasons.get(reason, 0) + 1
        B[row_idx, :] = row
        row_methods.append(method)
        method_counts[method] = method_counts.get(method, 0) + 1
    B[~np.isfinite(B)] = 0.0
    diagnostics: Dict[str, object] = {
        "row_count": int(len(internals)),
        "method_counts": method_counts,
        "fallback_reasons": fallback_reasons,
        "row_methods": row_methods,
    }
    return B, diagnostics


def _analytical_internal_coordinate_row(
    coords_A: np.ndarray,
    internal: InternalCoordinate,
    *,
    singular_tol: float,
    angle_sin_tol: float,
) -> tuple[np.ndarray | None, str, str]:
    if internal.source == "composed_coordinate":
        return None, "", "composed_coordinate"

    atoms = tuple(int(idx) for idx in internal.atoms0)
    if any(idx < 0 or idx >= coords_A.shape[0] for idx in atoms):
        raise ValueError("internal coordinate atom index is out of range for coords_A")

    kind = str(internal.kind).lower()
    if len(atoms) == 2 and _is_distance_like_kind(kind):
        row = _distance_analytical_row(coords_A, atoms[0], atoms[1], singular_tol=singular_tol)
        if row is None:
            return None, "", "zero_length_distance"
        return row, "analytical_distance", ""
    if len(atoms) == 3 and _is_regular_angle_kind(kind):
        row = _angle_analytical_row(
            coords_A,
            atoms[0],
            atoms[1],
            atoms[2],
            singular_tol=singular_tol,
            angle_sin_tol=angle_sin_tol,
        )
        if row is None:
            return None, "", "singular_or_near_linear_angle"
        return row, "analytical_angle", ""
    return None, "", "unsupported_coordinate_kind"


def _is_distance_like_kind(kind: str) -> bool:
    if "torsion" in kind or "angle" in kind or "bend" in kind:
        return False
    return (
        "stretch" in kind
        or kind in {"bond", "distance", "interfragment_distance"}
        or kind.endswith("_ha")
        or kind.endswith("_da")
    )


def _is_regular_angle_kind(kind: str) -> bool:
    if kind == "linear_bend_component":
        return False
    return kind == "bend" or "bend" in kind or "angle" in kind


def _distance_analytical_row(
    coords_A: np.ndarray,
    i: int,
    j: int,
    *,
    singular_tol: float,
) -> np.ndarray | None:
    delta = coords_A[i] - coords_A[j]
    distance = float(np.linalg.norm(delta))
    if distance <= singular_tol or not np.isfinite(distance):
        return None
    direction = delta / distance
    row = np.zeros(coords_A.size, dtype=float)
    row.reshape(coords_A.shape)[i, :] = direction
    row.reshape(coords_A.shape)[j, :] = -direction
    return row


def _angle_analytical_row(
    coords_A: np.ndarray,
    i: int,
    j: int,
    k: int,
    *,
    singular_tol: float,
    angle_sin_tol: float,
) -> np.ndarray | None:
    u = coords_A[i] - coords_A[j]
    v = coords_A[k] - coords_A[j]
    u_norm = float(np.linalg.norm(u))
    v_norm = float(np.linalg.norm(v))
    if u_norm <= singular_tol or v_norm <= singular_tol:
        return None
    cosine = float(np.dot(u, v) / (u_norm * v_norm))
    cosine = max(-1.0, min(1.0, cosine))
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    if sine <= max(singular_tol, angle_sin_tol):
        return None

    dtheta_du = -(
        v / (u_norm * v_norm)
        - cosine * u / (u_norm * u_norm)
    ) / sine
    dtheta_dv = -(
        u / (u_norm * v_norm)
        - cosine * v / (v_norm * v_norm)
    ) / sine
    scale = 180.0 / math.pi
    row3 = np.zeros_like(coords_A, dtype=float)
    row3[i, :] = dtheta_du * scale
    row3[k, :] = dtheta_dv * scale
    row3[j, :] = -(row3[i, :] + row3[k, :])
    return row3.reshape(-1)


def compose_b_row(B: np.ndarray, components: Sequence[ComposedCoordinateTerm | tuple[int, float]]) -> np.ndarray:
    """
    Build a composed-coordinate B row from already-computed primitive rows.

    The composed row is the coefficient-weighted sum of primitive rows. This is
    an additive PED-basis helper only; it does not change finite-difference
    evaluation, rank selection, or Stage 3D assignment behavior.
    """
    B_arr = np.asarray(B, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if not components:
        raise ValueError("components must contain at least one primitive coordinate")

    row = np.zeros(B_arr.shape[1], dtype=float)
    for component in components:
        if isinstance(component, ComposedCoordinateTerm):
            idx = int(component.coordinate_index)
            coefficient = float(component.coefficient)
        else:
            idx = int(component[0])
            coefficient = float(component[1])
        if idx < 0 or idx >= B_arr.shape[0]:
            raise ValueError("component coordinate_index is out of range for B")
        if not np.isfinite(coefficient):
            raise ValueError("component coefficient must be finite")
        row += coefficient * B_arr[idx, :]
    row[~np.isfinite(row)] = 0.0
    return row


def build_composed_candidate_b_matrix(
    B_primitive: np.ndarray,
    primitive_internals: Sequence[InternalCoordinate],
    composed_internals: Sequence[InternalCoordinate],
) -> tuple[List[InternalCoordinate], np.ndarray, Dict[str, object]]:
    """
    Append composed-coordinate rows to a primitive B matrix for PED candidates.

    This helper is intentionally side-effect free: primitive internals remain at
    the start of the returned list, composed rows are appended, and no rank
    selection or pipeline output behavior changes here.
    """
    B_arr = np.asarray(B_primitive, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B_primitive must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(primitive_internals):
        raise ValueError("B_primitive row count must match primitive_internals")

    composed_rows: List[np.ndarray] = []
    generation_rule_counts: Dict[str, int] = {}
    for composed in composed_internals:
        if composed.source != "composed_coordinate":
            raise ValueError("composed_internals must have source='composed_coordinate'")
        row = compose_b_row(B_arr, composed.composition)
        composed_rows.append(row)
        rule = composed.generation_rule or "unknown"
        generation_rule_counts[rule] = generation_rule_counts.get(rule, 0) + 1

    if composed_rows:
        B_candidates = np.vstack([B_arr, np.vstack(composed_rows)])
    else:
        B_candidates = B_arr.copy()

    candidate_internals = list(primitive_internals) + list(composed_internals)
    diagnostics: Dict[str, object] = {
        "primitive_count": int(len(primitive_internals)),
        "composed_count": int(len(composed_internals)),
        "candidate_count": int(len(candidate_internals)),
        "generation_rule_counts": generation_rule_counts,
    }
    return candidate_internals, B_candidates, diagnostics


def svd_rank_condition(B: np.ndarray, tol_abs: float = 1.0e-6) -> Tuple[int, float, np.ndarray]:
    if B.size == 0:
        return 0, float("inf"), np.array([])
    s = np.linalg.svd(B, compute_uv=False)
    rank = int(np.sum(s > tol_abs))
    if rank == 0:
        cond = float("inf")
    else:
        s_nonzero = s[s > tol_abs]
        cond = float(s_nonzero[0] / s_nonzero[-1])
    return rank, cond, s


def select_independent_coordinates(
    B: np.ndarray,
    internals: Sequence[InternalCoordinate],
    target_rank: int,
    tol_abs: float = 1.0e-6,
):
    selected_idx: List[int] = []
    current = np.zeros((0, B.shape[1]))
    current_rank = 0
    ordered = sorted(range(len(internals)), key=lambda i: (internals[i].priority, internals[i].name))
    for idx in ordered:
        candidate = np.vstack([current, B[idx:idx + 1, :]])
        rank, _, _ = svd_rank_condition(candidate, tol_abs=tol_abs)
        if rank > current_rank:
            selected_idx.append(idx)
            current = candidate
            current_rank = rank
        if current_rank >= target_rank:
            break

    redundant_rank, _, _ = svd_rank_condition(B, tol_abs=tol_abs)
    recoverable_rank = min(int(target_rank), int(redundant_rank))
    if current_rank < recoverable_rank:
        selected_idx = _select_independent_coordinates_pivoted_cholesky(
            B,
            internals,
            recoverable_rank,
            tol_abs=tol_abs,
        )
        current = B[selected_idx, :] if selected_idx else np.zeros((0, B.shape[1]))

    rank, cond, s = svd_rank_condition(current, tol_abs=tol_abs)
    return selected_idx, rank, cond, s


def ped_basis_localization_metrics(
    B: np.ndarray,
    normal_modes: np.ndarray,
    basis_idx: Sequence[int],
    mode_indices: Sequence[int] | None = None,
    *,
    tol: float = 1.0e-12,
) -> Dict[str, float]:
    """
    Score how localized a nonredundant internal-coordinate basis is for PED.

    This is an EPM-like diagnostic inspired by VEDA's internal-coordinate
    optimization goal: maximize dominant PED matrix elements. It operates on
    the current ORCAVEDA normalized B-matrix projection model and is not a claim
    of VEDA equivalence.
    """
    B_arr = np.asarray(B, dtype=float)
    modes = np.asarray(normal_modes, dtype=float)
    idx = [int(i) for i in basis_idx]
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if modes.ndim != 2:
        raise ValueError(f"normal_modes must be a 2D matrix, got shape {modes.shape}")
    if B_arr.shape[1] != modes.shape[0]:
        raise ValueError(f"B column count {B_arr.shape[1]} does not match normal-mode row count {modes.shape[0]}")
    if any(i < 0 or i >= B_arr.shape[0] for i in idx):
        raise ValueError("basis_idx contains an out-of-range internal coordinate index")
    if not idx:
        return {
            "mode_count": 0.0,
            "mean_top_percent": 0.0,
            "median_top_percent": 0.0,
            "min_top_percent": 0.0,
            "diffuse_mode_fraction": 1.0,
            "localization_score": 0.0,
        }

    if mode_indices is None:
        mode_idx = [i for i in range(modes.shape[1]) if np.linalg.norm(modes[:, i]) > tol]
    else:
        mode_idx = [int(i) for i in mode_indices]
    mode_idx = [i for i in mode_idx if 0 <= i < modes.shape[1] and np.linalg.norm(modes[:, i]) > tol]
    if not mode_idx:
        return {
            "mode_count": 0.0,
            "mean_top_percent": 0.0,
            "median_top_percent": 0.0,
            "min_top_percent": 0.0,
            "diffuse_mode_fraction": 1.0,
            "localization_score": 0.0,
        }

    basis_B = B_arr[idx, :]
    row_norms = np.linalg.norm(basis_B, axis=1)
    valid_rows = np.isfinite(row_norms) & (row_norms > tol)
    B_unit = np.zeros_like(basis_B, dtype=float)
    B_unit[valid_rows, :] = basis_B[valid_rows, :] / row_norms[valid_rows, None]

    mode_matrix = modes[:, mode_idx].astype(float, copy=True)
    mode_norms = np.linalg.norm(mode_matrix, axis=0)
    valid_modes = np.isfinite(mode_norms) & (mode_norms > tol)
    if not np.any(valid_modes):
        return {
            "mode_count": 0.0,
            "mean_top_percent": 0.0,
            "median_top_percent": 0.0,
            "min_top_percent": 0.0,
            "diffuse_mode_fraction": 1.0,
            "localization_score": 0.0,
        }
    mode_unit = mode_matrix[:, valid_modes] / mode_norms[valid_modes][None, :]
    weights = (B_unit @ mode_unit) ** 2
    weights[~np.isfinite(weights)] = 0.0
    totals = np.sum(weights, axis=0)
    good = np.isfinite(totals) & (totals > tol)
    if not np.any(good):
        return {
            "mode_count": 0.0,
            "mean_top_percent": 0.0,
            "median_top_percent": 0.0,
            "min_top_percent": 0.0,
            "diffuse_mode_fraction": 1.0,
            "localization_score": 0.0,
        }

    pct = 100.0 * weights[:, good] / totals[good][None, :]
    top = np.max(pct, axis=0)
    mean_top = float(np.mean(top))
    median_top = float(np.median(top))
    min_top = float(np.min(top))
    diffuse_fraction = float(np.mean(top < 25.0))
    localization_score = mean_top + 0.25 * median_top - 10.0 * diffuse_fraction
    return {
        "mode_count": float(top.size),
        "mean_top_percent": mean_top,
        "median_top_percent": median_top,
        "min_top_percent": min_top,
        "diffuse_mode_fraction": diffuse_fraction,
        "localization_score": float(localization_score),
    }


def optimize_independent_coordinates_for_ped(
    B: np.ndarray,
    internals: Sequence[InternalCoordinate],
    selected_idx: Sequence[int],
    normal_modes: np.ndarray,
    mode_indices: Sequence[int] | None = None,
    *,
    target_rank: int | None = None,
    tol_abs: float = 1.0e-6,
    max_passes: int = 2,
    improvement_tol: float = 1.0e-6,
):
    """
    Greedily swap independent coordinates to improve PED localization.

    The returned basis always preserves the recoverable rank of the starting
    basis. This is a conservative EPM-like optimization layer for PED/Wilson
    diagnostics; it does not change Stage 3D assignment-audit labels.
    """
    B_arr = np.asarray(B, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    selected = [int(i) for i in selected_idx]
    if any(i < 0 or i >= len(internals) for i in selected):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")
    if not selected:
        metrics = ped_basis_localization_metrics(B_arr, normal_modes, [], mode_indices)
        return [], {"changed": False, "swaps": 0, **{f"initial_{k}": v for k, v in metrics.items()}, **{f"optimized_{k}": v for k, v in metrics.items()}}

    start_rank, start_cond, _ = svd_rank_condition(B_arr[selected, :], tol_abs=tol_abs)
    required_rank = int(target_rank) if target_rank is not None else start_rank
    required_rank = min(required_rank, start_rank)
    initial_metrics = ped_basis_localization_metrics(B_arr, normal_modes, selected, mode_indices)
    best_score = float(initial_metrics["localization_score"])
    swaps = 0

    candidate_pool = [
        idx for idx in range(len(internals))
        if idx not in set(selected) and np.linalg.norm(B_arr[idx, :]) > tol_abs
    ]
    candidate_pool = sorted(candidate_pool, key=lambda i: (internals[i].priority, internals[i].name))

    for _pass in range(max(0, int(max_passes))):
        improved = False
        for pos, old_idx in enumerate(list(selected)):
            best_replacement: tuple[float, int, Dict[str, float]] | None = None
            for candidate_idx in candidate_pool:
                if candidate_idx in selected:
                    continue
                trial = list(selected)
                trial[pos] = candidate_idx
                if len(set(trial)) != len(trial):
                    continue
                rank, cond, _ = svd_rank_condition(B_arr[trial, :], tol_abs=tol_abs)
                if rank < required_rank:
                    continue
                metrics = ped_basis_localization_metrics(B_arr, normal_modes, trial, mode_indices)
                score = float(metrics["localization_score"])
                if score <= best_score + improvement_tol:
                    continue
                tie_break = (internals[candidate_idx].priority, internals[candidate_idx].name)
                if best_replacement is None or score > best_replacement[0] + improvement_tol:
                    best_replacement = (score, candidate_idx, metrics)
                elif best_replacement is not None and abs(score - best_replacement[0]) <= improvement_tol:
                    current_tie = (internals[best_replacement[1]].priority, internals[best_replacement[1]].name)
                    if tie_break < current_tie:
                        best_replacement = (score, candidate_idx, metrics)
            if best_replacement is None:
                continue
            selected[pos] = best_replacement[1]
            best_score = best_replacement[0]
            swaps += 1
            improved = True
        if not improved:
            break

    optimized_rank, optimized_cond, _ = svd_rank_condition(B_arr[selected, :], tol_abs=tol_abs)
    optimized_metrics = ped_basis_localization_metrics(B_arr, normal_modes, selected, mode_indices)
    report = {
        "changed": bool(swaps > 0),
        "swaps": int(swaps),
        "rank": int(optimized_rank),
        "condition": float(optimized_cond),
        "initial_rank": int(start_rank),
        "initial_condition": float(start_cond),
    }
    report.update({f"initial_{key}": value for key, value in initial_metrics.items()})
    report.update({f"optimized_{key}": value for key, value in optimized_metrics.items()})
    return selected, report


def select_rank_preserving_composed_ped_basis(
    B_candidates: np.ndarray,
    candidate_internals: Sequence[InternalCoordinate],
    starting_idx: Sequence[int],
    normal_modes: np.ndarray,
    mode_indices: Sequence[int] | None = None,
    *,
    tol_abs: float = 1.0e-6,
    max_passes: int = 2,
    improvement_tol: float = 1.0e-6,
):
    """
    Optimize a PED candidate basis while preserving the starting basis rank.

    The input `starting_idx` should normally be the already accepted primitive
    independent-coordinate basis. Composed candidates may replace primitive
    rows only when the original rank is preserved. This helper is not connected
    to Stage 3D or report generation.
    """
    B_arr = np.asarray(B_candidates, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B_candidates must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(candidate_internals):
        raise ValueError("B_candidates row count must match candidate_internals")
    starting = [int(i) for i in starting_idx]
    if any(i < 0 or i >= len(candidate_internals) for i in starting):
        raise ValueError("starting_idx contains an out-of-range candidate index")
    if not starting:
        raise ValueError("starting_idx must contain at least one accepted primitive coordinate")

    starting_rank, starting_condition, _ = svd_rank_condition(B_arr[starting, :], tol_abs=tol_abs)
    if starting_rank <= 0:
        raise ValueError("starting_idx has zero recoverable rank")

    optimized_idx, optimizer_report = optimize_independent_coordinates_for_ped(
        B_arr,
        candidate_internals,
        starting,
        normal_modes,
        mode_indices,
        target_rank=starting_rank,
        tol_abs=tol_abs,
        max_passes=max_passes,
        improvement_tol=improvement_tol,
    )
    optimized_rank, optimized_condition, _ = svd_rank_condition(B_arr[optimized_idx, :], tol_abs=tol_abs)
    if optimized_rank < starting_rank:
        raise RuntimeError("rank-preserving composed PED basis optimization lost rank")

    selected_composed_indices = [
        idx for idx in optimized_idx
        if candidate_internals[idx].source == "composed_coordinate"
    ]
    report = dict(optimizer_report)
    report.update(
        {
            "rank_preserved": True,
            "required_rank": int(starting_rank),
            "starting_rank": int(starting_rank),
            "starting_condition": float(starting_condition),
            "optimized_rank": int(optimized_rank),
            "optimized_condition": float(optimized_condition),
            "composed_selected_count": int(len(selected_composed_indices)),
            "selected_composed_indices": selected_composed_indices,
        }
    )
    return optimized_idx, report


def _select_independent_coordinates_pivoted_cholesky(
    B: np.ndarray,
    internals: Sequence[InternalCoordinate],
    target_rank: int,
    tol_abs: float = 1.0e-6,
) -> List[int]:
    """
    Select a deterministic row basis when priority-ordered greedy SVD stalls.

    Aromatic coordinate pools can contain groups of nearly dependent rows where
    no single low-priority row clears the absolute SVD threshold against the
    current greedy basis, even though the full redundant B matrix has additional
    rank. Pivoted Cholesky on the row Gram matrix selects rows by residual power
    and recovers the rank reported for the redundant matrix.
    """
    if B.size == 0 or target_rank <= 0:
        return []

    ordered = sorted(range(len(internals)), key=lambda i: (internals[i].priority, internals[i].name))
    ordered_B = np.asarray(B[ordered, :], dtype=float)
    gram = ordered_B @ ordered_B.T
    nrows = gram.shape[0]
    max_rank = min(int(target_rank), nrows)
    diag = np.diag(gram).astype(float).copy()
    factors = np.zeros((nrows, max_rank), dtype=float)
    selected_positions: List[int] = []
    disabled = np.zeros(nrows, dtype=bool)

    for k in range(max_rank):
        active_diag = np.where(disabled, -np.inf, diag)
        pivot_pos = int(np.argmax(active_diag))
        pivot_residual = float(active_diag[pivot_pos])
        if not np.isfinite(pivot_residual) or pivot_residual <= tol_abs * tol_abs:
            break

        selected_positions.append(pivot_pos)
        pivot = float(np.sqrt(pivot_residual))
        if k:
            projection = factors[:, :k] @ factors[pivot_pos, :k]
        else:
            projection = 0.0
        factors[:, k] = (gram[:, pivot_pos] - projection) / pivot
        diag = diag - factors[:, k] ** 2
        disabled[pivot_pos] = True

    return [ordered[pos] for pos in selected_positions]
