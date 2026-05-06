from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

import numpy as np

from orcaveda_models import InternalCoordinate


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
