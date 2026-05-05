from __future__ import annotations

import math
from typing import List, Sequence, Tuple

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
