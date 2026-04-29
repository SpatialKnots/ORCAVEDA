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
    rank, cond, s = svd_rank_condition(current, tol_abs=tol_abs)
    return selected_idx, rank, cond, s
