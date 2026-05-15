from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from b_matrix import finite_difference_B
from mode_assignment import _assignment_family_from_internal, _compact_coord_label, _stage3d_coord_class
from orcaveda_models import HessData, InternalCoordinate
from orca_parser import BOHR_TO_ANGSTROM
from ped import build_wilson_g_matrix, wilson_coordinate_scales


WILSON_GF_VALIDATION_METHOD = (
    "Wilson GF diagonalization validation prototype using selected nonredundant "
    "internal-coordinate basis; diagnostic only, not VEDA-equivalent PED"
)
VEDA_LIKE_PED_METHOD = (
    "Comparable VEDA-like closed Wilson GF/PED audit using ORCAVEDA-selected "
    "nonredundant internal-coordinate basis; diagnostic only, does not reproduce original VEDA"
)
VEDA_LIKE_PED_MATRIX_ORIENTATION = "mode_rows_by_coordinate_columns_long_form"

HARTREE_J = 4.3597447222071e-18
AMU_KG = 1.66053906892e-27
ANGSTROM_M = 1.0e-10
LIGHT_SPEED_CM_S = 2.99792458e10
WILSON_GF_FIXED_CONVERSION_CM1 = (
    np.sqrt(HARTREE_J / (AMU_KG * ANGSTROM_M * ANGSTROM_M))
    / (2.0 * np.pi * LIGHT_SPEED_CM_S)
)


@dataclass
class WilsonGFResult:
    filename: str
    atom_count: int
    cartesian_size: int
    internal_basis_size: int
    expected_vibrational_rank: int
    basis_indices: Tuple[int, ...]
    g_rank: int
    g_condition: float
    f_rank: int
    f_condition: float
    gf_eigenvalues: np.ndarray
    gf_eigenvectors: np.ndarray
    orca_mode_indices: np.ndarray
    orca_frequencies_cm1: np.ndarray
    reconstructed_frequencies_cm1: np.ndarray
    empirical_ratios: np.ndarray
    max_relative_error: float
    empirical_ratio_median: float
    empirical_ratio_std: float
    orca_nonpositive_mode_count: int
    orca_min_nonpositive_frequency_cm1: float
    gf_nonpositive_eigenvalue_count: int
    gf_min_nonpositive_eigenvalue: float
    validation_status: str
    warnings: Tuple[str, ...]
    mapping_method: str
    conversion_method: str
    epm_optimized: bool = False
    epm_swaps: int = 0
    epm_initial_localization_score: float = 0.0
    epm_optimized_localization_score: float = 0.0
    epm_initial_mean_top_percent: float = 0.0
    epm_optimized_mean_top_percent: float = 0.0
    epm_initial_diffuse_mode_fraction: float = 1.0
    epm_optimized_diffuse_mode_fraction: float = 1.0
    validation_internals: Tuple[InternalCoordinate, ...] = ()
    validation_B: np.ndarray | None = None


def _empty_wilson_gf_epm_metrics() -> Dict[str, float]:
    return {
        "mode_count": 0.0,
        "mean_top_percent": 0.0,
        "median_top_percent": 0.0,
        "min_top_percent": 0.0,
        "diffuse_mode_fraction": 1.0,
        "localization_score": 0.0,
    }


def _rank_condition(matrix: np.ndarray, tol: float) -> Tuple[int, float]:
    if matrix.size == 0:
        return 0, float("inf")
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    finite = singular_values[np.isfinite(singular_values)]
    if finite.size == 0:
        return 0, float("inf")
    rank = int(np.sum(finite > tol))
    if rank == 0:
        return 0, float("inf")
    nonzero = finite[finite > tol]
    return rank, float(nonzero[0] / nonzero[-1])


def _wilson_g_rank_condition(
    B_arr: np.ndarray,
    internals: Sequence[InternalCoordinate],
    basis_idx: Sequence[int],
    masses: np.ndarray,
    tol: float,
) -> Tuple[int, float]:
    basis_internals = [internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = B_arr[list(basis_idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, masses)
    return _rank_condition(G, tol)


def _wilson_gf_rank_condition(
    B_arr: np.ndarray,
    internals: Sequence[InternalCoordinate],
    basis_idx: Sequence[int],
    masses: np.ndarray,
    cartesian_hessian: np.ndarray,
    tol: float,
) -> Tuple[int, float, int, float]:
    basis_internals = [internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = B_arr[list(basis_idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, masses)
    F_internal = reconstruct_wilson_gf_internal_force_matrix(B_internal, masses, cartesian_hessian, G)
    g_rank, g_condition = _rank_condition(G, tol)
    f_rank, f_condition = _rank_condition(F_internal, tol)
    return g_rank, g_condition, f_rank, f_condition


def _select_mass_weighted_pivot_basis(
    B_arr: np.ndarray,
    internals: Sequence[InternalCoordinate],
    masses: np.ndarray,
    expected_rank: int,
) -> Tuple[int, ...]:
    if expected_rank <= 0:
        return ()
    B_base = np.asarray(B_arr, dtype=float)
    if B_base.ndim != 2 or B_base.shape[0] < expected_rank:
        return ()
    mass_vec = np.repeat(np.asarray(masses, dtype=float), 3)
    if mass_vec.size != B_base.shape[1] or np.any(mass_vec <= 0.0) or not np.all(np.isfinite(mass_vec)):
        return ()

    scales = wilson_coordinate_scales(internals)
    weighted = (B_base * scales[:, None]) / np.sqrt(mass_vec)[None, :]
    weighted[~np.isfinite(weighted)] = 0.0
    remaining = set(range(weighted.shape[0]))
    selected: List[int] = []
    q_basis = np.zeros((weighted.shape[1], 0), dtype=float)
    for _ in range(expected_rank):
        best_idx = None
        best_residual_norm = 0.0
        best_residual = None
        for idx in remaining:
            residual = weighted[idx, :].astype(float, copy=True)
            if q_basis.shape[1]:
                residual = residual - q_basis @ (q_basis.T @ residual)
            residual_norm = float(np.linalg.norm(residual))
            if residual_norm > best_residual_norm:
                best_idx = int(idx)
                best_residual_norm = residual_norm
                best_residual = residual
        if best_idx is None or best_residual is None or best_residual_norm <= 0.0:
            return ()
        selected.append(best_idx)
        remaining.remove(best_idx)
        q_basis = np.column_stack([q_basis, best_residual / best_residual_norm])
    return tuple(selected)


def _select_conditioned_wilson_basis(
    B_arr: np.ndarray,
    internals: Sequence[InternalCoordinate],
    selected_idx: Sequence[int],
    masses: np.ndarray,
    expected_rank: int,
    tol: float,
    coords_A: np.ndarray | None = None,
    cartesian_hessian: np.ndarray | None = None,
) -> Tuple[int, ...]:
    """Small-system fallback for a Wilson-GF-conditioned validation basis."""
    basis_idx = tuple(int(idx) for idx in selected_idx)
    if expected_rank <= 0 or len(basis_idx) != expected_rank:
        return basis_idx

    g_rank, g_condition = _wilson_g_rank_condition(B_arr, internals, basis_idx, masses, tol)
    near_linear_selected = (
        coords_A is not None and _has_near_linear_bend(internals, basis_idx, np.asarray(coords_A, dtype=float))
    )
    f_rank = 0
    f_condition = float("inf")
    if cartesian_hessian is not None:
        _, _, f_rank, f_condition = _wilson_gf_rank_condition(
            B_arr,
            internals,
            basis_idx,
            masses,
            cartesian_hessian,
            tol,
        )
    current_f_acceptable = cartesian_hessian is None or (
        f_rank >= expected_rank and np.isfinite(f_condition) and f_condition <= 1.0e12
    )
    if (
        g_rank >= expected_rank
        and np.isfinite(g_condition)
        and g_condition <= 1.0e12
        and current_f_acceptable
        and not near_linear_selected
    ):
        return basis_idx

    candidate_count = len(internals)
    if candidate_count > 24 or comb(candidate_count, expected_rank) > 500_000:
        pivot_basis = _select_mass_weighted_pivot_basis(B_arr, internals, masses, expected_rank)
        if len(pivot_basis) == expected_rank:
            if cartesian_hessian is None:
                pivot_g_rank, pivot_g_condition = _wilson_g_rank_condition(
                    B_arr,
                    internals,
                    pivot_basis,
                    masses,
                    tol,
                )
                pivot_f_rank, pivot_f_condition = 0, float("inf")
            else:
                pivot_g_rank, pivot_g_condition, pivot_f_rank, pivot_f_condition = _wilson_gf_rank_condition(
                    B_arr,
                    internals,
                    pivot_basis,
                    masses,
                    cartesian_hessian,
                    tol,
                )
            improves_g = not np.isfinite(g_condition) or pivot_g_condition < g_condition
            improves_f = (
                cartesian_hessian is not None
                and pivot_f_rank >= expected_rank
                and np.isfinite(pivot_f_condition)
                and (not np.isfinite(f_condition) or pivot_f_condition < f_condition)
            )
            if (
                pivot_g_rank >= expected_rank
                and np.isfinite(pivot_g_condition)
                and (improves_g or improves_f)
            ):
                return pivot_basis
        return basis_idx

    best_basis = basis_idx
    best_g_condition = g_condition if np.isfinite(g_condition) else float("inf")
    for combo in combinations(range(candidate_count), expected_rank):
        b_rank, _ = _rank_condition(B_arr[list(combo), :], 1.0e-6)
        if b_rank < expected_rank:
            continue
        candidate_g_rank, candidate_g_condition = _wilson_g_rank_condition(
            B_arr,
            internals,
            combo,
            masses,
            tol,
        )
        if candidate_g_rank < expected_rank or not np.isfinite(candidate_g_condition):
            continue
        if candidate_g_condition < best_g_condition:
            best_basis = tuple(int(idx) for idx in combo)
            best_g_condition = candidate_g_condition
    return best_basis


def _orthonormal_perpendicular_frame(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    axis_arr = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis_arr)
    if norm <= 0.0:
        raise ValueError("linear bend reference axis has zero length")
    unit_axis = axis_arr / norm
    ref = min(np.eye(3), key=lambda candidate: abs(float(np.dot(candidate, unit_axis))))
    first = ref - np.dot(ref, unit_axis) * unit_axis
    first /= np.linalg.norm(first)
    second = np.cross(unit_axis, first)
    second /= np.linalg.norm(second)
    return first, second


def _linear_bend_component_fn(i: int, j: int, k: int, component_axis: np.ndarray):
    component = np.asarray(component_axis, dtype=float)

    def coordinate(xyz: np.ndarray) -> float:
        left = xyz[i] - xyz[j]
        right = xyz[k] - xyz[j]
        left_norm = np.linalg.norm(left)
        right_norm = np.linalg.norm(right)
        if left_norm <= 0.0 or right_norm <= 0.0:
            return float("nan")
        bend_vector = left / left_norm + right / right_norm
        return float(np.dot(bend_vector, component) * 180.0 / np.pi)

    return coordinate


def _augment_linear_bend_coordinates(
    internals: Sequence[InternalCoordinate],
    B_arr: np.ndarray,
    coords_A: np.ndarray,
    *,
    angle_tol_degrees: float = 5.0,
) -> Tuple[List[InternalCoordinate], np.ndarray]:
    augmented = list(internals)
    if any(internal.kind == "linear_bend_component" for internal in augmented):
        return augmented, np.asarray(B_arr, dtype=float)
    extra: List[InternalCoordinate] = []
    coords = np.asarray(coords_A, dtype=float)
    for internal in internals:
        if internal.kind != "bend" or len(internal.atoms0) != 3:
            continue
        i, j, k = internal.atoms0
        left = coords[i] - coords[j]
        right = coords[k] - coords[j]
        norm = np.linalg.norm(left) * np.linalg.norm(right)
        if norm <= 0.0:
            continue
        angle = float(np.degrees(np.arccos(np.clip(float(np.dot(left, right) / norm), -1.0, 1.0))))
        if angle_tol_degrees < angle < 180.0 - angle_tol_degrees:
            continue
        first, second = _orthonormal_perpendicular_frame(left)
        label = internal.name.replace("ang(", "linbend(")
        extra.append(
            InternalCoordinate(
                f"{label}:perp1)",
                "linear_bend_component",
                internal.atoms0,
                max(1, int(internal.priority) - 1),
                _linear_bend_component_fn(i, j, k, first),
                "wilson_gf_validation",
                generation_rule="near_linear_bend_perpendicular_components",
            )
        )
        extra.append(
            InternalCoordinate(
                f"{label}:perp2)",
                "linear_bend_component",
                internal.atoms0,
                max(1, int(internal.priority) - 1),
                _linear_bend_component_fn(i, j, k, second),
                "wilson_gf_validation",
                generation_rule="near_linear_bend_perpendicular_components",
            )
        )
    if not extra:
        return augmented, np.asarray(B_arr, dtype=float)
    B_extra = finite_difference_B(coords, extra)
    return augmented + extra, np.vstack([np.asarray(B_arr, dtype=float), B_extra])


def _has_near_linear_bend(
    internals: Sequence[InternalCoordinate],
    basis_idx: Sequence[int],
    coords_A: np.ndarray,
    *,
    angle_tol_degrees: float = 5.0,
) -> bool:
    coords = np.asarray(coords_A, dtype=float)
    for idx in basis_idx:
        internal = internals[int(idx)]
        if internal.kind != "bend" or len(internal.atoms0) != 3:
            continue
        a, b, c = internal.atoms0
        v1 = coords[a] - coords[b]
        v2 = coords[c] - coords[b]
        norm = np.linalg.norm(v1) * np.linalg.norm(v2)
        if norm <= 0.0:
            continue
        cos_angle = float(np.dot(v1, v2) / norm)
        angle = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
        if angle <= angle_tol_degrees or angle >= 180.0 - angle_tol_degrees:
            return True
    return False


def _uses_linear_bend_component(internals: Sequence[InternalCoordinate], basis_idx: Sequence[int]) -> bool:
    return any(internals[int(idx)].kind == "linear_bend_component" for idx in basis_idx)


def _join_warning_tokens(*parts: str) -> str:
    tokens: List[str] = []
    for part in parts:
        for token in str(part or "").split(";"):
            cleaned = token.strip()
            if cleaned and cleaned not in tokens:
                tokens.append(cleaned)
    return "; ".join(tokens)


def symmetric_sqrt_decomp(A: np.ndarray, tol: float = 1.0e-12) -> Tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(A, dtype=float)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError(f"A must be a square 2D matrix, got shape {arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError("A must contain only finite values")
    if not np.allclose(arr, arr.T, atol=tol, rtol=0.0):
        raise ValueError("A must be symmetric")

    eigvals, eigvecs = np.linalg.eigh(0.5 * (arr + arr.T))
    if np.any(eigvals < -tol):
        raise ValueError("A must be positive semidefinite within tolerance")
    clipped = np.clip(eigvals, 0.0, None)
    sqrt_arr = (eigvecs * np.sqrt(clipped)) @ eigvecs.T
    sqrt_arr = 0.5 * (sqrt_arr + sqrt_arr.T)
    return sqrt_arr, eigvals


def solve_symmetric_gf_eigenproblem(
    G: np.ndarray,
    F: np.ndarray,
    *,
    tol: float = 1.0e-12,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    G_arr = np.asarray(G, dtype=float)
    F_arr = np.asarray(F, dtype=float)
    if G_arr.shape != F_arr.shape:
        raise ValueError(f"G shape {G_arr.shape} does not match F shape {F_arr.shape}")
    if G_arr.ndim != 2 or G_arr.shape[0] != G_arr.shape[1]:
        raise ValueError(f"G and F must be square 2D matrices, got {G_arr.shape}")
    if not np.allclose(F_arr, F_arr.T, atol=tol, rtol=0.0):
        raise ValueError("F must be symmetric")

    G_sqrt, _ = symmetric_sqrt_decomp(G_arr, tol=tol)
    symmetric_gf = G_sqrt @ (0.5 * (F_arr + F_arr.T)) @ G_sqrt
    symmetric_gf = 0.5 * (symmetric_gf + symmetric_gf.T)
    eigenvalues, eigenvectors_orth = np.linalg.eigh(symmetric_gf)
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors_orth = eigenvectors_orth[:, order]
    internal_eigenvectors = G_sqrt @ eigenvectors_orth
    return eigenvalues, internal_eigenvectors, symmetric_gf


def reconstruct_wilson_gf_internal_force_matrix(
    B_internal: np.ndarray,
    masses: np.ndarray,
    cartesian_hessian: np.ndarray,
    G: np.ndarray,
) -> np.ndarray:
    B_arr = np.asarray(B_internal, dtype=float)
    mass_vec = np.repeat(np.asarray(masses, dtype=float), 3)
    F_cart = np.asarray(cartesian_hessian, dtype=float)
    G_arr = np.asarray(G, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B_internal must be a 2D matrix, got shape {B_arr.shape}")
    if F_cart.shape != (B_arr.shape[1], B_arr.shape[1]):
        raise ValueError(f"cartesian_hessian shape {F_cart.shape} does not match B column count {B_arr.shape[1]}")
    if mass_vec.size != B_arr.shape[1]:
        raise ValueError(f"mass vector length {mass_vec.size} does not match B column count {B_arr.shape[1]}")
    if G_arr.shape != (B_arr.shape[0], B_arr.shape[0]):
        raise ValueError(f"G shape {G_arr.shape} does not match internal basis size {B_arr.shape[0]}")
    if np.any(mass_vec <= 0.0) or not np.all(np.isfinite(mass_vec)):
        raise ValueError("masses must be finite and positive for Wilson GF internal force reconstruction")

    F_cart_A = 0.5 * (F_cart + F_cart.T) / (BOHR_TO_ANGSTROM ** 2)
    mass_inverse = np.diag(1.0 / mass_vec)
    G_inverse = np.linalg.pinv(G_arr)
    cartesian_from_internal = mass_inverse @ B_arr.T @ G_inverse
    F_internal = cartesian_from_internal.T @ F_cart_A @ cartesian_from_internal
    F_internal = 0.5 * (F_internal + F_internal.T)
    F_internal[~np.isfinite(F_internal)] = 0.0
    return F_internal


def wilson_gf_ped_localization_metrics(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    basis_idx: Sequence[int],
    *,
    tol: float = 1.0e-12,
) -> Dict[str, float]:
    """
    Score a nonredundant basis using closed Wilson GF/PED contributions.

    This is an opt-in EPM-like diagnostic objective for ORCAVEDA's VEDA-like
    audit layer. It maximizes dominant closed GF/PED terms and is not evidence
    of original VEDA reproduction.
    """
    if hess.cartesian_hessian is None:
        raise ValueError("Wilson GF EPM optimization requires HessData.cartesian_hessian parsed from ORCA $hessian")
    B_arr = np.asarray(B, dtype=float)
    idx = tuple(int(i) for i in basis_idx)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(internals):
        raise ValueError(f"B row count {B_arr.shape[0]} does not match internal coordinate count {len(internals)}")
    if any(i < 0 or i >= len(internals) for i in idx):
        raise ValueError("basis_idx contains an out-of-range internal coordinate index")
    if not idx:
        return _empty_wilson_gf_epm_metrics()

    basis_internals = [internals[i] for i in idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = B_arr[list(idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, hess.masses)
    F_internal = reconstruct_wilson_gf_internal_force_matrix(B_internal, hess.masses, hess.cartesian_hessian, G)
    eigenvalues, eigenvectors, _ = solve_symmetric_gf_eigenproblem(G, F_internal, tol=tol)
    positive_cols = np.where(np.isfinite(eigenvalues) & (eigenvalues > tol))[0]
    positive_freq_mask = np.isfinite(hess.frequencies_cm1) & (hess.frequencies_cm1 > tol)
    pair_count = min(len(positive_cols), int(np.sum(positive_freq_mask)))
    if pair_count <= 0:
        return _empty_wilson_gf_epm_metrics()

    top_values: List[float] = []
    for col in positive_cols[:pair_count]:
        vector = np.asarray(eigenvectors[:, col], dtype=float)
        force_response = F_internal @ vector
        signed_terms = vector * force_response
        signed_terms[~np.isfinite(signed_terms)] = 0.0
        weights = np.abs(signed_terms)
        total = float(np.sum(weights))
        if total <= tol or not np.isfinite(total):
            continue
        pct = 100.0 * weights / total
        top_values.append(float(np.max(pct)))
    if not top_values:
        return _empty_wilson_gf_epm_metrics()

    top = np.asarray(top_values, dtype=float)
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


def optimize_wilson_gf_basis_for_epm(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    target_rank: int | None = None,
    tol: float = 1.0e-12,
    max_passes: int = 2,
    improvement_tol: float = 1.0e-6,
) -> Tuple[Tuple[int, ...], Dict[str, float]]:
    """
    Greedily swap internal coordinates to improve closed Wilson GF/PED locality.

    Rank preservation is mandatory. This optimizer is opt-in and scoped to the
    Wilson GF/VEDA-like diagnostic basis; it does not change Stage 3D assignment
    audit labels.
    """
    B_arr = np.asarray(B, dtype=float)
    selected = [int(i) for i in selected_idx]
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(internals):
        raise ValueError(f"B row count {B_arr.shape[0]} does not match internal coordinate count {len(internals)}")
    if any(i < 0 or i >= len(internals) for i in selected):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")
    if not selected:
        metrics = _empty_wilson_gf_epm_metrics()
        report = {"changed": False, "swaps": 0, "rank": 0, "condition": float("inf"), "initial_rank": 0, "initial_condition": float("inf")}
        report.update({f"initial_{key}": value for key, value in metrics.items()})
        report.update({f"optimized_{key}": value for key, value in metrics.items()})
        return (), report

    start_rank, start_condition, start_f_rank, start_f_condition = _wilson_gf_rank_condition(
        B_arr,
        internals,
        selected,
        hess.masses,
        hess.cartesian_hessian,
        tol,
    )
    required_rank = int(target_rank) if target_rank is not None else int(start_rank)
    required_rank = min(required_rank, int(start_rank))
    initial_metrics = wilson_gf_ped_localization_metrics(hess, internals, B_arr, selected, tol=tol)
    best_score = float(initial_metrics["localization_score"])
    swaps = 0

    candidate_pool = [
        idx for idx in range(len(internals))
        if idx not in set(selected) and np.linalg.norm(B_arr[idx, :]) > tol
    ]
    candidate_pool = sorted(candidate_pool, key=lambda i: (internals[i].priority, internals[i].name))

    for _pass in range(max(0, int(max_passes))):
        improved = False
        for pos, _old_idx in enumerate(list(selected)):
            best_replacement: tuple[float, int, Dict[str, float]] | None = None
            for candidate_idx in candidate_pool:
                if candidate_idx in selected:
                    continue
                trial = list(selected)
                trial[pos] = candidate_idx
                if len(set(trial)) != len(trial):
                    continue
                rank, condition, f_rank, f_condition = _wilson_gf_rank_condition(
                    B_arr,
                    internals,
                    trial,
                    hess.masses,
                    hess.cartesian_hessian,
                    tol,
                )
                if (
                    rank < required_rank
                    or f_rank < required_rank
                    or not np.isfinite(condition)
                    or not np.isfinite(f_condition)
                ):
                    continue
                metrics = wilson_gf_ped_localization_metrics(hess, internals, B_arr, trial, tol=tol)
                score = float(metrics["localization_score"])
                if score <= best_score + improvement_tol:
                    continue
                if best_replacement is None or score > best_replacement[0] + improvement_tol:
                    best_replacement = (score, candidate_idx, metrics)
                elif abs(score - best_replacement[0]) <= improvement_tol:
                    tie_break = (internals[candidate_idx].priority, internals[candidate_idx].name)
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

    optimized_rank, optimized_condition, optimized_f_rank, optimized_f_condition = _wilson_gf_rank_condition(
        B_arr,
        internals,
        selected,
        hess.masses,
        hess.cartesian_hessian,
        tol,
    )
    optimized_metrics = wilson_gf_ped_localization_metrics(hess, internals, B_arr, selected, tol=tol)
    report = {
        "changed": bool(swaps > 0),
        "swaps": int(swaps),
        "rank": int(optimized_rank),
        "condition": float(optimized_condition),
        "f_rank": int(optimized_f_rank),
        "f_condition": float(optimized_f_condition),
        "initial_rank": int(start_rank),
        "initial_condition": float(start_condition),
        "initial_f_rank": int(start_f_rank),
        "initial_f_condition": float(start_f_condition),
    }
    report.update({f"initial_{key}": value for key, value in initial_metrics.items()})
    report.update({f"optimized_{key}": value for key, value in optimized_metrics.items()})
    return tuple(int(idx) for idx in selected), report


def wilson_gf_diagonalization(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    tol: float = 1.0e-12,
    frequency_tol_relative: float = 1.0e-4,
    epm_optimize: bool = False,
    epm_max_passes: int = 2,
    epm_improvement_tol: float = 1.0e-6,
) -> WilsonGFResult:
    warnings: List[str] = []
    if hess.cartesian_hessian is None:
        raise ValueError("Wilson GF validation requires HessData.cartesian_hessian parsed from ORCA $hessian")

    atom_count = len(hess.atoms)
    cartesian_size = 3 * atom_count
    B_arr = np.asarray(B, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(internals):
        raise ValueError(f"B row count {B_arr.shape[0]} does not match internal coordinate count {len(internals)}")
    if B_arr.shape[1] != cartesian_size:
        raise ValueError(f"B column count {B_arr.shape[1]} does not match 3N={cartesian_size}")

    validation_internals, validation_B = _augment_linear_bend_coordinates(internals, B_arr, hess.coords_A)
    basis_idx = tuple(int(idx) for idx in selected_idx)
    if any(idx < 0 or idx >= len(internals) for idx in basis_idx):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")
    if not basis_idx:
        raise ValueError("selected_idx must contain at least one internal coordinate")

    expected_vibrational_rank = max(cartesian_size - 6, 0)
    basis_idx = _select_conditioned_wilson_basis(
        validation_B,
        validation_internals,
        basis_idx,
        hess.masses,
        expected_vibrational_rank,
        tol,
        hess.coords_A,
        hess.cartesian_hessian,
    )
    epm_report: Dict[str, float | int | bool] = {}
    if epm_optimize:
        basis_idx, epm_report = optimize_wilson_gf_basis_for_epm(
            hess,
            validation_internals,
            validation_B,
            basis_idx,
            target_rank=expected_vibrational_rank,
            tol=tol,
            max_passes=epm_max_passes,
            improvement_tol=epm_improvement_tol,
        )
    basis_internals = [validation_internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = validation_B[list(basis_idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, hess.masses)
    F_internal = reconstruct_wilson_gf_internal_force_matrix(B_internal, hess.masses, hess.cartesian_hessian, G)
    g_rank, g_condition = _rank_condition(G, tol)
    f_rank, f_condition = _rank_condition(F_internal, tol)

    if len(basis_idx) != expected_vibrational_rank:
        warnings.append("basis_size_mismatch_expected_vibrational_rank")
    if _has_near_linear_bend(validation_internals, basis_idx, hess.coords_A):
        warnings.append("near_linear_bend_coordinate")
    if _uses_linear_bend_component(validation_internals, basis_idx):
        warnings.append("linear_bend_coordinate_used")
    if g_rank < min(len(basis_idx), expected_vibrational_rank):
        warnings.append("basis_rank_below_expected")
    if not np.isfinite(g_condition) or g_condition > 1.0e12:
        warnings.append("g_ill_conditioned")
    if not np.isfinite(f_condition) or f_condition > 1.0e12:
        warnings.append("f_ill_conditioned")

    eigenvalues, eigenvectors, _ = solve_symmetric_gf_eigenproblem(G, F_internal, tol=tol)
    positive_mask = np.isfinite(eigenvalues) & (eigenvalues > tol)
    nonpositive_gf = eigenvalues[np.isfinite(eigenvalues) & (eigenvalues <= tol)]
    positive_gf = eigenvalues[positive_mask]
    positive_eigenvectors = eigenvectors[:, positive_mask]
    positive_freq_mask = np.isfinite(hess.frequencies_cm1) & (hess.frequencies_cm1 > tol)
    nonpositive_orca = np.asarray(hess.frequencies_cm1[np.isfinite(hess.frequencies_cm1) & (hess.frequencies_cm1 <= tol)], dtype=float)
    orca_mode_indices = np.where(positive_freq_mask)[0].astype(int)
    orca_freqs = np.asarray(hess.frequencies_cm1[positive_freq_mask], dtype=float)
    gf_order = np.argsort(positive_gf)
    positive_gf = positive_gf[gf_order]
    positive_eigenvectors = positive_eigenvectors[:, gf_order]
    orca_order = np.argsort(orca_freqs)
    orca_mode_indices = orca_mode_indices[orca_order]
    orca_freqs = orca_freqs[orca_order]

    pair_count = min(len(positive_gf), len(orca_freqs))
    if len(orca_freqs) < expected_vibrational_rank:
        warnings.append("positive_orca_mode_count_below_expected_vibrational_rank")
        if len(nonpositive_orca) > 6:
            warnings.append("nonpositive_orca_modes_within_expected_vibrational_space")
    if len(positive_gf) < expected_vibrational_rank:
        warnings.append("positive_gf_eigenvalue_count_below_expected_vibrational_rank")
        if len(nonpositive_gf) > 0:
            warnings.append("nonpositive_gf_eigenvalues_within_expected_vibrational_space")
    if len(positive_gf) != len(orca_freqs):
        warnings.append("mode_count_mismatch")
    positive_gf = positive_gf[:pair_count]
    positive_eigenvectors = positive_eigenvectors[:, :pair_count]
    orca_mode_indices = orca_mode_indices[:pair_count]
    orca_freqs = orca_freqs[:pair_count]

    reconstructed = WILSON_GF_FIXED_CONVERSION_CM1 * np.sqrt(np.clip(positive_gf, 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_errors = np.abs(reconstructed - orca_freqs) / np.maximum(np.abs(orca_freqs), tol)
        empirical = orca_freqs / np.sqrt(np.clip(positive_gf, tol, None))
    finite_errors = rel_errors[np.isfinite(rel_errors)]
    finite_empirical = empirical[np.isfinite(empirical)]
    max_relative_error = float(np.max(finite_errors)) if finite_errors.size else float("inf")
    empirical_median = float(np.median(finite_empirical)) if finite_empirical.size else float("nan")
    empirical_std = float(np.std(finite_empirical)) if finite_empirical.size else float("nan")
    orca_min_nonpositive = float(np.min(nonpositive_orca)) if nonpositive_orca.size else float("nan")
    gf_min_nonpositive = float(np.min(nonpositive_gf)) if nonpositive_gf.size else float("nan")

    fixed_conversion_failed = not np.isfinite(max_relative_error) or max_relative_error > frequency_tol_relative
    if fixed_conversion_failed:
        warnings.append("fixed_conversion_failed")
        if "near_linear_bend_coordinate" in warnings or "linear_bend_coordinate_used" in warnings:
            warnings.append("linear_or_near_linear_fixed_conversion_review")
    if fixed_conversion_failed and finite_empirical.size:
        warnings.append("empirical_ratio_only")

    if pair_count == 0 or "basis_rank_below_expected" in warnings or "mode_count_mismatch" in warnings:
        status = "FAIL"
    elif max_relative_error <= frequency_tol_relative:
        status = "PASS"
    elif finite_empirical.size:
        status = "WARN"
    else:
        status = "FAIL"

    return WilsonGFResult(
        filename=hess.filename,
        atom_count=atom_count,
        cartesian_size=cartesian_size,
        internal_basis_size=len(basis_idx),
        expected_vibrational_rank=expected_vibrational_rank,
        basis_indices=basis_idx,
        g_rank=g_rank,
        g_condition=g_condition,
        f_rank=f_rank,
        f_condition=f_condition,
        gf_eigenvalues=positive_gf,
        gf_eigenvectors=positive_eigenvectors,
        orca_mode_indices=orca_mode_indices,
        orca_frequencies_cm1=orca_freqs,
        reconstructed_frequencies_cm1=reconstructed,
        empirical_ratios=empirical,
        max_relative_error=max_relative_error,
        empirical_ratio_median=empirical_median,
        empirical_ratio_std=empirical_std,
        orca_nonpositive_mode_count=int(nonpositive_orca.size),
        orca_min_nonpositive_frequency_cm1=orca_min_nonpositive,
        gf_nonpositive_eigenvalue_count=int(nonpositive_gf.size),
        gf_min_nonpositive_eigenvalue=gf_min_nonpositive,
        validation_status=status,
        warnings=tuple(dict.fromkeys(warnings)),
        mapping_method="sorted_positive_gf_eigenvalues_to_sorted_positive_orca_frequencies",
        conversion_method=f"fixed_SI_hartree_per_amu_angstrom2_to_cm-1:{WILSON_GF_FIXED_CONVERSION_CM1:.12g}",
        epm_optimized=bool(epm_report.get("changed", False)),
        epm_swaps=int(epm_report.get("swaps", 0)),
        epm_initial_localization_score=float(epm_report.get("initial_localization_score", 0.0)),
        epm_optimized_localization_score=float(epm_report.get("optimized_localization_score", 0.0)),
        epm_initial_mean_top_percent=float(epm_report.get("initial_mean_top_percent", 0.0)),
        epm_optimized_mean_top_percent=float(epm_report.get("optimized_mean_top_percent", 0.0)),
        epm_initial_diffuse_mode_fraction=float(epm_report.get("initial_diffuse_mode_fraction", 1.0)),
        epm_optimized_diffuse_mode_fraction=float(epm_report.get("optimized_diffuse_mode_fraction", 1.0)),
        validation_internals=tuple(validation_internals),
        validation_B=validation_B,
    )


def build_wilson_gf_validation_dataframe(result: WilsonGFResult) -> pd.DataFrame:
    warnings = "; ".join(result.warnings)
    rows = []
    for row_idx in range(len(result.orca_frequencies_cm1)):
        rows.append(
            {
                "Source": "",
                "Filename": result.filename,
                "mode_index": int(result.orca_mode_indices[row_idx]),
                "orca_frequency_cm-1": float(result.orca_frequencies_cm1[row_idx]),
                "gf_eigenvalue": float(result.gf_eigenvalues[row_idx]),
                "reconstructed_frequency_cm-1": float(result.reconstructed_frequencies_cm1[row_idx]),
                "fixed_conversion_relative_error": (
                    abs(float(result.reconstructed_frequencies_cm1[row_idx]) - float(result.orca_frequencies_cm1[row_idx]))
                    / max(abs(float(result.orca_frequencies_cm1[row_idx])), 1.0e-12)
                ),
                "empirical_ratio_frequency_cm1_per_sqrt_lambda": float(result.empirical_ratios[row_idx]),
                "mapping_method": result.mapping_method,
                "conversion_method": result.conversion_method,
                "validation_status": result.validation_status,
                "max_relative_error": result.max_relative_error,
                "empirical_ratio_median": result.empirical_ratio_median,
                "empirical_ratio_std": result.empirical_ratio_std,
                "orca_nonpositive_mode_count": result.orca_nonpositive_mode_count,
                "orca_min_nonpositive_frequency_cm-1": result.orca_min_nonpositive_frequency_cm1,
                "gf_nonpositive_eigenvalue_count": result.gf_nonpositive_eigenvalue_count,
                "gf_min_nonpositive_eigenvalue": result.gf_min_nonpositive_eigenvalue,
                "basis_size": result.internal_basis_size,
                "expected_vibrational_rank": result.expected_vibrational_rank,
                "g_rank": result.g_rank,
                "g_condition": result.g_condition,
                "f_rank": result.f_rank,
                "f_condition": result.f_condition,
                "warnings": warnings,
                "method": WILSON_GF_VALIDATION_METHOD,
            }
        )
    if rows:
        return pd.DataFrame(rows)
    return pd.DataFrame(
        [
            {
                "Source": "",
                "Filename": result.filename,
                "mode_index": "",
                "orca_frequency_cm-1": "",
                "gf_eigenvalue": "",
                "reconstructed_frequency_cm-1": "",
                "fixed_conversion_relative_error": "",
                "empirical_ratio_frequency_cm1_per_sqrt_lambda": "",
                "mapping_method": result.mapping_method,
                "conversion_method": result.conversion_method,
                "validation_status": result.validation_status,
                "max_relative_error": result.max_relative_error,
                "empirical_ratio_median": result.empirical_ratio_median,
                "empirical_ratio_std": result.empirical_ratio_std,
                "orca_nonpositive_mode_count": result.orca_nonpositive_mode_count,
                "orca_min_nonpositive_frequency_cm-1": result.orca_min_nonpositive_frequency_cm1,
                "gf_nonpositive_eigenvalue_count": result.gf_nonpositive_eigenvalue_count,
                "gf_min_nonpositive_eigenvalue": result.gf_min_nonpositive_eigenvalue,
                "basis_size": result.internal_basis_size,
                "expected_vibrational_rank": result.expected_vibrational_rank,
                "g_rank": result.g_rank,
                "g_condition": result.g_condition,
                "f_rank": result.f_rank,
                "f_condition": result.f_condition,
                "epm_optimized": result.epm_optimized,
                "epm_swaps": result.epm_swaps,
                "epm_initial_localization_score": result.epm_initial_localization_score,
                "epm_optimized_localization_score": result.epm_optimized_localization_score,
                "epm_initial_mean_top_percent": result.epm_initial_mean_top_percent,
                "epm_optimized_mean_top_percent": result.epm_optimized_mean_top_percent,
                "epm_initial_diffuse_mode_fraction": result.epm_initial_diffuse_mode_fraction,
                "epm_optimized_diffuse_mode_fraction": result.epm_optimized_diffuse_mode_fraction,
                "warnings": warnings,
                "method": WILSON_GF_VALIDATION_METHOD,
            }
        ]
    )


def build_wilson_gf_basis_diagnostics_dataframe(result: WilsonGFResult) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Source": "",
                "Filename": result.filename,
                "basis_size": result.internal_basis_size,
                "expected_vibrational_rank": result.expected_vibrational_rank,
                "selected_indices": ";".join(str(idx) for idx in result.basis_indices),
                "g_rank": result.g_rank,
                "g_condition": result.g_condition,
                "f_rank": result.f_rank,
                "f_condition": result.f_condition,
                "positive_orca_mode_count": int(len(result.orca_frequencies_cm1)),
                "positive_gf_eigenvalue_count": int(len(result.gf_eigenvalues)),
                "orca_nonpositive_mode_count": result.orca_nonpositive_mode_count,
                "orca_min_nonpositive_frequency_cm-1": result.orca_min_nonpositive_frequency_cm1,
                "gf_nonpositive_eigenvalue_count": result.gf_nonpositive_eigenvalue_count,
                "gf_min_nonpositive_eigenvalue": result.gf_min_nonpositive_eigenvalue,
                "epm_optimized": result.epm_optimized,
                "epm_swaps": result.epm_swaps,
                "epm_initial_localization_score": result.epm_initial_localization_score,
                "epm_optimized_localization_score": result.epm_optimized_localization_score,
                "epm_initial_mean_top_percent": result.epm_initial_mean_top_percent,
                "epm_optimized_mean_top_percent": result.epm_optimized_mean_top_percent,
                "epm_initial_diffuse_mode_fraction": result.epm_initial_diffuse_mode_fraction,
                "epm_optimized_diffuse_mode_fraction": result.epm_optimized_diffuse_mode_fraction,
                "warnings": "; ".join(result.warnings),
            }
        ]
    )


def wilson_gf_closed_ped(
    result: WilsonGFResult,
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    top_n: int = 8,
    tol: float = 1.0e-12,
) -> pd.DataFrame:
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    basis_idx = tuple(int(idx) for idx in selected_idx)
    if tuple(basis_idx) != result.basis_indices:
        raise ValueError("selected_idx must match result.basis_indices")
    basis_internals = [internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = np.asarray(B, dtype=float)[list(basis_idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, hess.masses)
    F_internal = reconstruct_wilson_gf_internal_force_matrix(B_internal, hess.masses, hess.cartesian_hessian, G)

    positive_cols = np.where(np.isfinite(result.gf_eigenvalues) & (result.gf_eigenvalues > tol))[0]
    rows = []
    for local_mode, col in enumerate(positive_cols):
        if local_mode >= len(result.orca_frequencies_cm1):
            break
        vector = np.asarray(result.gf_eigenvectors[:, col], dtype=float)
        force_response = F_internal @ vector
        signed_terms = vector * force_response
        signed_terms[~np.isfinite(signed_terms)] = 0.0
        weights = np.abs(signed_terms)
        total = float(np.sum(weights))
        if total > tol and np.isfinite(total):
            pct = 100.0 * weights / total
        else:
            pct = np.zeros_like(weights)
        order = np.argsort(pct)[::-1] if pct.size else np.array([], dtype=int)
        if not len(order) or not np.any(pct > 0.0):
            rows.append(
                {
                    "Source": "",
                    "Filename": result.filename,
                    "mode": int(result.orca_mode_indices[local_mode]),
                    "frequency_cm-1": float(result.orca_frequencies_cm1[local_mode]),
                    "gf_rank": 0,
                    "coord_index": "",
                    "internal_coordinate": "",
                    "coordinate_kind": "",
                    "coordinate_family": "",
                    "signed_ped_fraction": 0.0,
                    "contribution_percent": 0.0,
                    "normalization_sum_percent": 0.0,
                    "basis_size": result.internal_basis_size,
                    "validation_status": result.validation_status,
                    "max_relative_error": result.max_relative_error,
                    "warnings": "; ".join(result.warnings),
                    "method": WILSON_GF_VALIDATION_METHOD,
                }
            )
            continue
        for rank, local_idx in enumerate(order[:top_n], start=1):
            if pct[local_idx] <= 0.0:
                continue
            ic = basis_internals[int(local_idx)]
            rows.append(
                {
                    "Source": "",
                    "Filename": result.filename,
                    "mode": int(result.orca_mode_indices[local_mode]),
                    "frequency_cm-1": float(result.orca_frequencies_cm1[local_mode]),
                    "gf_rank": rank,
                    "coord_index": basis_idx[int(local_idx)],
                    "internal_coordinate": _compact_coord_label(ic.name),
                    "coordinate_kind": str(ic.kind),
                    "coordinate_family": _assignment_family_from_internal(ic),
                    "coordinate_class": _stage3d_coord_class(ic),
                    "signed_ped_fraction": float(signed_terms[local_idx] / total) if total > tol else 0.0,
                    "contribution_percent": float(pct[local_idx]),
                    "normalization_sum_percent": float(np.sum(pct)),
                    "basis_size": result.internal_basis_size,
                    "validation_status": result.validation_status,
                    "max_relative_error": result.max_relative_error,
                    "warnings": "; ".join(result.warnings),
                    "method": WILSON_GF_VALIDATION_METHOD,
                }
            )
    return pd.DataFrame(rows)


def _veda_like_signed_terms(
    result: WilsonGFResult,
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    tol: float = 1.0e-12,
) -> Tuple[Tuple[int, ...], List[InternalCoordinate], np.ndarray, np.ndarray]:
    basis_idx = tuple(int(idx) for idx in selected_idx)
    if tuple(basis_idx) != result.basis_indices:
        raise ValueError("selected_idx must match result.basis_indices")
    if hess.cartesian_hessian is None:
        raise ValueError("VEDA-like PED requires HessData.cartesian_hessian parsed from ORCA $hessian")
    basis_internals = [internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = np.asarray(B, dtype=float)[list(basis_idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, hess.masses)
    F_internal = reconstruct_wilson_gf_internal_force_matrix(B_internal, hess.masses, hess.cartesian_hessian, G)

    pair_count = min(len(result.orca_frequencies_cm1), result.gf_eigenvectors.shape[1])
    signed = np.zeros((pair_count, len(basis_idx)), dtype=float)
    pct = np.zeros_like(signed)
    for row_idx in range(pair_count):
        vector = np.asarray(result.gf_eigenvectors[:, row_idx], dtype=float)
        force_response = F_internal @ vector
        signed_terms = vector * force_response
        signed_terms[~np.isfinite(signed_terms)] = 0.0
        weights = np.abs(signed_terms)
        total = float(np.sum(weights))
        signed[row_idx, :] = signed_terms
        if total > tol and np.isfinite(total):
            pct[row_idx, :] = 100.0 * weights / total
    return basis_idx, basis_internals, signed, pct


def _is_xh_stretch_family(family: str) -> bool:
    text = str(family or "").lower()
    return ("c-h" in text or "n-h" in text or "o-h" in text) and "stretch" in text


def _is_hbond_or_interfragment_family(family: str) -> bool:
    text = str(family or "").lower()
    return "h-bond" in text or "hbond" in text or "intermolecular" in text or "interfragment" in text


def _veda_like_mode_warning(
    frequency_cm1: float,
    row_pct: np.ndarray,
    basis_internals: Sequence[InternalCoordinate],
    base_warnings: str,
    *,
    high_frequency_cm1: float = 2800.0,
) -> str:
    warnings = str(base_warnings or "")
    if not np.isfinite(float(frequency_cm1)) or float(frequency_cm1) <= high_frequency_cm1 or row_pct.size == 0:
        return warnings

    order = np.argsort(row_pct)[::-1]
    if not len(order) or float(row_pct[order[0]]) <= 0.0:
        return warnings
    top_family = _assignment_family_from_internal(basis_internals[int(order[0])])
    if not _is_hbond_or_interfragment_family(top_family):
        return warnings
    for local_idx in order[1:]:
        if float(row_pct[local_idx]) <= 0.0:
            continue
        family = _assignment_family_from_internal(basis_internals[int(local_idx)])
        if _is_xh_stretch_family(family):
            return _join_warning_tokens(warnings, "high_frequency_hbond_dominates_xh_stretch_secondary")
    return warnings


def build_veda_like_ped_audit_dataframe(
    result: WilsonGFResult,
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    source_label: str = "",
    top_n: int = 8,
    tol: float = 1.0e-12,
) -> pd.DataFrame:
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    basis_idx, basis_internals, signed, pct = _veda_like_signed_terms(
        result,
        hess,
        internals,
        B,
        selected_idx,
        tol=tol,
    )
    rows = []
    warnings = "; ".join(result.warnings)
    for row_idx in range(pct.shape[0]):
        row_warnings = _veda_like_mode_warning(
            float(result.orca_frequencies_cm1[row_idx]),
            pct[row_idx, :],
            basis_internals,
            warnings,
        )
        order = np.argsort(pct[row_idx, :])[::-1] if pct.shape[1] else np.array([], dtype=int)
        if not len(order) or not np.any(pct[row_idx, :] > 0.0):
            rows.append(
                {
                    "Source": source_label,
                    "Filename": result.filename,
                    "mode": int(result.orca_mode_indices[row_idx]),
                    "frequency_cm-1": float(result.orca_frequencies_cm1[row_idx]),
                    "gf_eigenvector_index": int(row_idx),
                    "veda_like_rank": 0,
                    "coord_index": "",
                    "internal_coordinate": "",
                    "coordinate_kind": "",
                    "coordinate_family": "",
                    "coordinate_class": "",
                    "signed_ped_fraction": 0.0,
                    "contribution_percent": 0.0,
                    "normalization_sum_percent": 0.0,
                    "matrix_orientation": VEDA_LIKE_PED_MATRIX_ORIENTATION,
                    "basis_size": result.internal_basis_size,
                    "validation_status": result.validation_status,
                    "max_relative_error": result.max_relative_error,
                    "warnings": row_warnings or "zero_veda_like_ped_distribution",
                    "method": VEDA_LIKE_PED_METHOD,
                }
            )
            continue
        normalization = float(np.sum(pct[row_idx, :]))
        for rank, local_idx in enumerate(order[:top_n], start=1):
            percent = float(pct[row_idx, local_idx])
            if percent <= 0.0:
                continue
            ic = basis_internals[int(local_idx)]
            rows.append(
                {
                    "Source": source_label,
                    "Filename": result.filename,
                    "mode": int(result.orca_mode_indices[row_idx]),
                    "frequency_cm-1": float(result.orca_frequencies_cm1[row_idx]),
                    "gf_eigenvector_index": int(row_idx),
                    "veda_like_rank": rank,
                    "coord_index": basis_idx[int(local_idx)],
                    "internal_coordinate": _compact_coord_label(ic.name),
                    "coordinate_kind": str(ic.kind),
                    "coordinate_family": _assignment_family_from_internal(ic),
                    "coordinate_class": _stage3d_coord_class(ic),
                    "signed_ped_fraction": float(signed[row_idx, local_idx] / np.sum(np.abs(signed[row_idx, :])))
                    if np.sum(np.abs(signed[row_idx, :])) > tol
                    else 0.0,
                    "contribution_percent": percent,
                    "normalization_sum_percent": normalization,
                    "matrix_orientation": VEDA_LIKE_PED_MATRIX_ORIENTATION,
                    "basis_size": result.internal_basis_size,
                    "validation_status": result.validation_status,
                    "max_relative_error": result.max_relative_error,
                    "warnings": row_warnings,
                    "method": VEDA_LIKE_PED_METHOD,
                }
            )
    return pd.DataFrame(rows)


def build_veda_like_ped_matrix_dataframe(
    result: WilsonGFResult,
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    source_label: str = "",
    tol: float = 1.0e-12,
) -> pd.DataFrame:
    basis_idx, basis_internals, signed, pct = _veda_like_signed_terms(
        result,
        hess,
        internals,
        B,
        selected_idx,
        tol=tol,
    )
    rows = []
    warnings = "; ".join(result.warnings)
    for row_idx in range(pct.shape[0]):
        row_warnings = _veda_like_mode_warning(
            float(result.orca_frequencies_cm1[row_idx]),
            pct[row_idx, :],
            basis_internals,
            warnings,
        )
        denominator = float(np.sum(np.abs(signed[row_idx, :])))
        for local_idx, ic in enumerate(basis_internals):
            rows.append(
                {
                    "Source": source_label,
                    "Filename": result.filename,
                    "mode": int(result.orca_mode_indices[row_idx]),
                    "frequency_cm-1": float(result.orca_frequencies_cm1[row_idx]),
                    "gf_eigenvector_index": int(row_idx),
                    "coord_index": basis_idx[local_idx],
                    "internal_coordinate": _compact_coord_label(ic.name),
                    "coordinate_kind": str(ic.kind),
                    "coordinate_family": _assignment_family_from_internal(ic),
                    "coordinate_class": _stage3d_coord_class(ic),
                    "signed_ped_fraction": float(signed[row_idx, local_idx] / denominator) if denominator > tol else 0.0,
                    "contribution_percent": float(pct[row_idx, local_idx]),
                    "normalization_sum_percent": float(np.sum(pct[row_idx, :])),
                    "matrix_orientation": VEDA_LIKE_PED_MATRIX_ORIENTATION,
                    "basis_size": result.internal_basis_size,
                    "validation_status": result.validation_status,
                    "warnings": row_warnings,
                    "method": VEDA_LIKE_PED_METHOD,
                }
            )
    return pd.DataFrame(rows)


def build_veda_like_mode_correspondence_dataframe(
    result: WilsonGFResult,
    *,
    source_label: str = "",
) -> pd.DataFrame:
    rows = []
    warnings = "; ".join(result.warnings)
    for row_idx in range(len(result.orca_frequencies_cm1)):
        rows.append(
            {
                "Source": source_label,
                "Filename": result.filename,
                "mode": int(result.orca_mode_indices[row_idx]),
                "orca_frequency_cm-1": float(result.orca_frequencies_cm1[row_idx]),
                "gf_eigenvector_index": int(row_idx),
                "gf_eigenvalue": float(result.gf_eigenvalues[row_idx]),
                "reconstructed_frequency_cm-1": float(result.reconstructed_frequencies_cm1[row_idx]),
                "fixed_conversion_relative_error": (
                    abs(float(result.reconstructed_frequencies_cm1[row_idx]) - float(result.orca_frequencies_cm1[row_idx]))
                    / max(abs(float(result.orca_frequencies_cm1[row_idx])), 1.0e-12)
                ),
                "mapping_method": result.mapping_method,
                "conversion_method": result.conversion_method,
                "validation_status": result.validation_status,
                "max_relative_error": result.max_relative_error,
                "warnings": warnings,
                "method": VEDA_LIKE_PED_METHOD,
            }
        )
    if rows:
        return pd.DataFrame(rows)
    return pd.DataFrame(
        [
            {
                "Source": source_label,
                "Filename": result.filename,
                "mode": "",
                "orca_frequency_cm-1": "",
                "gf_eigenvector_index": "",
                "gf_eigenvalue": "",
                "reconstructed_frequency_cm-1": "",
                "fixed_conversion_relative_error": "",
                "mapping_method": result.mapping_method,
                "conversion_method": result.conversion_method,
                "validation_status": result.validation_status,
                "max_relative_error": result.max_relative_error,
                "warnings": warnings,
                "method": VEDA_LIKE_PED_METHOD,
            }
        ]
    )


def build_veda_like_basis_diagnostics_dataframe(
    result: WilsonGFResult,
    *,
    source_label: str = "",
) -> pd.DataFrame:
    basis = build_wilson_gf_basis_diagnostics_dataframe(result).copy()
    if basis.empty:
        return basis
    basis["Source"] = source_label
    basis["basis_scope"] = "veda_like_closed_wilson_gf_ped"
    basis["matrix_orientation"] = VEDA_LIKE_PED_MATRIX_ORIENTATION
    basis["validation_status"] = result.validation_status
    basis["method"] = VEDA_LIKE_PED_METHOD
    basis["method_boundary"] = "not VEDA-equivalent; original VEDA reference outputs not compared"
    return basis


def build_veda_like_metadata(
    result: WilsonGFResult,
    *,
    source_label: str = "",
) -> dict:
    return {
        "Source": source_label,
        "Filename": result.filename,
        "method": VEDA_LIKE_PED_METHOD,
        "method_boundary": "comparable VEDA-like closed Wilson GF/PED audit; does not reproduce original VEDA",
        "forbidden_claims": ["VEDA-equivalent", "original VEDA reproduced", "strict VEDA PED"],
        "normal_mode_orientation_rule": "normal_modes[:, mode]",
        "matrix_orientation": VEDA_LIKE_PED_MATRIX_ORIENTATION,
        "coordinate_basis": "Wilson GF validation selected nonredundant internal-coordinate basis",
        "basis_size": result.internal_basis_size,
        "expected_vibrational_rank": result.expected_vibrational_rank,
        "selected_indices": list(result.basis_indices),
        "validation_status": result.validation_status,
        "warnings": list(result.warnings),
        "epm_optimized": result.epm_optimized,
        "epm_swaps": result.epm_swaps,
        "epm_initial_localization_score": result.epm_initial_localization_score,
        "epm_optimized_localization_score": result.epm_optimized_localization_score,
        "epm_initial_mean_top_percent": result.epm_initial_mean_top_percent,
        "epm_optimized_mean_top_percent": result.epm_optimized_mean_top_percent,
        "epm_initial_diffuse_mode_fraction": result.epm_initial_diffuse_mode_fraction,
        "epm_optimized_diffuse_mode_fraction": result.epm_optimized_diffuse_mode_fraction,
        "mapping_method": result.mapping_method,
        "conversion_method": result.conversion_method,
        "max_relative_error": result.max_relative_error,
        "units": {
            "frequency": "cm-1",
            "geometry": "Angstrom",
            "contribution_percent": "dimensionless percent normalized per mode",
        },
    }
