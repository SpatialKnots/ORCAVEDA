from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd

from mode_assignment import _assignment_family_from_internal, _compact_coord_label, _stage3d_coord_class
from orcaveda_models import HessData, InternalCoordinate
from orca_parser import BOHR_TO_ANGSTROM
from ped import build_wilson_g_matrix, wilson_coordinate_scales


WILSON_GF_VALIDATION_METHOD = (
    "Wilson GF diagonalization validation prototype using selected nonredundant "
    "internal-coordinate basis; diagnostic only, not VEDA-equivalent PED"
)

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
    validation_status: str
    warnings: Tuple[str, ...]
    mapping_method: str
    conversion_method: str


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


def wilson_gf_diagonalization(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    *,
    tol: float = 1.0e-12,
    frequency_tol_relative: float = 1.0e-4,
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

    basis_idx = tuple(int(idx) for idx in selected_idx)
    if any(idx < 0 or idx >= len(internals) for idx in basis_idx):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")
    if not basis_idx:
        raise ValueError("selected_idx must contain at least one internal coordinate")

    expected_vibrational_rank = max(cartesian_size - 6, 0)
    basis_internals = [internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = B_arr[list(basis_idx), :] * scales[:, None]
    G = build_wilson_g_matrix(B_internal, hess.masses)
    F_internal = reconstruct_wilson_gf_internal_force_matrix(B_internal, hess.masses, hess.cartesian_hessian, G)
    g_rank, g_condition = _rank_condition(G, tol)
    f_rank, f_condition = _rank_condition(F_internal, tol)

    if len(basis_idx) != expected_vibrational_rank:
        warnings.append("basis_size_mismatch_expected_vibrational_rank")
    if g_rank < min(len(basis_idx), expected_vibrational_rank):
        warnings.append("basis_rank_below_expected")
    if not np.isfinite(g_condition) or g_condition > 1.0e12:
        warnings.append("g_ill_conditioned")
    if not np.isfinite(f_condition) or f_condition > 1.0e12:
        warnings.append("f_ill_conditioned")

    eigenvalues, eigenvectors, _ = solve_symmetric_gf_eigenproblem(G, F_internal, tol=tol)
    positive_mask = np.isfinite(eigenvalues) & (eigenvalues > tol)
    positive_gf = eigenvalues[positive_mask]
    positive_eigenvectors = eigenvectors[:, positive_mask]
    positive_freq_mask = np.isfinite(hess.frequencies_cm1) & (hess.frequencies_cm1 > tol)
    orca_mode_indices = np.where(positive_freq_mask)[0].astype(int)
    orca_freqs = np.asarray(hess.frequencies_cm1[positive_freq_mask], dtype=float)
    gf_order = np.argsort(positive_gf)
    positive_gf = positive_gf[gf_order]
    positive_eigenvectors = positive_eigenvectors[:, gf_order]
    orca_order = np.argsort(orca_freqs)
    orca_mode_indices = orca_mode_indices[orca_order]
    orca_freqs = orca_freqs[orca_order]

    pair_count = min(len(positive_gf), len(orca_freqs))
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

    fixed_conversion_failed = not np.isfinite(max_relative_error) or max_relative_error > frequency_tol_relative
    if fixed_conversion_failed:
        warnings.append("fixed_conversion_failed")
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
        validation_status=status,
        warnings=tuple(dict.fromkeys(warnings)),
        mapping_method="sorted_positive_gf_eigenvalues_to_sorted_positive_orca_frequencies",
        conversion_method=f"fixed_SI_hartree_per_amu_angstrom2_to_cm-1:{WILSON_GF_FIXED_CONVERSION_CM1:.12g}",
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
                "basis_size": result.internal_basis_size,
                "expected_vibrational_rank": result.expected_vibrational_rank,
                "g_rank": result.g_rank,
                "g_condition": result.g_condition,
                "f_rank": result.f_rank,
                "f_condition": result.f_condition,
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
