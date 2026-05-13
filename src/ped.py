from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mode_assignment import _assignment_family_from_internal, _compact_coord_label, _stage3d_coord_class
from orcaveda_models import HessData, InternalCoordinate


PED_V1_METHOD = (
    "PED v1 normalized B-matrix internal-coordinate projection; "
    "not force-constant Wilson GF PED"
)

PED_V2_METHOD = (
    "PED v2 force-aware normalized B-matrix projection using ORCA $hessian; "
    "mode-specific force-coupled internal-coordinate energy diagnostic, not full Wilson GF PED"
)

WILSON_PED_METHOD = (
    "Wilson GF-style PED audit using selected nonredundant internal-coordinate B matrix, "
    "G = B M^-1 B^T, reconstructed internal F matrix, and mode-projected potential-energy terms"
)

BOHR_TO_ANGSTROM = 0.529177210903


@dataclass
class PEDContribution:
    coord_index: int
    internal_coordinate: str
    coordinate_kind: str
    coordinate_family: str
    coordinate_class: str
    atoms0: Tuple[int, ...]
    source: str
    generation_rule: str
    signed_projection: float
    weight: float
    percent: float


@dataclass
class PEDModeResult:
    source_label: str
    filename: str
    mode: int
    frequency_cm1: float
    ir_intensity: float
    basis_size: int
    valid_coordinate_count: int
    total_weight: float
    normalization_sum_percent: float
    top1_percent: float
    warnings: Tuple[str, ...]
    contributions: Tuple[PEDContribution, ...]


def _normalize_rows(matrix: np.ndarray, tol: float) -> Tuple[np.ndarray, np.ndarray]:
    row_norms = np.linalg.norm(matrix, axis=1)
    valid = np.isfinite(row_norms) & (row_norms > tol)
    normalized = np.zeros_like(matrix, dtype=float)
    normalized[valid, :] = matrix[valid, :] / row_norms[valid, None]
    return normalized, valid


def compute_ped(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Optional[Sequence[int]] = None,
    *,
    source_label: str = "",
    top_n: int = 8,
    tol: float = 1.0e-12,
) -> List[PEDModeResult]:
    """
    Compute PED v1 as normalized B-matrix projections onto normal modes.

    This deliberately uses the same normal-mode orientation rule as Stage 3D:
        normal_mode_vector = normal_modes[:, mode]

    The result is an unweighted internal-coordinate projection audit. It does
    not use force constants and must not be described as full Wilson GF PED.
    """
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    B_arr = np.asarray(B, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(internals):
        raise ValueError(f"B row count {B_arr.shape[0]} does not match internal coordinate count {len(internals)}")
    if B_arr.shape[1] != hess.normal_modes.shape[0]:
        raise ValueError(
            f"B column count {B_arr.shape[1]} does not match normal-mode vector length {hess.normal_modes.shape[0]}"
        )

    if selected_idx is None:
        basis_idx = list(range(len(internals)))
    else:
        basis_idx = [int(idx) for idx in selected_idx]
    if any(idx < 0 or idx >= len(internals) for idx in basis_idx):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")

    if basis_idx:
        basis_B = B_arr[basis_idx, :]
        basis_internals = [internals[idx] for idx in basis_idx]
    else:
        basis_B = np.zeros((0, B_arr.shape[1]), dtype=float)
        basis_internals = []

    B_unit, valid_rows = _normalize_rows(basis_B, tol)
    results: List[PEDModeResult] = []

    for mode, freq in enumerate(hess.frequencies_cm1):
        warnings: List[str] = []
        mode_vec = np.asarray(hess.normal_modes[:, mode], dtype=float)
        mode_norm = float(np.linalg.norm(mode_vec))
        if not np.isfinite(mode_norm) or mode_norm <= tol:
            mode_unit = np.zeros_like(mode_vec)
            warnings.append("zero_or_invalid_normal_mode_vector")
        else:
            mode_unit = mode_vec / mode_norm

        if not basis_idx:
            warnings.append("no_internal_coordinate_basis")
            projections = np.zeros(0, dtype=float)
            weights = np.zeros(0, dtype=float)
        else:
            if not np.any(valid_rows):
                warnings.append("no_valid_internal_coordinate_rows")
            projections = B_unit @ mode_unit
            weights = projections ** 2
            weights[~valid_rows] = 0.0
            weights[~np.isfinite(weights)] = 0.0

        total = float(np.sum(weights))
        if total > tol and np.isfinite(total):
            pct = 100.0 * weights / total
        else:
            pct = np.zeros_like(weights)
            if basis_idx:
                warnings.append("zero_ped_projection_weight")

        order = np.argsort(pct)[::-1] if pct.size else np.array([], dtype=int)
        contributions: List[PEDContribution] = []
        for local_idx in order[:top_n]:
            percent = float(pct[local_idx])
            if percent <= 0.0:
                continue
            ic = basis_internals[int(local_idx)]
            coord_index = basis_idx[int(local_idx)]
            contributions.append(
                PEDContribution(
                    coord_index=coord_index,
                    internal_coordinate=_compact_coord_label(ic.name),
                    coordinate_kind=str(ic.kind),
                    coordinate_family=_assignment_family_from_internal(ic),
                    coordinate_class=_stage3d_coord_class(ic),
                    atoms0=tuple(int(a) for a in ic.atoms0),
                    source=str(ic.source),
                    generation_rule=str(ic.generation_rule or ""),
                    signed_projection=float(projections[local_idx]),
                    weight=float(weights[local_idx]),
                    percent=percent,
                )
            )

        top1_percent = float(contributions[0].percent) if contributions else 0.0
        if contributions and top1_percent < 25.0:
            warnings.append("diffuse_ped_contributions")

        results.append(
            PEDModeResult(
                source_label=source_label,
                filename=hess.filename,
                mode=int(mode),
                frequency_cm1=float(freq),
                ir_intensity=float(hess.ir_intensities[mode]),
                basis_size=len(basis_idx),
                valid_coordinate_count=int(np.sum(valid_rows)),
                total_weight=total,
                normalization_sum_percent=float(np.sum(pct)),
                top1_percent=top1_percent,
                warnings=tuple(warnings),
                contributions=tuple(contributions),
            )
        )

    return results


def compute_ped_v2_force_aware(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Optional[Sequence[int]] = None,
    *,
    source_label: str = "",
    top_n: int = 8,
    tol: float = 1.0e-12,
) -> List[PEDModeResult]:
    """
    Compute PED v2 as a force-aware internal-coordinate diagnostic.

    This uses the parsed ORCA Cartesian Hessian and row-normalized B matrix:
        K_int = B_unit @ F_cart @ B_unit.T
        p = B_unit @ mode_unit
        signed_weight_i = p_i * (K_int @ p)_i

    Percentages are normalized from absolute signed weights so coupling terms
    remain visible without allowing cancellation to hide contributions.

    This is stricter than PED v1 because force constants affect the reported
    contributors. It is still not a full Wilson GF PED implementation: it does
    not construct a complete internal-coordinate G/F treatment and does not
    claim VEDA equivalence.
    """
    if hess.cartesian_hessian is None:
        raise ValueError("PED v2 requires HessData.cartesian_hessian parsed from ORCA $hessian")
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    B_arr = np.asarray(B, dtype=float)
    F_cart = np.asarray(hess.cartesian_hessian, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(internals):
        raise ValueError(f"B row count {B_arr.shape[0]} does not match internal coordinate count {len(internals)}")
    if B_arr.shape[1] != hess.normal_modes.shape[0]:
        raise ValueError(
            f"B column count {B_arr.shape[1]} does not match normal-mode vector length {hess.normal_modes.shape[0]}"
        )
    if F_cart.shape != (B_arr.shape[1], B_arr.shape[1]):
        raise ValueError(f"cartesian_hessian shape {F_cart.shape} does not match B column count {B_arr.shape[1]}")

    F_cart = 0.5 * (F_cart + F_cart.T)
    if selected_idx is None:
        basis_idx = list(range(len(internals)))
    else:
        basis_idx = [int(idx) for idx in selected_idx]
    if any(idx < 0 or idx >= len(internals) for idx in basis_idx):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")

    if basis_idx:
        basis_B = B_arr[basis_idx, :]
        basis_internals = [internals[idx] for idx in basis_idx]
    else:
        basis_B = np.zeros((0, B_arr.shape[1]), dtype=float)
        basis_internals = []

    B_unit, valid_rows = _normalize_rows(basis_B, tol)
    K_internal = B_unit @ F_cart @ B_unit.T if basis_idx else np.zeros((0, 0), dtype=float)
    K_internal[~np.isfinite(K_internal)] = 0.0
    results: List[PEDModeResult] = []

    for mode, freq in enumerate(hess.frequencies_cm1):
        warnings: List[str] = []
        mode_vec = np.asarray(hess.normal_modes[:, mode], dtype=float)
        mode_norm = float(np.linalg.norm(mode_vec))
        if not np.isfinite(mode_norm) or mode_norm <= tol:
            mode_unit = np.zeros_like(mode_vec)
            warnings.append("zero_or_invalid_normal_mode_vector")
        else:
            mode_unit = mode_vec / mode_norm

        if not basis_idx:
            warnings.append("no_internal_coordinate_basis")
            projections = np.zeros(0, dtype=float)
            signed_weights = np.zeros(0, dtype=float)
        else:
            if not np.any(valid_rows):
                warnings.append("no_valid_internal_coordinate_rows")
            projections = B_unit @ mode_unit
            force_response = K_internal @ projections
            signed_weights = projections * force_response
            signed_weights[~valid_rows] = 0.0
            signed_weights[~np.isfinite(signed_weights)] = 0.0

        weights = np.abs(signed_weights)
        total = float(np.sum(weights))
        if total > tol and np.isfinite(total):
            pct = 100.0 * weights / total
        else:
            pct = np.zeros_like(weights)
            if basis_idx:
                warnings.append("zero_ped_v2_force_projection_weight")

        order = np.argsort(pct)[::-1] if pct.size else np.array([], dtype=int)
        contributions: List[PEDContribution] = []
        for local_idx in order[:top_n]:
            percent = float(pct[local_idx])
            if percent <= 0.0:
                continue
            ic = basis_internals[int(local_idx)]
            coord_index = basis_idx[int(local_idx)]
            contributions.append(
                PEDContribution(
                    coord_index=coord_index,
                    internal_coordinate=_compact_coord_label(ic.name),
                    coordinate_kind=str(ic.kind),
                    coordinate_family=_assignment_family_from_internal(ic),
                    coordinate_class=_stage3d_coord_class(ic),
                    atoms0=tuple(int(a) for a in ic.atoms0),
                    source=str(ic.source),
                    generation_rule=str(ic.generation_rule or ""),
                    signed_projection=float(projections[local_idx]),
                    weight=float(signed_weights[local_idx]),
                    percent=percent,
                )
            )

        top1_percent = float(contributions[0].percent) if contributions else 0.0
        if contributions and top1_percent < 25.0:
            warnings.append("diffuse_ped_v2_force_contributions")

        results.append(
            PEDModeResult(
                source_label=source_label,
                filename=hess.filename,
                mode=int(mode),
                frequency_cm1=float(freq),
                ir_intensity=float(hess.ir_intensities[mode]),
                basis_size=len(basis_idx),
                valid_coordinate_count=int(np.sum(valid_rows)),
                total_weight=total,
                normalization_sum_percent=float(np.sum(pct)),
                top1_percent=top1_percent,
                warnings=tuple(warnings),
                contributions=tuple(contributions),
            )
        )

    return results


def ped_results_to_dataframe(results: Sequence[PEDModeResult], *, method: str = PED_V1_METHOD) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for result in results:
        if result.contributions:
            for rank, contribution in enumerate(result.contributions, start=1):
                rows.append({
                    "Source": result.source_label,
                    "Filename": result.filename,
                    "mode": result.mode,
                    "frequency_cm-1": result.frequency_cm1,
                    "IR_intensity": result.ir_intensity,
                    "ped_rank": rank,
                    "coord_index": contribution.coord_index,
                    "internal_coordinate": contribution.internal_coordinate,
                    "coordinate_kind": contribution.coordinate_kind,
                    "coordinate_family": contribution.coordinate_family,
                    "coordinate_class": contribution.coordinate_class,
                    "atoms_1based": "-".join(str(a + 1) for a in contribution.atoms0),
                    "source": contribution.source,
                    "generation_rule": contribution.generation_rule,
                    "signed_projection": round(float(contribution.signed_projection), 8),
                    "weight": round(float(contribution.weight), 12),
                    "contribution_percent": round(float(contribution.percent), 6),
                    "basis_size": result.basis_size,
                    "valid_coordinate_count": result.valid_coordinate_count,
                    "normalization_sum_percent": round(float(result.normalization_sum_percent), 6),
                    "top1_percent": round(float(result.top1_percent), 6),
                    "ped_warnings": "; ".join(result.warnings),
                    "ped_method": method,
                })
        else:
            rows.append({
                "Source": result.source_label,
                "Filename": result.filename,
                "mode": result.mode,
                "frequency_cm-1": result.frequency_cm1,
                "IR_intensity": result.ir_intensity,
                "ped_rank": 0,
                "coord_index": "",
                "internal_coordinate": "",
                "coordinate_kind": "",
                "coordinate_family": "",
                "coordinate_class": "",
                "atoms_1based": "",
                "source": "",
                "generation_rule": "",
                "signed_projection": 0.0,
                "weight": 0.0,
                "contribution_percent": 0.0,
                "basis_size": result.basis_size,
                "valid_coordinate_count": result.valid_coordinate_count,
                "normalization_sum_percent": round(float(result.normalization_sum_percent), 6),
                "top1_percent": round(float(result.top1_percent), 6),
                "ped_warnings": "; ".join(result.warnings),
                "ped_method": method,
            })
    return pd.DataFrame(rows)


def build_ped_audit_dataframe(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Optional[Sequence[int]] = None,
    *,
    source_label: str = "",
    top_n: int = 8,
) -> pd.DataFrame:
    return ped_results_to_dataframe(
        compute_ped(
            hess,
            internals,
            B,
            selected_idx,
            source_label=source_label,
            top_n=top_n,
        )
    )


def build_ped_v2_force_audit_dataframe(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Optional[Sequence[int]] = None,
    *,
    source_label: str = "",
    top_n: int = 8,
) -> pd.DataFrame:
    return ped_results_to_dataframe(
        compute_ped_v2_force_aware(
            hess,
            internals,
            B,
            selected_idx,
            source_label=source_label,
            top_n=top_n,
        ),
        method=PED_V2_METHOD,
    )


def wilson_coordinate_scales(internals: Sequence[InternalCoordinate]) -> np.ndarray:
    scales = np.ones(len(internals), dtype=float)
    for idx, ic in enumerate(internals):
        kind = str(ic.kind).lower()
        name = str(ic.name).lower()
        if "bend" in kind or "angle" in kind or name.startswith("ang("):
            scales[idx] = np.pi / 180.0
    return scales


def build_wilson_g_matrix(B_internal: np.ndarray, masses: np.ndarray) -> np.ndarray:
    B_arr = np.asarray(B_internal, dtype=float)
    mass_vec = np.repeat(np.asarray(masses, dtype=float), 3)
    if B_arr.ndim != 2:
        raise ValueError(f"B_internal must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[1] != mass_vec.size:
        raise ValueError(f"B column count {B_arr.shape[1]} does not match 3N mass vector length {mass_vec.size}")
    if np.any(mass_vec <= 0.0) or not np.all(np.isfinite(mass_vec)):
        raise ValueError("masses must be finite and positive for Wilson G matrix")
    weighted = B_arr / mass_vec[None, :]
    G = weighted @ B_arr.T
    G = 0.5 * (G + G.T)
    G[~np.isfinite(G)] = 0.0
    return G


def reconstruct_internal_force_matrix(B_internal: np.ndarray, cartesian_hessian: np.ndarray) -> np.ndarray:
    B_arr = np.asarray(B_internal, dtype=float)
    F_cart = np.asarray(cartesian_hessian, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B_internal must be a 2D matrix, got shape {B_arr.shape}")
    if F_cart.shape != (B_arr.shape[1], B_arr.shape[1]):
        raise ValueError(f"cartesian_hessian shape {F_cart.shape} does not match B column count {B_arr.shape[1]}")
    # ORCA .hess Cartesian force constants are in Bohr-based Cartesian units;
    # the finite-difference B matrix is built against Angstrom coordinates.
    F_cart_A = 0.5 * (F_cart + F_cart.T) / (BOHR_TO_ANGSTROM ** 2)
    B_pinv = np.linalg.pinv(B_arr)
    F_internal = B_pinv.T @ F_cart_A @ B_pinv
    F_internal = 0.5 * (F_internal + F_internal.T)
    F_internal[~np.isfinite(F_internal)] = 0.0
    return F_internal


def build_wilson_ped_audit_dataframe(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Optional[Sequence[int]] = None,
    *,
    source_label: str = "",
    top_n: int = 8,
    tol: float = 1.0e-12,
) -> pd.DataFrame:
    if hess.cartesian_hessian is None:
        raise ValueError("Wilson PED requires HessData.cartesian_hessian parsed from ORCA $hessian")
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    B_arr = np.asarray(B, dtype=float)
    if B_arr.ndim != 2:
        raise ValueError(f"B must be a 2D matrix, got shape {B_arr.shape}")
    if B_arr.shape[0] != len(internals):
        raise ValueError(f"B row count {B_arr.shape[0]} does not match internal coordinate count {len(internals)}")
    if B_arr.shape[1] != hess.normal_modes.shape[0]:
        raise ValueError(
            f"B column count {B_arr.shape[1]} does not match normal-mode vector length {hess.normal_modes.shape[0]}"
        )

    if selected_idx is None:
        basis_idx = list(range(len(internals)))
    else:
        basis_idx = [int(idx) for idx in selected_idx]
    if any(idx < 0 or idx >= len(internals) for idx in basis_idx):
        raise ValueError("selected_idx contains an out-of-range internal coordinate index")

    if not basis_idx:
        return pd.DataFrame()

    basis_internals = [internals[idx] for idx in basis_idx]
    scales = wilson_coordinate_scales(basis_internals)
    B_internal = B_arr[basis_idx, :] * scales[:, None]
    valid_rows = np.linalg.norm(B_internal, axis=1) > tol
    G = build_wilson_g_matrix(B_internal, hess.masses)
    F_internal = reconstruct_internal_force_matrix(B_internal, hess.cartesian_hessian)
    gf = G @ F_internal
    gf_eigs = np.linalg.eigvals(gf) if gf.size else np.array([])
    g_singular_values = np.linalg.svd(G, compute_uv=False) if G.size else np.array([])
    f_singular_values = np.linalg.svd(F_internal, compute_uv=False) if F_internal.size else np.array([])
    g_rank = int(np.sum(g_singular_values > 1.0e-8))
    f_rank = int(np.sum(f_singular_values > 1.0e-8))
    g_condition = (
        float(g_singular_values[0] / g_singular_values[g_rank - 1])
        if g_rank > 0 else float("inf")
    )
    f_condition = (
        float(f_singular_values[0] / f_singular_values[f_rank - 1])
        if f_rank > 0 else float("inf")
    )

    rows: List[Dict[str, object]] = []
    for mode, freq in enumerate(hess.frequencies_cm1):
        warnings: List[str] = []
        mode_vec = np.asarray(hess.normal_modes[:, mode], dtype=float)
        mode_norm = float(np.linalg.norm(mode_vec))
        if not np.isfinite(mode_norm) or mode_norm <= tol:
            mode_unit = np.zeros_like(mode_vec)
            warnings.append("zero_or_invalid_normal_mode_vector")
        else:
            mode_unit = mode_vec / mode_norm

        internal_displacement = B_internal @ mode_unit
        force_response = F_internal @ internal_displacement
        signed_terms = internal_displacement * force_response
        signed_terms[~valid_rows] = 0.0
        signed_terms[~np.isfinite(signed_terms)] = 0.0
        weights = np.abs(signed_terms)
        total = float(np.sum(weights))
        if total > tol and np.isfinite(total):
            pct = 100.0 * weights / total
        else:
            pct = np.zeros_like(weights)
            warnings.append("zero_wilson_ped_energy_distribution")

        order = np.argsort(pct)[::-1] if pct.size else np.array([], dtype=int)
        if len(order) and float(pct[order[0]]) < 25.0:
            warnings.append("diffuse_wilson_ped_contributions")

        for rank, local_idx in enumerate(order[:top_n], start=1):
            percent = float(pct[local_idx])
            if percent <= 0.0:
                continue
            ic = basis_internals[int(local_idx)]
            rows.append(
                {
                    "Source": source_label,
                    "Filename": hess.filename,
                    "mode": int(mode),
                    "frequency_cm-1": float(freq),
                    "IR_intensity": float(hess.ir_intensities[mode]),
                    "wilson_rank": rank,
                    "coord_index": basis_idx[int(local_idx)],
                    "internal_coordinate": _compact_coord_label(ic.name),
                    "coordinate_kind": str(ic.kind),
                    "coordinate_family": _assignment_family_from_internal(ic),
                    "coordinate_class": _stage3d_coord_class(ic),
                    "atoms_1based": "-".join(str(a + 1) for a in ic.atoms0),
                    "source": str(ic.source),
                    "generation_rule": str(ic.generation_rule or ""),
                    "internal_displacement": round(float(internal_displacement[local_idx]), 10),
                    "force_response": round(float(force_response[local_idx]), 10),
                    "signed_potential_term": round(float(signed_terms[local_idx]), 12),
                    "contribution_percent": round(percent, 6),
                    "normalization_sum_percent": round(float(np.sum(pct)), 6),
                    "basis_size": len(basis_idx),
                    "valid_coordinate_count": int(np.sum(valid_rows)),
                    "wilson_g_rank": g_rank,
                    "wilson_f_rank": f_rank,
                    "wilson_g_condition": g_condition,
                    "wilson_f_condition": f_condition,
                    "gf_eigenvalue_min": float(np.min(np.real(gf_eigs))) if gf_eigs.size else "",
                    "gf_eigenvalue_max": float(np.max(np.real(gf_eigs))) if gf_eigs.size else "",
                    "wilson_ped_warnings": "; ".join(warnings),
                    "wilson_ped_method": WILSON_PED_METHOD,
                }
            )

        if not len(order) or not np.any(pct > 0.0):
            rows.append(
                {
                    "Source": source_label,
                    "Filename": hess.filename,
                    "mode": int(mode),
                    "frequency_cm-1": float(freq),
                    "IR_intensity": float(hess.ir_intensities[mode]),
                    "wilson_rank": 0,
                    "coord_index": "",
                    "internal_coordinate": "",
                    "coordinate_kind": "",
                    "coordinate_family": "",
                    "coordinate_class": "",
                    "atoms_1based": "",
                    "source": "",
                    "generation_rule": "",
                    "internal_displacement": 0.0,
                    "force_response": 0.0,
                    "signed_potential_term": 0.0,
                    "contribution_percent": 0.0,
                    "normalization_sum_percent": round(float(np.sum(pct)), 6),
                    "basis_size": len(basis_idx),
                    "valid_coordinate_count": int(np.sum(valid_rows)),
                    "wilson_g_rank": g_rank,
                    "wilson_f_rank": f_rank,
                    "wilson_g_condition": g_condition,
                    "wilson_f_condition": f_condition,
                    "gf_eigenvalue_min": float(np.min(np.real(gf_eigs))) if gf_eigs.size else "",
                    "gf_eigenvalue_max": float(np.max(np.real(gf_eigs))) if gf_eigs.size else "",
                    "wilson_ped_warnings": "; ".join(warnings),
                    "wilson_ped_method": WILSON_PED_METHOD,
                }
            )

    return pd.DataFrame(rows)
