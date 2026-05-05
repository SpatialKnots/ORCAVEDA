from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from chemistry import build_connectivity, split_fragments
from orca_parser import read_orca_hess
from orcaveda_models import HessData
from reports import output_prefix_for_hess_paths


def mode_matrix_for_tracking(hess: HessData) -> np.ndarray:
    nm = np.asarray(hess.normal_modes, dtype=float)
    nat = len(hess.atoms)
    if nm.shape[0] != 3 * nat:
        raise ValueError(f"{hess.filename}: normal_modes first dimension {nm.shape[0]} != 3N {3 * nat}")
    return np.stack([nm[:, i].reshape(nat, 3) for i in range(nm.shape[1])], axis=0)


def normalize_tracking_vector(vec: np.ndarray, masses: Optional[np.ndarray] = None, mass_weighted: bool = False) -> np.ndarray:
    arr = np.asarray(vec, dtype=float).copy()
    if mass_weighted:
        if masses is None:
            raise ValueError("Masses are required for mass-weighted mode tracking.")
        w = np.sqrt(np.repeat(np.asarray(masses, dtype=float), 3)).reshape(arr.shape)
        arr = arr * w
    norm = float(np.linalg.norm(arr.ravel()))
    if norm <= 0.0:
        return arr
    return arr / norm


def kabsch_rotation(reference_coords_A: np.ndarray, target_coords_A: np.ndarray) -> np.ndarray:
    ref = np.asarray(reference_coords_A, dtype=float)
    tgt = np.asarray(target_coords_A, dtype=float)
    ref_c = ref - ref.mean(axis=0)
    tgt_c = tgt - tgt.mean(axis=0)
    cov = tgt_c.T @ ref_c
    u, _s, vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(u @ vt))
    corr = np.diag([1.0, 1.0, d])
    return u @ corr @ vt


def same_size_tracking_compatible(reference: HessData, target: HessData, require_same_symbols: bool = True) -> Tuple[bool, str]:
    if len(reference.atoms) != len(target.atoms):
        return False, f"Different atom counts: {len(reference.atoms)} vs {len(target.atoms)}"
    if require_same_symbols and list(reference.atoms) != list(target.atoms):
        return False, "Atom symbols differ or atom order differs."
    if reference.normal_modes.shape != target.normal_modes.shape:
        return False, f"Normal-mode matrix shapes differ: {reference.normal_modes.shape} vs {target.normal_modes.shape}"
    return True, "compatible"


def compute_mode_overlap_matrix(
    reference: HessData,
    target: HessData,
    *,
    mass_weighted: bool = True,
    align: bool = True,
    require_same_symbols: bool = True,
) -> np.ndarray:
    ok, reason = same_size_tracking_compatible(reference, target, require_same_symbols=require_same_symbols)
    if not ok:
        raise ValueError(f"Mode tracking incompatible: {reason}")

    ref_modes = mode_matrix_for_tracking(reference)
    tgt_modes = mode_matrix_for_tracking(target)
    if align:
        rot = kabsch_rotation(reference.coords_A, target.coords_A)
        tgt_modes = np.einsum("mni,ij->mnj", tgt_modes, rot)

    nref = ref_modes.shape[0]
    ntgt = tgt_modes.shape[0]
    overlap = np.zeros((nref, ntgt), dtype=float)
    for i in range(nref):
        qi = normalize_tracking_vector(ref_modes[i], reference.masses, mass_weighted).ravel()
        for j in range(ntgt):
            qj = normalize_tracking_vector(tgt_modes[j], target.masses, mass_weighted).ravel()
            overlap[i, j] = abs(float(np.dot(qi, qj)))
    return overlap


def mode_overlap_matrix_table(reference: HessData, target: HessData, overlap: np.ndarray, pair_label: str) -> pd.DataFrame:
    rows = []
    for i in range(overlap.shape[0]):
        for j in range(overlap.shape[1]):
            rows.append({
                "pair": pair_label,
                "reference_file": reference.filename,
                "target_file": target.filename,
                "reference_mode": i,
                "target_mode": j,
                "reference_frequency_cm-1": float(reference.frequencies_cm1[i]) if i < len(reference.frequencies_cm1) else np.nan,
                "target_frequency_cm-1": float(target.frequencies_cm1[j]) if j < len(target.frequencies_cm1) else np.nan,
                "overlap": float(overlap[i, j]),
            })
    return pd.DataFrame(rows)


def mode_tracking_table(reference: HessData, target: HessData, overlap: np.ndarray, pair_label: str) -> pd.DataFrame:
    best_target_for_ref = np.argmax(overlap, axis=1)
    best_ref_for_target = np.argmax(overlap, axis=0)
    rows = []
    for i, j in enumerate(best_target_for_ref):
        sorted_j = np.argsort(overlap[i])[::-1]
        second = sorted_j[1] if len(sorted_j) > 1 else j
        rows.append({
            "pair": pair_label,
            "reference_file": reference.filename,
            "target_file": target.filename,
            "reference_mode": i,
            "target_mode": int(j),
            "reference_frequency_cm-1": float(reference.frequencies_cm1[i]) if i < len(reference.frequencies_cm1) else np.nan,
            "target_frequency_cm-1": float(target.frequencies_cm1[j]) if j < len(target.frequencies_cm1) else np.nan,
            "frequency_shift_cm-1": float(target.frequencies_cm1[j] - reference.frequencies_cm1[i]) if i < len(reference.frequencies_cm1) and j < len(target.frequencies_cm1) else np.nan,
            "best_overlap": float(overlap[i, j]),
            "second_best_target_mode": int(second),
            "second_best_overlap": float(overlap[i, second]),
            "overlap_gap": float(overlap[i, j] - overlap[i, second]),
            "reciprocal_best_match": bool(best_ref_for_target[j] == i),
            "tracking_confidence": "high" if overlap[i, j] >= 0.75 and (overlap[i, j] - overlap[i, second]) >= 0.15 else ("medium" if overlap[i, j] >= 0.50 else "low"),
        })
    return pd.DataFrame(rows)


def mode_mixing_warnings_table(reference: HessData, target: HessData, overlap: np.ndarray, pair_label: str) -> pd.DataFrame:
    rows = []
    for i in range(overlap.shape[0]):
        candidates = np.where(overlap[i] >= 0.35)[0]
        if len(candidates) >= 2:
            top = candidates[np.argsort(overlap[i, candidates])[::-1]]
            rows.append({
                "pair": pair_label,
                "reference_file": reference.filename,
                "target_file": target.filename,
                "reference_mode": i,
                "reference_frequency_cm-1": float(reference.frequencies_cm1[i]) if i < len(reference.frequencies_cm1) else np.nan,
                "warning_type": "reference_mode_splits_or_mixes",
                "candidate_target_modes": ";".join(map(str, top.tolist())),
                "candidate_overlaps": ";".join(f"{overlap[i, j]:.4f}" for j in top),
            })
    for j in range(overlap.shape[1]):
        candidates = np.where(overlap[:, j] >= 0.35)[0]
        if len(candidates) >= 2:
            top = candidates[np.argsort(overlap[candidates, j])[::-1]]
            rows.append({
                "pair": pair_label,
                "reference_file": reference.filename,
                "target_file": target.filename,
                "target_mode": j,
                "target_frequency_cm-1": float(target.frequencies_cm1[j]) if j < len(target.frequencies_cm1) else np.nan,
                "warning_type": "target_mode_combines_reference_modes",
                "candidate_reference_modes": ";".join(map(str, top.tolist())),
                "candidate_overlaps": ";".join(f"{overlap[i, j]:.4f}" for i in top),
            })
    return pd.DataFrame(rows)


def fragment_same_symbol_order(reference: HessData, target: HessData, fragment: Sequence[int]) -> bool:
    if len(fragment) != len(reference.atoms):
        return False
    return [target.atoms[i] for i in fragment] == list(reference.atoms)


def fragment_projected_overlap_matrix(
    reference: HessData,
    target: HessData,
    fragment: Sequence[int],
    *,
    mass_weighted: bool = True,
    align: bool = True,
) -> np.ndarray:
    fragment = list(fragment)
    if not fragment_same_symbol_order(reference, target, fragment):
        raise ValueError("Target fragment is not compatible with reference atom count/symbol order.")

    ref_modes = mode_matrix_for_tracking(reference)
    tgt_modes_all = mode_matrix_for_tracking(target)
    tgt_modes = tgt_modes_all[:, fragment, :]

    ref_coords = reference.coords_A
    tgt_coords = target.coords_A[fragment, :]
    if align:
        rot = kabsch_rotation(ref_coords, tgt_coords)
        tgt_modes = np.einsum("mni,ij->mnj", tgt_modes, rot)

    target_masses_fragment = target.masses[fragment]
    overlap = np.zeros((ref_modes.shape[0], tgt_modes.shape[0]), dtype=float)
    for i in range(ref_modes.shape[0]):
        qi = normalize_tracking_vector(ref_modes[i], reference.masses, mass_weighted).ravel()
        for j in range(tgt_modes.shape[0]):
            qj = normalize_tracking_vector(tgt_modes[j], target_masses_fragment, mass_weighted).ravel()
            overlap[i, j] = abs(float(np.dot(qi, qj)))
    return overlap


def mode_tracking_outputs_for_hess_files(
    paths: Sequence[str | Path],
    outdir: str | Path,
    *,
    mass_weighted: bool = True,
    align: bool = True,
    include_overlap_matrices: bool = False,
) -> Dict[str, pd.DataFrame]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_prefix_for_hess_paths(paths)
    hess_list = [read_orca_hess(path) for path in paths]

    tracking_frames: List[pd.DataFrame] = []
    warning_frames: List[pd.DataFrame] = []
    overlap_frames: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, object]] = []

    for a_idx, ref in enumerate(hess_list):
        for b_idx, tgt in enumerate(hess_list):
            if a_idx == b_idx:
                continue
            pair_label = f"{Path(ref.filename).stem}__to__{Path(tgt.filename).stem}"
            ok, reason = same_size_tracking_compatible(ref, tgt, require_same_symbols=True)
            if ok and a_idx < b_idx:
                try:
                    ov = compute_mode_overlap_matrix(ref, tgt, mass_weighted=mass_weighted, align=align)
                    tr = mode_tracking_table(ref, tgt, ov, pair_label)
                    warn = mode_mixing_warnings_table(ref, tgt, ov, pair_label)
                    tracking_frames.append(tr.assign(tracking_type="same_size"))
                    if not warn.empty:
                        warning_frames.append(warn.assign(tracking_type="same_size"))
                    if include_overlap_matrices:
                        overlap_frames.append(mode_overlap_matrix_table(ref, tgt, ov, pair_label).assign(tracking_type="same_size"))
                    summary_rows.append({
                        "pair": pair_label,
                        "tracking_type": "same_size",
                        "reference_file": ref.filename,
                        "target_file": tgt.filename,
                        "n_reference_modes": ov.shape[0],
                        "n_target_modes": ov.shape[1],
                        "median_best_overlap": float(np.median(np.max(ov, axis=1))),
                        "low_confidence_count": int((tr["tracking_confidence"] == "low").sum()),
                        "medium_confidence_count": int((tr["tracking_confidence"] == "medium").sum()),
                        "high_confidence_count": int((tr["tracking_confidence"] == "high").sum()),
                        "status": "OK",
                        "notes": "",
                    })
                except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                    summary_rows.append({"pair": pair_label, "tracking_type": "same_size", "reference_file": ref.filename, "target_file": tgt.filename, "status": "ERROR", "notes": str(exc)})
            elif a_idx < b_idx:
                summary_rows.append({"pair": pair_label, "tracking_type": "same_size", "reference_file": ref.filename, "target_file": tgt.filename, "status": "SKIPPED", "notes": reason})

            if len(ref.atoms) < len(tgt.atoms):
                try:
                    bonds = build_connectivity(tgt.atoms, tgt.coords_A)
                    frags = split_fragments(len(tgt.atoms), bonds)
                    compatible_frags = [frag for frag in frags if fragment_same_symbol_order(ref, tgt, frag)]
                    if not compatible_frags:
                        summary_rows.append({"pair": pair_label, "tracking_type": "fragment_projected", "reference_file": ref.filename, "target_file": tgt.filename, "status": "SKIPPED", "notes": "No target fragment with same atom count/symbol order."})
                    for frag_idx, frag in enumerate(compatible_frags):
                        ov = fragment_projected_overlap_matrix(ref, tgt, frag, mass_weighted=mass_weighted, align=align)
                        frag_label = f"{pair_label}__fragment_{frag_idx}"
                        tr = mode_tracking_table(ref, tgt, ov, frag_label)
                        warn = mode_mixing_warnings_table(ref, tgt, ov, frag_label)
                        tracking_frames.append(tr.assign(tracking_type="fragment_projected", target_fragment_index=frag_idx, target_fragment_atoms=";".join(map(str, frag))))
                        if not warn.empty:
                            warning_frames.append(warn.assign(tracking_type="fragment_projected", target_fragment_index=frag_idx, target_fragment_atoms=";".join(map(str, frag))))
                        if include_overlap_matrices:
                            overlap_frames.append(mode_overlap_matrix_table(ref, tgt, ov, frag_label).assign(tracking_type="fragment_projected", target_fragment_index=frag_idx, target_fragment_atoms=";".join(map(str, frag))))
                        summary_rows.append({"pair": frag_label, "tracking_type": "fragment_projected", "reference_file": ref.filename, "target_file": tgt.filename, "target_fragment_index": frag_idx, "n_reference_modes": ov.shape[0], "n_target_modes": ov.shape[1], "median_best_overlap": float(np.median(np.max(ov, axis=1))), "status": "OK", "notes": ""})
                except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                    summary_rows.append({"pair": pair_label, "tracking_type": "fragment_projected", "reference_file": ref.filename, "target_file": tgt.filename, "status": "ERROR", "notes": str(exc)})

    tables: Dict[str, pd.DataFrame] = {
        "mode_tracking_summary": pd.DataFrame(summary_rows),
        "mode_tracking": pd.concat(tracking_frames, ignore_index=True) if tracking_frames else pd.DataFrame(),
        "mode_mixing_warnings": pd.concat(warning_frames, ignore_index=True) if warning_frames else pd.DataFrame(),
    }
    if include_overlap_matrices:
        tables["mode_overlap_matrix"] = pd.concat(overlap_frames, ignore_index=True) if overlap_frames else pd.DataFrame()

    for name, df in tables.items():
        df.to_csv(outdir / f"{output_prefix}__{name}.csv", index=False)
    return tables
