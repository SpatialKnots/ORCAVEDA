# @title ORCAVEDA (V.5.0 / Stage 3D v5.0)
#!/usr/bin/env python3
"""
orca_ped_general_engine.py

Integrated Stage 3A/3C development module: general organic molecule engine for ORCA .hess files.

Scope:
    - ORCA 6 native .hess parsing
    - connectivity and fragment/system classification
    - graph-based atom environment typing
    - functional group detection for common small organic molecules
    - universal local-coordinate template generation
    - validation/confidence diagnostics suitable for monomers, dimers and clusters

Critical ORCA .hess convention:
    normal mode vector = normal_modes[:, mode]
not:
    normal_modes[mode, :]

This module is conservative: unsupported chemistry is reported explicitly rather than guessed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple, Optional
import math
import re
import json

import numpy as np
import pandas as pd

from chemistry import (
    SUPPORTED_ELEMENTS,
    adjacency,
    annotate_chemical_system as chemistry_annotate_chemical_system,
    atom_environment_table as chemistry_atom_environment_table,
    build_connectivity as chemistry_build_connectivity,
    classify_system as chemistry_classify_system,
    detect_functional_groups as chemistry_detect_functional_groups,
    detect_interfragment_hbonds as chemistry_detect_interfragment_hbonds,
    expected_vibrational_rank as chemistry_expected_vibrational_rank,
    formula_string as chemistry_formula_string,
    get_active_backend_name as chemistry_get_active_backend_name,
    list_backends as chemistry_list_backends,
    set_active_backend as chemistry_set_active_backend,
    set_active_backend_from_env as chemistry_set_active_backend_from_env,
    get_supported_elements as chemistry_get_supported_elements,
    split_fragments as chemistry_split_fragments,
)
from b_matrix import (
    finite_difference_B as bmatrix_finite_difference_B,
    select_independent_coordinates as bmatrix_select_independent_coordinates,
    svd_rank_condition as bmatrix_svd_rank_condition,
)
from internal_coordinates import (
    angle_deg_from_vectors as internal_angle_deg_from_vectors,
    angle_fn as internal_angle_fn,
    build_internal_coordinates as internal_build_internal_coordinates,
    dihedral_rad as internal_dihedral_rad,
    distance_fn as internal_distance_fn,
    torsion_fn as internal_torsion_fn,
)
from mode_tracking import (
    mode_tracking_outputs_for_hess_files as tracking_mode_tracking_outputs_for_hess_files,
)
from mode_assignment import (
    build_stage3d_assignment_audit as mode_assignment_build_stage3d_assignment_audit,
)
from orca_parser import (
    parse_atoms as parser_parse_atoms,
    parse_block_matrix as parser_parse_block_matrix,
    parse_frequencies as parser_parse_frequencies,
    parse_ir_spectrum as parser_parse_ir_spectrum,
    parse_scalar_section as parser_parse_scalar_section,
    read_orca_hess as parser_read_orca_hess,
    split_orca_hess_sections as parser_split_orca_hess_sections,
)
from reports import (
    build_spectrum_payload as reports_build_spectrum_payload,
    normalize_sheet_name as reports_normalize_sheet_name,
    output_prefix_for_hess_paths as reports_output_prefix_for_hess_paths,
    safe_output_stem as reports_safe_output_stem,
    write_interactive_spectrum_viewer as reports_write_interactive_spectrum_viewer,
    write_xlsx_report as reports_write_xlsx_report,
)
from orcaveda_cli import (
    cli_main as external_cli_main,
    colab_upload_and_run as external_colab_upload_and_run,
    is_google_colab as external_is_google_colab,
)
from orcaveda_models import ChemicalSystemAnnotation, FunctionalGroup, HessData, InternalCoordinate


BOHR_TO_ANGSTROM = 0.529177210903
EPS_FD_A = 1.0e-4


def parse_orca_out_metadata(path: str | Path) -> Dict[str, object]:
    p = Path(path)
    if not p.exists():
        return {}
    text = p.read_text(errors="ignore")
    meta: Dict[str, object] = {"out_file": p.name}
    meta["has_orca_termination_normal"] = "ORCA TERMINATED NORMALLY" in text
    meta["has_freq_job"] = bool(re.search(r"!\s*.*\bFreq\b", text, flags=re.IGNORECASE))
    m = re.search(r"^\s*!\s*(.+)$", text, flags=re.MULTILINE)
    if m:
        meta["route_line"] = m.group(1).strip()
    m = re.search(r"epsilon\s+([0-9]+(?:\.[0-9]+)?)", text, flags=re.IGNORECASE)
    if m:
        meta["cpcm_epsilon"] = float(m.group(1))
    m = re.search(r"Temp\s+([0-9]+(?:\.[0-9]+)?)", text, flags=re.IGNORECASE)
    if m:
        meta["temperature_K_out"] = float(m.group(1))
    return meta


def build_connectivity(atoms: Sequence[str], coords_A: np.ndarray, scale: float = 1.25, extra_A: float = 0.15):
    bonds = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if atoms[i] == atoms[j] == "H":
                continue
            ri = COVALENT_RADII_A.get(atoms[i], 0.75)
            rj = COVALENT_RADII_A.get(atoms[j], 0.75)
            d = float(np.linalg.norm(coords_A[i] - coords_A[j]))
            cutoff = scale * (ri + rj) + extra_A
            if d <= cutoff:
                bonds.append((i, j, d))
    return bonds


def adjacency(natoms: int, bonds) -> Dict[int, set]:
    adj = {i: set() for i in range(natoms)}
    for i, j, _ in bonds:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def bond_distance(bonds, i, j) -> Optional[float]:
    a, b = sorted((i, j))
    for x, y, d in bonds:
        if x == a and y == b:
            return d
    return None


def split_fragments(natoms: int, bonds):
    adj = adjacency(natoms, bonds)
    seen, fragments = set(), []
    for i in range(natoms):
        if i in seen:
            continue
        stack, fragment = [i], []
        seen.add(i)
        while stack:
            u = stack.pop()
            fragment.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        fragments.append(sorted(fragment))
    return fragments


def formula_string(atoms: Sequence[str]) -> str:
    counts = {a: atoms.count(a) for a in sorted(set(atoms))}
    order = ["C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    parts = []
    for el in order:
        if counts.get(el, 0):
            n = counts.pop(el)
            parts.append(f"{el}{'' if n == 1 else n}")
    for el, n in sorted(counts.items()):
        parts.append(f"{el}{'' if n == 1 else n}")
    return "".join(parts)


def classify_system(fragments: Sequence[Sequence[int]]) -> str:
    n = len(fragments)
    sizes = sorted(len(f) for f in fragments)
    if n == 1:
        return "monomer"
    if n == 2 and sizes[0] == sizes[1]:
        return "homodimer"
    if n == 2:
        return "heterodimer"
    return f"cluster_{n}_fragments"


def expected_vibrational_rank(natoms: int, linear: bool = False) -> int:
    return 3 * natoms - (5 if linear else 6)


def _angle_deg_from_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return float("nan")
    c = float(np.dot(v1, v2) / denom)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def detect_interfragment_hbonds(atoms, coords_A, bonds, fragments, h_a_max_A=2.70, d_a_max_A=3.50, angle_min_deg=120.0):
    frag_id = {atom: k for k, frag in enumerate(fragments) for atom in frag}
    bonded = {tuple(sorted((i, j))) for i, j, _ in bonds}
    donors = []
    for i, j, _ in bonds:
        if atoms[i] in ("O", "N") and atoms[j] == "H":
            donors.append((i, j))
        if atoms[j] in ("O", "N") and atoms[i] == "H":
            donors.append((j, i))
    acceptors = [i for i, atom in enumerate(atoms) if atom in ("O", "N")]
    hbonds = []
    for D, H in donors:
        for A in acceptors:
            if A == D or frag_id.get(D) == frag_id.get(A) or tuple(sorted((H, A))) in bonded:
                continue
            r_HA = float(np.linalg.norm(coords_A[H] - coords_A[A]))
            r_DA = float(np.linalg.norm(coords_A[D] - coords_A[A]))
            angle = _angle_deg_from_vectors(coords_A[D] - coords_A[H], coords_A[A] - coords_A[H])
            if r_HA <= h_a_max_A and r_DA <= d_a_max_A and angle >= angle_min_deg:
                hbonds.append(
                    {
                        "D0": D, "H0": H, "A0": A,
                        "D": D + 1, "H": H + 1, "A": A + 1,
                        "type": f"{atoms[D]}-H···{atoms[A]}",
                        "rHA_A": r_HA, "rDA_A": r_DA, "angle_deg": angle,
                    }
                )
    return sorted(hbonds, key=lambda x: x["rHA_A"])


def atom_environment_table(atoms, coords_A, bonds) -> pd.DataFrame:
    adj = adjacency(len(atoms), bonds)
    rows = []
    for i, a in enumerate(atoms):
        neigh = sorted(adj[i])
        neigh_elems = [atoms[j] for j in neigh]
        heavy_neigh = [j for j in neigh if atoms[j] != "H"]
        h_count = sum(1 for j in neigh if atoms[j] == "H")
        env = f"{a};deg{len(neigh)};H{h_count};heavy{len(heavy_neigh)};nbrs={','.join(neigh_elems)}"
        rows.append({
            "atom": i+1, "element": a,
            "degree": len(neigh),
            "H_neighbors": h_count,
            "heavy_neighbors": len(heavy_neigh),
            "neighbor_elements": ",".join(neigh_elems),
            "environment_label": env,
        })
    return pd.DataFrame(rows)


def detect_rings(natoms: int, bonds, max_size: int = 8) -> List[Tuple[int, ...]]:
    # Small simple-cycle finder adequate for organic ring flags; deduplicated by atom set.
    adj = adjacency(natoms, bonds)
    rings = set()
    def dfs(start, current, path):
        if len(path) > max_size:
            return
        for nb in adj[current]:
            if nb == start and len(path) >= 3:
                rings.add(tuple(sorted(path)))
            elif nb > start and nb not in path:
                dfs(start, nb, path + [nb])
    for start in range(natoms):
        dfs(start, start, [start])
    # remove supersets where smaller ring has same core? keep unique minimal sets
    unique = sorted(rings, key=lambda x: (len(x), x))
    minimal = []
    for r in unique:
        sr = set(r)
        if not any(set(m).issubset(sr) for m in minimal):
            minimal.append(r)
    return minimal


def detect_functional_groups(atoms, coords_A, bonds) -> List[FunctionalGroup]:
    adj = adjacency(len(atoms), bonds)
    groups: List[FunctionalGroup] = []
    rings = detect_rings(len(atoms), bonds)

    def add(group, atoms0, desc, conf, evidence):
        key = (group, tuple(sorted(atoms0)))
        if key not in {(g.group, tuple(sorted(g.atoms0))) for g in groups}:
            groups.append(FunctionalGroup(group, tuple(atoms0), desc, conf, evidence))

    # Alcohol / hydroxyl, ether
    for o, el in enumerate(atoms):
        if el != "O":
            continue
        ns = sorted(adj[o])
        h_ns = [n for n in ns if atoms[n] == "H"]
        c_ns = [n for n in ns if atoms[n] == "C"]
        s_ns = [n for n in ns if atoms[n] == "S"]
        if h_ns and c_ns:
            add("alcohol", (o, h_ns[0], c_ns[0]), "O-H bonded to carbon-bearing O", "high", "O bonded to H and C")
        if len(c_ns) >= 2:
            add("ether", (o, c_ns[0], c_ns[1]), "C-O-C oxygen", "high", "O bonded to two carbons")
        if s_ns:
            d = bond_distance(bonds, o, s_ns[0])
            if d is not None and d < 1.62:
                add("sulfoxide_S=O", (s_ns[0], o), "short S-O bond assigned as sulfoxide S=O candidate", "medium", f"S-O distance={d:.3f} A")

    # Carbonyls, amides/lactams, ketones/aldehydes/acids/esters
    carbonyl_C = []
    for c, el in enumerate(atoms):
        if el != "C":
            continue
        o_short = []
        for n in adj[c]:
            if atoms[n] == "O":
                d = bond_distance(bonds, c, n)
                if d is not None and d < 1.35:
                    o_short.append(n)
        for o in o_short:
            carbonyl_C.append(c)
            n_ns = [n for n in adj[c] if atoms[n] == "N"]
            c_ns = [n for n in adj[c] if atoms[n] == "C"]
            h_ns = [n for n in adj[c] if atoms[n] == "H"]
            o_single = [n for n in adj[c] if atoms[n] == "O" and n != o]
            add("carbonyl_C=O", (c, o), "short C-O bond assigned as carbonyl candidate", "high", f"C-O distance={bond_distance(bonds,c,o):.3f} A")
            if n_ns:
                ring_flag = any(c in r and any(n in r for n in n_ns) for r in rings)
                add("lactam_amide" if ring_flag else "amide", (c, o, n_ns[0]), "carbonyl C bonded to N", "high", "C=O carbon bonded to N")
            elif len(c_ns) >= 2:
                add("ketone", (c, o, c_ns[0], c_ns[1]), "carbonyl C bonded to two carbons", "high", "C=O carbon bonded to two C atoms")
            elif h_ns and c_ns:
                add("aldehyde", (c, o, h_ns[0], c_ns[0]), "carbonyl C bonded to H and C", "high", "C=O carbon bonded to H and C")
            elif o_single:
                o2 = o_single[0]
                if any(atoms[n] == "H" for n in adj[o2]):
                    add("carboxylic_acid", (c, o, o2), "carbonyl with hydroxyl oxygen", "high", "C=O and C-OH in same carbon")
                elif any(atoms[n] == "C" for n in adj[o2]):
                    add("ester", (c, o, o2), "carbonyl with alkoxy oxygen", "high", "C=O and C-OR in same carbon")

    # Nitrile: short C-N with carbon attached to exactly one non-H, N terminal.
    for c, el in enumerate(atoms):
        if el != "C":
            continue
        for n in adj[c]:
            if atoms[n] == "N":
                d = bond_distance(bonds, c, n)
                if d is not None and d < 1.22:
                    add("nitrile_C≡N", (c, n), "short terminal C-N bond assigned as nitrile", "high", f"C-N distance={d:.3f} A")

    # Amines/amides around nitrogen.
    for n, el in enumerate(atoms):
        if el != "N":
            continue
        ns = sorted(adj[n])
        h_ns = [x for x in ns if atoms[x] == "H"]
        c_ns = [x for x in ns if atoms[x] == "C"]
        if any(c in carbonyl_C for c in c_ns):
            # amide/lactam already added from carbonyl side
            continue
        if h_ns and c_ns:
            if len(h_ns) == 2:
                add("primary_amine", (n, *h_ns, c_ns[0]), "N bonded to carbon and two hydrogens", "high", "N-H2 and N-C")
            elif len(h_ns) == 1:
                add("secondary_amine", (n, h_ns[0], *c_ns[:2]), "N bonded to carbon(s) and one hydrogen", "high", "N-H and N-C")
        elif len(c_ns) >= 3:
            add("tertiary_amine", (n, *c_ns[:3]), "N bonded to three carbons", "high", "N-C3")
        elif len(c_ns) >= 2:
            add("dialkyl_amide_N_or_amine_N", (n, *c_ns[:2]), "N bonded to at least two carbons", "medium", "N-C2; carbonyl context checked separately")

    # Sulfoxide S center
    for s, el in enumerate(atoms):
        if el != "S":
            continue
        o_short = [n for n in adj[s] if atoms[n] == "O" and (bond_distance(bonds, s, n) or 9) < 1.62]
        c_ns = [n for n in adj[s] if atoms[n] == "C"]
        if o_short and len(c_ns) >= 2:
            add("sulfoxide", (s, o_short[0], *c_ns[:2]), "S=O with two carbon substituents", "high", "short S-O and two S-C bonds")

    # Alkyl CH groups
    for c, el in enumerate(atoms):
        if el != "C":
            continue
        h_ns = [n for n in adj[c] if atoms[n] == "H"]
        if len(h_ns) == 3:
            add("methyl", (c, *h_ns), "CH3 group", "high", "C bonded to three H")
        elif len(h_ns) == 2:
            add("methylene", (c, *h_ns), "CH2 group", "high", "C bonded to two H")
        elif len(h_ns) == 1:
            add("methine", (c, h_ns[0]), "CH group", "high", "C bonded to one H")

    # Ring flags
    for r in rings:
        add("ring", r, f"{len(r)}-membered ring candidate", "medium", "simple cycle detected in covalent graph")

    return groups


# Extracted parser and chemistry entry points are rebound here so the rest of the
# Stage 3D engine uses the new modules without changing downstream call sites.
split_orca_hess_sections = parser_split_orca_hess_sections
parse_atoms = parser_parse_atoms
parse_frequencies = parser_parse_frequencies
parse_scalar_section = parser_parse_scalar_section
parse_ir_spectrum = parser_parse_ir_spectrum
parse_block_matrix = parser_parse_block_matrix
read_orca_hess = parser_read_orca_hess
build_connectivity = chemistry_build_connectivity
split_fragments = chemistry_split_fragments
formula_string = chemistry_formula_string
classify_system = chemistry_classify_system
expected_vibrational_rank = chemistry_expected_vibrational_rank
detect_interfragment_hbonds = chemistry_detect_interfragment_hbonds
atom_environment_table = chemistry_atom_environment_table
detect_functional_groups = chemistry_detect_functional_groups
annotate_chemical_system = chemistry_annotate_chemical_system
distance_fn = internal_distance_fn
angle_fn = internal_angle_fn
dihedral_rad = internal_dihedral_rad
torsion_fn = internal_torsion_fn
build_internal_coordinates = internal_build_internal_coordinates
finite_difference_B = bmatrix_finite_difference_B
svd_rank_condition = bmatrix_svd_rank_condition
select_independent_coordinates = bmatrix_select_independent_coordinates
safe_output_stem = reports_safe_output_stem
output_prefix_for_hess_paths = reports_output_prefix_for_hess_paths
normalize_sheet_name = reports_normalize_sheet_name
write_xlsx_report = reports_write_xlsx_report
mode_tracking_outputs_for_hess_files = tracking_mode_tracking_outputs_for_hess_files


def distance_fn(i: int, j: int) -> Callable[[np.ndarray], float]:
    return lambda xyz: float(np.linalg.norm(xyz[i] - xyz[j]))


def angle_fn(i: int, j: int, k: int) -> Callable[[np.ndarray], float]:
    return lambda xyz: _angle_deg_from_vectors(xyz[i] - xyz[j], xyz[k] - xyz[j])


def dihedral_rad(p0, p1, p2, p3) -> float:
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    n = np.linalg.norm(b1)
    if n == 0:
        return float("nan")
    b1 = b1 / n
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(math.atan2(y, x))


def torsion_fn(i: int, j: int, k: int, l: int) -> Callable[[np.ndarray], float]:
    return lambda xyz: dihedral_rad(xyz[i], xyz[j], xyz[k], xyz[l])


def build_internal_coordinates(atoms, coords_A, bonds, fragments, hbonds, groups: Optional[List[FunctionalGroup]] = None) -> List[InternalCoordinate]:
    adj = adjacency(len(atoms), bonds)
    coords: List[InternalCoordinate] = []

    # Bond stretches.
    for i, j, _ in bonds:
        label = f"r({atoms[i]}{i+1}-{atoms[j]}{j+1})"
        pri = 10 if "H" not in (atoms[i], atoms[j]) else 25
        coords.append(InternalCoordinate(label, "stretch", (i, j), pri, distance_fn(i, j), "primitive"))

    # Angle bends.
    seen_angles = set()
    for j in range(len(atoms)):
        neigh = sorted(adj[j])
        for a in range(len(neigh)):
            for b in range(a + 1, len(neigh)):
                i, k = neigh[a], neigh[b]
                key = (i, j, k)
                if key in seen_angles:
                    continue
                seen_angles.add(key)
                pri = 35 if atoms[j] != "H" else 65
                label = f"ang({atoms[i]}{i+1}-{atoms[j]}{j+1}-{atoms[k]}{k+1})"
                coords.append(InternalCoordinate(label, "bend", key, pri, angle_fn(i, j, k), "primitive"))

    # Proper torsions around bonds.
    seen_torsions = set()
    for j, k, _ in bonds:
        for i in adj[j] - {k}:
            for l in adj[k] - {j}:
                if i == l:
                    continue
                key = (i, j, k, l)
                rev = (l, k, j, i)
                if key in seen_torsions or rev in seen_torsions:
                    continue
                seen_torsions.add(key)
                label = f"tor({atoms[i]}{i+1}-{atoms[j]}{j+1}-{atoms[k]}{k+1}-{atoms[l]}{l+1})"
                coords.append(InternalCoordinate(label, "torsion", key, 55, torsion_fn(i, j, k, l), "primitive"))

    # Functional-group local coordinate templates: implemented as named coordinate aliases
    # and prioritized above generic primitives where chemically diagnostic.
    if groups:
        for g in groups:
            if g.group == "alcohol":
                o, h, c = g.atoms0[:3]
                coords.append(InternalCoordinate(f"FG_alcohol_OH_stretch({atoms[o]}{o+1}-{atoms[h]}{h+1})", "fg_OH_stretch", (o, h), 5, distance_fn(o, h), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_alcohol_CO_stretch({atoms[c]}{c+1}-{atoms[o]}{o+1})", "fg_CO_stretch", (c, o), 6, distance_fn(c, o), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_alcohol_COH_bend({atoms[c]}{c+1}-{atoms[o]}{o+1}-{atoms[h]}{h+1})", "fg_COH_bend", (c, o, h), 20, angle_fn(c, o, h), "functional_group_template"))
            elif g.group in ("carbonyl_C=O", "ketone", "amide", "lactam_amide"):
                c, o = g.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_carbonyl_CO_stretch({atoms[c]}{c+1}={atoms[o]}{o+1})", "fg_carbonyl_stretch", (c, o), 4, distance_fn(c, o), "functional_group_template"))
            elif g.group == "nitrile_C≡N":
                c, n = g.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_nitrile_CN_stretch({atoms[c]}{c+1}≡{atoms[n]}{n+1})", "fg_nitrile_stretch", (c, n), 4, distance_fn(c, n), "functional_group_template"))
            elif g.group in ("sulfoxide", "sulfoxide_S=O"):
                s, o = g.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_sulfoxide_SO_stretch({atoms[s]}{s+1}={atoms[o]}{o+1})", "fg_sulfoxide_stretch", (s, o), 4, distance_fn(s, o), "functional_group_template"))
            elif g.group in ("methyl", "methylene", "methine"):
                c = g.atoms0[0]
                hs = [x for x in g.atoms0[1:] if atoms[x] == "H"]
                for h in hs:
                    coords.append(InternalCoordinate(f"FG_{g.group}_CH_stretch({atoms[c]}{c+1}-{atoms[h]}{h+1})", "fg_CH_stretch", (c, h), 15, distance_fn(c, h), "functional_group_template"))

    # Interfragment hydrogen-bond coordinates.
    for h in hbonds:
        D, H, A = h["D0"], h["H0"], h["A0"]
        coords.append(InternalCoordinate(f"Hbond_rHA({atoms[H]}{H+1}···{atoms[A]}{A+1})", "hbond_HA", (H, A), 5, distance_fn(H, A), "cluster_template"))
        coords.append(InternalCoordinate(f"Hbond_rDA({atoms[D]}{D+1}···{atoms[A]}{A+1})", "hbond_DA", (D, A), 6, distance_fn(D, A), "cluster_template"))
        coords.append(InternalCoordinate(f"Hbond_ang({atoms[D]}{D+1}-{atoms[H]}{H+1}···{atoms[A]}{A+1})", "hbond_angle", (D, H, A), 7, angle_fn(D, H, A), "cluster_template"))

    # Interfragment heavy-atom distances for cluster relative motion.
    if len(fragments) >= 2:
        frag_id = {atom_index: frag_index for frag_index, frag in enumerate(fragments) for atom_index in frag}
        heavy = [i for i, a in enumerate(atoms) if a != "H"]
        for pos, i in enumerate(heavy):
            for j in heavy[pos + 1:]:
                if frag_id.get(i) != frag_id.get(j):
                    label = f"interfrag_R({atoms[i]}{i+1}···{atoms[j]}{j+1})"
                    coords.append(InternalCoordinate(label, "interfragment_distance", (i, j), 40, distance_fn(i, j), "cluster_template"))

    return sorted(coords, key=lambda c: (c.priority, c.name))


def finite_difference_B(coords_A: np.ndarray, internals: Sequence[InternalCoordinate], eps: float = EPS_FD_A) -> np.ndarray:
    n = coords_A.size
    B = np.zeros((len(internals), n), dtype=float)
    flat0 = coords_A.reshape(-1)
    for r, ic in enumerate(internals):
        for c in range(n):
            plus = flat0.copy(); minus = flat0.copy()
            plus[c] += eps; minus[c] -= eps
            xp = plus.reshape(coords_A.shape); xm = minus.reshape(coords_A.shape)
            vp = ic.fn(xp); vm = ic.fn(xm)
            if ic.kind == "torsion":
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


def select_independent_coordinates(B: np.ndarray, internals: Sequence[InternalCoordinate], target_rank: int, tol_abs: float = 1.0e-6):
    selected_idx: List[int] = []
    current = np.zeros((0, B.shape[1]))
    current_rank = 0
    ordered = sorted(range(len(internals)), key=lambda i: (internals[i].priority, internals[i].name))
    for idx in ordered:
        candidate = np.vstack([current, B[idx:idx+1, :]])
        rank, _, _ = svd_rank_condition(candidate, tol_abs=tol_abs)
        if rank > current_rank:
            selected_idx.append(idx)
            current = candidate
            current_rank = rank
        if current_rank >= target_rank:
            break
    rank, cond, s = svd_rank_condition(current, tol_abs=tol_abs)
    return selected_idx, rank, cond, s


def classify_mode_region(freq: float) -> str:
    if freq < 0:
        return "imaginary"
    if freq < 100:
        return "very_low"
    if freq < 300:
        return "low_frequency"
    if freq < 1800:
        return "fingerprint"
    if freq < 2300:
        return "triple_bond_or_combination_region"
    if freq < 2800:
        return "mid_silent_gap"
    if freq < 3000:
        return "CH_stretch_edge"
    if freq < 3600:
        return "CH_NH_stretch_region"
    return "OH_NH_high_stretch_region"


def fragment_motion_metrics(hess: HessData, fragments, mode: int):
    disp = hess.normal_modes[:, mode].reshape(len(hess.atoms), 3)
    total = float(np.sum(disp**2))
    translation, internal, means = 0.0, 0.0, []
    for fragment in fragments:
        D = disp[fragment]
        mean = D.mean(axis=0)
        means.append(mean)
        translation += len(fragment) * float(np.sum(mean**2))
        internal += float(np.sum((D - mean) ** 2))
    return {
        "fragment_translation_percent": 100.0 * translation / total if total else 0.0,
        "fragment_internal_deformation_percent": 100.0 * internal / total if total else 0.0,
        "fragment_mean_relative_amplitude": float(np.linalg.norm(means[0] - means[1])) if len(means) >= 2 else 0.0,
    }



# =============================================================================
# Stage 3D v2: user-facing assignment audit
# =============================================================================

def _compact_coord_label(name: str, max_len: Optional[int] = None) -> str:
    """
    Return a readable coordinate label for report tables.

    Stage 3D v4.3 rule:
        no truncation in audit fields. Long labels are audit evidence and must
        remain reproducible in CSV/XLSX outputs.
    """
    return str(name).replace("FG_", "").replace("Hbond_", "Hbond:")


def _assignment_family_from_internal(ic: InternalCoordinate) -> str:
    """Map an internal coordinate to a chemically readable assignment family.

    Stage 3D v3 fix:
    generic primitive labels such as r(N7-H10) and ang(H10-N7-H11)
    must be recognized chemically. The v2 logic recognized functional-group
    CH/OH labels but could mislabel primitive N-H stretches when CH terms were
    also present in the top list.
    """
    kind = ic.kind.lower()
    name = ic.name.lower()

    if "hbond" in kind or "hbond" in name:
        return "H-bond / intermolecular"
    if "interfragment" in kind or "interfrag" in name:
        return "intermolecular / cluster motion"

    # Scissor / local angle labels from primitive angle names.
    if re.search(r"ang\(h\d+-n\d+-h\d+\)", name):
        return "NH2 scissor"
    if re.search(r"ang\(h\d+-c\d+-h\d+\)", name):
        return "CH2 scissor"

    # Diagnostic stretches first. Handle both functional-template names
    # and primitive r(Xn-Hm) names.
    if (
        "oh_stretch" in kind
        or "alcohol_oh_stretch" in name
        or re.search(r"r\(o\d+-h\d+\)", name)
        or re.search(r"r\(h\d+-o\d+\)", name)
        or ("stretch" in kind and re.search(r"o\d+-h\d+|h\d+-o\d+", name))
    ):
        return "O-H stretch"
    if (
        "nh_stretch" in kind
        or re.search(r"r\(n\d+-h\d+\)", name)
        or re.search(r"r\(h\d+-n\d+\)", name)
        or ("stretch" in kind and re.search(r"n\d+-h\d+|h\d+-n\d+", name))
    ):
        return "N-H stretch"
    if (
        "ch_stretch" in kind
        or re.search(r"r\(c\d+-h\d+\)", name)
        or re.search(r"r\(h\d+-c\d+\)", name)
        or ("stretch" in kind and re.search(r"c\d+-h\d+|h\d+-c\d+", name))
    ):
        return "C-H stretch"

    if "carbonyl" in kind or "carbonyl" in name:
        return "C=O stretch"
    if "nitrile" in kind or "nitrile" in name:
        return "C≡N stretch"
    if "sulfoxide" in kind or "sulfoxide" in name:
        return "S=O stretch"
    if "co_stretch" in kind:
        return "C-O stretch"
    if "stretch" in kind:
        return "bond stretch"

    if "bend" in kind or "angle" in kind:
        if "hbond" in name:
            return "H-bond angle bend"
        return "angle bend"
    if "torsion" in kind or "tor(" in name:
        return "torsion"
    return kind or "internal coordinate"


def _stage3d_coord_class(ic: InternalCoordinate) -> str:
    """Coarse coordinate class used for totals and confidence diagnostics."""
    kind = ic.kind.lower()
    name = ic.name.lower()
    if "hbond" in kind or "hbond" in name:
        return "hbond"
    if "interfragment" in kind or "interfrag" in name:
        return "interfragment"
    if "stretch" in kind or kind == "stretch":
        return "stretch"
    if "bend" in kind or "angle" in kind:
        return "bend"
    if "torsion" in kind or "tor(" in name:
        return "torsion"
    return "other"


def _stage3d_xh_stretch_info(ic: InternalCoordinate, atoms: Sequence[str]) -> Optional[Tuple[int, str, int]]:
    """Return (heavy_atom_index, heavy_element, hydrogen_index) for X-H stretch coordinates."""
    kind = ic.kind.lower()
    name = ic.name.lower()
    if "stretch" not in kind and "stretch" not in name and not name.startswith("r("):
        return None
    if len(ic.atoms0) != 2:
        return None
    a, b = ic.atoms0
    ea, eb = atoms[a], atoms[b]
    if ea == "H" and eb in ("C", "N", "O"):
        return b, eb, a
    if eb == "H" and ea in ("C", "N", "O"):
        return a, ea, b
    return None


def _stage3d_xh_center_h_count(
    internals: Sequence[InternalCoordinate],
    atoms: Sequence[str],
    heavy: int,
    elem: Optional[str] = None,
) -> int:
    """
    Count unique H atoms attached to an X-H stretch center as represented in
    the redundant/local coordinate pool.

    Stage 3D v4.8 guard:
        CH2 symmetric/asymmetric labels are allowed only for carbon centers
        with exactly two H stretch partners. Methyl/methine centers must remain
        CH3/generic C-H, not CH2.
    """
    hs = set()
    for ic in internals:
        info = _stage3d_xh_stretch_info(ic, atoms)
        if info is None:
            continue
        hvy, el, h = info
        if int(hvy) == int(heavy) and (elem is None or str(el) == str(elem)):
            hs.add(int(h))
    return len(hs)


def _stage3d_local_xh2_symmetry_result(
    internals: Sequence[InternalCoordinate],
    atoms: Sequence[str],
    projections: np.ndarray,
    pct: np.ndarray,
    freq: float,
) -> Tuple[Optional[str], Dict[str, object]]:
    """
    Assign local XH2 symmetric/asymmetric stretch and report balance diagnostics.

    The label is based on the signs of the two strongest X-H stretch projections
    on the same heavy atom. Confidence is penalized later when the two reported
    percentage contributions are strongly imbalanced.
    """
    empty = {
        "xh2_center_atom": "",
        "xh2_center_element": "",
        "xh2_pair_coords": "",
        "xh2_pair_percent_1": 0.0,
        "xh2_pair_percent_2": 0.0,
        "xh2_balance_ratio": 0.0,
        "xh2_imbalance_warning": "",
    }
    if freq < 2800.0:
        return None, empty

    centers: Dict[Tuple[int, str], List[Tuple[int, int, float, float, str]]] = {}
    for idx, ic in enumerate(internals):
        info = _stage3d_xh_stretch_info(ic, atoms)
        if info is None:
            continue
        heavy, elem, h = info
        centers.setdefault((heavy, elem), []).append(
            (idx, h, float(pct[idx]), float(projections[idx]), _compact_coord_label(ic.name))
        )

    best = None
    best_total = 0.0
    for (heavy, elem), vals in centers.items():
        vals_sorted = sorted(vals, key=lambda x: x[2], reverse=True)
        if len(vals_sorted) < 2:
            continue
        total = vals_sorted[0][2] + vals_sorted[1][2]
        if total > best_total:
            best_total = total
            best = (heavy, elem, vals_sorted[:2])

    if best is None or best_total < 35.0:
        return None, empty

    heavy, elem, vals2 = best
    p1 = vals2[0][3]
    p2 = vals2[1][3]
    pct1 = float(vals2[0][2])
    pct2 = float(vals2[1][2])
    hi = max(pct1, pct2)
    lo = min(pct1, pct2)
    balance_ratio = lo / hi if hi > 0.0 else 0.0
    same_sign = (p1 == 0.0 or p2 == 0.0) or (p1 * p2 > 0.0)

    label = None
    if elem == "N":
        label = "NH2 symmetric stretch" if same_sign else "NH2 asymmetric stretch"
    elif elem == "C":
        c_h_count = _stage3d_xh_center_h_count(internals, atoms, heavy, elem)
        label = (
            "CH2 symmetric stretch" if same_sign else "CH2 asymmetric stretch"
        ) if c_h_count == 2 else "C-H stretch"
    elif elem == "O":
        label = "O-H stretch"

    diag = {
        "xh2_center_atom": heavy + 1,
        "xh2_center_element": elem,
        "xh2_pair_coords": " | ".join([vals2[0][4], vals2[1][4]]),
        "xh2_pair_percent_1": round(pct1, 3),
        "xh2_pair_percent_2": round(pct2, 3),
        "xh2_balance_ratio": round(float(balance_ratio), 3),
        "xh2_imbalance_warning": (
            "imbalanced_XH2_pair_contributions"
            if label and elem in ("C", "N") and balance_ratio < 0.35
            else ""
        ),
    }
    return label, diag



def _stage3d_protected_xh_stretch_audit(
    internals: Sequence[InternalCoordinate],
    atoms: Sequence[str],
    B: np.ndarray,
    mode_unit: np.ndarray,
    freq: float,
) -> Dict[str, object]:
    """
    Protected X-H stretch audit over the full redundant/local coordinate pool.

    Stage 3D v4.3 rationale:
        independent-coordinate selection is rank-driven and may exclude a
        chemically diagnostic X-H coordinate from the selected basis. High
        frequency modes must therefore be checked against all primitive and
        functional-group X-H stretch coordinates before being left unassigned.

    This remains an assignment-audit layer only. It does not modify frequencies,
    normal modes, the B matrix, or the independent-basis diagnostics.
    """
    empty = {
        "protected_xh_assignment": "",
        "protected_xh_top_coordinates": "",
        "protected_xh_top1_coord": "",
        "protected_xh_top1_percent": 0.0,
        "protected_xh_total_percent": 0.0,
        "protected_xh_center_atom": "",
        "protected_xh_center_element": "",
        "protected_xh_pair_coords": "",
        "protected_xh_pair_percent_1": 0.0,
        "protected_xh_pair_percent_2": 0.0,
        "protected_xh_balance_ratio": 0.0,
        "protected_xh_used": False,
    }
    if freq < 2800.0 or B.size == 0:
        return empty.copy()

    candidates = []
    for idx, ic in enumerate(internals):
        info = _stage3d_xh_stretch_info(ic, atoms)
        if info is None:
            continue
        row = np.asarray(B[idx], dtype=float)
        row_norm = float(np.linalg.norm(row))
        if row_norm <= 0.0 or not np.isfinite(row_norm):
            continue
        proj = float(np.dot(row / row_norm, mode_unit))
        weight = (proj ** 2) * _stage3d_frequency_region_priority(float(freq), ic)
        if not np.isfinite(weight) or weight <= 0.0:
            continue
        heavy, elem, h = info
        candidates.append({
            "idx": idx,
            "ic": ic,
            "heavy": heavy,
            "elem": elem,
            "h": h,
            "projection": proj,
            "weight": float(weight),
            "coord": _compact_coord_label(ic.name),
            "family": _assignment_family_from_internal(ic),
        })

    if not candidates:
        return empty.copy()

    # Deduplicate exact atom-pair aliases by keeping the strongest coordinate
    # for each heavy-H pair. This prevents primitive and functional-group aliases
    # from double-counting the same physical stretch.
    pair_best: Dict[Tuple[int, int], Dict[str, object]] = {}
    for c in candidates:
        key = tuple(sorted((int(c["heavy"]), int(c["h"]))))
        old = pair_best.get(key)
        if old is None or float(c["weight"]) > float(old["weight"]):
            pair_best[key] = c
    candidates = list(pair_best.values())

    total_weight = float(sum(float(c["weight"]) for c in candidates))
    if total_weight <= 0.0 or not np.isfinite(total_weight):
        return empty.copy()
    for c in candidates:
        c["percent"] = 100.0 * float(c["weight"]) / total_weight

    top_candidates = sorted(candidates, key=lambda x: float(x["percent"]), reverse=True)
    top_terms = "; ".join(f"{c['coord']}={float(c['percent']):.1f}%" for c in top_candidates[:8])
    top1 = top_candidates[0]

    # Prefer the strongest local XH2 center when two H atoms are present;
    # otherwise use the strongest individual X-H stretch.
    centers: Dict[Tuple[int, str], List[Dict[str, object]]] = {}
    for c in candidates:
        centers.setdefault((int(c["heavy"]), str(c["elem"])), []).append(c)

    best_center = None
    best_center_total = 0.0
    for (heavy, elem), vals in centers.items():
        vals_sorted = sorted(vals, key=lambda x: float(x["percent"]), reverse=True)
        center_total = sum(float(v["percent"]) for v in vals_sorted[:2])
        if center_total > best_center_total:
            best_center_total = center_total
            best_center = (heavy, elem, vals_sorted[:2])

    assignment = ""
    pair_coords = ""
    p1 = p2 = 0.0
    balance_ratio = 0.0
    center_atom = ""
    center_element = ""

    if best_center is not None:
        heavy, elem, vals = best_center
        center_atom = heavy + 1
        center_element = elem
        if elem == "O":
            assignment = "O-H stretch"
        elif elem == "N" and len(vals) >= 2:
            same_sign = (
                float(vals[0]["projection"]) == 0.0
                or float(vals[1]["projection"]) == 0.0
                or float(vals[0]["projection"]) * float(vals[1]["projection"]) > 0.0
            )
            assignment = "NH2 symmetric stretch" if same_sign else "NH2 asymmetric stretch"
        elif elem == "N":
            assignment = "N-H stretch"
        elif elem == "C" and len(vals) >= 2:
            same_sign = (
                float(vals[0]["projection"]) == 0.0
                or float(vals[1]["projection"]) == 0.0
                or float(vals[0]["projection"]) * float(vals[1]["projection"]) > 0.0
            )
            c_h_count = _stage3d_xh_center_h_count(internals, atoms, heavy, elem)
            assignment = (
                "CH2 symmetric stretch" if same_sign else "CH2 asymmetric stretch"
            ) if c_h_count == 2 else "C-H stretch"
        elif elem == "C":
            assignment = "C-H stretch"

        if len(vals) >= 2:
            p1 = float(vals[0]["percent"])
            p2 = float(vals[1]["percent"])
            hi = max(p1, p2)
            lo = min(p1, p2)
            balance_ratio = lo / hi if hi > 0.0 else 0.0
            pair_coords = " | ".join([str(vals[0]["coord"]), str(vals[1]["coord"])])

    out = empty.copy()
    out.update({
        "protected_xh_assignment": assignment,
        "protected_xh_top_coordinates": top_terms,
        "protected_xh_top1_coord": str(top1["coord"]),
        "protected_xh_top1_percent": round(float(top1["percent"]), 3),
        "protected_xh_total_percent": 100.0,
        "protected_xh_center_atom": center_atom,
        "protected_xh_center_element": center_element,
        "protected_xh_pair_coords": pair_coords,
        "protected_xh_pair_percent_1": round(p1, 3),
        "protected_xh_pair_percent_2": round(p2, 3),
        "protected_xh_balance_ratio": round(float(balance_ratio), 3),
        "protected_xh_used": False,
    })
    return out


def _stage3d_local_xh_symmetry_assignment(
    internals: Sequence[InternalCoordinate],
    atoms: Sequence[str],
    projections: np.ndarray,
    pct: np.ndarray,
    freq: float,
) -> Optional[str]:
    """Backward-compatible wrapper retained for external callers."""
    label, _diag = _stage3d_local_xh2_symmetry_result(internals, atoms, projections, pct, freq)
    return label



def _stage3d_frequency_region_priority(freq: float, ic: InternalCoordinate) -> float:
    """
    Empirical audit weighting for chemically interpretable assignments.

    This is not strict VEDA/PED. It prevents angular coordinates from dominating
    X-H stretching regions when the topological coordinate basis contains many
    angle rows with large directional overlap.
    """
    cls = _stage3d_coord_class(ic)
    fam = _assignment_family_from_internal(ic)

    weight = 1.0
    if cls == "stretch":
        weight *= 1.35
    elif cls == "bend":
        weight *= 0.45
    elif cls == "torsion":
        weight *= 0.30
    elif cls in ("hbond", "interfragment"):
        weight *= 0.60

    # High-frequency region: X-H stretches must be prioritized for assignment audit.
    if freq >= 2800.0:
        if cls == "stretch":
            weight *= 4.0
        elif cls == "bend":
            weight *= 0.08
        elif cls == "torsion":
            weight *= 0.04

        if fam in ("O-H stretch", "N-H stretch", "C-H stretch"):
            weight *= 2.0

    # Mid-frequency heteroatom diagnostic stretches.
    if 1500.0 <= freq <= 2300.0 and fam in ("C=O stretch", "C≡N stretch", "S=O stretch"):
        weight *= 3.0

    # Low-frequency cluster modes should not be over-interpreted as local bends.
    if freq < 300.0 and cls in ("hbond", "interfragment"):
        weight *= 1.8

    return float(weight)


def _stage3d_assignment_from_weighted_terms(
    top_terms: Sequence[Tuple[InternalCoordinate, float]],
    totals: Dict[str, float],
    freq: float,
) -> str:
    """Create the final readable assignment from weighted dominant terms."""
    if not top_terms:
        return "unassigned"

    family_totals: Dict[str, float] = {}
    for ic, pct in top_terms:
        fam = _assignment_family_from_internal(ic)
        family_totals[fam] = family_totals.get(fam, 0.0) + float(pct)

    ordered = sorted(family_totals.items(), key=lambda x: x[1], reverse=True)
    primary, primary_pct = ordered[0]

    # Stage 3D v4.8: mid-frequency diagnostic stretch priority.
    #
    # If a chemically diagnostic functional-group stretch is the strongest
    # individual coordinate, do not let the aggregate angle-bend family become
    # the leading human-readable label solely because several bend coordinates
    # sum together. This affects the assignment label only; percentages,
    # projections, B-matrix diagnostics, and confidence scoring are unchanged.
    top1_ic, top1_pct = top_terms[0]
    top1_family = _assignment_family_from_internal(top1_ic)
    diagnostic_stretches = {"C-O stretch", "C=O stretch", "C≡N stretch", "S=O stretch"}
    if (
        800.0 <= float(freq) <= 2400.0
        and top1_family in diagnostic_stretches
        and float(top1_pct) >= 30.0
        and primary == "angle bend"
    ):
        secondary = primary
        return f"{top1_family} mixed with {secondary}"

    # Conservative high-frequency override: do not call a 3000+ cm-1 mode a bend
    # when any chemically diagnostic X-H stretch survives the weighted audit.
    if freq >= 2800.0:
        xh = [(fam, val) for fam, val in ordered if fam in ("O-H stretch", "N-H stretch", "C-H stretch")]
        if xh and totals.get("stretch", 0.0) >= 8.0:
            xh = sorted(xh, key=lambda x: x[1], reverse=True)
            primary, primary_pct = xh[0]

    if len(ordered) == 1 or ordered[1][1] < 10.0:
        return primary

    secondary, secondary_pct = ordered[1]
    if secondary == primary:
        return primary
    return f"{primary} mixed with {secondary}"


def _stage3d_assignment_confidence(
    top_percent: float,
    second_percent: float,
    totals: Dict[str, float],
    warnings: str,
    severity: float,
    freq: float,
    assignment: str,
    xh2_balance_ratio: float = 1.0,
) -> Tuple[str, float]:
    """Score assignment confidence for the audit table."""
    gap = max(0.0, float(top_percent) - float(second_percent))
    score = 35.0 + 0.45 * float(top_percent) + 0.25 * gap

    dominant_total = max(totals.values()) if totals else 0.0
    score += 0.15 * dominant_total

    if "mixed with" in str(assignment):
        score -= 8.0
    if "diffuse_internal_coordinate_contributions" in str(warnings):
        score -= 14.0
    if "near_degenerate" in str(warnings):
        score -= 6.0
    if "imbalanced_XH2_pair_contributions" in str(warnings):
        # Stage 3D v4: symmetric/asymmetric XH2 labels are allowed, but a very
        # one-sided pair contribution cannot be assigned high confidence.
        score -= 28.0 if float(xh2_balance_ratio) < 0.20 else 18.0
    score -= 3.0 * float(severity)

    # Guard against false high confidence in the X-H region.
    if freq >= 2800.0 and "stretch" not in str(assignment):
        score -= 25.0

    score = max(0.0, min(100.0, round(score, 1)))
    if score >= 75.0:
        return "high", score
    if score >= 50.0:
        return "medium", score
    return "low", score



def build_stage3d_assignment_audit(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    mode_df: pd.DataFrame,
    source_label: str,
    top_n: int = 8,
) -> pd.DataFrame:
    """
    Build a user-facing per-mode assignment audit table.

    Stage 3D v2 method:
        1. Use only selected independent internal coordinates.
        2. Normalize each B row before projection to remove raw unit dominance.
        3. Apply transparent audit weights by coordinate class and frequency region.
        4. Report class totals and top coordinates separately.

    This is an assignment audit layer. It is not strict Wilson GF and not
    VEDA-equivalent PED.
    """
    rows: List[Dict[str, object]] = []

    selected_idx = list(selected_idx)
    if not selected_idx:
        for mode, freq in enumerate(hess.frequencies_cm1):
            rows.append({
                "Source": source_label,
                "Filename": hess.filename,
                "mode": mode,
                "frequency_cm-1": float(freq),
                "IR_intensity": float(hess.ir_intensities[mode]),
                "region": classify_mode_region(float(freq)),
                "top_internal_coordinates": "",
                "top1_coord": "",
                "top1_kind": "",
                "top1_percent": 0.0,
                "top2_coord": "",
                "top2_kind": "",
                "top2_percent": 0.0,
                "stretch_total_percent": 0.0,
                "bend_total_percent": 0.0,
                "torsion_total_percent": 0.0,
                "hbond_total_percent": 0.0,
                "interfragment_total_percent": 0.0,
                "xh2_center_atom": "",
                "xh2_center_element": "",
                "xh2_pair_coords": "",
                "xh2_pair_percent_1": 0.0,
                "xh2_pair_percent_2": 0.0,
                "xh2_balance_ratio": 0.0,
                "protected_xh_assignment": "",
                "protected_xh_used": False,
                "protected_xh_top_coordinates": "",
                "protected_xh_top1_coord": "",
                "protected_xh_top1_percent": 0.0,
                "protected_xh_center_atom": "",
                "protected_xh_center_element": "",
                "protected_xh_pair_coords": "",
                "protected_xh_pair_percent_1": 0.0,
                "protected_xh_pair_percent_2": 0.0,
                "protected_xh_balance_ratio": 0.0,
                "functional_group_assignment": "unassigned",
                "assignment_confidence": "low",
                "assignment_confidence_score": 0.0,
                "warnings": "no_independent_internal_coordinates_selected",
                "assignment_method": "Stage 3D v4.3 weighted independent-coordinate audit with protected X-H stretch audit",
            })
        return pd.DataFrame(rows)

    selected_B = np.asarray(B[selected_idx, :], dtype=float)
    selected_internals = [internals[i] for i in selected_idx]

    row_norms = np.linalg.norm(selected_B, axis=1)
    valid = row_norms > 1.0e-12
    selected_B_unit = np.zeros_like(selected_B)
    selected_B_unit[valid, :] = selected_B[valid, :] / row_norms[valid, None]

    warning_lookup = {}
    if mode_df is not None and not mode_df.empty and "mode" in mode_df.columns:
        warning_lookup = {int(r["mode"]): r for _, r in mode_df.iterrows()}

    for mode, freq in enumerate(hess.frequencies_cm1):
        mode_vec = np.asarray(hess.normal_modes[:, mode], dtype=float)
        mode_norm = float(np.linalg.norm(mode_vec))
        if mode_norm > 0.0:
            mode_unit = mode_vec / mode_norm
        else:
            mode_unit = mode_vec

        projections = selected_B_unit @ mode_unit
        weights = projections ** 2

        # Transparent chemical audit weights. This intentionally does not alter
        # frequencies or normal modes; it only affects the assignment label.
        audit_weights = np.array([
            _stage3d_frequency_region_priority(float(freq), ic)
            for ic in selected_internals
        ], dtype=float)
        weights = weights * audit_weights
        weights[~np.isfinite(weights)] = 0.0

        total = float(np.sum(weights))
        if total > 0.0 and np.isfinite(total):
            pct = 100.0 * weights / total
        else:
            pct = np.zeros_like(weights)

        class_totals = {
            "stretch": 0.0,
            "bend": 0.0,
            "torsion": 0.0,
            "hbond": 0.0,
            "interfragment": 0.0,
            "other": 0.0,
        }
        for ic, value in zip(selected_internals, pct):
            cls = _stage3d_coord_class(ic)
            class_totals[cls] = class_totals.get(cls, 0.0) + float(value)

        order = np.argsort(pct)[::-1]
        top = [(selected_internals[i], float(pct[i])) for i in order[:top_n] if float(pct[i]) > 0.0]
        top_terms_str = "; ".join(f"{_compact_coord_label(ic.name)}={value:.1f}%" for ic, value in top)

        top1_ic, top1_pct = top[0] if top else (None, 0.0)
        top2_ic, top2_pct = top[1] if len(top) > 1 else (None, 0.0)

        assignment = _stage3d_assignment_from_weighted_terms(top, class_totals, float(freq))
        local_xh_assignment, xh2_diag = _stage3d_local_xh2_symmetry_result(
            selected_internals,
            hess.atoms,
            projections,
            pct,
            float(freq),
        )
        if local_xh_assignment:
            assignment = local_xh_assignment

        protected_xh_diag = _stage3d_protected_xh_stretch_audit(
            internals,
            hess.atoms,
            B,
            mode_unit,
            float(freq),
        )
        protected_assignment = str(protected_xh_diag.get("protected_xh_assignment", "") or "")
        protected_used = False

        # Protected X-H override:
        # If independent-basis audit failed to expose a clear stretch in the
        # high-frequency region, use the full redundant/local X-H stretch audit
        # rather than leaving the mode unassigned.
        if (
            float(freq) >= 2800.0
            and protected_assignment
            and (
                assignment == "unassigned"
                or class_totals.get("stretch", 0.0) < 8.0
                or "stretch" not in str(assignment)
            )
        ):
            assignment = protected_assignment
            protected_used = True
            protected_xh_diag["protected_xh_used"] = True
            class_totals["stretch"] = max(
                float(class_totals.get("stretch", 0.0)),
                float(protected_xh_diag.get("protected_xh_total_percent", 0.0) or 0.0),
            )
            # Preserve independent-basis evidence when present, but fill empty
            # top fields with protected X-H evidence for audit readability.
            if not top_terms_str:
                top_terms_str = str(protected_xh_diag.get("protected_xh_top_coordinates", "") or "")
            if top1_ic is None:
                top1_pct = float(protected_xh_diag.get("protected_xh_top1_percent", 0.0) or 0.0)

        warn_row = warning_lookup.get(mode, {})
        warnings = str(warn_row.get("warnings", "") or "")
        warning_count = int(warn_row.get("warning_count", 0) or 0)
        severity = float(warn_row.get("severity_score", 0.0) or 0.0)

        extra_warnings = []
        if top1_pct < 25.0:
            extra_warnings.append("diffuse_internal_coordinate_contributions")
        if freq >= 2800.0 and class_totals.get("stretch", 0.0) < 8.0 and not protected_used:
            extra_warnings.append("high_frequency_mode_without_clear_stretch_coordinate")
        if protected_used:
            extra_warnings.append("protected_XH_stretch_assignment_used")
        if xh2_diag.get("xh2_imbalance_warning"):
            extra_warnings.append(str(xh2_diag["xh2_imbalance_warning"]))
        if warning_count:
            extra_warnings.append("mode_level_warning_present")

        all_warnings = "; ".join([w for w in [warnings, "; ".join(extra_warnings)] if w])

        conf_label, conf_score = _stage3d_assignment_confidence(
            top1_pct,
            top2_pct,
            class_totals,
            all_warnings,
            severity,
            float(freq),
            assignment,
            float(xh2_diag.get("xh2_balance_ratio", 1.0) or 0.0),
        )

        rows.append({
            "Source": source_label,
            "Filename": hess.filename,
            "mode": mode,
            "frequency_cm-1": float(freq),
            "IR_intensity": float(hess.ir_intensities[mode]),
            "region": classify_mode_region(float(freq)),
            "top_internal_coordinates": top_terms_str,
            "top1_coord": (
                _compact_coord_label(top1_ic.name)
                if top1_ic else str(protected_xh_diag.get("protected_xh_top1_coord", "") or "")
            ),
            "top1_kind": top1_ic.kind if top1_ic else ("protected_xh_stretch" if protected_used else ""),
            "top1_percent": round(float(top1_pct), 3),
            "top2_coord": _compact_coord_label(top2_ic.name) if top2_ic else "",
            "top2_kind": top2_ic.kind if top2_ic else "",
            "top2_percent": round(float(top2_pct), 3),
            "stretch_total_percent": round(float(class_totals.get("stretch", 0.0)), 3),
            "bend_total_percent": round(float(class_totals.get("bend", 0.0)), 3),
            "torsion_total_percent": round(float(class_totals.get("torsion", 0.0)), 3),
            "hbond_total_percent": round(float(class_totals.get("hbond", 0.0)), 3),
            "interfragment_total_percent": round(float(class_totals.get("interfragment", 0.0)), 3),
            "xh2_center_atom": xh2_diag.get("xh2_center_atom", ""),
            "xh2_center_element": xh2_diag.get("xh2_center_element", ""),
            "xh2_pair_coords": xh2_diag.get("xh2_pair_coords", ""),
            "xh2_pair_percent_1": xh2_diag.get("xh2_pair_percent_1", 0.0),
            "xh2_pair_percent_2": xh2_diag.get("xh2_pair_percent_2", 0.0),
            "xh2_balance_ratio": xh2_diag.get("xh2_balance_ratio", 0.0),
            "protected_xh_assignment": protected_xh_diag.get("protected_xh_assignment", ""),
            "protected_xh_used": bool(protected_xh_diag.get("protected_xh_used", False)),
            "protected_xh_top_coordinates": protected_xh_diag.get("protected_xh_top_coordinates", ""),
            "protected_xh_top1_coord": protected_xh_diag.get("protected_xh_top1_coord", ""),
            "protected_xh_top1_percent": protected_xh_diag.get("protected_xh_top1_percent", 0.0),
            "protected_xh_center_atom": protected_xh_diag.get("protected_xh_center_atom", ""),
            "protected_xh_center_element": protected_xh_diag.get("protected_xh_center_element", ""),
            "protected_xh_pair_coords": protected_xh_diag.get("protected_xh_pair_coords", ""),
            "protected_xh_pair_percent_1": protected_xh_diag.get("protected_xh_pair_percent_1", 0.0),
            "protected_xh_pair_percent_2": protected_xh_diag.get("protected_xh_pair_percent_2", 0.0),
            "protected_xh_balance_ratio": protected_xh_diag.get("protected_xh_balance_ratio", 0.0),
            "functional_group_assignment": assignment,
            "assignment_confidence": conf_label,
            "assignment_confidence_score": conf_score,
            "warnings": all_warnings,
            "assignment_method": "Stage 3D v4.3 weighted independent-coordinate audit with protected X-H stretch audit",
        })

    return pd.DataFrame(rows)


def mode_level_warnings(hess: HessData, fragments, hbonds) -> pd.DataFrame:
    rows = []
    freqs = hess.frequencies_cm1
    for mode, freq in enumerate(freqs):
        warnings = []
        severity = 0
        if freq < -1.0:
            warnings.append("imaginary_or_negative_frequency")
            severity += 4
        if mode >= 6:
            neighbors = [abs(freq - freqs[j]) for j in range(6, len(freqs)) if j != mode]
            if neighbors and min(neighbors) < 10.0:
                warnings.append("near_degenerate_within_10_cm-1")
                severity += 1
        if mode >= 6 and freq < 300.0:
            if len(fragments) >= 2:
                fm_tmp = fragment_motion_metrics(hess, fragments, mode)
                if fm_tmp["fragment_translation_percent"] > 50.0:
                    warnings.append("collective_fragment_motion_candidate")
                else:
                    warnings.append("low_frequency_internal_mixing_candidate")
            else:
                warnings.append("low_frequency_intramolecular_mode")
            severity += 1
        if freq > 4200.0:
            warnings.append("unusually_high_stretch_frequency_check_scaling/units")
            severity += 2
        fm = fragment_motion_metrics(hess, fragments, mode) if len(fragments) >= 2 else {
            "fragment_translation_percent": 0.0,
            "fragment_internal_deformation_percent": 0.0,
            "fragment_mean_relative_amplitude": 0.0,
        }
        rows.append({
            "mode": mode,
            "frequency_cm-1": float(freq),
            "IR_intensity": float(hess.ir_intensities[mode]),
            "region": classify_mode_region(float(freq)),
            "warning_count": len(warnings),
            "severity_score": severity,
            "warnings": "; ".join(warnings),
            **fm,
        })
    return pd.DataFrame(rows)


def confidence_from_general(system_flags: List[str], mode_df: pd.DataFrame, rank_ok: bool, cond_independent: float) -> float:
    score = 100.0
    score -= 18.0 * (not rank_ok)
    if np.isfinite(cond_independent):
        if cond_independent > 1.0e4:
            score -= 12.0
        elif cond_independent > 1.0e3:
            score -= 7.0
        elif cond_independent > 1.0e2:
            score -= 3.0
    else:
        score -= 15.0
    score -= min(20.0, float(mode_df["severity_score"].sum()) * 0.30)
    score -= min(10.0, 2.0 * len(system_flags))
    return max(0.0, min(100.0, round(score, 1)))


def build_sanity_check_monoethanolamine_monomer(
    hess: HessData,
    assignment_df: pd.DataFrame,
    source_label: str,
) -> pd.DataFrame:
    """
    Stage 3D v4.3 sanity benchmark for monoethanolamine monomer.

    This check is intentionally narrow. It only evaluates the validated
    monoethanolamine monomer assignment window when the source has formula
    C2H7NO and 33 normal modes.
    """
    rows: List[Dict[str, object]] = []
    formula = formula_string(hess.atoms)
    is_candidate = formula == "C2H7NO" and len(hess.atoms) == 11 and len(hess.frequencies_cm1) == 33

    expected = {
        25: ("NH2 scissor",),
        26: ("CH2 symmetric stretch",),
        27: ("CH2 symmetric stretch",),
        28: ("CH2 asymmetric stretch",),
        29: ("CH2 asymmetric stretch",),
        30: ("NH2 symmetric stretch",),
        31: ("NH2 asymmetric stretch",),
        32: ("O-H stretch",),
    }

    for mode, accepted in expected.items():
        row = assignment_df.loc[assignment_df["mode"].astype(int) == mode]
        if not is_candidate:
            actual = ""
            freq = np.nan
            conf = ""
            score = np.nan
            verdict = "SKIPPED"
            notes = f"not_monoethanolamine_monomer_candidate; formula={formula}; natoms={len(hess.atoms)}; nmodes={len(hess.frequencies_cm1)}"
        elif row.empty:
            actual = ""
            freq = np.nan
            conf = ""
            score = np.nan
            verdict = "FAIL"
            notes = "mode_missing_from_assignment_audit"
        else:
            r = row.iloc[0]
            actual = str(r.get("functional_group_assignment", ""))
            freq = float(r.get("frequency_cm-1", np.nan))
            conf = str(r.get("assignment_confidence", ""))
            score = float(r.get("assignment_confidence_score", np.nan))
            verdict = "PASS" if any(token in actual for token in accepted) else "FAIL"
            notes = "" if verdict == "PASS" else "assignment_mismatch"

        rows.append({
            "Source": source_label,
            "Filename": hess.filename,
            "sanity_check": "monoethanolamine_monomer_modes_25_32",
            "mode": mode,
            "frequency_cm-1": freq,
            "expected_assignment": " | ".join(accepted),
            "actual_assignment": actual,
            "assignment_confidence": conf,
            "assignment_confidence_score": score,
            "verdict": verdict,
            "notes": notes,
        })

    if rows:
        n_pass = sum(1 for r in rows if r["verdict"] == "PASS")
        n_fail = sum(1 for r in rows if r["verdict"] == "FAIL")
        n_skip = sum(1 for r in rows if r["verdict"] == "SKIPPED")
        overall = "PASS" if is_candidate and n_fail == 0 else ("SKIPPED" if not is_candidate else "FAIL")
        rows.append({
            "Source": source_label,
            "Filename": hess.filename,
            "sanity_check": "monoethanolamine_monomer_modes_25_32_summary",
            "mode": "",
            "frequency_cm-1": "",
            "expected_assignment": "8 benchmark assignments",
            "actual_assignment": f"PASS={n_pass}; FAIL={n_fail}; SKIPPED={n_skip}",
            "assignment_confidence": "",
            "assignment_confidence_score": "",
            "verdict": overall,
            "notes": (
                "benchmark_applicable"
                if is_candidate
                else f"benchmark_not_applicable; formula={formula}; natoms={len(hess.atoms)}; nmodes={len(hess.frequencies_cm1)}"
            ),
        })

    return pd.DataFrame(rows)




def safe_output_stem(name: str) -> str:
    """
    Convert a .hess/.out filename stem into a filesystem-safe output prefix.

    The output prefix intentionally preserves the current molecule/cluster name
    from the input file instead of using a hard-coded molecule name.
    """
    stem = Path(str(name)).name
    for suffix in ("_freq.hess", ".hess", "_freq.out", ".out"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    stem = re.sub(r"[^A-Za-z0-9._+-]+", "_", stem).strip("_")
    return stem or "ORCAVEDA_output"


def output_prefix_for_hess_paths(paths: Sequence[str | Path]) -> str:
    """
    Prefix used for CSV/XLSX/JSON outputs.

    Single input:
        <hess-stem>__assignment_audit.csv

    Multiple inputs:
        <stem1>__<stem2>...__multi_file_N__assignment_audit.csv
    """
    stems = [safe_output_stem(str(p)) for p in paths]
    if not stems:
        return "ORCAVEDA_output"
    if len(stems) == 1:
        return stems[0]
    joined = "__".join(stems[:3])
    if len(stems) > 3:
        joined += f"__plus_{len(stems) - 3}_files"
    return f"{joined}__multi_file_{len(stems)}"

def analyze_general_hess_files(hess_paths: Sequence[str | Path], outdir: str | Path, out_paths: Optional[Sequence[str | Path]] = None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_prefix_for_hess_paths(hess_paths)

    out_map = {}
    if out_paths:
        for p in out_paths:
            p = Path(p)
            stem = p.name.replace("_freq.out", "").replace(".out", "")
            out_map[stem] = p

    source_rows, summary_rows, group_rows, atom_rows = [], [], [], []
    hbond_rows, mode_frames, assignment_frames, basis_rows, selected_rows = [], [], [], [], []
    sanity_rows = []

    for source_index, hpath in enumerate(hess_paths, start=1):
        hpath = Path(hpath)
        hess = read_orca_hess(hpath)
        stem = hpath.name.replace("_freq.hess", "").replace(".hess", "")
        outp = out_map.get(stem, hpath.with_suffix(".out"))
        out_meta = parse_orca_out_metadata(outp) if outp.exists() else {}

        natoms = len(hess.atoms)
        n3 = 3 * natoms
        chemical_annotation: ChemicalSystemAnnotation = annotate_chemical_system(hess.atoms, hess.coords_A)
        bonds = list(chemical_annotation.bonds)
        fragments = [list(fragment) for fragment in chemical_annotation.fragments]
        system_type = chemical_annotation.system_type
        expected_rank = expected_vibrational_rank(natoms, linear=False)
        hbonds = list(chemical_annotation.interfragment_hbonds)
        groups = list(chemical_annotation.functional_groups)

        internals = build_internal_coordinates(hess.atoms, hess.coords_A, bonds, fragments, hbonds, groups)
        B = finite_difference_B(hess.coords_A, internals)
        rank_red, cond_red, _ = svd_rank_condition(B)
        selected_idx, rank_ind, cond_ind, _ = select_independent_coordinates(B, internals, expected_rank)

        mode_df = mode_level_warnings(hess, fragments, hbonds)
        mode_df.insert(0, "Filename", hess.filename)
        mode_df.insert(0, "Source", f"[{source_index}]")
        mode_frames.append(mode_df)

        assignment_df = build_stage3d_assignment_audit(
            hess,
            internals,
            B,
            selected_idx,
            mode_df,
            source_label=f"[{source_index}]",
            top_n=6,
        )
        assignment_frames.append(assignment_df)

        sanity_df = build_sanity_check_monoethanolamine_monomer(
            hess,
            assignment_df,
            f"[{source_index}]",
        )
        if not sanity_df.empty:
            sanity_rows.extend(sanity_df.to_dict("records"))

        unsupported = sorted(set(hess.atoms) - chemistry_get_supported_elements())
        flags = []
        if unsupported:
            flags.append("unsupported_elements:" + ",".join(unsupported))
        if rank_red < expected_rank:
            flags.append("redundant_B_rank_below_expected")
        if rank_ind != expected_rank:
            flags.append("independent_basis_rank_not_expected")
        if np.any(hess.frequencies_cm1[6:] < -1.0):
            flags.append("negative_vibrational_frequency_after_first_6")
        if not np.all(np.isfinite(hess.ir_intensities)):
            flags.append("nonfinite_IR_intensity")
        if len(groups) == 0:
            flags.append("no_functional_group_detected")
        if system_type != "monomer" and len(hbonds) == 0:
            flags.append("cluster_without_interfragment_hbond_detected")

        confidence = confidence_from_general(flags, mode_df, rank_ind == expected_rank, cond_ind)

        source_rows.append({
            "Source": f"[{source_index}]",
            "Filename": hess.filename,
            "Paired_out_file": outp.name if outp.exists() else "",
            "Language": "ORCA .hess numeric text",
            "File_type": "DATASET/HESSIAN",
            "Completeness": "complete_required_sections" if True else "not_checked",
            "natoms": natoms,
            "formula": chemical_annotation.formula,
            "system_type": system_type,
            "route_line_out": out_meta.get("route_line", ""),
            "cpcm_epsilon_out": out_meta.get("cpcm_epsilon", ""),
            "temperature_K_hess": hess.temperature_K if hess.temperature_K is not None else "",
            "temperature_K_out": out_meta.get("temperature_K_out", ""),
        })

        summary_rows.append({
            "Source": f"[{source_index}]",
            "Filename": hess.filename,
            "formula": chemical_annotation.formula,
            "natoms": natoms,
            "3N": n3,
            "system_type": system_type,
            "fragments": len(fragments),
            "fragment_sizes": ";".join(str(size) for size in chemical_annotation.fragment_sizes),
            "bonds_detected": len(bonds),
            "functional_groups_detected": "; ".join(chemical_annotation.functional_group_labels),
            "functional_group_count": len(groups),
            "interfragment_hbond_count": len(hbonds),
            "internal_coordinates_redundant": len(internals),
            "expected_rank_3N_minus_6": expected_rank,
            "rank_B_redundant": rank_red,
            "condition_B_redundant": cond_red,
            "rank_B_independent": rank_ind,
            "condition_B_independent": cond_ind,
            "selected_independent_coordinates": len(selected_idx),
            "negative_freq_count_after_first_6": int(np.sum(hess.frequencies_cm1[6:] < -1.0)),
            "near_degenerate_modes_count": int(mode_df["warnings"].str.contains("near_degenerate", regex=False).sum()),
            "confidence_score_0_100": confidence,
            "system_flags": "; ".join(flags),
            "normal_mode_orientation_rule": "PASS: uses normal_modes[:, mode].reshape(natoms, 3)",
        })

        for g in groups:
            group_rows.append({
                "Source": f"[{source_index}]",
                "Filename": hess.filename,
                "group": g.group,
                "atoms_1based": "-".join(str(i+1) for i in g.atoms0),
                "description": g.description,
                "confidence": g.confidence,
                "evidence": g.evidence,
            })

        env_df = atom_environment_table(hess.atoms, hess.coords_A, bonds)
        env_df.insert(0, "Filename", hess.filename)
        env_df.insert(0, "Source", f"[{source_index}]")
        atom_rows.extend(env_df.to_dict("records"))

        for h in hbonds:
            hbond_rows.append({
                "Source": f"[{source_index}]", "Filename": hess.filename,
                "type": h["type"], "D": h["D"], "H": h["H"], "A": h["A"],
                "rHA_A": h["rHA_A"], "rDA_A": h["rDA_A"], "angle_deg": h["angle_deg"],
                "chem_type": h.get("chem_type", ""),
                "context_label": h.get("context_label", ""),
                "donor_group": h.get("donor_group", ""),
                "acceptor_group": h.get("acceptor_group", ""),
                "fragment_pair": h.get("fragment_pair", ""),
            })

        for i, ic in enumerate(internals):
            basis_rows.append({
                "Source": f"[{source_index}]", "Filename": hess.filename, "coord_index": i,
                "name": ic.name, "kind": ic.kind, "atoms_1based": "-".join(str(a+1) for a in ic.atoms0),
                "priority": ic.priority, "source": ic.source,
            })
        for order, idx in enumerate(selected_idx, start=1):
            ic = internals[idx]
            selected_rows.append({
                "Source": f"[{source_index}]", "Filename": hess.filename, "selection_order": order,
                "coord_index_redundant": idx, "name": ic.name, "kind": ic.kind,
                "atoms_1based": "-".join(str(a+1) for a in ic.atoms0),
                "priority": ic.priority, "source": ic.source,
            })

    tables = {
        "source_map": pd.DataFrame(source_rows),
        "general_summary": pd.DataFrame(summary_rows),
        "functional_groups": pd.DataFrame(group_rows),
        "atom_environments": pd.DataFrame(atom_rows),
        "interfragment_hbonds": pd.DataFrame(hbond_rows),
        "mode_warnings": pd.concat(mode_frames, ignore_index=True) if mode_frames else pd.DataFrame(),
        "assignment_audit": pd.concat(assignment_frames, ignore_index=True) if assignment_frames else pd.DataFrame(),
        "sanity_check": pd.DataFrame(sanity_rows),
        "redundant_basis": pd.DataFrame(basis_rows),
        "independent_basis": pd.DataFrame(selected_rows),
    }

    for name, df in tables.items():
        df.to_csv(outdir / f"{output_prefix}__{name}.csv", index=False)

    (outdir / f"{output_prefix}__general_engine_manifest.json").write_text(
        json.dumps({
            "status": "Stage 3A general organic engine + Stage 3D assignment audit completed; output filenames are prefixed by input .hess stem",
            "tables": list(tables),
            "chemistry_backend": chemistry_get_active_backend_name(),
            "normal_mode_orientation_rule": "normal_modes[:, mode].reshape(natoms, 3)",
        }, indent=2),
        encoding="utf-8"
    )
    return tables



def normalize_sheet_name(name: str) -> str:
    """Excel sheet names are limited to 31 characters and cannot contain []:*?/\\."""
    bad = set('[]:*?/\\')
    clean = ''.join('_' if c in bad else c for c in name)
    return clean[:31] if len(clean) > 31 else clean




def write_xlsx_report(report_tables: Dict[str, pd.DataFrame], xlsx_path: str | Path) -> Path:
    """
    Write all report tables into one XLSX workbook.

    Uses xlsxwriter when available because it supports formatting.
    Falls back to openpyxl when xlsxwriter is absent.
    If neither package is available, CSV files are still written by upstream functions
    and XLSX export is skipped without aborting the analysis.
    """
    xlsx_path = Path(xlsx_path)
    xlsx_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        import xlsxwriter  # noqa: F401
        engine = "xlsxwriter"
    except ModuleNotFoundError:
        try:
            import openpyxl  # noqa: F401
            engine = "openpyxl"
        except ModuleNotFoundError:
            print("WARNING: neither xlsxwriter nor openpyxl is installed; XLSX report skipped.")
            return xlsx_path

    with pd.ExcelWriter(xlsx_path, engine=engine) as writer:
        if engine == "xlsxwriter":
            workbook = writer.book
            header_fmt = workbook.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1})
            warning_fmt = workbook.add_format({"bg_color": "#FFF2CC"})
            critical_fmt = workbook.add_format({"bg_color": "#F4CCCC"})
            number_fmt = workbook.add_format({"num_format": "0.000"})
            sci_fmt = workbook.add_format({"num_format": "0.00E+00"})
        else:
            header_fmt = warning_fmt = critical_fmt = number_fmt = sci_fmt = None

        for raw_name, df in report_tables.items():
            if df is None:
                continue
            if not isinstance(df, pd.DataFrame):
                df = pd.DataFrame(df)

            sheet = normalize_sheet_name(str(raw_name).replace(".csv", ""))
            df.to_excel(writer, index=False, sheet_name=sheet)

            # Formatting is xlsxwriter-specific. openpyxl fallback writes plain sheets.
            if engine != "xlsxwriter":
                continue

            ws = writer.sheets[sheet]

            for col_idx, col_name in enumerate(df.columns):
                ws.write(0, col_idx, col_name, header_fmt)
            ws.freeze_panes(1, 0)
            ws.autofilter(0, 0, max(len(df), 1), max(len(df.columns) - 1, 0))

            for col_idx, col_name in enumerate(df.columns):
                series = df[col_name].astype(str) if len(df) else pd.Series([str(col_name)])
                width = min(max(len(str(col_name)), int(series.str.len().quantile(0.90)) if len(series) else 10) + 2, 60)
                ws.set_column(col_idx, col_idx, width)

            for col_idx, col_name in enumerate(df.columns):
                low = str(col_name).lower()
                if "condition" in low:
                    ws.set_column(col_idx, col_idx, 14, sci_fmt)
                elif any(key in low for key in ["score", "freq", "intensity", "angle", "rha", "rda"]):
                    ws.set_column(col_idx, col_idx, 14, number_fmt)

            if "warnings" in df.columns and len(df):
                wcol = df.columns.get_loc("warnings")
                ws.conditional_format(1, wcol, len(df), wcol, {
                    "type": "text", "criteria": "containing", "value": "negative", "format": critical_fmt
                })
                ws.conditional_format(1, wcol, len(df), wcol, {
                    "type": "text", "criteria": "containing", "value": "near_degenerate", "format": warning_fmt
                })

            if "system_flags" in df.columns and len(df):
                fcol = df.columns.get_loc("system_flags")
                ws.conditional_format(1, fcol, len(df), fcol, {
                    "type": "text", "criteria": "not containing", "value": "", "format": warning_fmt
                })

    return xlsx_path




# =============================================================================
# Stage 3C: Mode tracking integration
# =============================================================================

def _mode_matrix_for_tracking(hess: HessData) -> np.ndarray:
    """
    Return normal modes as shape (n_modes, natoms, 3).

    Critical ORCA .hess convention:
        normal_modes[:, mode] is the Cartesian displacement vector.
    """
    nm = np.asarray(hess.normal_modes, dtype=float)
    nat = len(hess.atoms)
    if nm.shape[0] != 3 * nat:
        raise ValueError(f"{hess.filename}: normal_modes first dimension {nm.shape[0]} != 3N {3*nat}")
    return np.stack([nm[:, i].reshape(nat, 3) for i in range(nm.shape[1])], axis=0)


def _normalize_tracking_vector(vec: np.ndarray, masses: Optional[np.ndarray] = None, mass_weighted: bool = False) -> np.ndarray:
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
    """
    Rotation matrix aligning target coordinates to reference coordinates.
    Translation is handled separately by centroid subtraction.
    """
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
    """
    Compute |Q_i(ref)^T Q_j(target)| for same-size systems.

    If align=True, target displacement vectors are rotated using Kabsch alignment
    of target coordinates onto reference coordinates.
    """
    ok, reason = same_size_tracking_compatible(reference, target, require_same_symbols=require_same_symbols)
    if not ok:
        raise ValueError(f"Mode tracking incompatible: {reason}")

    ref_modes = _mode_matrix_for_tracking(reference)
    tgt_modes = _mode_matrix_for_tracking(target)

    if align:
        rot = kabsch_rotation(reference.coords_A, target.coords_A)
        tgt_modes = np.einsum("mni,ij->mnj", tgt_modes, rot)

    nref = ref_modes.shape[0]
    ntgt = tgt_modes.shape[0]
    s = np.zeros((nref, ntgt), dtype=float)

    for i in range(nref):
        qi = _normalize_tracking_vector(ref_modes[i], reference.masses, mass_weighted).ravel()
        for j in range(ntgt):
            qj = _normalize_tracking_vector(tgt_modes[j], target.masses, mass_weighted).ravel()
            s[i, j] = abs(float(np.dot(qi, qj)))
    return s


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
            "frequency_shift_cm-1": (
                float(target.frequencies_cm1[j] - reference.frequencies_cm1[i])
                if i < len(reference.frequencies_cm1) and j < len(target.frequencies_cm1) else np.nan
            ),
            "best_overlap": float(overlap[i, j]),
            "second_best_target_mode": int(second),
            "second_best_overlap": float(overlap[i, second]),
            "overlap_gap": float(overlap[i, j] - overlap[i, second]),
            "reciprocal_best_match": bool(best_ref_for_target[j] == i),
            "tracking_confidence": (
                "high" if overlap[i, j] >= 0.75 and (overlap[i, j] - overlap[i, second]) >= 0.15 else
                "medium" if overlap[i, j] >= 0.50 else
                "low"
            ),
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


def _fragment_same_symbol_order(reference: HessData, target: HessData, fragment: Sequence[int]) -> bool:
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
    """
    Project reference modes onto one same-size fragment of a larger target system.

    This is intended for monomer -> dimer/cluster mode inheritance.
    It requires identical atom-symbol order within the selected target fragment.
    """
    fragment = list(fragment)
    if not _fragment_same_symbol_order(reference, target, fragment):
        raise ValueError("Target fragment is not compatible with reference atom count/symbol order.")

    ref_modes = _mode_matrix_for_tracking(reference)
    tgt_modes_all = _mode_matrix_for_tracking(target)
    tgt_modes = tgt_modes_all[:, fragment, :]

    ref_coords = reference.coords_A
    tgt_coords = target.coords_A[fragment, :]
    if align:
        rot = kabsch_rotation(ref_coords, tgt_coords)
        tgt_modes = np.einsum("mni,ij->mnj", tgt_modes, rot)

    target_masses_fragment = target.masses[fragment]
    s = np.zeros((ref_modes.shape[0], tgt_modes.shape[0]), dtype=float)
    for i in range(ref_modes.shape[0]):
        qi = _normalize_tracking_vector(ref_modes[i], reference.masses, mass_weighted).ravel()
        for j in range(tgt_modes.shape[0]):
            qj = _normalize_tracking_vector(tgt_modes[j], target_masses_fragment, mass_weighted).ravel()
            s[i, j] = abs(float(np.dot(qi, qj)))
    return s


def mode_tracking_outputs_for_hess_files(
    paths: Sequence[str | Path],
    outdir: str | Path,
    *,
    mass_weighted: bool = True,
    align: bool = True,
    include_overlap_matrices: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Generate mode-tracking outputs for all useful file pairs.

    Same-size tracking is attempted for pairs with identical atom count/symbol order.
    Fragment-projected tracking is attempted from smaller systems to larger systems
    when a target fragment has the same atom count and atom-symbol order.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_prefix_for_hess_paths(paths)
    hess_list = [read_orca_hess(p) for p in paths]

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
                except Exception as exc:
                    summary_rows.append({
                        "pair": pair_label,
                        "tracking_type": "same_size",
                        "reference_file": ref.filename,
                        "target_file": tgt.filename,
                        "status": "ERROR",
                        "notes": str(exc),
                    })
            elif a_idx < b_idx:
                summary_rows.append({
                    "pair": pair_label,
                    "tracking_type": "same_size",
                    "reference_file": ref.filename,
                    "target_file": tgt.filename,
                    "status": "SKIPPED",
                    "notes": reason,
                })

            # Smaller reference -> larger target fragment-projected tracking.
            if len(ref.atoms) < len(tgt.atoms):
                try:
                    bonds = build_connectivity(tgt.atoms, tgt.coords_A)
                    frags = split_fragments(len(tgt.atoms), bonds)
                    compatible_frags = [frag for frag in frags if _fragment_same_symbol_order(ref, tgt, frag)]
                    if not compatible_frags:
                        summary_rows.append({
                            "pair": pair_label,
                            "tracking_type": "fragment_projected",
                            "reference_file": ref.filename,
                            "target_file": tgt.filename,
                            "status": "SKIPPED",
                            "notes": "No target fragment with same atom count/symbol order.",
                        })
                    for frag_idx, frag in enumerate(compatible_frags):
                        ov = fragment_projected_overlap_matrix(ref, tgt, frag, mass_weighted=mass_weighted, align=align)
                        frag_label = f"{pair_label}__fragment_{frag_idx}"
                        tr = mode_tracking_table(ref, tgt, ov, frag_label)
                        warn = mode_mixing_warnings_table(ref, tgt, ov, frag_label)
                        tracking_frames.append(tr.assign(
                            tracking_type="fragment_projected",
                            target_fragment_index=frag_idx,
                            target_fragment_atoms=";".join(map(str, frag)),
                        ))
                        if not warn.empty:
                            warning_frames.append(warn.assign(
                                tracking_type="fragment_projected",
                                target_fragment_index=frag_idx,
                                target_fragment_atoms=";".join(map(str, frag)),
                            ))
                        if include_overlap_matrices:
                            overlap_frames.append(mode_overlap_matrix_table(ref, tgt, ov, frag_label).assign(
                                tracking_type="fragment_projected",
                                target_fragment_index=frag_idx,
                                target_fragment_atoms=";".join(map(str, frag)),
                            ))
                        summary_rows.append({
                            "pair": frag_label,
                            "tracking_type": "fragment_projected",
                            "reference_file": ref.filename,
                            "target_file": tgt.filename,
                            "target_fragment_index": frag_idx,
                            "n_reference_modes": ov.shape[0],
                            "n_target_modes": ov.shape[1],
                            "median_best_overlap": float(np.median(np.max(ov, axis=1))),
                            "status": "OK",
                            "notes": "",
                        })
                except Exception as exc:
                    summary_rows.append({
                        "pair": pair_label,
                        "tracking_type": "fragment_projected",
                        "reference_file": ref.filename,
                        "target_file": tgt.filename,
                        "status": "ERROR",
                        "notes": str(exc),
                    })

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


# Rebind extracted modules after legacy in-file definitions so downstream code
# and external imports use the dedicated modules without changing call sites.
distance_fn = internal_distance_fn
angle_fn = internal_angle_fn
dihedral_rad = internal_dihedral_rad
torsion_fn = internal_torsion_fn
build_internal_coordinates = internal_build_internal_coordinates
finite_difference_B = bmatrix_finite_difference_B
svd_rank_condition = bmatrix_svd_rank_condition
select_independent_coordinates = bmatrix_select_independent_coordinates
safe_output_stem = reports_safe_output_stem
output_prefix_for_hess_paths = reports_output_prefix_for_hess_paths
normalize_sheet_name = reports_normalize_sheet_name
write_xlsx_report = reports_write_xlsx_report
mode_tracking_outputs_for_hess_files = tracking_mode_tracking_outputs_for_hess_files

def general_outputs_for_hess_files(paths: Sequence[str | Path], outdir: str | Path, out_paths: Optional[Sequence[str | Path]] = None) -> Dict[str, pd.DataFrame]:
    """
    Pipeline hook for the main PED workflow.

    This computes general organic recognition/validation tables and writes CSV outputs.
    Use the returned dictionary together with PED tables:
        all_tables = {**ped_tables, **general_outputs_for_hess_files(hess_paths, outdir, out_paths)}
        write_xlsx_report(all_tables, "orca_ped_like_report.xlsx")
    """
    return analyze_general_hess_files(paths, outdir, out_paths)


def analyze_orca_ped_like(paths: Sequence[str | Path], outdir: str | Path, out_paths: Optional[Sequence[str | Path]] = None) -> Dict[str, pd.DataFrame]:
    """
    Integrated entry point for the current development version.

    Current integrated layers:
        - ORCA .hess parser
        - connectivity / fragment / system classification
        - functional-group recognition
        - automatic functional-group local coordinate templates
        - redundant and independent B-matrix diagnostics
        - confidence diagnostics
        - XLSX export hook

    In the complete PED workflow, call this after or alongside true PED matrix
    generation so assignments can use functional-group template coordinates.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_prefix_for_hess_paths(paths)
    tables = general_outputs_for_hess_files(paths, outdir, out_paths)

    # Stage 3C mode tracking is meaningful only when two or more .hess files are supplied.
    if len(paths) >= 2:
        tracking_tables = mode_tracking_outputs_for_hess_files(
            paths,
            outdir,
            mass_weighted=True,
            align=True,
            include_overlap_matrices=False,
        )
        tables.update(tracking_tables)

    xlsx_path = write_xlsx_report(tables, outdir / f"{output_prefix}__orca_ped_like_stage3C_integrated_report.xlsx")
    hess_list = [read_orca_hess(p) for p in paths]
    spectrum_payload = reports_build_spectrum_payload(hess_list, tables.get("assignment_audit"))
    spectrum_json_path = outdir / f"{output_prefix}__spectrum_data.json"
    spectrum_html_path = reports_write_interactive_spectrum_viewer(
        spectrum_payload,
        outdir / f"{output_prefix}__interactive_spectrum.html",
        json_path=spectrum_json_path,
    )
    manifest = {
        "integration_status": "Stage 3A general organic engine + Stage 3C mode tracking + Stage 3D assignment audit integrated into ORCAVEDA",
        "xlsx_report": str(xlsx_path),
        "interactive_spectrum_html": str(spectrum_html_path),
        "interactive_spectrum_data_json": str(spectrum_json_path),
        "tables": list(tables.keys()),
        "chemistry_backend": chemistry_get_active_backend_name(),
        "functional_group_templates": "automatic in build_internal_coordinates(..., groups=detect_functional_groups(...))",
        "mode_tracking": "automatic when len(hess_paths) >= 2; same-size and fragment-projected tracking attempted",
        "stage3d_assignment_audit": "prefixed assignment_audit table generated; filenames use current .hess stem",
        "normal_mode_orientation_rule": "normal_modes[:, mode].reshape(natoms, 3)",
    }
    (outdir / f"{output_prefix}__integration_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return tables


# Backward-compatible aliases.
run_orca_ped_like = analyze_orca_ped_like
run_orca_ped_like_with_general_engine = analyze_orca_ped_like


def cli_main():
    chemistry_set_active_backend_from_env()
    external_cli_main(
        run_orca_ped_like,
        set_chem_backend=chemistry_set_active_backend,
        list_chem_backends=chemistry_list_backends,
        default_chem_backend=chemistry_get_active_backend_name,
    )


def colab_upload_and_run():
    external_colab_upload_and_run(run_orca_ped_like)


def is_google_colab():
    return external_is_google_colab()




# =============================================================================
# Stage 3D v4.4 hotfix: direct Cartesian protected X-H stretch fallback
# =============================================================================

_stage3d_v43_build_stage3d_assignment_audit = build_stage3d_assignment_audit


def _stage3d_direct_xh_stretch_fallback(hess: HessData, internals: Sequence[InternalCoordinate], mode: int) -> Dict[str, object]:
    """
    Direct Cartesian X-H stretch audit independent of selected internal-coordinate basis.

    This is a protected high-frequency fallback for cases where rank-selected internal
    coordinates give diffuse/zero contributions although an X-H stretching normal-mode
    displacement is present.

    Projection for each X-H pair:
        p_XH = ((dH - dX) dot u_XH)^2
    where u_XH is the equilibrium unit vector from heavy atom X to H.
    """
    empty = {
        "assignment": "",
        "used": False,
        "top_coordinates": "",
        "top1_coord": "",
        "top1_percent": 0.0,
        "center_atom": "",
        "center_element": "",
        "pair_coords": "",
        "pair_percent_1": 0.0,
        "pair_percent_2": 0.0,
        "balance_ratio": 0.0,
        "score": 0.0,
    }
    try:
        disp = hess.normal_modes[:, mode].reshape(len(hess.atoms), 3)
    except Exception:
        return empty.copy()

    pair_best: Dict[Tuple[int, int], Dict[str, object]] = {}
    for ic in internals:
        if len(ic.atoms0) != 2:
            continue
        a, b = ic.atoms0
        ea, eb = hess.atoms[a], hess.atoms[b]
        if ea == "H" and eb in ("C", "N", "O"):
            H, X = a, b
        elif eb == "H" and ea in ("C", "N", "O"):
            H, X = b, a
        else:
            continue

        vec = hess.coords_A[H] - hess.coords_A[X]
        norm = float(np.linalg.norm(vec))
        if norm <= 0.0 or not np.isfinite(norm):
            continue
        unit = vec / norm
        signed = float(np.dot(disp[H] - disp[X], unit))
        weight = signed * signed
        if not np.isfinite(weight) or weight <= 0.0:
            continue

        key = tuple(sorted((X, H)))
        label = _compact_coord_label(ic.name)
        old = pair_best.get(key)
        if old is None or weight > float(old["weight"]):
            pair_best[key] = {
                "heavy": X,
                "h": H,
                "element": hess.atoms[X],
                "coord": label,
                "signed": signed,
                "weight": weight,
            }

    candidates = list(pair_best.values())
    total = float(sum(float(c["weight"]) for c in candidates))
    if total <= 0.0 or not np.isfinite(total):
        return empty.copy()

    for c in candidates:
        c["percent"] = 100.0 * float(c["weight"]) / total
    top = sorted(candidates, key=lambda c: float(c["percent"]), reverse=True)
    top1 = top[0]
    top_terms = "; ".join(f"{c['coord']}={float(c['percent']):.1f}%" for c in top[:8])

    assignment = {
        "O": "O-H stretch",
        "N": "N-H stretch",
        "C": "C-H stretch",
    }.get(str(top1["element"]), "X-H stretch")

    # XH2 symmetry label when the two strongest X-H coordinates share the same C/N center.
    pair_coords = ""
    p1 = float(top1["percent"])
    p2 = 0.0
    balance = 0.0
    if len(top) >= 2:
        p2 = float(top[1]["percent"])
        pair_coords = f"{top[0]['coord']} | {top[1]['coord']}"
        balance = min(p1, p2) / max(p1, p2) if max(p1, p2) > 0 else 0.0
        if int(top[0]["heavy"]) == int(top[1]["heavy"]) and str(top[0]["element"]) in ("C", "N"):
            product = float(top[0]["signed"]) * float(top[1]["signed"])
            if str(top[0]["element"]) == "C":
                c_h_count = _stage3d_xh_center_h_count(internals, hess.atoms, int(top[0]["heavy"]), "C")
                assignment = (
                    "CH2 symmetric stretch" if product > 0 else "CH2 asymmetric stretch"
                ) if c_h_count == 2 else "C-H stretch"
            elif str(top[0]["element"]) == "N":
                assignment = "NH2 symmetric stretch" if product > 0 else "NH2 asymmetric stretch"
    else:
        pair_coords = str(top1["coord"])

    return {
        "assignment": assignment,
        "used": True,
        "top_coordinates": top_terms,
        "top1_coord": str(top1["coord"]),
        "top1_percent": round(float(top1["percent"]), 3),
        "center_atom": int(top1["heavy"]) + 1,
        "center_element": str(top1["element"]),
        "pair_coords": pair_coords,
        "pair_percent_1": round(float(p1), 3),
        "pair_percent_2": round(float(p2), 3),
        "balance_ratio": round(float(balance), 3),
        "score": round(float(top1["percent"]), 1),
    }


def build_stage3d_assignment_audit(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    mode_df: pd.DataFrame,
    source_label: str,
    top_n: int = 8,
) -> pd.DataFrame:
    """
    Stage 3D v4.4 wrapper.

    Starts from v4.3 independent-coordinate audit, then repairs high-frequency
    unassigned X-H modes using direct Cartesian X-H stretch projections over the
    redundant/local coordinate pool.
    """
    df = _stage3d_v43_build_stage3d_assignment_audit(
        hess=hess,
        internals=internals,
        B=B,
        selected_idx=selected_idx,
        mode_df=mode_df,
        source_label=source_label,
        top_n=top_n,
    )

    for ridx, row in df.iterrows():
        freq = float(row.get("frequency_cm-1", 0.0))
        assignment = str(row.get("functional_group_assignment", ""))
        warnings = "" if pd.isna(row.get("warnings", "")) else str(row.get("warnings", ""))
        needs_repair = (
            freq >= 2800.0
            and (
                assignment == "unassigned"
                or "high_frequency_mode_without_clear_stretch_coordinate" in warnings
                or bool(row.get("protected_xh_used", False)) is False and float(row.get("top1_percent", 0.0) or 0.0) <= 0.0
            )
        )
        if not needs_repair:
            continue

        mode = int(row["mode"])
        fb = _stage3d_direct_xh_stretch_fallback(hess, internals, mode)
        if not fb["used"]:
            continue

        df.at[ridx, "protected_xh_assignment"] = fb["assignment"]
        df.at[ridx, "protected_xh_used"] = True
        df.at[ridx, "protected_xh_top_coordinates"] = fb["top_coordinates"]
        df.at[ridx, "protected_xh_top1_coord"] = fb["top1_coord"]
        df.at[ridx, "protected_xh_top1_percent"] = fb["top1_percent"]
        df.at[ridx, "protected_xh_center_atom"] = fb["center_atom"]
        df.at[ridx, "protected_xh_center_element"] = fb["center_element"]
        df.at[ridx, "protected_xh_pair_coords"] = fb["pair_coords"]
        df.at[ridx, "protected_xh_pair_percent_1"] = fb["pair_percent_1"]
        df.at[ridx, "protected_xh_pair_percent_2"] = fb["pair_percent_2"]
        df.at[ridx, "protected_xh_balance_ratio"] = fb["balance_ratio"]

        df.at[ridx, "functional_group_assignment"] = fb["assignment"]
        df.at[ridx, "assignment_confidence"] = "high" if float(fb["score"]) >= 75.0 else ("medium" if float(fb["score"]) >= 50.0 else "low")
        df.at[ridx, "assignment_confidence_score"] = fb["score"]

        new_warnings = [w.strip() for w in warnings.split(";") if w.strip() and w.strip() != "high_frequency_mode_without_clear_stretch_coordinate"]
        if "direct_cartesian_protected_XH_stretch_assignment_used" not in new_warnings:
            new_warnings.append("direct_cartesian_protected_XH_stretch_assignment_used")
        df.at[ridx, "warnings"] = "; ".join(new_warnings)

    df["assignment_method"] = "Stage 3D v4.4 weighted independent-coordinate audit with direct Cartesian protected X-H stretch fallback"
    return df




# =============================================================================
# Stage 3D v4.6 hotfix: topology-based direct Cartesian X-H fallback
# =============================================================================
#
# v4.4 attempted direct Cartesian X-H repair through the internal-coordinate list.
# v4.5 removes that dependency: every covalent X-H bond detected from topology is
# evaluated directly from Cartesian normal-mode displacement.
#
# This remains an assignment-audit fallback, not strict PED.

_stage3d_v44_build_stage3d_assignment_audit = build_stage3d_assignment_audit


def _stage3d_topology_direct_xh_fallback(hess: HessData, mode: int) -> Dict[str, object]:
    empty = {
        "assignment": "",
        "used": False,
        "top_coordinates": "",
        "top1_coord": "",
        "top1_percent": 0.0,
        "center_atom": "",
        "center_element": "",
        "pair_coords": "",
        "pair_percent_1": 0.0,
        "pair_percent_2": 0.0,
        "balance_ratio": 0.0,
        "score": 0.0,
        "normal_mode_norm": 0.0,
        "xh_total_power": 0.0,
    }

    try:
        disp = hess.normal_modes[:, mode].reshape(len(hess.atoms), 3)
    except Exception:
        return empty.copy()

    mode_norm = float(np.linalg.norm(disp))
    if not np.isfinite(mode_norm) or mode_norm <= 0.0:
        out = empty.copy()
        out["normal_mode_norm"] = mode_norm
        return out

    try:
        bonds = build_connectivity(hess.atoms, hess.coords_A)
    except Exception:
        bonds = []

    rows = []
    for a, b, _dist in bonds:
        ea, eb = hess.atoms[a], hess.atoms[b]
        if ea == "H" and eb in ("C", "N", "O"):
            H, X = a, b
        elif eb == "H" and ea in ("C", "N", "O"):
            H, X = b, a
        else:
            continue

        bond_vec = hess.coords_A[H] - hess.coords_A[X]
        bond_norm = float(np.linalg.norm(bond_vec))
        if not np.isfinite(bond_norm) or bond_norm <= 0.0:
            continue

        unit = bond_vec / bond_norm
        rel_disp = disp[H] - disp[X]
        projection = float(np.dot(rel_disp, unit))
        power = projection * projection
        coord = f"direct_cart_XH_stretch({hess.atoms[X]}{X+1}-{hess.atoms[H]}{H+1})"

        rows.append({
            "coord": coord,
            "center_atom": X + 1,
            "center_element": hess.atoms[X],
            "H_atom": H + 1,
            "power": power,
        })

    total = float(sum(r["power"] for r in rows))
    if total <= 0.0 or not np.isfinite(total):
        out = empty.copy()
        out["normal_mode_norm"] = mode_norm
        out["xh_total_power"] = total
        return out

    rows = sorted(rows, key=lambda r: r["power"], reverse=True)
    for r in rows:
        r["percent"] = 100.0 * float(r["power"]) / total

    top1 = rows[0]
    assignment = {
        "O": "O-H stretch",
        "N": "N-H stretch",
        "C": "C-H stretch",
    }.get(str(top1["center_element"]), "X-H stretch")

    # Pair diagnostics for C/N XH2 centers.
    same_center = [r for r in rows if r["center_atom"] == top1["center_atom"]]
    same_center = sorted(same_center, key=lambda r: r["percent"], reverse=True)

    p1 = float(top1["percent"])
    p2 = float(same_center[1]["percent"]) if len(same_center) >= 2 else 0.0
    balance = min(p1, p2) / max(p1, p2) if max(p1, p2) > 0.0 else 0.0
    pair_coords = ""
    if len(same_center) >= 2:
        pair_coords = f"{same_center[0]['coord']} | {same_center[1]['coord']}"

    if str(top1["center_element"]) in ("C", "N") and len(same_center) >= 2:
        # Determine symmetric/asymmetric from signs of the two strongest direct projections.
        # Recompute signed projections for the pair.
        signs = []
        for r in same_center[:2]:
            # parse atom indices from stored center/H fields
            X = int(r["center_atom"]) - 1
            H = int(r["H_atom"]) - 1
            unit = (hess.coords_A[H] - hess.coords_A[X])
            unit = unit / float(np.linalg.norm(unit))
            signs.append(float(np.dot(disp[H] - disp[X], unit)))
        symmetry = "symmetric" if signs[0] * signs[1] >= 0.0 else "asymmetric"
        assignment = f"{'CH2' if top1['center_element'] == 'C' else 'NH2'} {symmetry} stretch"

    return {
        "assignment": assignment,
        "used": True,
        "top_coordinates": "; ".join(f"{r['coord']}={r['percent']:.1f}%" for r in rows[:8]),
        "top1_coord": str(top1["coord"]),
        "top1_percent": round(p1, 3),
        "center_atom": int(top1["center_atom"]),
        "center_element": str(top1["center_element"]),
        "pair_coords": pair_coords,
        "pair_percent_1": round(p1, 3),
        "pair_percent_2": round(p2, 3),
        "balance_ratio": round(float(balance), 3),
        "score": round(p1, 1),
        "normal_mode_norm": mode_norm,
        "xh_total_power": total,
    }


def build_stage3d_assignment_audit(
    hess: HessData,
    internals: Sequence[InternalCoordinate],
    B: np.ndarray,
    selected_idx: Sequence[int],
    mode_df: pd.DataFrame,
    source_label: str,
    top_n: int = 8,
) -> pd.DataFrame:
    df = _stage3d_v44_build_stage3d_assignment_audit(
        hess=hess,
        internals=internals,
        B=B,
        selected_idx=selected_idx,
        mode_df=mode_df,
        source_label=source_label,
        top_n=top_n,
    )

    for extra_col, default in [
        ("protected_xh_normal_mode_norm", pd.NA),
        ("protected_xh_total_power", pd.NA),
    ]:
        if extra_col not in df.columns:
            df[extra_col] = default

    for ridx, row in df.iterrows():
        freq = float(row.get("frequency_cm-1", 0.0))
        assignment = "" if pd.isna(row.get("functional_group_assignment", "")) else str(row.get("functional_group_assignment", ""))
        warnings = "" if pd.isna(row.get("warnings", "")) else str(row.get("warnings", ""))

        needs_repair = (
            freq >= 2800.0
            and (
                assignment == "unassigned"
                or "high_frequency_mode_without_clear_stretch_coordinate" in warnings
                or (not bool(row.get("protected_xh_used", False)) and float(row.get("top1_percent", 0.0) or 0.0) <= 0.0)
            )
        )
        if not needs_repair:
            continue

        mode = int(row["mode"])
        fb = _stage3d_topology_direct_xh_fallback(hess, mode)

        df.at[ridx, "protected_xh_normal_mode_norm"] = fb.get("normal_mode_norm", 0.0)
        df.at[ridx, "protected_xh_total_power"] = fb.get("xh_total_power", 0.0)

        if not fb["used"]:
            new_warnings = [w.strip() for w in warnings.split(";") if w.strip()]
            if fb.get("normal_mode_norm", 0.0) <= 0.0 and "normal_mode_vector_zero_or_unreadable" not in new_warnings:
                new_warnings.append("normal_mode_vector_zero_or_unreadable")
            if "topology_direct_XH_fallback_failed" not in new_warnings:
                new_warnings.append("topology_direct_XH_fallback_failed")
            df.at[ridx, "warnings"] = "; ".join(new_warnings)
            continue

        df.at[ridx, "protected_xh_assignment"] = fb["assignment"]
        df.at[ridx, "protected_xh_used"] = True
        df.at[ridx, "protected_xh_top_coordinates"] = fb["top_coordinates"]
        df.at[ridx, "protected_xh_top1_coord"] = fb["top1_coord"]
        df.at[ridx, "protected_xh_top1_percent"] = fb["top1_percent"]
        df.at[ridx, "protected_xh_center_atom"] = fb["center_atom"]
        df.at[ridx, "protected_xh_center_element"] = fb["center_element"]
        df.at[ridx, "protected_xh_pair_coords"] = fb["pair_coords"]
        df.at[ridx, "protected_xh_pair_percent_1"] = fb["pair_percent_1"]
        df.at[ridx, "protected_xh_pair_percent_2"] = fb["pair_percent_2"]
        df.at[ridx, "protected_xh_balance_ratio"] = fb["balance_ratio"]

        df.at[ridx, "functional_group_assignment"] = fb["assignment"]
        df.at[ridx, "assignment_confidence"] = "high" if float(fb["score"]) >= 75.0 else ("medium" if float(fb["score"]) >= 50.0 else "low")
        df.at[ridx, "assignment_confidence_score"] = fb["score"]

        new_warnings = [
            w.strip()
            for w in warnings.split(";")
            if w.strip()
            and w.strip() != "high_frequency_mode_without_clear_stretch_coordinate"
            and w.strip() != "topology_direct_XH_fallback_failed"
        ]
        if "topology_direct_cartesian_XH_stretch_assignment_used" not in new_warnings:
            new_warnings.append("topology_direct_cartesian_XH_stretch_assignment_used")
        df.at[ridx, "warnings"] = "; ".join(new_warnings)

    df["assignment_method"] = "Stage 3D v4.6 weighted independent-coordinate audit with single-column ORCA block parsing fix and topology-direct Cartesian X-H fallback"

    # Stage 3D v5.0 diagnostic-output cleanup:
    # topology-direct protected X-H fallback diagnostics are only meaningful
    # for rows where the fallback was actually used. Leave assignment evidence
    # columns untouched, but blank the fallback power/norm fields otherwise.
    if "protected_xh_used" in df.columns:
        unused_protected_mask = ~df["protected_xh_used"].fillna(False).astype(bool)
        for _col in ("protected_xh_normal_mode_norm", "protected_xh_total_power"):
            if _col in df.columns:
                df.loc[unused_protected_mask, _col] = pd.NA

    return df


# Final assignment-layer rebind: use the extracted module implementation,
# which preserves the same v4.3/v4.4/v4.6 wrapper chain out of the monolith.
build_stage3d_assignment_audit = mode_assignment_build_stage3d_assignment_audit



if __name__ == "__main__":
    if is_google_colab():
        print("Google Colab detected -> starting upload mode")
        colab_upload_and_run()
    else:
        print("Terminal mode detected -> using CLI")
        cli_main()
