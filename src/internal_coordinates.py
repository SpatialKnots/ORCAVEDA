from __future__ import annotations

import math
from typing import Callable, List, Optional, Sequence

import numpy as np

from chemistry import adjacency
from orcaveda_models import FunctionalGroup, InternalCoordinate


def angle_deg_from_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return float("nan")
    cosine = float(np.dot(v1, v2) / denom)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def distance_fn(i: int, j: int) -> Callable[[np.ndarray], float]:
    return lambda xyz: float(np.linalg.norm(xyz[i] - xyz[j]))


def angle_fn(i: int, j: int, k: int) -> Callable[[np.ndarray], float]:
    return lambda xyz: angle_deg_from_vectors(xyz[i] - xyz[j], xyz[k] - xyz[j])


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


def build_internal_coordinates(
    atoms,
    coords_A,
    bonds,
    fragments,
    hbonds,
    groups: Optional[List[FunctionalGroup]] = None,
) -> List[InternalCoordinate]:
    adj = adjacency(len(atoms), bonds)
    coords: List[InternalCoordinate] = []

    for i, j, _ in bonds:
        label = f"r({atoms[i]}{i+1}-{atoms[j]}{j+1})"
        priority = 10 if "H" not in (atoms[i], atoms[j]) else 25
        coords.append(InternalCoordinate(label, "stretch", (i, j), priority, distance_fn(i, j), "primitive"))

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
                priority = 35 if atoms[j] != "H" else 65
                label = f"ang({atoms[i]}{i+1}-{atoms[j]}{j+1}-{atoms[k]}{k+1})"
                coords.append(InternalCoordinate(label, "bend", key, priority, angle_fn(i, j, k), "primitive"))

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

    if groups:
        for group in groups:
            if group.group == "alcohol":
                o, h, c = group.atoms0[:3]
                coords.append(InternalCoordinate(f"FG_alcohol_OH_stretch({atoms[o]}{o+1}-{atoms[h]}{h+1})", "fg_OH_stretch", (o, h), 5, distance_fn(o, h), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_alcohol_CO_stretch({atoms[c]}{c+1}-{atoms[o]}{o+1})", "fg_CO_stretch", (c, o), 6, distance_fn(c, o), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_alcohol_COH_bend({atoms[c]}{c+1}-{atoms[o]}{o+1}-{atoms[h]}{h+1})", "fg_COH_bend", (c, o, h), 20, angle_fn(c, o, h), "functional_group_template"))
            elif group.group == "phenol":
                o, h, c = group.atoms0[:3]
                coords.append(InternalCoordinate(f"FG_phenol_OH_stretch({atoms[o]}{o+1}-{atoms[h]}{h+1})", "fg_phenolic_OH_stretch", (o, h), 4, distance_fn(o, h), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_phenol_CO_stretch({atoms[c]}{c+1}-{atoms[o]}{o+1})", "fg_phenolic_CO_stretch", (c, o), 6, distance_fn(c, o), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_phenol_COH_bend({atoms[c]}{c+1}-{atoms[o]}{o+1}-{atoms[h]}{h+1})", "fg_phenolic_COH_bend", (c, o, h), 18, angle_fn(c, o, h), "functional_group_template"))
            elif group.group == "aniline":
                n = group.atoms0[0]
                hs = [x for x in group.atoms0[1:] if atoms[x] == "H"]
                c_atoms = [x for x in group.atoms0[1:] if atoms[x] == "C"]
                for h in hs:
                    coords.append(InternalCoordinate(f"FG_aniline_NH_stretch({atoms[n]}{n+1}-{atoms[h]}{h+1})", "fg_aryl_amine_NH_stretch", (n, h), 5, distance_fn(n, h), "functional_group_template"))
                if c_atoms:
                    coords.append(InternalCoordinate(f"FG_aniline_CN_stretch({atoms[c_atoms[0]]}{c_atoms[0]+1}-{atoms[n]}{n+1})", "fg_aryl_amine_CN_stretch", (c_atoms[0], n), 10, distance_fn(c_atoms[0], n), "functional_group_template"))
            elif group.group == "secondary_aryl_amine":
                n = group.atoms0[0]
                hs = [x for x in group.atoms0[1:] if atoms[x] == "H"]
                c_atoms = [x for x in group.atoms0[1:] if atoms[x] == "C"]
                for h in hs:
                    coords.append(InternalCoordinate(f"FG_secondary_aryl_amine_NH_stretch({atoms[n]}{n+1}-{atoms[h]}{h+1})", "fg_secondary_aryl_amine_NH_stretch", (n, h), 5, distance_fn(n, h), "functional_group_template"))
                if c_atoms:
                    coords.append(InternalCoordinate(f"FG_secondary_aryl_amine_CN_stretch({atoms[c_atoms[0]]}{c_atoms[0]+1}-{atoms[n]}{n+1})", "fg_secondary_aryl_amine_CN_stretch", (c_atoms[0], n), 10, distance_fn(c_atoms[0], n), "functional_group_template"))
            elif group.group == "aryl_amide":
                n, c = group.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_aryl_amide_CN_stretch({atoms[c]}{c+1}-{atoms[n]}{n+1})", "fg_aryl_amide_CN_stretch", (c, n), 9, distance_fn(c, n), "functional_group_template"))
            elif group.group == "aryl_ether":
                c_ar, o, c_other = group.atoms0[:3]
                coords.append(InternalCoordinate(f"FG_aryl_ether_CO_stretch({atoms[c_ar]}{c_ar+1}-{atoms[o]}{o+1})", "fg_aryl_ether_CO_stretch", (c_ar, o), 7, distance_fn(c_ar, o), "functional_group_template"))
                coords.append(InternalCoordinate(f"FG_aryl_ether_OC_stretch({atoms[o]}{o+1}-{atoms[c_other]}{c_other+1})", "fg_aryl_ether_OC_stretch", (o, c_other), 8, distance_fn(o, c_other), "functional_group_template"))
            elif group.group == "heteroaromatic_N":
                n, c = group.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_heteroaromatic_CN_stretch({atoms[n]}{n+1}-{atoms[c]}{c+1})", "fg_heteroaromatic_CN_stretch", (n, c), 9, distance_fn(n, c), "functional_group_template"))
            elif group.group in ("carbonyl_C=O", "ketone", "amide", "lactam_amide"):
                c, o = group.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_carbonyl_CO_stretch({atoms[c]}{c+1}={atoms[o]}{o+1})", "fg_carbonyl_stretch", (c, o), 4, distance_fn(c, o), "functional_group_template"))
            elif group.group == "nitrile_C≡N":
                c, n = group.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_nitrile_CN_stretch({atoms[c]}{c+1}≡{atoms[n]}{n+1})", "fg_nitrile_stretch", (c, n), 4, distance_fn(c, n), "functional_group_template"))
            elif group.group in ("sulfoxide", "sulfoxide_S=O"):
                s, o = group.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_sulfoxide_SO_stretch({atoms[s]}{s+1}={atoms[o]}{o+1})", "fg_sulfoxide_stretch", (s, o), 4, distance_fn(s, o), "functional_group_template"))
            elif group.group in ("methyl", "methylene", "methine"):
                c = group.atoms0[0]
                hs = [x for x in group.atoms0[1:] if atoms[x] == "H"]
                for h in hs:
                    coords.append(InternalCoordinate(f"FG_{group.group}_CH_stretch({atoms[c]}{c+1}-{atoms[h]}{h+1})", "fg_CH_stretch", (c, h), 15, distance_fn(c, h), "functional_group_template"))
            elif group.group == "aromatic_CH":
                c, h = group.atoms0[:2]
                coords.append(InternalCoordinate(f"FG_aromatic_CH_stretch({atoms[c]}{c+1}-{atoms[h]}{h+1})", "fg_aromatic_CH_stretch", (c, h), 8, distance_fn(c, h), "functional_group_template"))
                for nbr in sorted(x for x in adj[c] if x != h and atoms[x] != "H"):
                    coords.append(InternalCoordinate(f"FG_aromatic_CH_bend({atoms[nbr]}{nbr+1}-{atoms[c]}{c+1}-{atoms[h]}{h+1})", "fg_aromatic_CH_bend", (nbr, c, h), 22, angle_fn(nbr, c, h), "functional_group_template"))
            elif group.group == "aromatic_ring":
                ring_atoms = list(group.atoms0)
                ring_size = len(ring_atoms)
                if ring_size >= 3:
                    for pos, a in enumerate(ring_atoms):
                        b = ring_atoms[(pos + 1) % ring_size]
                        if b in adj[a]:
                            coords.append(InternalCoordinate(f"FG_aromatic_ring_CC_stretch({atoms[a]}{a+1}-{atoms[b]}{b+1})", "fg_aromatic_ring_stretch", (a, b), 12, distance_fn(a, b), "functional_group_template"))

    for hbond in hbonds:
        d, h, a = hbond["D0"], hbond["H0"], hbond["A0"]
        coords.append(InternalCoordinate(f"Hbond_rHA({atoms[h]}{h+1}···{atoms[a]}{a+1})", "hbond_HA", (h, a), 5, distance_fn(h, a), "cluster_template"))
        coords.append(InternalCoordinate(f"Hbond_rDA({atoms[d]}{d+1}···{atoms[a]}{a+1})", "hbond_DA", (d, a), 6, distance_fn(d, a), "cluster_template"))
        coords.append(InternalCoordinate(f"Hbond_ang({atoms[d]}{d+1}-{atoms[h]}{h+1}···{atoms[a]}{a+1})", "hbond_angle", (d, h, a), 7, angle_fn(d, h, a), "cluster_template"))

    if len(fragments) >= 2:
        frag_id = {atom_index: frag_index for frag_index, frag in enumerate(fragments) for atom_index in frag}
        heavy = [i for i, atom in enumerate(atoms) if atom != "H"]
        for pos, i in enumerate(heavy):
            for j in heavy[pos + 1:]:
                if frag_id.get(i) != frag_id.get(j):
                    label = f"interfrag_R({atoms[i]}{i+1}···{atoms[j]}{j+1})"
                    coords.append(InternalCoordinate(label, "interfragment_distance", (i, j), 40, distance_fn(i, j), "cluster_template"))

    return sorted(coords, key=lambda coord: (coord.priority, coord.name))
