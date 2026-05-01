from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import re

from chemistry import build_connectivity
from orcaveda_models import HessData, InternalCoordinate


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

def _compact_coord_label(name: str, max_len: Optional[int] = None) -> str:
    """
    Return a readable coordinate label for report tables.

    Stage 3D v4.3 rule:
        no truncation in audit fields. Long labels are audit evidence and must
        remain reproducible in CSV/XLSX outputs.
    """
    return str(name).replace("FG_", "").replace("Hbond_", "Hbond:")


def _normalize_element_symbol(text: str) -> str:
    token = str(text or "").strip()
    if not token:
        return token
    return token[0].upper() + token[1:].lower()


def _primitive_bond_label(name: str) -> Optional[str]:
    match = re.search(r"r\(([A-Za-z]{1,2})\d+-([A-Za-z]{1,2})\d+\)", str(name), flags=re.IGNORECASE)
    if not match:
        return None
    left = _normalize_element_symbol(match.group(1))
    right = _normalize_element_symbol(match.group(2))
    return f"{left}-{right} stretch"


def _primitive_angle_label(name: str) -> Optional[str]:
    match = re.search(
        r"ang\(([A-Za-z]{1,2})\d+-([A-Za-z]{1,2})\d+-([A-Za-z]{1,2})\d+\)",
        str(name),
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    left = _normalize_element_symbol(match.group(1))
    center = _normalize_element_symbol(match.group(2))
    right = _normalize_element_symbol(match.group(3))
    return f"{left}-{center}-{right} bend"


def _primitive_angle_tokens(name: str) -> Optional[Tuple[str, str, str]]:
    match = re.search(
        r"ang\(([A-Za-z]{1,2})\d+-([A-Za-z]{1,2})\d+-([A-Za-z]{1,2})\d+\)",
        str(name),
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    return (
        _normalize_element_symbol(match.group(1)),
        _normalize_element_symbol(match.group(2)),
        _normalize_element_symbol(match.group(3)),
    )


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

    if "acid_dimer" in kind or "acid_dimer" in name:
        return "acid dimer H-bond / intermolecular"
    if "carboxylate_hbond" in kind or "carboxylate_hbond" in name:
        return "carboxylate H-bond / intermolecular"
    if "carboxylic_acid_hbond" in kind or "carboxylic_acid_hbond" in name:
        return "carboxylic acid H-bond / intermolecular"
    if "intermolecular_oh_o_hbond" in kind or "intermolecular_oh_o_hbond" in name:
        return "intermolecular O-H···O H-bond"
    if "intermolecular_nh_o_hbond" in kind or "intermolecular_nh_o_hbond" in name:
        return "intermolecular N-H···O H-bond"
    if "intermolecular_oh_n_hbond" in kind or "intermolecular_oh_n_hbond" in name:
        return "intermolecular O-H···N H-bond"
    if "phenol_oh_stretch" in name or "phenolic_oh_stretch" in kind or "phenolic_oh_stretch" in name:
        return "phenolic O-H stretch"
    if "phenol_co_stretch" in name or "phenolic_co_stretch" in kind or "phenolic_co_stretch" in name:
        return "phenolic C-O stretch"
    if "phenol_coh_bend" in name or "phenolic_coh_bend" in kind or "phenolic_coh_bend" in name:
        return "phenolic O-H bend"
    if "carboxylic_oh_stretch" in name or "carboxylic_oh_stretch" in kind:
        return "carboxylic O-H stretch"
    if "carboxylate_co_stretch" in name or "carboxylate_co_stretch" in kind:
        return "carboxylate C-O stretch"
    if "carboxylate_oco_bend" in name or "carboxylate_oco_bend" in kind:
        return "carboxylate O-C-O bend"
    if "carboxylic_co_single_stretch" in name or "carboxylic_co_stretch" in kind:
        return "carboxylic C-O stretch"
    if "carboxylic_coh_bend" in name or "carboxylic_coh_bend" in kind:
        return "carboxylic O-H bend"
    if "carboxylic_carbonyl_stretch" in name or "carboxylic_carbonyl_stretch" in kind:
        return "carboxylic C=O stretch"
    if "secondary_aryl_amine_nh_stretch" in name or "secondary_aryl_amine_nh_stretch" in kind:
        return "secondary aryl amine N-H stretch"
    if "secondary_aryl_amine_cn_stretch" in name or "secondary_aryl_amine_cn_stretch" in kind:
        return "aryl C-N stretch"
    if "aniline_nh_stretch" in name or "aryl_amine_nh_stretch" in kind or "aryl_amine_nh_stretch" in name:
        return "aryl amine N-H stretch"
    if "aniline_cn_stretch" in name or "aryl_amine_cn_stretch" in kind or "aryl_amine_cn_stretch" in name:
        return "aryl C-N stretch"
    if "aryl_amide_nh_stretch" in name or "aryl_amide_nh_stretch" in kind:
        return "aryl amide N-H stretch"
    if "aryl_amide_cn_stretch" in name or "aryl_amide_cn_stretch" in kind:
        return "aryl amide C-N stretch"
    if "aryl_amide_nc_stretch" in name or "aryl_amide_nc_stretch" in kind:
        return "aryl amide C-N stretch"
    if "aryl_amide_co_stretch" in name or "aryl_amide_co_stretch" in kind:
        return "aryl amide C=O stretch"
    if "aryl_ketone_co_stretch" in name or "aryl_ketone_co_stretch" in kind:
        return "aryl-conjugated C=O stretch"
    if "heteroaromatic_cn_stretch" in name or "heteroaromatic_cn_stretch" in kind:
        return "heteroaromatic C-N stretch"
    if "aryl_ether_co_stretch" in name or "aryl_ether_co_stretch" in kind:
        return "aryl ether C-O stretch"
    if "aryl_ether_oc_stretch" in name or "aryl_ether_oc_stretch" in kind:
        return "aryl ether O-C stretch"
    if "peroxide_oo_stretch" in name or "peroxide_oo_stretch" in kind:
        return "O-O stretch"
    if "peroxide_ooh_bend" in name or "peroxide_ooh_bend" in kind:
        return "H-O-O bend"
    if "aromatic_ch_stretch" in kind or "aromatic_ch_stretch" in name:
        return "aromatic C-H stretch"
    if "aromatic_ch_bend" in kind or "aromatic_ch_bend" in name:
        return "aromatic C-H bend"
    if "aromatic_ring_stretch" in kind or "aromatic_ring_stretch" in name:
        return "aromatic ring stretch"
    if "aromatic_ring_deformation" in kind or "aromatic_ring_deformation" in name:
        return "aromatic ring deformation"
    if "lactam_ocn_bend" in kind or "lactam_ocn_bend" in name:
        return "amide-adjacent C-N bend"
    if "lactam_ring_deformation" in kind or "lactam_ring_deformation" in name:
        return "lactam ring deformation"

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
        primitive = _primitive_bond_label(ic.name)
        if primitive:
            return primitive
        return "bond stretch"

    if "bend" in kind or "angle" in kind:
        if "hbond" in name:
            return "H-bond angle bend"
        primitive_tokens = _primitive_angle_tokens(ic.name)
        if primitive_tokens == ("O", "C", "N") or primitive_tokens == ("N", "C", "O"):
            return "amide-adjacent C-N bend"
        primitive = _primitive_angle_label(ic.name)
        if primitive:
            return primitive
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


def _stage3d_contextual_xh_label(
    elem: str,
    coord_names: Sequence[str],
    *,
    same_sign: Optional[bool] = None,
    c_h_count: Optional[int] = None,
) -> str:
    names = [str(name).lower() for name in coord_names]
    if elem == "O":
        if any("carboxylic_oh_stretch" in name for name in names):
            return "carboxylic O-H stretch"
        if any("phenol_oh_stretch" in name or "phenolic_oh_stretch" in name for name in names):
            return "phenolic O-H stretch"
        return "O-H stretch"
    if elem == "N":
        if any("aryl_amide_nh_stretch" in name for name in names):
            return "aryl amide N-H stretch"
        if any("secondary_aryl_amine_nh_stretch" in name for name in names):
            return "secondary aryl amine N-H stretch"
        if any("aniline_nh_stretch" in name or "aryl_amine_nh_stretch" in name for name in names):
            if same_sign is None:
                return "aryl amine N-H stretch"
            return "aryl amine NH2 symmetric stretch" if same_sign else "aryl amine NH2 asymmetric stretch"
        if same_sign is None:
            return "N-H stretch"
        return "NH2 symmetric stretch" if same_sign else "NH2 asymmetric stretch"
    if elem == "C":
        if any("aromatic_ch_stretch" in name for name in names):
            return "aromatic C-H stretch"
        if same_sign is None:
            return "C-H stretch"
        if c_h_count == 2:
            return "CH2 symmetric stretch" if same_sign else "CH2 asymmetric stretch"
        return "C-H stretch"
    return "X-H stretch"


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
        label = _stage3d_contextual_xh_label(elem, [vals2[0][4], vals2[1][4]], same_sign=same_sign)
    elif elem == "C":
        c_h_count = _stage3d_xh_center_h_count(internals, atoms, heavy, elem)
        label = _stage3d_contextual_xh_label(elem, [vals2[0][4], vals2[1][4]], same_sign=same_sign, c_h_count=c_h_count)
    elif elem == "O":
        label = _stage3d_contextual_xh_label(elem, [vals2[0][4], vals2[1][4]], same_sign=same_sign)

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
        coord_names = [str(v["coord"]) for v in vals]
        if elem == "O":
            assignment = _stage3d_contextual_xh_label(elem, coord_names)
        elif elem == "N" and len(vals) >= 2:
            same_sign = (
                float(vals[0]["projection"]) == 0.0
                or float(vals[1]["projection"]) == 0.0
                or float(vals[0]["projection"]) * float(vals[1]["projection"]) > 0.0
            )
            assignment = _stage3d_contextual_xh_label(elem, coord_names, same_sign=same_sign)
        elif elem == "N":
            assignment = _stage3d_contextual_xh_label(elem, coord_names)
        elif elem == "C" and len(vals) >= 2:
            same_sign = (
                float(vals[0]["projection"]) == 0.0
                or float(vals[1]["projection"]) == 0.0
                or float(vals[0]["projection"]) * float(vals[1]["projection"]) > 0.0
            )
            c_h_count = _stage3d_xh_center_h_count(internals, atoms, heavy, elem)
            assignment = _stage3d_contextual_xh_label(elem, coord_names, same_sign=same_sign, c_h_count=c_h_count)
        elif elem == "C":
            assignment = _stage3d_contextual_xh_label(elem, coord_names)

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
        xh_families = {
            "O-H stretch",
            "N-H stretch",
            "C-H stretch",
            "carboxylic O-H stretch",
            "phenolic O-H stretch",
            "aryl amide N-H stretch",
            "aryl amine N-H stretch",
            "secondary aryl amine N-H stretch",
            "aromatic C-H stretch",
        }
        xh = [(fam, val) for fam, val in ordered if fam in xh_families]
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

    assignment = _stage3d_contextual_xh_label(str(top1["element"]), [str(top1["coord"])])

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
                assignment = _stage3d_contextual_xh_label(
                    "C",
                    [str(top[0]["coord"]), str(top[1]["coord"])],
                    same_sign=(product > 0),
                    c_h_count=c_h_count,
                )
            elif str(top[0]["element"]) == "N":
                assignment = _stage3d_contextual_xh_label(
                    "N",
                    [str(top[0]["coord"]), str(top[1]["coord"])],
                    same_sign=(product > 0),
                )
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
    assignment = _stage3d_contextual_xh_label(str(top1["center_element"]), [str(top1["coord"])])

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
        if str(top1["center_element"]) == "C":
            c_h_count = _stage3d_xh_center_h_count(internals=[], atoms=hess.atoms, heavy=int(top1["center_atom"]) - 1, elem="C")
            assignment = _stage3d_contextual_xh_label("C", [str(same_center[0]["coord"]), str(same_center[1]["coord"])], same_sign=(symmetry == "symmetric"), c_h_count=c_h_count)
        else:
            assignment = _stage3d_contextual_xh_label("N", [str(same_center[0]["coord"]), str(same_center[1]["coord"])], same_sign=(symmetry == "symmetric"))

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

