from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Protocol, Sequence, Tuple

import numpy as np
import pandas as pd

from orcaveda_models import (
    AtomEnvironmentAnnotation,
    ChemicalSystemAnnotation,
    FunctionalGroup,
)


COVALENT_RADII_A = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05,
    "F": 0.57, "Cl": 1.02, "Br": 1.20, "I": 1.39, "P": 1.07,
}
SUPPORTED_ELEMENTS = set(COVALENT_RADII_A)


def build_connectivity(atoms: Sequence[str], coords_A: np.ndarray, scale: float = 1.25, extra_A: float = 0.15):
    bonds = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if atoms[i] == atoms[j] == "H":
                continue
            ri = COVALENT_RADII_A.get(atoms[i], 0.75)
            rj = COVALENT_RADII_A.get(atoms[j], 0.75)
            distance = float(np.linalg.norm(coords_A[i] - coords_A[j]))
            cutoff = scale * (ri + rj) + extra_A
            if distance <= cutoff:
                bonds.append((i, j, distance))
    return bonds


def adjacency(natoms: int, bonds) -> Dict[int, set]:
    adj = {i: set() for i in range(natoms)}
    for i, j, _ in bonds:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def bond_distance(bonds, i, j) -> Optional[float]:
    a, b = sorted((i, j))
    for x, y, distance in bonds:
        if x == a and y == b:
            return distance
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
            atom = stack.pop()
            fragment.append(atom)
            for neighbor in adj[atom]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        fragments.append(sorted(fragment))
    return fragments


def formula_string(atoms: Sequence[str]) -> str:
    counts = {atom: atoms.count(atom) for atom in sorted(set(atoms))}
    order = ["C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    parts = []
    for element in order:
        if counts.get(element, 0):
            count = counts.pop(element)
            parts.append(f"{element}{'' if count == 1 else count}")
    for element, count in sorted(counts.items()):
        parts.append(f"{element}{'' if count == 1 else count}")
    return "".join(parts)


def classify_system(fragments: Sequence[Sequence[int]]) -> str:
    n = len(fragments)
    sizes = sorted(len(fragment) for fragment in fragments)
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
    cosine = float(np.dot(v1, v2) / denom)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


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
    for donor, hydrogen in donors:
        for acceptor in acceptors:
            if acceptor == donor or frag_id.get(donor) == frag_id.get(acceptor) or tuple(sorted((hydrogen, acceptor))) in bonded:
                continue
            r_HA = float(np.linalg.norm(coords_A[hydrogen] - coords_A[acceptor]))
            r_DA = float(np.linalg.norm(coords_A[donor] - coords_A[acceptor]))
            angle = _angle_deg_from_vectors(coords_A[donor] - coords_A[hydrogen], coords_A[acceptor] - coords_A[hydrogen])
            if r_HA <= h_a_max_A and r_DA <= d_a_max_A and angle >= angle_min_deg:
                hbonds.append(
                    {
                        "D0": donor,
                        "H0": hydrogen,
                        "A0": acceptor,
                        "D": donor + 1,
                        "H": hydrogen + 1,
                        "A": acceptor + 1,
                        "type": f"{atoms[donor]}-HВ·В·В·{atoms[acceptor]}",
                        "rHA_A": r_HA,
                        "rDA_A": r_DA,
                        "angle_deg": angle,
                    }
                )
    return sorted(hbonds, key=lambda row: row["rHA_A"])


def atom_environment_table(atoms, coords_A, bonds) -> pd.DataFrame:
    adj = adjacency(len(atoms), bonds)
    rows: List[AtomEnvironmentAnnotation] = []
    for i, element in enumerate(atoms):
        neighbors = sorted(adj[i])
        neighbor_elements = [atoms[j] for j in neighbors]
        heavy_neighbors = [j for j in neighbors if atoms[j] != "H"]
        h_count = sum(1 for j in neighbors if atoms[j] == "H")
        rows.append(
            AtomEnvironmentAnnotation(
                atom=i + 1,
                element=element,
                degree=len(neighbors),
                h_neighbors=h_count,
                heavy_neighbors=len(heavy_neighbors),
                neighbor_elements=",".join(neighbor_elements),
                environment_label=f"{element};deg{len(neighbors)};H{h_count};heavy{len(heavy_neighbors)};nbrs={','.join(neighbor_elements)}",
            )
        )
    return pd.DataFrame(
        {
            "atom": row.atom,
            "element": row.element,
            "degree": row.degree,
            "H_neighbors": row.h_neighbors,
            "heavy_neighbors": row.heavy_neighbors,
            "neighbor_elements": row.neighbor_elements,
            "environment_label": row.environment_label,
        }
        for row in rows
    )


def detect_rings(natoms: int, bonds, max_size: int = 8) -> List[Tuple[int, ...]]:
    adj = adjacency(natoms, bonds)
    rings = set()

    def dfs(start, current, path):
        if len(path) > max_size:
            return
        for neighbor in adj[current]:
            if neighbor == start and len(path) >= 3:
                rings.add(tuple(sorted(path)))
            elif neighbor > start and neighbor not in path:
                dfs(start, neighbor, path + [neighbor])

    for start in range(natoms):
        dfs(start, start, [start])

    unique = sorted(rings, key=lambda ring: (len(ring), ring))
    minimal = []
    for ring in unique:
        set_ring = set(ring)
        if not any(set(existing).issubset(set_ring) for existing in minimal):
            minimal.append(ring)
    return minimal


def detect_functional_groups(atoms, coords_A, bonds) -> List[FunctionalGroup]:
    adj = adjacency(len(atoms), bonds)
    groups: List[FunctionalGroup] = []
    rings = detect_rings(len(atoms), bonds)

    def add(group, atoms0, desc, conf, evidence):
        key = (group, tuple(sorted(atoms0)))
        if key not in {(existing.group, tuple(sorted(existing.atoms0))) for existing in groups}:
            groups.append(FunctionalGroup(group, tuple(atoms0), desc, conf, evidence))

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
            distance = bond_distance(bonds, o, s_ns[0])
            if distance is not None and distance < 1.62:
                add("sulfoxide_S=O", (s_ns[0], o), "short S-O bond assigned as sulfoxide S=O candidate", "medium", f"S-O distance={distance:.3f} A")

    carbonyl_C = []
    for c, el in enumerate(atoms):
        if el != "C":
            continue
        o_short = []
        for n in adj[c]:
            if atoms[n] == "O":
                distance = bond_distance(bonds, c, n)
                if distance is not None and distance < 1.35:
                    o_short.append(n)
        for o in o_short:
            carbonyl_C.append(c)
            n_ns = [n for n in adj[c] if atoms[n] == "N"]
            c_ns = [n for n in adj[c] if atoms[n] == "C"]
            h_ns = [n for n in adj[c] if atoms[n] == "H"]
            o_single = [n for n in adj[c] if atoms[n] == "O" and n != o]
            add("carbonyl_C=O", (c, o), "short C-O bond assigned as carbonyl candidate", "high", f"C-O distance={bond_distance(bonds, c, o):.3f} A")
            if n_ns:
                ring_flag = any(c in ring and any(n in ring for n in n_ns) for ring in rings)
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

    for c, el in enumerate(atoms):
        if el != "C":
            continue
        for n in adj[c]:
            if atoms[n] == "N":
                distance = bond_distance(bonds, c, n)
                if distance is not None and distance < 1.22:
                    add("nitrile_C≡N", (c, n), "short terminal C-N bond assigned as nitrile", "high", f"C-N distance={distance:.3f} A")

    for n, el in enumerate(atoms):
        if el != "N":
            continue
        ns = sorted(adj[n])
        h_ns = [x for x in ns if atoms[x] == "H"]
        c_ns = [x for x in ns if atoms[x] == "C"]
        if any(c in carbonyl_C for c in c_ns):
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

    for s, el in enumerate(atoms):
        if el != "S":
            continue
        o_short = [n for n in adj[s] if atoms[n] == "O" and (bond_distance(bonds, s, n) or 9) < 1.62]
        c_ns = [n for n in adj[s] if atoms[n] == "C"]
        if o_short and len(c_ns) >= 2:
            add("sulfoxide", (s, o_short[0], *c_ns[:2]), "S=O with two carbon substituents", "high", "short S-O and two S-C bonds")

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

    for ring in rings:
        add("ring", ring, f"{len(ring)}-membered ring candidate", "medium", "simple cycle detected in covalent graph")

    return groups


def annotate_chemical_system(atoms, coords_A) -> ChemicalSystemAnnotation:
    bonds = tuple(build_connectivity(atoms, coords_A))
    fragments = tuple(tuple(fragment) for fragment in split_fragments(len(atoms), bonds))
    hbonds = tuple(detect_interfragment_hbonds(atoms, coords_A, bonds, fragments))
    groups = tuple(detect_functional_groups(atoms, coords_A, bonds))
    return ChemicalSystemAnnotation(
        formula=formula_string(atoms),
        system_type=classify_system(fragments),
        bonds=bonds,
        fragments=fragments,
        functional_groups=groups,
        interfragment_hbonds=hbonds,
    )


_legacy_build_connectivity = build_connectivity
_legacy_detect_interfragment_hbonds = detect_interfragment_hbonds
_legacy_atom_environment_table = atom_environment_table
_legacy_detect_functional_groups = detect_functional_groups
_legacy_annotate_chemical_system = annotate_chemical_system


class ChemistryBackend(Protocol):
    name: str
    supported_elements: set[str]

    def build_connectivity(self, atoms: Sequence[str], coords_A: np.ndarray, scale: float = 1.25, extra_A: float = 0.15):
        ...

    def detect_interfragment_hbonds(
        self,
        atoms,
        coords_A,
        bonds,
        fragments,
        h_a_max_A: float = 2.70,
        d_a_max_A: float = 3.50,
        angle_min_deg: float = 120.0,
    ):
        ...

    def atom_environment_table(self, atoms, coords_A, bonds) -> pd.DataFrame:
        ...

    def detect_functional_groups(self, atoms, coords_A, bonds) -> List[FunctionalGroup]:
        ...

    def annotate_chemical_system(self, atoms, coords_A) -> ChemicalSystemAnnotation:
        ...


@dataclass(frozen=True)
class LegacyChemistryBackend:
    name: str = "legacy"
    supported_elements: set[str] = field(default_factory=lambda: set(SUPPORTED_ELEMENTS))

    def build_connectivity(self, atoms: Sequence[str], coords_A: np.ndarray, scale: float = 1.25, extra_A: float = 0.15):
        return _legacy_build_connectivity(atoms, coords_A, scale=scale, extra_A=extra_A)

    def detect_interfragment_hbonds(
        self,
        atoms,
        coords_A,
        bonds,
        fragments,
        h_a_max_A: float = 2.70,
        d_a_max_A: float = 3.50,
        angle_min_deg: float = 120.0,
    ):
        return _legacy_detect_interfragment_hbonds(
            atoms,
            coords_A,
            bonds,
            fragments,
            h_a_max_A=h_a_max_A,
            d_a_max_A=d_a_max_A,
            angle_min_deg=angle_min_deg,
        )

    def atom_environment_table(self, atoms, coords_A, bonds) -> pd.DataFrame:
        return _legacy_atom_environment_table(atoms, coords_A, bonds)

    def detect_functional_groups(self, atoms, coords_A, bonds) -> List[FunctionalGroup]:
        return _legacy_detect_functional_groups(atoms, coords_A, bonds)

    def annotate_chemical_system(self, atoms, coords_A) -> ChemicalSystemAnnotation:
        return _legacy_annotate_chemical_system(atoms, coords_A)


_BACKENDS: Dict[str, ChemistryBackend] = {
    "legacy": LegacyChemistryBackend(),
}
_ACTIVE_BACKEND_NAME = "legacy"


def register_backend(backend: ChemistryBackend) -> None:
    _BACKENDS[str(backend.name)] = backend


def get_backend(name: str) -> ChemistryBackend:
    if name not in _BACKENDS:
        raise KeyError(f"Unknown chemistry backend: {name}")
    return _BACKENDS[name]


def get_active_backend() -> ChemistryBackend:
    return get_backend(_ACTIVE_BACKEND_NAME)


def get_active_backend_name() -> str:
    return _ACTIVE_BACKEND_NAME


def set_active_backend(name: str) -> ChemistryBackend:
    global _ACTIVE_BACKEND_NAME
    backend = get_backend(name)
    _ACTIVE_BACKEND_NAME = str(backend.name)
    return backend


def get_supported_elements() -> set[str]:
    return set(get_active_backend().supported_elements)


def build_connectivity(atoms: Sequence[str], coords_A: np.ndarray, scale: float = 1.25, extra_A: float = 0.15):
    return get_active_backend().build_connectivity(atoms, coords_A, scale=scale, extra_A=extra_A)


def detect_interfragment_hbonds(atoms, coords_A, bonds, fragments, h_a_max_A=2.70, d_a_max_A=3.50, angle_min_deg=120.0):
    return get_active_backend().detect_interfragment_hbonds(
        atoms,
        coords_A,
        bonds,
        fragments,
        h_a_max_A=h_a_max_A,
        d_a_max_A=d_a_max_A,
        angle_min_deg=angle_min_deg,
    )


def atom_environment_table(atoms, coords_A, bonds) -> pd.DataFrame:
    return get_active_backend().atom_environment_table(atoms, coords_A, bonds)


def detect_functional_groups(atoms, coords_A, bonds) -> List[FunctionalGroup]:
    return get_active_backend().detect_functional_groups(atoms, coords_A, bonds)


def annotate_chemical_system(atoms, coords_A) -> ChemicalSystemAnnotation:
    return get_active_backend().annotate_chemical_system(atoms, coords_A)
