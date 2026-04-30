from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

from chemistry import SUPPORTED_ELEMENTS, formula_string
from orcaveda_models import ChemicalSystemAnnotation, FunctionalGroup


@dataclass(frozen=True)
class RDKitChemistryBackend:
    name: str = "rdkit"
    supported_elements: set[str] = field(default_factory=lambda: set(SUPPORTED_ELEMENTS))
    default_charge: int = 0
    cov_factor: float = 1.25

    def __post_init__(self):
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import rdDetermineBonds  # noqa: F401
        from rdkit.Geometry import Point3D  # noqa: F401

    def _build_mol(self, atoms: Sequence[str], coords_A: np.ndarray):
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Geometry import Point3D

        mol = Chem.RWMol()
        for atom_symbol in atoms:
            mol.AddAtom(Chem.Atom(str(atom_symbol)))

        conf = Chem.Conformer(len(atoms))
        for idx, (x, y, z) in enumerate(np.asarray(coords_A, dtype=float)):
            conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
        mol.AddConformer(conf)

        rdDetermineBonds.DetermineConnectivity(
            mol,
            useHueckel=False,
            charge=self.default_charge,
            covFactor=self.cov_factor,
            useVdw=True,
        )
        try:
            rdDetermineBonds.DetermineBondOrders(
                mol,
                charge=self.default_charge,
                allowChargedFragments=True,
                embedChiral=True,
                useAtomMap=False,
                maxIterations=0,
            )
        except Exception:
            # Connectivity alone is still useful even if bond-order perception
            # cannot be completed from coordinates only.
            pass
        return mol.GetMol()

    def _bond_length(self, mol, i: int, j: int) -> float:
        conf = mol.GetConformer()
        pi = conf.GetAtomPosition(int(i))
        pj = conf.GetAtomPosition(int(j))
        return float(np.linalg.norm(np.array([pi.x, pi.y, pi.z]) - np.array([pj.x, pj.y, pj.z])))

    def build_connectivity(
        self,
        atoms: Sequence[str],
        coords_A: np.ndarray,
        scale: float = 1.25,
        extra_A: float = 0.15,
    ):
        del scale, extra_A
        mol = self._build_mol(atoms, coords_A)
        conf = mol.GetConformer()
        bonds = []
        for bond in mol.GetBonds():
            i = int(bond.GetBeginAtomIdx())
            j = int(bond.GetEndAtomIdx())
            pi = conf.GetAtomPosition(i)
            pj = conf.GetAtomPosition(j)
            distance = float(np.linalg.norm(np.array([pi.x, pi.y, pi.z]) - np.array([pj.x, pj.y, pj.z])))
            bonds.append((min(i, j), max(i, j), distance))
        return sorted(bonds)

    def _fragments_from_mol(self, mol) -> Tuple[Tuple[int, ...], ...]:
        from rdkit import Chem

        return tuple(tuple(int(i) for i in frag) for frag in Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False))

    def _add_group(
        self,
        groups: List[FunctionalGroup],
        group: str,
        atoms0: Iterable[int],
        description: str,
        confidence: str,
        evidence: str,
    ) -> None:
        atoms_tuple = tuple(int(a) for a in atoms0)
        key = (group, tuple(sorted(atoms_tuple)))
        if key in {(existing.group, tuple(sorted(existing.atoms0))) for existing in groups}:
            return
        groups.append(
            FunctionalGroup(
                group=group,
                atoms0=atoms_tuple,
                description=description,
                confidence=confidence,
                evidence=evidence,
            )
        )

    def _functional_groups_from_mol(self, mol, atoms: Sequence[str]) -> List[FunctionalGroup]:
        del atoms
        groups: List[FunctionalGroup] = []
        ring_info = mol.GetRingInfo()

        aromatic_rings = [
            tuple(int(atom_idx) for atom_idx in ring)
            for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(int(atom_idx)).GetIsAromatic() for atom_idx in ring)
        ]
        aromatic_atom_indices = {int(atom.GetIdx()) for atom in mol.GetAtoms() if atom.GetIsAromatic()}
        aromatic_carbon_indices = {
            idx for idx in aromatic_atom_indices
            if mol.GetAtomWithIdx(idx).GetSymbol() == "C"
        }
        aromatic_nitrogen_indices = {
            idx for idx in aromatic_atom_indices
            if mol.GetAtomWithIdx(idx).GetSymbol() == "N"
        }

        for bond in mol.GetBonds():
            begin = bond.GetBeginAtom()
            end = bond.GetEndAtom()
            symbols = {begin.GetSymbol(), end.GetSymbol()}

            if symbols == {"C", "O"} and bond.GetBondTypeAsDouble() >= 1.5:
                c_atom = begin if begin.GetSymbol() == "C" else end
                o_atom = end if c_atom is begin else begin
                c_idx = int(c_atom.GetIdx())
                o_idx = int(o_atom.GetIdx())
                self._add_group(groups, "carbonyl_C=O", (c_idx, o_idx), "RDKit carbonyl bond perception", "high", "RDKit C=O bond order")

                n_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetIdx() != o_idx and nbr.GetSymbol() == "N"]
                c_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetIdx() != o_idx and nbr.GetSymbol() == "C"]
                h_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetIdx() != o_idx and nbr.GetSymbol() == "H"]
                o_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetIdx() != o_idx and nbr.GetSymbol() == "O"]

                if n_neighbors:
                    n_idx = int(n_neighbors[0].GetIdx())
                    same_ring = any(c_idx in ring and n_idx in ring for ring in ring_info.AtomRings())
                    self._add_group(groups, "lactam_amide" if same_ring else "amide", (c_idx, o_idx, n_idx), "RDKit amide context", "high", "C=O carbon bonded to N")
                elif len(c_neighbors) >= 2:
                    self._add_group(groups, "ketone", (c_idx, o_idx, int(c_neighbors[0].GetIdx()), int(c_neighbors[1].GetIdx())), "RDKit ketone context", "high", "C=O carbon bonded to two carbons")
                    aromatic_neighbors = [nbr for nbr in c_neighbors if int(nbr.GetIdx()) in aromatic_carbon_indices]
                    if aromatic_neighbors:
                        aromatic_c = int(aromatic_neighbors[0].GetIdx())
                        alkyl_c = next(int(nbr.GetIdx()) for nbr in c_neighbors if int(nbr.GetIdx()) != aromatic_c)
                        self._add_group(groups, "aryl_ketone", (c_idx, o_idx, aromatic_c, alkyl_c), "RDKit aryl ketone context", "high", "Ketone carbonyl conjugated to aromatic carbon")
                elif h_neighbors and c_neighbors:
                    self._add_group(groups, "aldehyde", (c_idx, o_idx, int(h_neighbors[0].GetIdx()), int(c_neighbors[0].GetIdx())), "RDKit aldehyde context", "high", "C=O carbon bonded to H and C")
                elif o_neighbors:
                    o2 = o_neighbors[0]
                    if any(nbr.GetSymbol() == "H" for nbr in o2.GetNeighbors()):
                        self._add_group(groups, "carboxylic_acid", (c_idx, o_idx, int(o2.GetIdx())), "RDKit carboxylic acid context", "high", "C=O and C-OH in same carbon")
                    elif any(nbr.GetSymbol() == "C" and int(nbr.GetIdx()) != c_idx for nbr in o2.GetNeighbors()):
                        self._add_group(groups, "ester", (c_idx, o_idx, int(o2.GetIdx())), "RDKit ester context", "high", "C=O and C-OR in same carbon")

            if symbols == {"C", "N"} and bond.GetBondTypeAsDouble() >= 2.5:
                c_atom = begin if begin.GetSymbol() == "C" else end
                n_atom = end if c_atom is begin else begin
                self._add_group(groups, "nitrile_C≡N", (int(c_atom.GetIdx()), int(n_atom.GetIdx())), "RDKit nitrile bond perception", "high", "RDKit C#N bond order")

            if symbols == {"S", "O"} and bond.GetBondTypeAsDouble() >= 1.5:
                s_atom = begin if begin.GetSymbol() == "S" else end
                o_atom = end if s_atom is begin else begin
                self._add_group(groups, "sulfoxide_S=O", (int(s_atom.GetIdx()), int(o_atom.GetIdx())), "RDKit S=O bond perception", "high", "RDKit S=O bond order")

            if symbols == {"O", "O"} and 0.5 <= bond.GetBondTypeAsDouble() <= 1.5:
                o1_idx = int(begin.GetIdx())
                o2_idx = int(end.GetIdx())
                self._add_group(groups, "peroxide", (o1_idx, o2_idx), "RDKit peroxide O-O context", "high", "Single O-O bond")

        for atom in mol.GetAtoms():
            if atom.GetSymbol() != "C":
                continue
            c_idx = int(atom.GetIdx())
            o_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
            if len(o_neighbors) < 2:
                continue
            o_pair = sorted(int(nbr.GetIdx()) for nbr in o_neighbors[:2])
            o_has_h = any(any(nbr.GetSymbol() == "H" for nbr in mol.GetAtomWithIdx(o_idx).GetNeighbors()) for o_idx in o_pair)
            d1 = self._bond_length(mol, c_idx, o_pair[0])
            d2 = self._bond_length(mol, c_idx, o_pair[1])
            if not o_has_h and max(d1, d2) <= 1.32 and abs(d1 - d2) <= 0.08:
                self._add_group(groups, "carboxylate", (c_idx, o_pair[0], o_pair[1]), "RDKit carboxylate resonance context", "high", f"Two short C-O bonds with |Δ|={abs(d1-d2):.3f} A")

        for ring in aromatic_rings:
            self._add_group(groups, "aromatic_ring", ring, f"{len(ring)}-membered aromatic ring from RDKit aromaticity", "high", "RDKit aromatic ring perception")

        for atom in mol.GetAtoms():
            idx = int(atom.GetIdx())
            symbol = atom.GetSymbol()
            neighbors = list(atom.GetNeighbors())
            hydrogen_neighbors = [nbr for nbr in neighbors if nbr.GetSymbol() == "H"]
            carbon_neighbors = [nbr for nbr in neighbors if nbr.GetSymbol() == "C"]
            sulfur_neighbors = [nbr for nbr in neighbors if nbr.GetSymbol() == "S"]

            if symbol == "O":
                if hydrogen_neighbors and carbon_neighbors:
                    c_idx = int(carbon_neighbors[0].GetIdx())
                    h_idx = int(hydrogen_neighbors[0].GetIdx())
                    c_atom = mol.GetAtomWithIdx(c_idx)
                    is_carboxylic_oh = any(
                        nbr.GetSymbol() == "O"
                        and int(nbr.GetIdx()) != idx
                        and mol.GetBondBetweenAtoms(c_idx, int(nbr.GetIdx())) is not None
                        and mol.GetBondBetweenAtoms(c_idx, int(nbr.GetIdx())).GetBondTypeAsDouble() >= 1.5
                        for nbr in c_atom.GetNeighbors()
                    )
                    if not is_carboxylic_oh:
                        self._add_group(groups, "alcohol", (idx, h_idx, c_idx), "RDKit alcohol context", "high", "O bonded to H and C")
                    if any(int(cnbr.GetIdx()) in aromatic_carbon_indices for cnbr in carbon_neighbors):
                        aromatic_c = next(int(cnbr.GetIdx()) for cnbr in carbon_neighbors if int(cnbr.GetIdx()) in aromatic_carbon_indices)
                        self._add_group(groups, "phenol", (idx, h_idx, aromatic_c), "RDKit phenol SMARTS context", "high", "O-H bonded to aromatic carbon")
                if len(carbon_neighbors) >= 2:
                    self._add_group(groups, "ether", (idx, int(carbon_neighbors[0].GetIdx()), int(carbon_neighbors[1].GetIdx())), "RDKit ether context", "high", "O bonded to two carbons")
                if len(carbon_neighbors) >= 2 and any(int(cnbr.GetIdx()) in aromatic_carbon_indices for cnbr in carbon_neighbors):
                    aromatic_c = next(int(cnbr.GetIdx()) for cnbr in carbon_neighbors if int(cnbr.GetIdx()) in aromatic_carbon_indices)
                    other_c = next(int(cnbr.GetIdx()) for cnbr in carbon_neighbors if int(cnbr.GetIdx()) != aromatic_c)
                    self._add_group(groups, "aryl_ether", (aromatic_c, idx, other_c), "RDKit aryl ether context", "high", "O bonded to aromatic carbon and carbon substituent")
                if sulfur_neighbors and any(group.group == "sulfoxide_S=O" and idx in group.atoms0 for group in groups):
                    carbon_substituents = [int(nbr.GetIdx()) for nbr in sulfur_neighbors[0].GetNeighbors() if nbr.GetSymbol() == "C"][:2]
                    if len(carbon_substituents) == 2:
                        self._add_group(groups, "sulfoxide", (int(sulfur_neighbors[0].GetIdx()), idx, carbon_substituents[0], carbon_substituents[1]), "RDKit sulfoxide context", "high", "S=O with carbon substituents")

            if symbol == "N":
                carbon_neighbors = [nbr for nbr in neighbors if nbr.GetSymbol() == "C"]
                aromatic_carbon_neighbors = [nbr for nbr in carbon_neighbors if int(nbr.GetIdx()) in aromatic_carbon_indices]
                is_amide_n = any(group.group in {"amide", "lactam_amide"} and idx in group.atoms0 for group in groups)
                if idx in aromatic_nitrogen_indices:
                    for nbr in carbon_neighbors[:2]:
                        self._add_group(groups, "heteroaromatic_N", (idx, int(nbr.GetIdx())), "RDKit heteroaromatic nitrogen context", "high", "Aromatic ring nitrogen")
                    continue
                if is_amide_n:
                    if aromatic_carbon_neighbors:
                        self._add_group(groups, "aryl_amide", (idx, int(aromatic_carbon_neighbors[0].GetIdx())), "RDKit aryl amide context", "high", "Amide nitrogen bonded to aromatic carbon")
                    continue
                if len(hydrogen_neighbors) == 2 and carbon_neighbors:
                    c_idx = int(carbon_neighbors[0].GetIdx())
                    h1 = int(hydrogen_neighbors[0].GetIdx())
                    h2 = int(hydrogen_neighbors[1].GetIdx())
                    self._add_group(groups, "primary_amine", (idx, h1, h2, c_idx), "RDKit primary amine context", "high", "N-H2 and N-C")
                    if aromatic_carbon_neighbors:
                        aromatic_c = int(aromatic_carbon_neighbors[0].GetIdx())
                        self._add_group(groups, "aniline", (idx, h1, h2, aromatic_c), "RDKit aniline SMARTS context", "high", "NH2 bonded to aromatic carbon")
                elif len(hydrogen_neighbors) == 1 and len(carbon_neighbors) >= 1:
                    self._add_group(groups, "secondary_amine", (idx, int(hydrogen_neighbors[0].GetIdx()), *[int(nbr.GetIdx()) for nbr in carbon_neighbors[:2]]), "RDKit secondary amine context", "high", "N-H and N-C")
                    if aromatic_carbon_neighbors:
                        aromatic_c = int(aromatic_carbon_neighbors[0].GetIdx())
                        self._add_group(groups, "secondary_aryl_amine", (idx, int(hydrogen_neighbors[0].GetIdx()), aromatic_c), "RDKit secondary aryl amine context", "high", "N-H bonded to aromatic carbon")
                elif len(carbon_neighbors) >= 3:
                    self._add_group(groups, "tertiary_amine", (idx, int(carbon_neighbors[0].GetIdx()), int(carbon_neighbors[1].GetIdx()), int(carbon_neighbors[2].GetIdx())), "RDKit tertiary amine context", "high", "N-C3")
                elif len(carbon_neighbors) >= 2 and not aromatic_carbon_neighbors:
                    self._add_group(groups, "dialkyl_amide_N_or_amine_N", (idx, int(carbon_neighbors[0].GetIdx()), int(carbon_neighbors[1].GetIdx())), "RDKit dialkyl N context", "medium", "N-C2")

            if symbol == "C":
                if len(hydrogen_neighbors) == 3:
                    self._add_group(groups, "methyl", (idx, *[int(nbr.GetIdx()) for nbr in hydrogen_neighbors]), "RDKit CH3 context", "high", "C bonded to three H")
                elif len(hydrogen_neighbors) == 2:
                    self._add_group(groups, "methylene", (idx, *[int(nbr.GetIdx()) for nbr in hydrogen_neighbors]), "RDKit CH2 context", "high", "C bonded to two H")
                elif len(hydrogen_neighbors) == 1:
                    h_idx = int(hydrogen_neighbors[0].GetIdx())
                    self._add_group(groups, "methine", (idx, h_idx), "RDKit CH context", "high", "C bonded to one H")
                    if atom.GetIsAromatic():
                        self._add_group(groups, "aromatic_CH", (idx, h_idx), "RDKit aromatic C-H context", "high", "Aromatic carbon bonded to H")

        for ring in ring_info.AtomRings():
            self._add_group(groups, "ring", ring, f"{len(ring)}-membered ring from RDKit ring perception", "high", "RDKit ring perception")

        return groups

    def _carboxylic_roles_for_atom(
        self,
        atom_idx: int,
        groups: Sequence[FunctionalGroup],
    ) -> List[Tuple[FunctionalGroup, str]]:
        roles: List[Tuple[FunctionalGroup, str]] = []
        for group in groups:
            if group.group != "carboxylic_acid":
                continue
            c_idx, o_carbonyl, o_hydroxyl = (int(x) for x in group.atoms0[:3])
            if int(atom_idx) == o_hydroxyl:
                roles.append((group, "hydroxyl_O"))
            elif int(atom_idx) == o_carbonyl:
                roles.append((group, "carbonyl_O"))
            elif int(atom_idx) == c_idx:
                roles.append((group, "carbonyl_C"))
        return roles

    def _carboxylate_roles_for_atom(
        self,
        atom_idx: int,
        groups: Sequence[FunctionalGroup],
    ) -> List[Tuple[FunctionalGroup, str]]:
        roles: List[Tuple[FunctionalGroup, str]] = []
        for group in groups:
            if group.group != "carboxylate":
                continue
            c_idx, o1, o2 = (int(x) for x in group.atoms0[:3])
            if int(atom_idx) in {o1, o2}:
                roles.append((group, "carboxylate_O"))
            elif int(atom_idx) == c_idx:
                roles.append((group, "carboxylate_C"))
        return roles

    def _annotate_hbond_contexts(
        self,
        hbonds: Sequence[dict],
        groups: Sequence[FunctionalGroup],
        atoms: Sequence[str],
        fragments: Sequence[Sequence[int]],
    ) -> List[dict]:
        frag_id = {int(atom): int(frag_idx) for frag_idx, frag in enumerate(fragments) for atom in frag}
        annotated: List[dict] = []
        for hbond in hbonds:
            row = dict(hbond)
            donor = int(row["D0"])
            acceptor = int(row["A0"])
            donor_roles = self._carboxylic_roles_for_atom(donor, groups)
            acceptor_roles = self._carboxylic_roles_for_atom(acceptor, groups)
            acceptor_carboxylate_roles = self._carboxylate_roles_for_atom(acceptor, groups)

            chem_type = "intermolecular_hbond"
            context_label = f"intermolecular {atoms[donor]}-H···{atoms[acceptor]} H-bond"
            donor_group = ""
            acceptor_group = ""
            donor_group_atoms = ()
            acceptor_group_atoms = ()

            donor_hydroxyl = [(group, role) for group, role in donor_roles if role == "hydroxyl_O"]
            acceptor_carbonyl = [(group, role) for group, role in acceptor_roles if role == "carbonyl_O"]

            if donor_hydroxyl and acceptor_carbonyl:
                chem_type = "carboxylic_acid_hbond"
                context_label = "intermolecular carboxylic O-H···O=C H-bond"
                donor_group = "carboxylic_acid"
                acceptor_group = "carboxylic_acid"
                donor_group_atoms = tuple(int(x) for x in donor_hydroxyl[0][0].atoms0)
                acceptor_group_atoms = tuple(int(x) for x in acceptor_carbonyl[0][0].atoms0)
            elif donor_hydroxyl and acceptor_carboxylate_roles:
                chem_type = "carboxylate_hbond"
                context_label = "intermolecular carboxylic O-H···O(carboxylate) H-bond"
                donor_group = "carboxylic_acid"
                acceptor_group = "carboxylate"
                donor_group_atoms = tuple(int(x) for x in donor_hydroxyl[0][0].atoms0)
                acceptor_group_atoms = tuple(int(x) for x in acceptor_carboxylate_roles[0][0].atoms0)
            elif atoms[donor] == "O" and atoms[acceptor] == "O":
                chem_type = "intermolecular_OH_O_hbond"
                context_label = "intermolecular O-H···O H-bond"
            elif atoms[donor] == "N" and atoms[acceptor] == "O":
                chem_type = "intermolecular_NH_O_hbond"
                context_label = "intermolecular N-H···O H-bond"
            elif atoms[donor] == "O" and atoms[acceptor] == "N":
                chem_type = "intermolecular_OH_N_hbond"
                context_label = "intermolecular O-H···N H-bond"

            row["chem_type"] = chem_type
            row["context_label"] = context_label
            row["donor_group"] = donor_group
            row["acceptor_group"] = acceptor_group
            row["donor_group_atoms"] = "-".join(str(int(x) + 1) for x in donor_group_atoms) if donor_group_atoms else ""
            row["acceptor_group_atoms"] = "-".join(str(int(x) + 1) for x in acceptor_group_atoms) if acceptor_group_atoms else ""
            row["fragment_pair"] = f"{frag_id.get(donor, -1)+1}-{frag_id.get(acceptor, -1)+1}"
            annotated.append(row)
        return annotated

    def _augment_supramolecular_groups(
        self,
        groups: Sequence[FunctionalGroup],
        hbonds: Sequence[dict],
    ) -> List[FunctionalGroup]:
        augmented = list(groups)
        seen = {(group.group, tuple(sorted(int(a) for a in group.atoms0))) for group in augmented}
        by_pair = {}
        for row in hbonds:
            if row.get("chem_type") != "carboxylic_acid_hbond":
                continue
            donor_atoms = tuple(sorted(int(x) for x in str(row.get("donor_group_atoms", "")).split("-") if x))
            acceptor_atoms = tuple(sorted(int(x) for x in str(row.get("acceptor_group_atoms", "")).split("-") if x))
            if not donor_atoms or not acceptor_atoms:
                continue
            key = frozenset((donor_atoms, acceptor_atoms))
            by_pair.setdefault(key, set()).add((donor_atoms, acceptor_atoms))

        for key, directions in by_pair.items():
            if len(directions) < 2:
                continue
            atom_union = sorted({atom - 1 for atoms1 in key for atom in atoms1})
            group_key = ("carboxylic_acid_dimer", tuple(atom_union))
            if group_key in seen:
                continue
            seen.add(group_key)
            augmented.append(
                FunctionalGroup(
                    group="carboxylic_acid_dimer",
                    atoms0=tuple(atom_union),
                    description="Two carboxylic acid groups linked by reciprocal intermolecular H-bonds",
                    confidence="high",
                    evidence="paired carboxylic O-H···O=C H-bonds across fragments",
                )
            )
            for row in hbonds:
                donor_atoms = tuple(sorted(int(x) for x in str(row.get("donor_group_atoms", "")).split("-") if x))
                acceptor_atoms = tuple(sorted(int(x) for x in str(row.get("acceptor_group_atoms", "")).split("-") if x))
                if frozenset((donor_atoms, acceptor_atoms)) == key and row.get("chem_type") == "carboxylic_acid_hbond":
                    row["chem_type"] = "acid_dimer_hbond"
                    row["context_label"] = "carboxylic acid dimer H-bond"
        return augmented

    def _filter_legacy_conflicts(
        self,
        primary: Sequence[FunctionalGroup],
        fallback: Sequence[FunctionalGroup],
    ) -> List[FunctionalGroup]:
        primary_groups = {group.group for group in primary}
        carboxyl_atoms = {
            int(atom)
            for group in primary
            if group.group == "carboxylic_acid"
            for atom in group.atoms0
        }
        heteroaromatic_n_atoms = {
            int(atom)
            for group in primary
            if group.group == "heteroaromatic_N"
            for atom in group.atoms0
        }

        filtered: List[FunctionalGroup] = []
        for group in fallback:
            atom_set = {int(atom) for atom in group.atoms0}

            if "carboxylic_acid" in primary_groups and group.group in {"alcohol", "ester", "carbonyl_C=O"}:
                if atom_set & carboxyl_atoms:
                    continue

            if "heteroaromatic_N" in primary_groups and group.group == "dialkyl_amide_N_or_amine_N":
                if atom_set & heteroaromatic_n_atoms:
                    continue

            filtered.append(group)
        return filtered

    def _merge_functional_groups(
        self,
        primary: Sequence[FunctionalGroup],
        fallback: Sequence[FunctionalGroup],
    ) -> List[FunctionalGroup]:
        merged: List[FunctionalGroup] = []
        seen = set()
        filtered_fallback = self._filter_legacy_conflicts(primary, fallback)
        for group in list(primary) + list(filtered_fallback):
            key = (group.group, tuple(sorted(int(a) for a in group.atoms0)))
            if key in seen:
                continue
            seen.add(key)
            merged.append(group)
        return merged

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
        from chemistry import _legacy_detect_interfragment_hbonds
        hbonds = _legacy_detect_interfragment_hbonds(
            atoms,
            coords_A,
            bonds,
            fragments,
            h_a_max_A=h_a_max_A,
            d_a_max_A=d_a_max_A,
            angle_min_deg=angle_min_deg,
        )
        groups = self.detect_functional_groups(atoms, coords_A, bonds)
        return self._annotate_hbond_contexts(hbonds, groups, atoms, fragments)

    def atom_environment_table(self, atoms, coords_A, bonds) -> pd.DataFrame:
        from chemistry import _legacy_atom_environment_table

        return _legacy_atom_environment_table(atoms, coords_A, bonds)

    def detect_functional_groups(self, atoms, coords_A, bonds) -> List[FunctionalGroup]:
        from chemistry import _legacy_detect_functional_groups

        mol = self._build_mol(atoms, coords_A)
        rdkit_groups = self._functional_groups_from_mol(mol, atoms)
        legacy_groups = _legacy_detect_functional_groups(atoms, coords_A, bonds)
        return self._merge_functional_groups(rdkit_groups, legacy_groups)

    def annotate_chemical_system(self, atoms, coords_A) -> ChemicalSystemAnnotation:
        from chemistry import classify_system

        mol = self._build_mol(atoms, coords_A)
        bonds = tuple(self.build_connectivity(atoms, coords_A))
        fragments = self._fragments_from_mol(mol)
        groups = list(self.detect_functional_groups(atoms, coords_A, bonds))
        hbonds = list(self.detect_interfragment_hbonds(atoms, coords_A, bonds, fragments))
        groups = self._augment_supramolecular_groups(groups, hbonds)
        return ChemicalSystemAnnotation(
            formula=formula_string(atoms),
            system_type=classify_system(fragments),
            bonds=bonds,
            fragments=fragments,
            functional_groups=tuple(groups),
            interfragment_hbonds=tuple(hbonds),
        )
