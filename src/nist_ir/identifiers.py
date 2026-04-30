from __future__ import annotations

from typing import Dict

import numpy as np
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D


def structure_to_mol(atoms, coords_A, charge: int = 0):
    atoms = [str(atom) for atom in atoms]
    coords = np.asarray(coords_A, dtype=float)
    if coords.shape != (len(atoms), 3):
        raise ValueError(f"coords_A shape mismatch: expected {(len(atoms), 3)}, got {coords.shape}")

    mol = Chem.RWMol()
    for atom_symbol in atoms:
        mol.AddAtom(Chem.Atom(atom_symbol))

    conf = Chem.Conformer(len(atoms))
    for idx, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf)

    mol = mol.GetMol()
    rdDetermineBonds.DetermineBonds(mol, charge=int(charge))
    Chem.SanitizeMol(mol)
    return mol


def structure_to_identifiers(atoms, coords_A, charge: int = 0) -> Dict[str, str]:
    mol = structure_to_mol(atoms, coords_A, charge=charge)
    collapsed = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
    target = collapsed
    selected_fragment_strategy = "full_structure"

    frags = list(Chem.GetMolFrags(collapsed, asMols=True, sanitizeFrags=True))
    if len(frags) > 1:
        frag_smiles = [Chem.MolToSmiles(frag, canonical=True) for frag in frags]
        if len(set(frag_smiles)) == 1:
            target = frags[0]
            selected_fragment_strategy = "first_identical_fragment"

    canonical_smiles = Chem.MolToSmiles(target, canonical=True)
    inchi_str = inchi.MolToInchi(target)
    inchikey = inchi.MolToInchiKey(target)

    return {
        "input_smiles": canonical_smiles,
        "canonical_smiles": canonical_smiles,
        "inchi": inchi_str,
        "inchikey": inchikey,
        "fragment_count": str(len(frags)),
        "selected_fragment_strategy": selected_fragment_strategy,
    }


def smiles_to_identifiers(smiles: str) -> Dict[str, str]:
    raw = str(smiles or "").strip()
    if not raw:
        raise ValueError("SMILES is empty")

    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {raw}")

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    inchi_str = inchi.MolToInchi(mol)
    inchikey = inchi.MolToInchiKey(mol)

    return {
        "input_smiles": raw,
        "canonical_smiles": canonical_smiles,
        "inchi": inchi_str,
        "inchikey": inchikey,
    }


def hess_to_identifiers(hess, charge: int = 0) -> Dict[str, str]:
    if not hasattr(hess, "atoms") or not hasattr(hess, "coords_A"):
        raise TypeError("hess object must expose .atoms and .coords_A")
    return structure_to_identifiers(hess.atoms, hess.coords_A, charge=charge)
