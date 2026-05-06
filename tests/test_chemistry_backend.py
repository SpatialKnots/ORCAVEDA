from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from chemistry import annotate_chemical_system, get_active_backend_name, list_backends, set_active_backend  # noqa: E402
from chemistry_rdkit_backend import RDKitChemistryBackend  # noqa: E402
from mode_assignment import _assignment_family_from_internal  # noqa: E402
from orcaveda_models import FunctionalGroup, InternalCoordinate  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402


def test_default_chemistry_backend_is_legacy():
    assert "legacy" in list_backends()
    assert get_active_backend_name() == "legacy"
    set_active_backend("legacy")
    assert get_active_backend_name() == "legacy"
    set_active_backend("LeGaCy")
    assert get_active_backend_name() == "legacy"


def test_legacy_backend_annotation_matches_current_regression_case():
    set_active_backend("legacy")
    hess = read_orca_hess(ROOT / "data" / "hess" / "CH3CN_freq.hess")
    annotation = annotate_chemical_system(hess.atoms, hess.coords_A)

    assert annotation.formula == "C2H3N"
    assert annotation.system_type == "monomer"
    assert set(annotation.functional_group_labels) == {"methyl", "nitrile_C≡N"}


def test_rdkit_backend_is_registered_when_available():
    try:
        import rdkit  # noqa: F401
    except ModuleNotFoundError:
        return

    assert "rdkit" in list_backends()


def test_rdkit_detects_extended_golden_semantic_groups():
    Chem = pytest.importorskip("rdkit.Chem")
    AllChem = pytest.importorskip("rdkit.Chem.AllChem")

    backend = RDKitChemistryBackend()

    def groups_for(smiles: str) -> set[str]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        return {group.group for group in backend._functional_groups_from_mol(mol, [atom.GetSymbol() for atom in mol.GetAtoms()])}

    assert "thiol" in groups_for("CS")
    assert "nitro" in groups_for("C[N+](=O)[O-]")
    assert {"thiocarbonyl_C=S", "thioamide"} <= groups_for("CC(N)=S")
    assert {"alkyne_C#C", "terminal_alkyne_C#C-H"} <= groups_for("CC#C")
    assert {"alkene_C=C", "vinylic_C-H"} <= groups_for("C=C")
    sulfone = groups_for("CS(=O)(=O)C")
    assert "sulfone" in sulfone
    assert "sulfoxide" not in sulfone
    assert "sulfoxide_S=O" not in sulfone
    assert "thioether" in groups_for("CSC")
    assert "acyl_chloride" in groups_for("CC(=O)Cl")
    assert "acid_anhydride" in groups_for("CC(=O)OC(C)=O")
    assert "carbonate_ester" in groups_for("COC(=O)OC")
    assert "epoxide" in groups_for("C1CO1")
    assert "cyclic_acetal" in groups_for("C1OCO1")
    assert "oxime" in groups_for("c1ccccc1C=NO")

    benzaldimine = groups_for("c1ccccc1C=N")
    assert "imine_C=N" in benzaldimine
    assert "alkene_C=C" not in benzaldimine
    assert "vinylic_C-H" not in benzaldimine
    assert "secondary_amine" not in benzaldimine

    phenyl_isocyanate = groups_for("O=C=Nc1ccccc1")
    assert "isocyanate_NCO" in phenyl_isocyanate
    assert "amide" not in phenyl_isocyanate
    assert "carbonyl_C=O" not in phenyl_isocyanate
    assert not any(group.startswith("nitrile") for group in phenyl_isocyanate)


def test_rdkit_merge_filters_legacy_isocyanate_nitrile_overlap():
    pytest.importorskip("rdkit")
    backend = RDKitChemistryBackend()
    primary = [
        FunctionalGroup("nitrile_C≡N", (12, 1), "rdkit", "high", "RDKit C#N"),
    ]
    fallback = [
        FunctionalGroup("isocyanate_NCO", (1, 12, 13), "legacy", "high", "N=C=O"),
        FunctionalGroup("nitrile_C≡N", (12, 1), "legacy", "high", "short C-N"),
    ]
    merged = backend._merge_functional_groups(primary, fallback)
    labels = {group.group for group in merged}
    assert "isocyanate_NCO" in labels
    assert not any(group.startswith("nitrile") for group in labels)


def test_rdkit_merge_filters_legacy_sulfone_and_imine_overlaps():
    pytest.importorskip("rdkit")
    backend = RDKitChemistryBackend()
    merged = backend._merge_functional_groups(
        [
            FunctionalGroup("sulfone", (0, 1, 2, 3, 4), "rdkit", "high", "SO2"),
            FunctionalGroup("imine_C=N", (5, 6), "rdkit", "high", "C=N"),
        ],
        [
            FunctionalGroup("sulfoxide", (0, 1, 3, 4), "legacy", "high", "S=O"),
            FunctionalGroup("sulfoxide_S=O", (0, 2), "legacy", "high", "S=O"),
            FunctionalGroup("secondary_amine", (6, 7, 5), "legacy", "high", "N-H and N-C"),
        ],
    )
    labels = {group.group for group in merged}
    assert "sulfone" in labels
    assert "imine_C=N" in labels
    assert "sulfoxide" not in labels
    assert "sulfoxide_S=O" not in labels
    assert "secondary_amine" not in labels


def test_extended_semantic_internal_coordinate_labels():
    def dummy(_xyz):
        return 0.0
    cases = {
        "fg_nitro_NO_stretch": "nitro N-O stretch",
        "fg_thiol_SH_stretch": "thiol S-H stretch",
        "fg_thiocarbonyl_CS_stretch": "thiocarbonyl C=S stretch",
        "fg_isocyanate_NC_stretch": "isocyanate N=C stretch",
        "fg_alkyne_CC_stretch": "alkyne C#C stretch",
        "fg_terminal_alkyne_CH_stretch": "terminal alkyne C-H stretch",
        "fg_alkene_CC_stretch": "alkene C=C stretch",
        "fg_vinylic_CH_stretch": "vinylic C-H stretch",
        "fg_sulfone_SO_stretch": "sulfone S=O stretch",
        "fg_imine_CN_stretch": "imine C=N stretch",
        "fg_oxime_NO_stretch": "oxime N-O stretch",
        "fg_acyl_chloride_CCl_stretch": "acyl chloride C-Cl stretch",
        "fg_anhydride_CO_stretch": "anhydride C=O stretch",
        "fg_epoxide_ring_bend": "epoxide ring bend",
    }
    for kind, expected in cases.items():
        coord = InternalCoordinate(f"FG_test({kind})", kind, (0, 1), 1, dummy, "functional_group_template")
        assert _assignment_family_from_internal(coord) == expected
