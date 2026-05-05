from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

pytest.importorskip("rdkit")

from chemistry_rdkit_backend import RDKitChemistryBackend  # noqa: E402
from internal_coordinates import build_internal_coordinates  # noqa: E402
from mode_assignment import _assignment_family_from_internal  # noqa: E402
from orcaveda_models import FunctionalGroup, InternalCoordinate  # noqa: E402


def test_rdkit_detects_synthetic_carboxylate():
    backend = RDKitChemistryBackend()
    atoms = ["C", "O", "O", "C", "H", "H", "H"]
    coords = np.array(
        [
            [0.00, 0.00, 0.00],
            [1.26, 0.00, 0.00],
            [-0.63, 1.09, 0.00],
            [-1.48, -0.18, 0.00],
            [-1.88, 0.37, 0.89],
            [-1.86, -1.23, -0.05],
            [-1.86, 0.45, -0.84],
        ],
        dtype=float,
    )

    bonds = backend.build_connectivity(atoms, coords)
    groups = backend.detect_functional_groups(atoms, coords, bonds)
    labels = {group.group for group in groups}

    assert "carboxylate" in labels


def test_acid_dimer_hbond_context_and_group_augmentation():
    backend = RDKitChemistryBackend()
    atoms = ["C", "O", "O", "H", "C", "O", "O", "H"]
    fragments = ((0, 1, 2, 3), (4, 5, 6, 7))
    groups = [
        FunctionalGroup("carboxylic_acid", (0, 1, 2), "", "high", ""),
        FunctionalGroup("carboxylic_acid", (4, 5, 6), "", "high", ""),
    ]
    hbonds = [
        {"D0": 2, "H0": 3, "A0": 5, "D": 3, "H": 4, "A": 6, "type": "O-H···O", "rHA_A": 1.70, "rDA_A": 2.65, "angle_deg": 176.0},
        {"D0": 6, "H0": 7, "A0": 1, "D": 7, "H": 8, "A": 2, "type": "O-H···O", "rHA_A": 1.71, "rDA_A": 2.66, "angle_deg": 175.0},
    ]

    annotated = backend._annotate_hbond_contexts(hbonds, groups, atoms, fragments)
    assert all(row["chem_type"] == "carboxylic_acid_hbond" for row in annotated)

    augmented = backend._augment_supramolecular_groups(groups, annotated)
    assert "carboxylic_acid_dimer" in {group.group for group in augmented}
    assert all(row["chem_type"] == "acid_dimer_hbond" for row in annotated)


def test_supramolecular_internal_coordinate_labels():
    atoms = ["C", "O", "O", "H", "C", "O", "O", "H"]
    coords = np.zeros((8, 3), dtype=float)
    bonds = (
        (0, 1, 1.23),
        (0, 2, 1.34),
        (2, 3, 0.98),
        (4, 5, 1.23),
        (4, 6, 1.34),
        (6, 7, 0.98),
    )
    fragments = ((0, 1, 2, 3), (4, 5, 6, 7))
    groups = [
        FunctionalGroup("carboxylic_acid", (0, 1, 2), "", "high", ""),
        FunctionalGroup("carboxylic_acid", (4, 5, 6), "", "high", ""),
        FunctionalGroup("carboxylic_acid_dimer", (0, 1, 2, 4, 5, 6), "", "high", ""),
        FunctionalGroup("carboxylate", (0, 1, 2), "", "high", ""),
    ]
    hbonds = [
        {
            "D0": 2,
            "H0": 3,
            "A0": 5,
            "D": 3,
            "H": 4,
            "A": 6,
            "type": "O-H···O",
            "rHA_A": 1.70,
            "rDA_A": 2.65,
            "angle_deg": 176.0,
            "chem_type": "acid_dimer_hbond",
        }
    ]

    internals = build_internal_coordinates(atoms, coords, bonds, fragments, hbonds, groups)
    labels = {_assignment_family_from_internal(ic) for ic in internals if "acid_dimer" in ic.name.lower() or "carboxylate" in ic.name.lower()}

    assert "acid dimer H-bond / intermolecular" in labels
    assert "carboxylate C-O stretch" in labels
    assert "carboxylate O-C-O bend" in labels


def test_generic_primitive_coordinate_labels_are_specific():
    dummy_fn = lambda coords: 0.0
    cn = InternalCoordinate("r(C1-N2)", "stretch", (0, 1), 10, dummy_fn)
    ccn = InternalCoordinate("ang(C1-C2-N3)", "angle", (0, 1, 2), 10, dummy_fn)
    ocn = InternalCoordinate("ang(O1-C2-N3)", "angle", (0, 1, 2), 10, dummy_fn)
    cnh = InternalCoordinate("ang(C1-N2-H3)", "angle", (0, 1, 2), 10, dummy_fn)
    lactam = InternalCoordinate("FG_lactam_ring_CNC_bend(C1-N2-C3)", "fg_lactam_ring_deformation", (0, 1, 2), 10, dummy_fn)
    aromatic_ring = InternalCoordinate("FG_aromatic_ring_CCC_bend(C1-C2-C3)", "fg_aromatic_ring_deformation", (0, 1, 2), 10, dummy_fn)

    assert _assignment_family_from_internal(cn) == "C-N stretch"
    assert _assignment_family_from_internal(ccn) == "C-C-N bend"
    assert _assignment_family_from_internal(ocn) == "amide carbonyl-adjacent C-N bend"
    assert _assignment_family_from_internal(cnh) == "N-H bend"
    assert _assignment_family_from_internal(lactam) == "lactam ring deformation"
    assert _assignment_family_from_internal(aromatic_ring) == "aromatic ring deformation"
