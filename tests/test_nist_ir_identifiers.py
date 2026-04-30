from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir.identifiers import hess_to_identifiers, smiles_to_identifiers  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402


def test_smiles_to_identifiers_acetophenone():
    ids = smiles_to_identifiers("CC(=O)c1ccccc1")
    assert ids["canonical_smiles"] == "CC(=O)c1ccccc1"
    assert ids["inchi"].startswith("InChI=1S/C8H8O")
    assert ids["inchikey"] == "KWOLFJPFCHCOCG-UHFFFAOYSA-N"


def test_hess_to_identifiers_acetophenone():
    hess = read_orca_hess(ROOT / "data" / "hess" / "acetophenone.hess")
    ids = hess_to_identifiers(hess)
    assert ids["canonical_smiles"] == "CC(=O)c1ccccc1"
    assert ids["inchi"].startswith("InChI=1S/C8H8O")
    assert ids["inchikey"] == "KWOLFJPFCHCOCG-UHFFFAOYSA-N"


def test_hess_to_identifiers_identical_dimer_collapses_to_monomer():
    hess = read_orca_hess(ROOT / "data" / "hess" / "monoethanolamine_dimer_OH_to_N_DFT.hess")
    ids = hess_to_identifiers(hess)
    assert ids["canonical_smiles"] == "NCCO"
    assert ids["inchi"].startswith("InChI=1S/C2H7NO")
    assert ids["inchikey"] == "HZAXFHJVJLSVMW-UHFFFAOYSA-N"
    assert ids["fragment_count"] == "2"
    assert ids["selected_fragment_strategy"] == "first_identical_fragment"
