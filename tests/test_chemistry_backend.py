from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from chemistry import annotate_chemical_system, get_active_backend_name, set_active_backend  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402


def test_default_chemistry_backend_is_legacy():
    assert get_active_backend_name() == "legacy"
    set_active_backend("legacy")
    assert get_active_backend_name() == "legacy"


def test_legacy_backend_annotation_matches_current_regression_case():
    set_active_backend("legacy")
    hess = read_orca_hess(ROOT / "data" / "hess" / "CH3CN_freq.hess")
    annotation = annotate_chemical_system(hess.atoms, hess.coords_A)

    assert annotation.formula == "C2H3N"
    assert annotation.system_type == "monomer"
    assert set(annotation.functional_group_labels) == {"methyl", "nitrile_C≡N"}
