from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

pytest.importorskip("rdkit")

from ORCAVEDA_patched_stage3D_v5_0 import analyze_orca_ped_like  # noqa: E402
from chemistry import get_active_backend_name, set_active_backend  # noqa: E402


TARGET_HESS = [
    "benzoic_acid.hess",
    "pyridine.hess",
    "anisole.hess",
    "N-methylaniline.hess",
]


@pytest.fixture(scope="module")
def extended_rdkit_outputs():
    outdir = ROOT / "outputs" / "pytest_extended_golden_rdkit_outputs"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    hess_paths = [ROOT / "data" / "hess" / name for name in TARGET_HESS]

    previous_backend = get_active_backend_name()
    set_active_backend("rdkit")
    try:
        analyze_orca_ped_like(hess_paths, outdir)
    finally:
        set_active_backend(previous_backend)

    functional_groups = pd.read_csv(next(outdir.glob("*__functional_groups.csv")))
    assignment_audit = pd.read_csv(next(outdir.glob("*__assignment_audit.csv")))
    return {
        "functional_groups": functional_groups,
        "assignment_audit": assignment_audit,
    }


def _groups_for(df: pd.DataFrame, filename: str) -> set[str]:
    rows = df[df["Filename"].astype(str) == filename]
    if rows.empty:
        return set()
    return set(rows["group"].astype(str))


def _assignments_for(df: pd.DataFrame, filename: str) -> pd.Series:
    rows = df[df["Filename"].astype(str) == filename]
    return rows["functional_group_assignment"].fillna("").astype(str)


def test_rdkit_removes_known_false_positives(extended_rdkit_outputs):
    groups = extended_rdkit_outputs["functional_groups"]

    benzoic = _groups_for(groups, "benzoic_acid.hess")
    assert "carboxylic_acid" in benzoic
    assert "alcohol" not in benzoic
    assert "ester" not in benzoic

    pyridine = _groups_for(groups, "pyridine.hess")
    assert "heteroaromatic_N" in pyridine
    assert "dialkyl_amide_N_or_amine_N" not in pyridine


def test_rdkit_extended_semantic_groups_present(extended_rdkit_outputs):
    groups = extended_rdkit_outputs["functional_groups"]

    assert "aryl_ether" in _groups_for(groups, "anisole.hess")
    assert "secondary_aryl_amine" in _groups_for(groups, "N-methylaniline.hess")


def test_rdkit_extended_semantic_assignments_present(extended_rdkit_outputs):
    audit = extended_rdkit_outputs["assignment_audit"]

    anisole_assignments = _assignments_for(audit, "anisole.hess")
    assert anisole_assignments.str.contains("aryl ether C-O stretch", regex=False).any()

    pyridine_assignments = _assignments_for(audit, "pyridine.hess")
    assert pyridine_assignments.str.contains("heteroaromatic C-N stretch", regex=False).any()

    n_methylaniline_assignments = _assignments_for(audit, "N-methylaniline.hess")
    assert n_methylaniline_assignments.str.contains("secondary aryl amine N-H stretch", regex=False).any()
