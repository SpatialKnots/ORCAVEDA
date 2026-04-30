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
    "acetanilide.hess",
    "acetophenone.hess",
    "H2O2_freq.hess",
]


@pytest.fixture(scope="module")
def evening_outputs():
    outdir = ROOT / "outputs" / "pytest_evening_semantics"
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


def test_evening_semantic_groups(evening_outputs):
    groups = evening_outputs["functional_groups"]

    assert "carboxylic_acid" in _groups_for(groups, "benzoic_acid.hess")
    assert "aryl_amide" in _groups_for(groups, "acetanilide.hess")
    assert "aryl_ketone" in _groups_for(groups, "acetophenone.hess")
    assert "peroxide" in _groups_for(groups, "H2O2_freq.hess")


def test_evening_semantic_assignments(evening_outputs):
    audit = evening_outputs["assignment_audit"]

    benzoic = _assignments_for(audit, "benzoic_acid.hess")
    assert benzoic.str.contains("carboxylic O-H stretch", regex=False).any()
    assert benzoic.str.contains("carboxylic O-H bend", regex=False).any()
    assert benzoic.str.contains("carboxylic C-O stretch", regex=False).any()
    assert benzoic.str.contains("carboxylic C=O stretch", regex=False).any()

    acetanilide = _assignments_for(audit, "acetanilide.hess")
    assert acetanilide.str.contains("aryl amide N-H stretch", regex=False).any()
    assert acetanilide.str.contains("aryl amide C-N stretch", regex=False).any()
    assert acetanilide.str.contains("aryl amide C=O stretch", regex=False).any()

    acetophenone = _assignments_for(audit, "acetophenone.hess")
    assert acetophenone.str.contains("aryl-conjugated C=O stretch", regex=False).any()

    h2o2 = _assignments_for(audit, "H2O2_freq.hess")
    assert h2o2.str.contains("O-O stretch", regex=False).any()
    assert h2o2.str.contains("H-O-O bend", regex=False).any()
    assert h2o2.str.contains("O-H stretch", regex=False).any()
    assert h2o2.str.contains("torsion", regex=False).any()
