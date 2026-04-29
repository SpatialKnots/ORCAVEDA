from __future__ import annotations

import sys
import shutil
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


GOLDEN_HESS = [
    "H2O_freq.hess",
    "NH3.hess",
    "acetaldehyde.hess",
    "acetamide.hess",
    "aniline.hess",
    "phenol.hess",
    "benzene.hess",
]

EXPECTED_FUNCTIONAL_GROUPS = {
    "H2O_freq.hess": set(),
    "NH3.hess": set(),
    "acetaldehyde.hess": {"aldehyde", "carbonyl_C=O", "methine", "methyl"},
    "acetamide.hess": {"amide", "carbonyl_C=O", "methyl"},
    "aniline.hess": {"aniline", "aromatic_CH", "aromatic_ring", "primary_amine", "ring"},
    "phenol.hess": {"phenol", "aromatic_CH", "aromatic_ring", "alcohol", "ring"},
    "benzene.hess": {"aromatic_CH", "aromatic_ring", "ring"},
}


@pytest.fixture(scope="module")
def golden_rdkit_outputs():
    outdir = ROOT / "outputs" / "pytest_golden_rdkit_outputs"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    hess_paths = [ROOT / "data" / "hess" / name for name in GOLDEN_HESS]

    previous_backend = get_active_backend_name()
    set_active_backend("rdkit")
    try:
        tables = analyze_orca_ped_like(hess_paths, outdir)
    finally:
        set_active_backend(previous_backend)

    summary = pd.read_csv(next(outdir.glob("*__general_summary.csv")))
    functional_groups = pd.read_csv(next(outdir.glob("*__functional_groups.csv")))
    assignment_audit = pd.read_csv(next(outdir.glob("*__assignment_audit.csv")))
    return {
        "outdir": outdir,
        "tables": tables,
        "summary": summary,
        "functional_groups": functional_groups,
        "assignment_audit": assignment_audit,
    }


def _assignment_rows(df: pd.DataFrame, filename: str) -> pd.DataFrame:
    return df[df["Filename"].astype(str) == filename].copy()


def _assert_has_assignment(df: pd.DataFrame, filename: str, text: str) -> None:
    rows = _assignment_rows(df, filename)
    assignments = rows["functional_group_assignment"].fillna("").astype(str)
    assert assignments.str.contains(text, regex=False).any(), f"{filename}: missing assignment containing {text!r}"


def test_golden_rdkit_summary_and_functional_groups(golden_rdkit_outputs):
    summary = golden_rdkit_outputs["summary"]
    functional_groups = golden_rdkit_outputs["functional_groups"]

    observed_files = set(summary["Filename"].astype(str))
    assert observed_files == set(GOLDEN_HESS)

    for filename, expected_groups in EXPECTED_FUNCTIONAL_GROUPS.items():
        fg_rows = functional_groups[functional_groups["Filename"].astype(str) == filename]
        observed_groups = set(fg_rows["group"].astype(str)) if not fg_rows.empty else set()
        assert expected_groups.issubset(observed_groups), (
            f"{filename}: expected groups {sorted(expected_groups)} not subset of observed {sorted(observed_groups)}"
        )


def test_golden_rdkit_key_assignments(golden_rdkit_outputs):
    audit = golden_rdkit_outputs["assignment_audit"]

    _assert_has_assignment(audit, "H2O_freq.hess", "O-H stretch")
    _assert_has_assignment(audit, "NH3.hess", "NH2 symmetric stretch")
    _assert_has_assignment(audit, "NH3.hess", "NH2 asymmetric stretch")
    _assert_has_assignment(audit, "NH3.hess", "NH2 scissor")
    _assert_has_assignment(audit, "acetaldehyde.hess", "C=O stretch")
    _assert_has_assignment(audit, "acetamide.hess", "C=O stretch")
    _assert_has_assignment(audit, "acetamide.hess", "NH2 symmetric stretch")
    _assert_has_assignment(audit, "acetamide.hess", "NH2 asymmetric stretch")
    _assert_has_assignment(audit, "aniline.hess", "aromatic C-H stretch")
    _assert_has_assignment(audit, "aniline.hess", "aryl amine NH2 symmetric stretch")
    _assert_has_assignment(audit, "aniline.hess", "aryl amine NH2 asymmetric stretch")
    _assert_has_assignment(audit, "phenol.hess", "phenolic O-H stretch")
    _assert_has_assignment(audit, "phenol.hess", "aromatic C-H stretch")
    _assert_has_assignment(audit, "benzene.hess", "aromatic C-H stretch")


def test_golden_rdkit_no_unassigned_high_frequency_modes(golden_rdkit_outputs):
    audit = golden_rdkit_outputs["assignment_audit"]
    high = audit[pd.to_numeric(audit["frequency_cm-1"], errors="coerce") >= 2800.0].copy()
    unassigned = high[high["functional_group_assignment"].fillna("").astype(str).str.lower().eq("unassigned")]
    assert unassigned.empty, unassigned[["Filename", "mode", "frequency_cm-1"]].to_string(index=False)


def test_golden_rank_diagnostics_document_current_aromatic_deficits(golden_rdkit_outputs):
    summary = golden_rdkit_outputs["summary"]
    deficits = {
        str(row["Filename"]): int(row["expected_rank_3N_minus_6"]) - int(row["rank_B_independent"])
        for _, row in summary.iterrows()
        if int(row["expected_rank_3N_minus_6"]) - int(row["rank_B_independent"]) > 0
    }

    # This is a documented current limitation for aromatic systems rather than
    # a random failure. We lock it down so future chemistry work can improve it
    # intentionally instead of drifting silently.
    assert deficits == {
        "aniline.hess": 2,
        "phenol.hess": 2,
        "benzene.hess": 1,
    }
