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
    "acetic_anhydride.hess",
    "acetaldehyde.hess",
    "acetamide.hess",
    "aniline.hess",
    "benzaldimine.hess",
    "benzaldoxime.hess",
    "phenol.hess",
    "benzene.hess",
    "cyclohexane_chair.hess",
    "dimethyl_carbonate.hess",
    "dimethyl_sulfide.hess",
    "dimethyl_sulfone.hess",
    "dimethylamine.hess",
    "ethene.hess",
    "ethylene_oxide.hess",
    "ethyne.hess",
    "gamma-butyrolactone.hess",
    "methanethiol.hess",
    "methyl_acetate.hess",
    "methylamine.hess",
    "nitrobenzene.hess",
    "nitromethane.hess",
    "phenyl_isocyanate.hess",
    "piperidine.hess",
    "propene.hess",
    "propyne.hess",
    "pyrrole.hess",
    "succinimide.hess",
    "tetrahydrofuran.hess",
    "trimethylamine.hess",
]

NEW_CLEAN_GOLDEN_HESS = [
    "acetic_anhydride.hess",
    "benzaldimine.hess",
    "benzaldoxime.hess",
    "cyclohexane_chair.hess",
    "dimethyl_carbonate.hess",
    "dimethyl_sulfide.hess",
    "dimethyl_sulfone.hess",
    "dimethylamine.hess",
    "ethene.hess",
    "ethylene_oxide.hess",
    "ethyne.hess",
    "gamma-butyrolactone.hess",
    "methanethiol.hess",
    "methyl_acetate.hess",
    "methylamine.hess",
    "nitrobenzene.hess",
    "nitromethane.hess",
    "phenyl_isocyanate.hess",
    "piperidine.hess",
    "propene.hess",
    "propyne.hess",
    "pyrrole.hess",
    "succinimide.hess",
    "tetrahydrofuran.hess",
    "trimethylamine.hess",
]

EXPECTED_FUNCTIONAL_GROUPS = {
    "H2O_freq.hess": set(),
    "NH3.hess": set(),
    "acetic_anhydride.hess": {"acid_anhydride", "carbonyl_C=O", "ester", "ether", "methyl"},
    "acetaldehyde.hess": {"aldehyde", "carbonyl_C=O", "methine", "methyl"},
    "acetamide.hess": {"amide", "carbonyl_C=O", "methyl"},
    "aniline.hess": {"aniline", "aromatic_CH", "aromatic_ring", "primary_amine", "ring"},
    "benzaldimine.hess": {"aromatic_CH", "aromatic_ring", "imine_C=N", "methine", "ring"},
    "benzaldoxime.hess": {"aromatic_CH", "aromatic_ring", "imine_C=N", "methine", "oxime", "ring"},
    "phenol.hess": {"phenol", "aromatic_CH", "aromatic_ring", "alcohol", "ring"},
    "benzene.hess": {"aromatic_CH", "aromatic_ring", "ring"},
    "cyclohexane_chair.hess": {"methylene", "ring"},
    "dimethyl_carbonate.hess": {"carbonate_ester", "carbonyl_C=O", "ester", "ether", "methyl"},
    "dimethyl_sulfide.hess": {"methyl", "thioether"},
    "dimethyl_sulfone.hess": {"methyl", "sulfone"},
    "dimethylamine.hess": {"methyl", "secondary_amine"},
    "ethene.hess": {"alkene_C=C", "methylene", "vinylic_C-H"},
    "ethylene_oxide.hess": {"epoxide", "ether", "methylene", "ring"},
    "ethyne.hess": {"alkyne_C#C", "methine", "terminal_alkyne_C#C-H"},
    "gamma-butyrolactone.hess": {"carbonyl_C=O", "ester", "ether", "lactone", "methylene", "ring"},
    "methanethiol.hess": {"methyl", "thiol"},
    "methyl_acetate.hess": {"carbonyl_C=O", "ester", "ether", "methyl"},
    "methylamine.hess": {"methyl", "primary_amine"},
    "nitrobenzene.hess": {"methine", "nitro", "ring"},
    "nitromethane.hess": {"methyl", "nitro"},
    "phenyl_isocyanate.hess": {"aromatic_CH", "aromatic_ring", "isocyanate_NCO", "methine", "ring"},
    "piperidine.hess": {"methylene", "ring", "secondary_amine"},
    "propene.hess": {"alkene_C=C", "methine", "methyl", "methylene", "vinylic_C-H"},
    "propyne.hess": {"alkyne_C#C", "methine", "methyl", "terminal_alkyne_C#C-H"},
    "pyrrole.hess": {"aromatic_CH", "aromatic_ring", "methine", "ring"},
    "succinimide.hess": {"carbonyl_C=O", "lactam_amide", "methylene", "ring"},
    "tetrahydrofuran.hess": {"ether", "methylene", "ring"},
    "trimethylamine.hess": {"methyl", "tertiary_amine"},
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


def test_added_clean_golden_hess_have_no_negative_vibrational_modes_after_first_six(golden_rdkit_outputs):
    summary = golden_rdkit_outputs["summary"]
    added = summary[summary["Filename"].astype(str).isin(NEW_CLEAN_GOLDEN_HESS)]
    flagged = added[
        added["system_flags"].fillna("").astype(str).str.contains("negative_vibrational_frequency_after_first_6", regex=False)
    ]
    assert flagged.empty, flagged[["Filename", "negative_freq_count_after_first_6", "system_flags"]].to_string(index=False)


def test_golden_rank_diagnostics_have_no_aromatic_deficits(golden_rdkit_outputs):
    summary = golden_rdkit_outputs["summary"]
    deficits = {
        str(row["Filename"]): int(row["expected_rank_3N_minus_6"]) - int(row["rank_B_independent"])
        for _, row in summary.iterrows()
        if int(row["expected_rank_3N_minus_6"]) - int(row["rank_B_independent"]) > 0
    }

    assert deficits == {}
