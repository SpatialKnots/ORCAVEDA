from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ORCAVEDA_patched_stage3D_v5_0 import analyze_orca_ped_like  # noqa: E402
from run_regression_tests import run_checks  # noqa: E402


REGRESSION_HESS = [
    "Acetone_freq.hess",
    "CH3CN_freq.hess",
    "DMF_freq.hess",
    "DMSO_freq.hess",
    "EtOH_freq.hess",
    "MeOH_freq.hess",
    "NMP_freq.hess",
    "iPrOH_freq.hess",
]

EXPECTED_SUMMARY = {
    "Acetone_freq.hess": {
        "formula": "C3H6O",
        "rank_B_independent": 24,
        "expected_rank_3N_minus_6": 24,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "carbonyl_C=O; ketone; methyl",
    },
    "CH3CN_freq.hess": {
        "formula": "C2H3N",
        "rank_B_independent": 12,
        "expected_rank_3N_minus_6": 12,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "methyl; nitrile_C≡N",
    },
    "DMF_freq.hess": {
        "formula": "C3H7NO",
        "rank_B_independent": 30,
        "expected_rank_3N_minus_6": 30,
        "negative_freq_count_after_first_6": 2,
        "functional_groups_detected": "amide; carbonyl_C=O; methine; methyl",
    },
    "DMSO_freq.hess": {
        "formula": "C2H6OS",
        "rank_B_independent": 24,
        "expected_rank_3N_minus_6": 24,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "methyl; sulfoxide; sulfoxide_S=O",
    },
    "EtOH_freq.hess": {
        "formula": "C2H6O",
        "rank_B_independent": 21,
        "expected_rank_3N_minus_6": 21,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "alcohol; methyl; methylene",
    },
    "MeOH_freq.hess": {
        "formula": "CH4O",
        "rank_B_independent": 12,
        "expected_rank_3N_minus_6": 12,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "alcohol; methyl",
    },
    "NMP_freq.hess": {
        "formula": "C5H9NO",
        "rank_B_independent": 42,
        "expected_rank_3N_minus_6": 42,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "carbonyl_C=O; lactam_amide; methyl; methylene; ring",
    },
    "iPrOH_freq.hess": {
        "formula": "C3H8O",
        "rank_B_independent": 30,
        "expected_rank_3N_minus_6": 30,
        "negative_freq_count_after_first_6": 0,
        "functional_groups_detected": "alcohol; methine; methyl",
    },
}


def test_regression_baseline_outputs():
    outdir = ROOT / "outputs" / "pytest_regression_outputs"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    hess_paths = [ROOT / "data" / "hess" / name for name in REGRESSION_HESS]

    tables = analyze_orca_ped_like(hess_paths, outdir)

    assert "general_summary" in tables
    assert "assignment_audit" in tables

    expectations = json.loads((ROOT / "expectations" / "regression_expectations_stage3D_v5_0.json").read_text(encoding="utf-8"))
    results = run_checks(outdir, expectations)
    failed = results[results["status"] == "FAIL"]
    assert failed.empty, failed.to_string(index=False)

    summary_path = next(outdir.glob("*__general_summary.csv"))
    summary = pd.read_csv(summary_path)

    observed = {
        row["Filename"]: {
            "formula": row["formula"],
            "rank_B_independent": int(row["rank_B_independent"]),
            "expected_rank_3N_minus_6": int(row["expected_rank_3N_minus_6"]),
            "negative_freq_count_after_first_6": int(row["negative_freq_count_after_first_6"]),
            "functional_groups_detected": row["functional_groups_detected"],
        }
        for _, row in summary.iterrows()
    }

    assert observed == EXPECTED_SUMMARY
