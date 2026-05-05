"""
Pytest-compatible wrapper for ORCAVEDA Stage 3D v5.0 regression outputs.

Usage:
  pytest tests/test_stage3d_outputs.py --outdir <ORCAVEDA_OUTPUT_DIR> --expectations <EXPECTATIONS_JSON>

If pytest option injection is not configured, use run_regression_tests.py directly.
"""
from pathlib import Path
import json
import shutil
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

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


def test_stage3d_regression_outputs():
    outdir = ROOT / "outputs" / "pytest_stage3d_outputs"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    hess_paths = [ROOT / "data" / "hess" / name for name in REGRESSION_HESS]
    analyze_orca_ped_like(hess_paths, outdir)

    expectations = ROOT / "expectations" / "regression_expectations_stage3D_v5_0.json"
    results = run_checks(outdir, json.loads(expectations.read_text(encoding="utf-8")))
    failed = results[results["status"] == "FAIL"]
    assert failed.empty, failed.to_string(index=False)
