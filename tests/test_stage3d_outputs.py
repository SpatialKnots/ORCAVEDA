"""
Pytest-compatible wrapper for ORCAVEDA Stage 3D v5.0 regression outputs.

Usage:
  pytest tests/test_stage3d_outputs.py --outdir <ORCAVEDA_OUTPUT_DIR> --expectations <EXPECTATIONS_JSON>

If pytest option injection is not configured, use run_regression_tests.py directly.
"""
from pathlib import Path
import json
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from run_regression_tests import run_checks  # noqa: E402


def test_stage3d_regression_outputs():
    # Default paths are useful when the repository is used with copied reports.
    outdir = Path.cwd()
    expectations = ROOT / "expectations" / "regression_expectations_stage3D_v5_0.json"
    results = run_checks(outdir, json.loads(expectations.read_text(encoding="utf-8")))
    failed = results[results["status"] == "FAIL"]
    assert failed.empty, failed.to_string(index=False)
