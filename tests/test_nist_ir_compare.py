from __future__ import annotations

import sys
from pathlib import Path
import shutil

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir.compare import (  # noqa: E402
    load_orcaveda_assignments,
    match_reference_to_orcaveda,
    pick_reference_peaks,
)


def test_pick_reference_peaks_simple():
    spectrum = pd.DataFrame(
        {
            "wavenumber_cm-1": [1000, 1004, 1008, 1012, 1016, 1020, 1024],
            "intensity": [0.1, 0.5, 0.2, 0.8, 0.3, 0.6, 0.1],
        }
    )
    peaks = pick_reference_peaks(spectrum, top_n=2, min_separation_cm1=4.0)
    assert len(peaks) == 2
    assert set(round(v, 1) for v in peaks["wavenumber_cm-1"]) == {1012.0, 1020.0}


def test_match_reference_to_orcaveda_simple():
    tmp_path = ROOT / "outputs" / "pytest_nist_ir_compare"
    if tmp_path.exists():
        shutil.rmtree(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    audit_path = tmp_path / "audit.csv"
    pd.DataFrame(
        [
            {
                "mode": 10,
                "frequency_cm-1": 1710.0,
                "IR_intensity": 300.0,
                "functional_group_assignment": "C=O stretch",
            },
            {
                "mode": 11,
                "frequency_cm-1": 1260.0,
                "IR_intensity": 120.0,
                "functional_group_assignment": "C-O stretch",
            },
        ]
    ).to_csv(audit_path, index=False, encoding="utf-8")

    audit = load_orcaveda_assignments(audit_path, scale_factor=0.96)
    reference = pd.DataFrame(
        [
            {"wavenumber_cm-1": 1640.0, "intensity": 0.8},
            {"wavenumber_cm-1": 1210.0, "intensity": 0.4},
        ]
    )
    matched = match_reference_to_orcaveda(reference, audit)

    assert list(matched["orcaveda_mode"]) == [10, 11]
    assert list(matched["orcaveda_assignment"]) == ["C=O stretch", "C-O stretch"]
