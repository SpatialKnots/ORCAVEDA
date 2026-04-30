from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orca_parser import read_orca_hess  # noqa: E402
from reports import build_spectrum_payload, write_interactive_spectrum_viewer  # noqa: E402


def test_interactive_spectrum_viewer_artifacts():
    hess = read_orca_hess(ROOT / "data" / "hess" / "H2O_freq.hess")
    positive_mode = next(idx for idx, freq in enumerate(hess.frequencies_cm1) if float(freq) > 0.0)
    assignment_audit = pd.DataFrame(
        [
            {
                "Filename": hess.filename,
                "mode": positive_mode,
                "frequency_cm-1": float(hess.frequencies_cm1[positive_mode]),
                "IR_intensity": float(hess.ir_intensities[positive_mode]),
                "functional_group_assignment": "O-H stretch",
                "top_internal_coordinates": "r(O1-H2)=50.0%; r(O1-H3)=50.0%",
                "warnings": "",
            }
        ]
    )

    payload = build_spectrum_payload([hess], assignment_audit)
    outdir = ROOT / "outputs" / "pytest_interactive_spectrum_viewer"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / "viewer.html"
    json_path = outdir / "viewer.json"
    write_interactive_spectrum_viewer(payload, html_path, json_path=json_path)

    assert html_path.exists()
    assert json_path.exists()

    html_text = html_path.read_text(encoding="utf-8")
    assert "Interactive IR Spectrum" in html_text
    assert "3D Molecule Viewer" in html_text
    assert "moleculeViewer" in html_text
    assert "3Dmol-min.js" in html_text
    assert "molStyle" in html_text

    json_text = json_path.read_text(encoding="utf-8")
    assert "frequency_cm1" in json_text
    assert hess.filename in json_text
    assert "O-H stretch" in json_text
    assert "\"geometry\"" in json_text
    assert "\"atoms\"" in json_text
    assert "\"bonds\"" in json_text
