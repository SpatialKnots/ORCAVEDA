from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd
import json

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orca_parser import read_orca_hess  # noqa: E402
from reports import build_spectrum_payload, write_interactive_spectrum_viewer  # noqa: E402
from reports import attach_nist_reference_set  # noqa: E402


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


def test_interactive_spectrum_viewer_with_nist_reference():
    payload = {
        "viewer_title": "ORCAVEDA Interactive IR Spectrum",
        "default_scale_factor": 1.0,
        "default_lorentz_hwhm": 12.0,
        "files": [
            {
                "filename": "acetophenone.hess",
                "title": "acetophenone",
                "summary": {},
                "modes": [{"mode": 1, "frequency_cm1": 1700.0, "intensity": 100.0, "assignment": "C=O stretch", "top_internal_coordinates": "", "warnings": ""}],
                "geometry": {"atoms": [], "bonds": []},
            }
        ],
    }
    outdir = ROOT / "outputs" / "pytest_interactive_spectrum_viewer_nist"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_csv = outdir / "ref.csv"
    pd.DataFrame({"wavenumber_cm-1": [1600.0, 1700.0, 1800.0], "intensity": [0.1, 0.5, 0.2]}).to_csv(ref_csv, index=False)
    ref_meta = outdir / "ref_meta.json"
    ref_meta.write_text(json.dumps({"jcamp_metadata": {"YUNITS": "ABSORBANCE", "STATE": "gas"}}), encoding="utf-8")
    manifest_path = outdir / "refset.json"
    manifest_path.write_text(
        json.dumps(
            {
                "inchikey": "TEST",
                "canonical_smiles": "CC",
                "reference_spectra": [
                    {
                        "csv": str(ref_csv),
                        "meta_json": str(ref_meta),
                        "index": "0",
                        "phase_tag": "gas",
                        "phase_label": "gas",
                        "selection_priority": 100,
                        "description": "gas",
                    }
                ],
                "preferred_reference": {"index": "0", "phase_tag": "gas", "phase_label": "gas"},
            }
        ),
        encoding="utf-8",
    )

    payload = attach_nist_reference_set(payload, manifest_path)
    ref_payload = payload["nist_reference_sets"]["acetophenone"]["reference_spectra"][0]["scale_engine_payload"]
    assert "engine_table" in ref_payload
    assert "engine_fits" in ref_payload
    assert "global_ls" in ref_payload["engine_fits"]

    html_path = outdir / "viewer.html"
    json_path = outdir / "viewer.json"
    write_interactive_spectrum_viewer(payload, html_path, json_path=json_path)

    html_text = html_path.read_text(encoding="utf-8")
    assert "NIST Reference" in html_text
    assert "Scale Engine" in html_text
    assert "scaleEngine" in html_text
    assert "engineTableBody" in html_text
    assert "Mean %Δ" in html_text
    assert "RMS %Δ" in html_text
    assert "Max %Δ" in html_text
    assert "nistReference" in html_text
    assert "Auto-fit scale" in html_text
    assert "fitSummary" in html_text
    json_text = json_path.read_text(encoding="utf-8")
    assert "reference_spectra" in json_text
    assert "scale_engine_payload" in json_text
    assert "global_weighted_ls" in json_text
    assert "mean_percent_deviation" in json_text
