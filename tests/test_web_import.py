from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from web_import import import_hess_files_for_web  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402
from reports import build_spectrum_payload, write_interactive_spectrum_viewer  # noqa: E402


def test_web_hess_import_copies_inputs_and_writes_manifest():
    outroot = ROOT / "outputs" / "pytest_web_import_fake"
    if outroot.exists():
        shutil.rmtree(outroot)

    def fake_pipeline(paths, outdir):
        outdir.mkdir(parents=True, exist_ok=True)
        manifest = {
            "xlsx_report": str(outdir / "fake.xlsx"),
            "interactive_spectrum_html": str(outdir / "fake.html"),
            "interactive_spectrum_data_json": str(outdir / "fake.json"),
        }
        (outdir / "sample__integration_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return {"assignment_audit": pd.DataFrame()}

    result = import_hess_files_for_web(
        [ROOT / "data" / "hess" / "H2O_freq.hess"],
        import_root=outroot,
        run_id="unit_run",
        pipeline_runner=fake_pipeline,
    )

    assert result.status == "completed"
    assert result.input_files == ("H2O_freq.hess",)
    assert Path(result.input_paths[0]).exists()
    assert result.artifacts["interactive_spectrum_html"].endswith("fake.html")
    manifest_path = outroot / "unit_run" / "web_import_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["run_id"] == "unit_run"
    assert "accepted_hess_upload:H2O_freq.hess" in manifest["diagnostics"]


def test_web_hess_import_rejects_non_hess_file():
    outroot = ROOT / "outputs" / "pytest_web_import_reject"
    if outroot.exists():
        shutil.rmtree(outroot)

    with pytest.raises(ValueError, match="expected .hess"):
        import_hess_files_for_web([ROOT / "README.md"], import_root=outroot, run_id="bad_upload")


def test_web_hess_import_rejects_invalid_run_id():
    outroot = ROOT / "outputs" / "pytest_web_import_bad_run_id"
    if outroot.exists():
        shutil.rmtree(outroot)

    with pytest.raises(ValueError, match="Invalid run_id"):
        import_hess_files_for_web(
            [ROOT / "data" / "hess" / "H2O_freq.hess"],
            import_root=outroot,
            run_id="../bad",
            pipeline_runner=lambda _paths, _outdir: {"assignment_audit": pd.DataFrame()},
        )


def test_web_hess_import_runs_existing_pipeline_smoke():
    outroot = ROOT / "outputs" / "pytest_web_import_smoke"
    if outroot.exists():
        shutil.rmtree(outroot)

    result = import_hess_files_for_web(
        [ROOT / "data" / "hess" / "H2O_freq.hess"],
        import_root=outroot,
        run_id="h2o_smoke",
    )

    assert result.status == "completed"
    assert Path(result.artifacts["interactive_spectrum_html"]).exists()
    assert Path(result.artifacts["interactive_spectrum_data_json"]).exists()
    assert Path(result.artifacts["run_manifest_json"]).exists()
    assert any(item.startswith("pipeline_tables:") for item in result.diagnostics)


def test_web_hess_import_can_attach_nist_reference_set():
    outroot = ROOT / "outputs" / "pytest_web_import_nist"
    if outroot.exists():
        shutil.rmtree(outroot)

    def fake_pipeline(paths, outdir):
        outdir.mkdir(parents=True, exist_ok=True)
        hess = read_orca_hess(paths[0])
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
        html_path = outdir / "H2O__interactive_spectrum.html"
        json_path = outdir / "H2O__spectrum_data.json"
        write_interactive_spectrum_viewer(payload, html_path, json_path=json_path)
        manifest = {
            "xlsx_report": str(outdir / "H2O.xlsx"),
            "interactive_spectrum_html": str(html_path),
            "interactive_spectrum_data_json": str(json_path),
        }
        (outdir / "H2O__integration_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return {"assignment_audit": assignment_audit}

    def fake_nist_runner(_hess_path, nist_dir):
        nist_dir.mkdir(parents=True, exist_ok=True)
        csv_path = nist_dir / "ref.csv"
        meta_path = nist_dir / "ref_meta.json"
        pd.DataFrame(
            {
                "wavenumber_cm-1": [1500.0, 1600.0, 3700.0, 3800.0],
                "intensity": [0.1, 0.2, 1.0, 0.3],
            }
        ).to_csv(csv_path, index=False)
        meta_path.write_text(json.dumps({"jcamp_metadata": {"YUNITS": "ABSORBANCE", "STATE": "gas"}}), encoding="utf-8")
        reference_manifest = {
            "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            "canonical_smiles": "O",
            "nist_page_url": "https://webbook.nist.gov/cgi/inchi?ID=C7732185",
            "reference_spectra": [
                {
                    "csv": str(csv_path),
                    "meta_json": str(meta_path),
                    "jcamp_url": "https://webbook.nist.gov/cgi/inchi?JCAMP=C7732185&Index=0&Type=IR",
                    "nist_id": "C7732185",
                    "index": "0",
                    "phase_tag": "gas",
                    "phase_label": "gas",
                    "selection_priority": 100,
                    "description": "gas",
                }
            ],
            "preferred_reference": {},
        }
        (nist_dir / "XLYOFNOQVPJJNP-UHFFFAOYSA-N_reference_set.json").write_text(
            json.dumps(reference_manifest),
            encoding="utf-8",
        )
        return reference_manifest["reference_spectra"]

    result = import_hess_files_for_web(
        [ROOT / "data" / "hess" / "H2O_freq.hess"],
        import_root=outroot,
        run_id="nist_attach",
        pipeline_runner=fake_pipeline,
        include_nist_ir=True,
        nist_ir_runner=fake_nist_runner,
    )

    spectrum_json = Path(result.artifacts["interactive_spectrum_data_json"])
    payload = json.loads(spectrum_json.read_text(encoding="utf-8"))
    ref_sets = payload["nist_reference_sets"]
    assert "H2O" in ref_sets
    ref_item = ref_sets["H2O"]["reference_spectra"][0]
    assert ref_item["suitable_for_matching"] is True
    assert "scale_engine_payload" in ref_item
    assert Path(result.artifacts["nist_reference_set_H2O"]).exists()
    assert "nist_ir_attached:H2O_freq.hess:references=1" in result.diagnostics
