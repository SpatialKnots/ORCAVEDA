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
