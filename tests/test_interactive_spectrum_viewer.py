from __future__ import annotations

import shutil
import sys
import subprocess
from pathlib import Path

import pandas as pd
import json
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orca_parser import read_orca_hess  # noqa: E402
from reports import build_ped_driven_final_assignment_table, build_ped_stage3d_agreement_table, build_spectrum_payload, classify_composed_ped_diagnostic_policy, classify_composed_ped_evidence_origin, classify_composed_ped_warning, triage_composed_ped_diagnostic_hint, write_interactive_spectrum_viewer  # noqa: E402
from reports import attach_nist_reference_set, classify_reference_suitability  # noqa: E402


CHROME_EXE = Path("C:/Program Files/Google/Chrome/Application/chrome.exe")


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
    wilson_ped = pd.DataFrame(
        [
            {
                "Filename": hess.filename,
                "mode": positive_mode,
                "frequency_cm-1": float(hess.frequencies_cm1[positive_mode]),
                "wilson_rank": 1,
                "coordinate_family": "H-O-H bend",
                "internal_coordinate": "ang(H2-O1-H3)",
                "coordinate_class": "bend",
                "contribution_percent": 99.9,
                "wilson_ped_method": "Wilson GF-style PED audit test method",
            }
        ]
    )
    composed_wilson_ped = pd.DataFrame(
        [
            {
                "Filename": hess.filename,
                "mode": positive_mode,
                "frequency_cm-1": float(hess.frequencies_cm1[positive_mode]),
                "wilson_rank": 1,
                "coordinate_family": "composed O-H stretch",
                "internal_coordinate": "composed_symmetric_XH_stretch(O1:H2,H3)",
                "coordinate_class": "stretch",
                "contribution_percent": 99.8,
                "source": "composed_coordinate",
                "wilson_ped_method": "Composed Wilson GF-style PED audit test method",
            }
        ]
    )

    payload = build_spectrum_payload(
        [hess],
        assignment_audit,
        wilson_ped_audit=wilson_ped,
        composed_wilson_ped_audit=composed_wilson_ped,
    )
    target_mode = next(mode for mode in payload["files"][0]["modes"] if mode["mode"] == positive_mode)
    assert target_mode["assignment"] == "O-H stretch"
    assert target_mode["final_assignment"] == "O-H stretch"
    assert target_mode["final_assignment_source"] == "ORCAVEDA assignment audit"
    assert target_mode["final_assignment_policy"] == "stage3d_fallback_due_to_ped_disagreement"
    assert target_mode["stage3d_assignment"] == "O-H stretch"
    assert target_mode["ped_source"] == "Wilson GF-style PED audit"
    assert target_mode["ped_top_percent"] == 99.9
    assert target_mode["composed_ped_source"] == "Composed-coordinate Wilson GF-style PED audit"
    assert target_mode["composed_ped_assignment"] == "composed O-H stretch"
    assert target_mode["composed_ped_top_percent"] == 99.8
    assert "composed_symmetric_XH_stretch" in target_mode["composed_ped_top_contributors"]
    assert target_mode["composed_ped_policy_hint"] == "diagnostic_hint_composed_differs_from_baseline"
    assert target_mode["composed_ped_triage_category"] == "baseline_preferred_composed_lower_localization"
    assert target_mode["composed_ped_evidence_origin"] == "composed_coordinate_top"
    assert target_mode["composed_ped_warning"] == ""
    assert target_mode["composed_ped_warning_reason"] == ""
    assert target_mode["composed_ped_localization_delta_percent"] == -0.1
    assert target_mode["composed_ped_semantic_status"] == "FAIL"
    assert target_mode["composed_ped_semantic_reason"] == "motion_family_mismatch"
    assert target_mode["ped_agreement_status"] == "disagrees"
    assert "ped_stage3d_semantic_disagreement" in target_mode["ped_policy_warning"]
    fallback_mode = next(mode for mode in payload["files"][0]["modes"] if mode["mode"] != positive_mode)
    assert fallback_mode["ped_source"] == ""
    assert fallback_mode["composed_ped_source"] == ""
    assert fallback_mode["stage3d_assignment"] == ""
    assert fallback_mode["ped_agreement_status"] == "not_available"

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
    assert 'data-section="summary-first"' in html_text
    assert "3D Molecule Viewer" in html_text
    assert "Summary" in html_text
    assert "Evidence" in html_text
    assert "NIST / Scaling" in html_text
    assert "Raw diagnostics" in html_text
    assert '<details id="advancedDiagnostics" class="advanced-diagnostics">' in html_text
    assert "Final Assignment" in html_text
    assert "Final Assignment Policy" in html_text
    assert "PED Contribution" in html_text
    assert "PED Diagnostic Interpretation" in html_text
    assert "PED Agreement Status" in html_text
    assert "PED Policy Warning" in html_text
    assert "ORCAVEDA Assignment" in html_text
    assert "PED Contributors" in html_text
    assert "Evidence Layer" in html_text
    assert "Composed Hint" in html_text
    assert "composedHintFilter" in html_text
    assert "All modes" in html_text
    assert "Composed hints" in html_text
    assert "Better localization" in html_text
    assert "Differs from baseline" in html_text
    assert "No modes match the selected composed hint filter." in html_text
    assert 'appendEmptyRow(peakTable, 6' in html_text
    assert "Selected Evidence Interpretation" in html_text
    assert "Composed PED Interpretation" in html_text
    assert "Composed PED Contributors" in html_text
    assert "Composed PED Policy Hint" in html_text
    assert "Composed PED Triage" in html_text
    assert "Composed PED Evidence Origin" in html_text
    assert "Composed PED Warning" in html_text
    assert "Composed PED Localization Delta" in html_text
    assert "Composed PED Semantic Status" in html_text
    assert "moleculeViewer" in html_text
    assert "3Dmol-min.js" in html_text
    assert "molStyle" in html_text

    json_text = json_path.read_text(encoding="utf-8")
    assert "frequency_cm1" in json_text
    assert hess.filename in json_text
    assert "final_assignment" in json_text
    assert "stage3d_fallback_due_to_ped_disagreement" in json_text
    assert "stage3d_assignment" in json_text
    assert "ped_agreement_status" in json_text
    assert "ped_policy_warning" in json_text
    assert "composed_ped_assignment" in json_text
    assert "composed_ped_policy_hint" in json_text
    assert "composed_ped_triage_category" in json_text
    assert "composed_ped_evidence_origin" in json_text
    assert "composed_ped_warning" in json_text
    assert "composed_ped_localization_delta_percent" in json_text
    assert "composed_ped_semantic_status" in json_text
    assert "Composed-coordinate Wilson GF-style PED audit" in json_text
    assert "composed_symmetric_XH_stretch" in json_text
    assert "O-H stretch" in json_text
    assert "\"geometry\"" in json_text
    assert "\"atoms\"" in json_text
    assert "\"bonds\"" in json_text


def test_composed_ped_diagnostic_policy_is_viewer_only_and_conservative():
    hint, delta, status, reason = classify_composed_ped_diagnostic_policy(
        "C=O stretch",
        45.9,
        "C=O stretch",
        63.0,
    )
    assert hint == "composed_confirms_with_better_localization"
    assert round(delta, 1) == 17.1
    assert status == "PASS"
    assert reason == "baseline_ped_semantic_overlap"

    hint, _delta, status, reason = classify_composed_ped_diagnostic_policy(
        "H-O-H bend",
        99.9,
        "O-H stretch",
        99.8,
    )
    assert hint == "diagnostic_hint_composed_differs_from_baseline"
    assert status == "FAIL"
    assert reason == "motion_family_mismatch"


def test_composed_ped_hint_triage_separates_noise_from_review_targets():
    category, recommendation = triage_composed_ped_diagnostic_hint(
        "diagnostic_hint_composed_differs_from_baseline",
        "motion_family_mismatch",
        -4.1,
        3196.0,
    )
    assert category == "baseline_preferred_composed_lower_localization"
    assert recommendation == "do_not_promote_composed_evidence"

    category, recommendation = triage_composed_ped_diagnostic_hint(
        "diagnostic_hint_composed_differs_from_baseline",
        "motion_family_mismatch",
        53.9,
        3073.0,
        "C-H torsion mixed with CH2 scissor",
        "C-H stretch",
    )
    assert category == "high_frequency_xh_stretch_recovery"
    assert recommendation == "keep_composed_xh_stretch_as_diagnostic_evidence"

    category, recommendation = triage_composed_ped_diagnostic_hint(
        "diagnostic_hint_composed_differs_from_baseline",
        "motion_family_mismatch",
        53.9,
        3073.0,
        "C-C-C bend",
        "C-C stretch",
    )
    assert category == "high_frequency_motion_family_review"
    assert recommendation == "inspect_xh_stretch_vs_bend_or_torsion_coordinate_generation"

    category, recommendation = triage_composed_ped_diagnostic_hint(
        "diagnostic_hint_composed_differs_from_baseline",
        "motion_family_mismatch",
        28.5,
        289.7,
        "C-N-C bend",
        "C-H torsion",
        "primitive",
        False,
    )
    assert category == "primitive_row_optimizer_substitution"
    assert recommendation == "inspect_optimizer_substitution_before_coordinate_generation"

    category, recommendation = triage_composed_ped_diagnostic_hint(
        "composed_confirms_with_better_localization",
        "baseline_ped_semantic_overlap",
        17.1,
        1683.0,
    )
    assert category == "confirmation_candidate"
    assert recommendation == "keep_as_diagnostic_confirmation"


def test_composed_ped_evidence_origin_uses_top_provenance():
    assert (
        classify_composed_ped_evidence_origin(
            "C-H stretch",
            "composed_coordinate",
            True,
        )
        == "composed_coordinate_top"
    )
    assert (
        classify_composed_ped_evidence_origin(
            "C-H torsion",
            "primitive",
            False,
        )
        == "primitive_substitution_top"
    )
    assert (
        classify_composed_ped_evidence_origin(
            "",
            "",
            False,
        )
        == "baseline_or_no_composed_top"
    )


def test_composed_ped_warning_requires_origin_and_triage():
    assert classify_composed_ped_warning(
        "primitive_substitution_top",
        "primitive_row_optimizer_substitution",
    ) == (
        "primitive_row_optimizer_substitution_warning",
        "composed_ped_top_contributor_is_primitive_after_optimizer_substitution",
    )
    assert classify_composed_ped_warning(
        "primitive_substitution_top",
        "viewer_evidence_only",
    ) == ("", "")
    assert classify_composed_ped_warning(
        "composed_coordinate_top",
        "primitive_row_optimizer_substitution",
    ) == ("", "")


def test_ped_stage3d_agreement_table_policy_statuses():
    assignment_audit = pd.DataFrame(
        [
            {
                "Source": "[1]",
                "Filename": "synthetic.hess",
                "mode": 1,
                "frequency_cm-1": 1700.0,
                "functional_group_assignment": "C=O stretch",
            },
            {
                "Source": "[1]",
                "Filename": "synthetic.hess",
                "mode": 2,
                "frequency_cm-1": 1000.0,
                "functional_group_assignment": "O-H stretch",
            },
            {
                "Source": "[1]",
                "Filename": "synthetic.hess",
                "mode": 3,
                "frequency_cm-1": 900.0,
                "functional_group_assignment": "C-C-H bend",
            },
        ]
    )
    wilson_ped = pd.DataFrame(
        [
            {
                "Filename": "synthetic.hess",
                "mode": 1,
                "frequency_cm-1": 1700.0,
                "wilson_rank": 1,
                "coordinate_family": "C=O stretch",
                "internal_coordinate": "carbonyl_CO_stretch(C1=O2)",
                "contribution_percent": 72.0,
            },
            {
                "Filename": "synthetic.hess",
                "mode": 2,
                "frequency_cm-1": 1000.0,
                "wilson_rank": 1,
                "coordinate_family": "C-C bend",
                "internal_coordinate": "ang(C1-C2-C3)",
                "contribution_percent": 60.0,
            },
            {
                "Filename": "synthetic.hess",
                "mode": 3,
                "frequency_cm-1": 900.0,
                "wilson_rank": 1,
                "coordinate_family": "C-C-H bend",
                "internal_coordinate": "ang(C1-C2-H3)",
                "contribution_percent": 18.0,
            },
        ]
    )

    agreement = build_ped_stage3d_agreement_table(assignment_audit, wilson_ped_audit=wilson_ped)
    final = build_ped_driven_final_assignment_table(agreement)

    by_mode = {int(row["mode"]): row for _, row in agreement.iterrows()}
    assert by_mode[1]["ped_agreement_status"] == "confirms"
    assert by_mode[2]["ped_agreement_status"] == "disagrees"
    assert "ped_stage3d_semantic_disagreement" in by_mode[2]["ped_policy_warning"]
    assert by_mode[3]["ped_agreement_status"] == "adds_context"
    assert "diffuse_ped_contributions" in by_mode[3]["ped_policy_warning"]
    final_by_mode = {int(row["mode"]): row for _, row in final.iterrows()}
    assert final_by_mode[1]["final_assignment"] == "C=O stretch"
    assert final_by_mode[1]["final_assignment_policy"] == "ped_confirms_stage3d"
    assert final_by_mode[2]["final_assignment"] == "O-H stretch"
    assert final_by_mode[2]["final_assignment_policy"] == "stage3d_fallback_due_to_ped_disagreement"
    assert final_by_mode[3]["final_assignment"] == "C-C-H bend"
    assert final_by_mode[3]["final_assignment_policy"] == "ped_adds_context"


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
                    "nist_page_url": "https://webbook.nist.gov/cgi/inchi?ID=C98862",
                    "inchikey": "TEST",
                    "canonical_smiles": "CC",
                    "reference_spectra": [
                        {
                            "csv": str(ref_csv),
                            "meta_json": str(ref_meta),
                            "nist_id": "C98862",
                            "jcamp_url": "https://webbook.nist.gov/cgi/inchi?JCAMP=C98862&Index=0&Type=IR",
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
    assert "matching_layers" in ref_payload
    assert "high_confidence" in ref_payload["matching_layers"]
    assert "extended" in ref_payload["matching_layers"]
    assert "global_ls" in ref_payload["engine_fits"]
    assert payload["nist_reference_sets"]["acetophenone"]["reference_spectra"][0]["nist_spectrum_url"].startswith("https://webbook.nist.gov/")

    html_path = outdir / "viewer.html"
    json_path = outdir / "viewer.json"
    write_interactive_spectrum_viewer(payload, html_path, json_path=json_path)

    html_text = html_path.read_text(encoding="utf-8")
    assert "NIST Reference" in html_text
    assert "Open on NIST" in html_text
    assert "nistReferenceLink" in html_text
    assert "Index ${item.index} - ${phase}${units}${description}${suitability}" in html_text
    assert "Scale Engine" in html_text
    assert "scaleEngine" in html_text
    assert "Matching Layer" in html_text
    assert "matchingLayer" in html_text
    assert "engineTableBody" in html_text
    assert "matchingLayerSummary" in html_text
    assert "matchingLayerTableBody" in html_text
    assert "engineLayerMatrixBody" in html_text
    assert "Active layer:" in html_text
    assert "Coverage" in html_text
    assert "Nearest %Δ" in html_text
    assert "High-Conf %Δ" in html_text
    assert "Extended %Δ" in html_text
    assert "Mean %Δ" in html_text
    assert "RMS %Δ" in html_text
    assert "Max %Δ" in html_text
    assert "nistReference" in html_text
    assert "Auto-fit scale" in html_text
    assert "fitSummary" in html_text
    json_text = json_path.read_text(encoding="utf-8")
    assert "reference_spectra" in json_text
    assert "nist_spectrum_url" in json_text
    assert "scale_engine_payload" in json_text
    assert "matching_layers" in json_text
    assert "matching_layer_overview" in json_text
    assert "engine_layer_matrix" in json_text
    assert "high_confidence_matched_pairs" in json_text
    assert "nearest_matched_pairs" in json_text
    assert "global_weighted_ls" in json_text
    assert "mean_percent_deviation" in json_text


def test_unsuitable_reference_is_kept_but_skipped_for_matching():
    assert classify_reference_suitability(y_units="ABSORPTION INDEX", phase_tag="liquid", description="liquid") == (
        False,
        "absorption_index_reference",
    )

    payload = {
        "viewer_title": "ORCAVEDA Interactive IR Spectrum",
        "default_scale_factor": 1.0,
        "default_lorentz_hwhm": 12.0,
        "files": [
            {
                "filename": "MeOH_freq.hess",
                "title": "MeOH",
                "summary": {},
                "modes": [{"mode": 1, "frequency_cm1": 1030.0, "intensity": 100.0, "assignment": "C-O stretch", "top_internal_coordinates": "", "warnings": ""}],
                "geometry": {"atoms": [], "bonds": []},
            }
        ],
    }
    outdir = ROOT / "outputs" / "pytest_interactive_spectrum_viewer_unsuitable"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_csv = outdir / "ref.csv"
    pd.DataFrame({"wavenumber_cm-1": [1000.0, 1010.0, 1020.0], "intensity": [1.0, 2.0, 3.0]}).to_csv(ref_csv, index=False)
    ref_meta = outdir / "ref_meta.json"
    ref_meta.write_text(json.dumps({"jcamp_metadata": {"YUNITS": "ABSORPTION INDEX", "STATE": "liquid"}}), encoding="utf-8")
    manifest_path = outdir / "refset.json"
    manifest_path.write_text(
        json.dumps(
            {
                "inchikey": "TEST",
                "canonical_smiles": "CO",
                "reference_spectra": [
                    {
                        "csv": str(ref_csv),
                        "meta_json": str(ref_meta),
                        "index": "28",
                        "phase_tag": "liquid",
                        "phase_label": "liquid",
                        "selection_priority": 10,
                        "description": "absorption index",
                    }
                ],
                "preferred_reference": {"index": "28", "phase_tag": "liquid", "phase_label": "liquid"},
            }
        ),
        encoding="utf-8",
    )
    payload = attach_nist_reference_set(payload, manifest_path)
    ref_item = payload["nist_reference_sets"]["MeOH"]["reference_spectra"][0]
    assert ref_item["suitable_for_matching"] is False
    assert ref_item["suitability_reason"] == "absorption_index_reference"
    assert ref_item["scale_engine_payload"] is None
    assert ref_item["matching_status"] == "skipped_unsuitable_reference"


def test_interactive_spectrum_viewer_can_inline_local_3dmol_asset():
    payload = {
        "viewer_title": "ORCAVEDA Interactive IR Spectrum",
        "default_scale_factor": 1.0,
        "default_lorentz_hwhm": 12.0,
        "files": [
            {
                "filename": "H2O_freq.hess",
                "title": "H2O",
                "summary": {},
                "modes": [{"mode": 1, "frequency_cm1": 3700.0, "intensity": 1.0, "assignment": "O-H stretch"}],
                "geometry": {"atoms": [], "bonds": []},
            }
        ],
    }
    outdir = ROOT / "outputs" / "pytest_interactive_spectrum_viewer_local_3dmol"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    asset_path = outdir / "3Dmol-min.js"
    asset_path.write_text(
        "window.$3Dmol={createViewer:function(){return {clear:function(){},addModel:function(){},setStyle:function(){},zoomTo:function(){},resize:function(){},render:function(){}};}}; // </script safety",
        encoding="utf-8",
    )
    html_path = outdir / "viewer.html"

    write_interactive_spectrum_viewer(payload, html_path, three_dmol_js_path=asset_path)
    html_text = html_path.read_text(encoding="utf-8")

    assert "cdn.jsdelivr.net/npm/3dmol" not in html_text
    assert "Inlined local 3Dmol.js asset: 3Dmol-min.js" in html_text
    assert "<\\/script safety" in html_text


def test_interactive_spectrum_viewer_escapes_payload_script_breakout_text():
    payload = {
        "viewer_title": "ORCAVEDA Interactive IR Spectrum",
        "default_scale_factor": 1.0,
        "default_lorentz_hwhm": 12.0,
        "files": [
            {
                "filename": "synthetic.hess",
                "title": "synthetic",
                "summary": {},
                "modes": [
                    {
                        "mode": 1,
                        "frequency_cm1": 1200.0,
                        "intensity": 1.0,
                        "assignment": '</script><div id="breakout">bad</div>',
                    }
                ],
                "geometry": {"atoms": [], "bonds": []},
            }
        ],
    }
    outdir = ROOT / "outputs" / "pytest_interactive_spectrum_viewer_script_escape"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / "viewer.html"

    write_interactive_spectrum_viewer(payload, html_path)
    html_text = html_path.read_text(encoding="utf-8")

    assert '</script><div id="breakout">bad</div>' not in html_text
    assert "\\u003c/script\\u003e" in html_text
    assert "\\u003cdiv id=" in html_text
    assert "\\u003c/div\\u003e" in html_text


def test_interactive_spectrum_viewer_headless_chrome_smoke_without_cdn():
    if not CHROME_EXE.exists():
        pytest.skip("Chrome executable is not available for viewer smoke test.")

    hess = read_orca_hess(ROOT / "data" / "hess" / "H2O_freq.hess")
    positive_modes = [idx for idx, freq in enumerate(hess.frequencies_cm1) if float(freq) > 0.0]
    assignment_audit = pd.DataFrame(
        [
            {
                "Filename": hess.filename,
                "mode": mode,
                "frequency_cm-1": float(hess.frequencies_cm1[mode]),
                "IR_intensity": float(hess.ir_intensities[mode]),
                "functional_group_assignment": "O-H stretch",
                "top_internal_coordinates": "r(O1-H2)=50.0%; r(O1-H3)=50.0%",
                "warnings": "",
            }
            for mode in positive_modes
        ]
    )
    payload = build_spectrum_payload([hess], assignment_audit)
    outdir = ROOT / "outputs" / "pytest_interactive_spectrum_viewer_chrome"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / "viewer.html"
    screenshot_path = outdir / "viewer.png"
    write_interactive_spectrum_viewer(payload, html_path)

    chrome_base = [
        str(CHROME_EXE),
        "--headless=new",
        "--disable-gpu",
        "--disable-crash-reporter",
        "--disable-crashpad",
        "--no-first-run",
        "--allow-file-access-from-files",
        "--host-resolver-rules=MAP cdn.jsdelivr.net 0.0.0.0",
        "--window-size=1400,1100",
        "--virtual-time-budget=3500",
    ]
    dom_result = subprocess.run(
        [*chrome_base, "--dump-dom", html_path.as_uri()],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=30,
    )
    if dom_result.returncode != 0 and "crashpad" in dom_result.stderr.lower() and "0x5" in dom_result.stderr:
        pytest.skip(f"Headless Chrome is blocked by local crashpad permissions: {dom_result.stderr.strip()}")
    if dom_result.returncode != 0 and "gpu process isn't usable" in dom_result.stderr.lower():
        pytest.skip(f"Headless Chrome GPU process is unavailable in this local environment: {dom_result.stderr.strip()}")
    assert dom_result.returncode == 0, dom_result.stderr
    dom = dom_result.stdout
    assert 'class="panel panel-spectrum"' in dom
    assert 'data-section="summary-first"' in dom
    assert 'id="advancedDiagnostics" class="advanced-diagnostics"' in dom
    assert "Final Assignment" in dom
    assert 'id="summaryGrid"' in dom
    assert 'id="moleculeViewer"' in dom
    assert 'id="peakTable"' in dom
    assert 'data-renderer="native-fallback"' in dom
    assert "Native 2D projection used because 3Dmol.js is unavailable." in dom
    assert "O-H stretch" in dom

    screenshot_result = subprocess.run(
        [*chrome_base, f"--screenshot={screenshot_path}", html_path.as_uri()],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=30,
    )
    if screenshot_result.returncode != 0 and "gpu process isn't usable" in screenshot_result.stderr.lower():
        pytest.skip(f"Headless Chrome GPU process is unavailable in this local environment: {screenshot_result.stderr.strip()}")
    assert screenshot_result.returncode == 0, screenshot_result.stderr
    assert screenshot_path.exists()
    image = Image.open(screenshot_path).convert("RGB")
    colors = image.getcolors(maxcolors=1_000_000)
    assert colors is not None and len(colors) > 100
