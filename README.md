# ORCAVEDA Stage 3D v5.0 Regression Test Structure

## Purpose

Stage 3D v5.0 formalizes regression testing as a project component.

Chemistry / assignment logic is inherited from Stage 3D v4.9. The v5.0 change is organizational:

- machine-readable expectations
- reusable regression runner
- pytest-compatible wrapper
- CSV and JSON test outputs
- explicit PASS/FAIL accounting

## Scope

The current regression set contains:

- Acetone
- CH3CN
- DMF
- DMSO
- EtOH
- MeOH
- NMP
- iPrOH

The harness checks:

1. expected molecules present
2. `rank_B_independent == expected_rank_3N_minus_6`
3. normal-mode orientation rule reported
4. expected functional groups detected
5. mid-frequency diagnostic assignments detected
6. forbidden CH2 labels absent where applicable
7. no high-frequency unassigned modes
8. monoethanolamine sanity checks skipped for external molecules
9. protected X-H fallback cleanup is preserved
10. expected DMF warning is tolerated

## Run

From this directory:

```bash
python run_regression_tests.py \
  --outdir /path/to/orcaveda_outputs \
  --expectations expectations/regression_expectations_stage3D_v5_0.json
```

The runner writes:

- `regression_harness_results_stage3D_v5_0.csv`
- `regression_harness_summary_stage3D_v5_0.json`

## Fresh Setup

From a clean checkout on Windows PowerShell:

```powershell
py -3.12 -m venv .venv312
.\.venv312\Scripts\python.exe -m pip install --upgrade pip
.\.venv312\Scripts\python.exe -m pip install -r requirements.txt
.\.venv312\Scripts\python.exe -m pytest -q
```

The root pytest command is configured by `pytest.ini` to collect only `tests/`.

To generate fresh Stage 3D outputs and validate them against the checked-in
expectations:

```powershell
.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py `
  data\hess\Acetone_freq.hess data\hess\CH3CN_freq.hess data\hess\DMF_freq.hess data\hess\DMSO_freq.hess `
  data\hess\EtOH_freq.hess data\hess\MeOH_freq.hess data\hess\NMP_freq.hess data\hess\iPrOH_freq.hess `
  --outdir outputs\regression_live
.\.venv312\Scripts\python.exe run_regression_tests.py --outdir outputs\regression_live --expectations expectations\regression_expectations_stage3D_v5_0.json
```

For the local web import UI:

```powershell
.\.venv312\Scripts\python.exe src\web_app.py
```

Then open `http://127.0.0.1:8765/` and upload one or more ORCA `.hess` files.

## Offline Viewer Assets

Generated interactive HTML reports use the public 3Dmol.js CDN by default and
fall back to a native 2D molecule projection if 3Dmol.js cannot load. For fully
self-contained offline reports, call `write_interactive_spectrum_viewer(...,
three_dmol_js_path="path\\to\\3Dmol-min.js")`; the local asset is inlined into
the generated HTML.

## Epistemic Status

This remains a PED-like / assignment-audit regression harness. It is not a strict VEDA-equivalent PED validation suite.
