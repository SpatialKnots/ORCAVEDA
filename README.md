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

## Epistemic Status

This remains a PED-like / assignment-audit regression harness. It is not a strict VEDA-equivalent PED validation suite.
