# ORCAVEDA Progress Checkpoint — Stage 3D v5.0

## Current Baseline

`ORCAVEDA_patched_stage3D_v5_0.py`

## Status

Stage 3D v5.0 is a regression-infrastructure baseline.

Chemistry and assignment behavior are inherited from Stage 3D v4.9. The v5.0 change is formalization of regression testing.

## What v5.0 Adds

- `expectations/regression_expectations_stage3D_v5_0.json`
- `run_regression_tests.py`
- `tests/test_stage3d_outputs.py`
- regression result CSV
- regression summary JSON
- project README for test execution

## Regression Set

- Acetone
- CH3CN
- DMF
- DMSO
- EtOH
- MeOH
- NMP
- iPrOH

## Checks

- expected molecule outputs present
- independent B-rank equals expected vibrational rank
- normal-mode orientation rule is reported
- expected functional groups detected
- mid-frequency diagnostic assignments present
- forbidden CH2 labels absent where no methylene is expected
- no high-frequency unassigned modes
- monoethanolamine sanity checks skipped for external molecules
- protected X-H unused norm/power cleanup preserved
- DMF negative-frequency warning tolerated as expected input-quality warning

## Result

Regression harness status: PASS

## Remaining Warning

DMF retains an expected input-quality warning: negative vibrational frequencies after the first six modes.

## Method Boundary

Correct label:

`PED-like / assignment-audit layer`

Do not claim:

- strict VEDA-equivalent PED
- full Wilson GF formalism
- publication-grade universal benchmark validation
