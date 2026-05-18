# VEDA Reference Comparison Harness

This directory contains the inactive-by-default comparison harness for GAP 3.
It compares ORCAVEDA opt-in `veda_like_*` artifacts against checked-in original
VEDA reference CSV rows when those references are available.

The harness must not be used to claim original VEDA reproduction unless real
reference outputs are present and the comparison summary reports `PASS`.

## Reference Schema

Place reference CSV files in a separate checked-in reference directory, for
example `benchmarks/veda_compare/references/<molecule_or_set>/`. Do not put
generated ORCAVEDA outputs in that directory.

The required matrix file is:

- `veda_reference_ped_matrix.csv`

Required columns:

- `Filename`
- `mode`
- `internal_coordinate`
- `contribution_percent`

Rows are matched to ORCAVEDA `*__veda_like_ped_matrix.csv` by the exact key:

- `Filename`
- `mode`
- `internal_coordinate`

`contribution_percent` is compared with an absolute percent tolerance. The
default tolerance is 5.0 percentage points because this first harness is a
diagnostic gate, not a bit-for-bit reproduction claim.

Optional dominant-contributor file:

- `veda_reference_dominant_assignments.csv`

Required columns when present:

- `Filename`
- `mode`
- `internal_coordinate`

Rows are matched to ORCAVEDA rank-1 `*__veda_like_ped_audit.csv` rows by:

- `Filename`
- `mode`

The dominant `internal_coordinate` must match exactly.

## Status Semantics

- `PASS`: checked-in reference rows exist and all implemented comparisons are
  within the declared tolerance.
- `FAIL`: checked-in reference rows exist, but required ORCAVEDA artifacts are
  absent, reference rows are missing from ORCAVEDA outputs, PED percentages are
  outside tolerance, or dominant contributors mismatch.
- `SKIP`: original VEDA reference rows are absent. This is not evidence of
  agreement and must not be reported as validation success.

## Command

From the repository root:

    .\.venv312\Scripts\python.exe benchmarks\veda_compare\compare_veda_outputs.py --orcaveda outputs\veda_like_full_sweep_live --reference <checked-in-veda-reference-dir> --out outputs\veda_like_reference_comparison_live

If the reference directory or required reference matrix is absent, the command
writes a `SKIP` summary. It does not fabricate reference rows or report `PASS`.
