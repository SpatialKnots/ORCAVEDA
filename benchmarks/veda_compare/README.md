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

## Reference Ingest

Use `convert_veda_reference.py` to normalize checked-in VEDA reference CSV rows
before comparison. This tool accepts either an already-normalized
`veda_reference_ped_matrix.csv` file or an explicit CSV column mapping. It does
not parse unknown native VEDA formats by inference.

Normalized pass-through:

    .\.venv312\Scripts\python.exe benchmarks\veda_compare\convert_veda_reference.py --raw-reference <raw-reference-dir> --out benchmarks\veda_compare\references\<set-name>

Explicit CSV mapping:

    .\.venv312\Scripts\python.exe benchmarks\veda_compare\convert_veda_reference.py --raw-reference <raw-reference-dir> --matrix-csv <raw-reference-dir>\matrix.csv --filename-column file --mode-column mode_no --coordinate-column veda_coord --percent-column ped_percent --out benchmarks\veda_compare\references\<set-name>

If no convertible reference matrix exists, the converter writes
`veda_reference_ingest_summary.json` with `conversion_status=SKIP`. If required
columns or numeric values are invalid, it writes `conversion_status=FAIL`.
