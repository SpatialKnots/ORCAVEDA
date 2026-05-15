# Add VEDA Reference Validation Harness

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

GAP 3 adds a reproducible harness for comparing ORCAVEDA opt-in `veda_like_*` artifacts with original VEDA reference outputs when those outputs are checked in. The immediate user-visible behavior is conservative: if no VEDA references exist, the tool reports `SKIP`, not `PASS`.

## Progress

- [x] (2026-05-15) Initial GAP 3 plan written.
- [x] (2026-05-15) Added skip-safe comparison harness under `benchmarks/veda_compare`.
- [x] (2026-05-15) Added focused synthetic tests for missing-reference SKIP, matching PASS, percent-delta FAIL, and dominant-contributor FAIL.
- [x] (2026-05-15) Focused tests and missing-reference CLI probe completed.
- [x] (2026-05-15) Documented the reference fixture contract and status semantics.

## Surprises & Discoveries

- Observation: The broader VEDA-like ExecPlan already named a future comparison harness but left it unimplemented.
  Evidence: `docs/full_veda_implementation_execplan.md`.

- Observation: No checked-in original VEDA reference outputs were found during this patch.
  Evidence: repository file listing and absence of `veda_reference_ped_matrix.csv`.

## Decision Log

- Decision: Missing VEDA references produce `comparison_status=SKIP` and `acceptance_status=SKIP`.
  Rationale: Absence of reference data is not evidence of agreement.
  Date/Author: 2026-05-15 / Codex

- Decision: Use a simple CSV reference schema for the first harness.
  Rationale: The smallest useful comparison is matrix rows plus optional dominant contributors; richer VEDA exports can be adapted later without claiming equivalence now.
  Date/Author: 2026-05-15 / Codex

- Decision: Match matrix rows by exact `Filename`, `mode`, and `internal_coordinate`.
  Rationale: This makes the first harness deterministic and prevents fuzzy relabeling from hiding coordinate-basis disagreements.
  Date/Author: 2026-05-15 / Codex

## Outcomes & Retrospective

Validation outcome:

- `.\.venv312\Scripts\python.exe -m py_compile benchmarks\veda_compare\compare_veda_outputs.py tests\test_veda_reference_compare.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_veda_reference_compare.py tests\test_validate_veda_like_outputs.py -q` returned `9 passed`.
- `.\.venv312\Scripts\python.exe benchmarks\veda_compare\compare_veda_outputs.py --orcaveda outputs\pytest_veda_like_epm_opt_in_h2o --reference data\veda_reference_missing --out outputs\veda_reference_compare_skip_probe` returned `comparison_status=SKIP`, `acceptance_status=SKIP`, and `reason=veda_reference_directory_missing`.

## Context and Orientation

Current ORCAVEDA VEDA-like outputs are generated only when `--veda-like-ped` is enabled. The relevant artifacts are `*__veda_like_ped_matrix.csv` and `*__veda_like_ped_audit.csv`.

The reference harness lives in `benchmarks/veda_compare/compare_veda_outputs.py`. It reads ORCAVEDA artifacts from an output directory and original VEDA reference CSVs from a separate reference directory.

## Plan of Work

First, keep the harness inactive by default and make absence of reference outputs explicit. Second, compare matrix rows by `Filename`, `mode`, and `internal_coordinate`, with an absolute PED percent tolerance. Third, optionally compare rank-1 dominant internal coordinates when `veda_reference_dominant_assignments.csv` exists. Fourth, write JSON and CSV diagnostics for review.

## Concrete Steps

Run focused tests:

    .\.venv312\Scripts\python.exe -m pytest tests\test_veda_reference_compare.py -q

Run the CLI against any existing ORCAVEDA `veda_like_*` output directory and an absent reference directory:

    .\.venv312\Scripts\python.exe benchmarks\veda_compare\compare_veda_outputs.py --orcaveda outputs\veda_like_full_sweep_live --reference data\veda_reference_missing --out outputs\veda_like_reference_comparison_live

Expected result without references: JSON summary with `comparison_status` set to `SKIP`.

## Validation and Acceptance

Accepted when missing references skip clearly, synthetic matching rows pass, synthetic out-of-tolerance PED percentages fail, synthetic dominant-coordinate mismatches fail, and no default ORCAVEDA output schemas are changed.

## Idempotence and Recovery

The harness is read-only for ORCAVEDA outputs and VEDA references. It writes only to the requested comparison output directory. Re-running the command overwrites comparison summaries.

## Artifacts and Notes

Added artifacts:

- `benchmarks/veda_compare/compare_veda_outputs.py`
- `benchmarks/veda_compare/README.md`
- `benchmarks/veda_compare/examples/missing_reference_probe/veda_reference_comparison_summary.json`
- `tests/test_veda_reference_compare.py`

## Interfaces and Dependencies

Reference matrix file:

- `veda_reference_ped_matrix.csv`

Required columns:

- `Filename`
- `mode`
- `internal_coordinate`
- `contribution_percent`

Optional dominant reference file:

- `veda_reference_dominant_assignments.csv`

Required columns:

- `Filename`
- `mode`
- `internal_coordinate`
