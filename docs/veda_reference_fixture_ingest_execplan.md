# Add VEDA Reference Fixture Ingest

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA now has a skip-safe GAP 3 comparison harness for `veda_like_*` outputs, but it still needs a reproducible way to prepare original VEDA reference rows. This milestone adds a conservative ingest tool that normalizes checked-in reference CSV files into the comparison schema. The tool must not infer unknown native VEDA formats or fabricate reference rows.

## Progress

- [x] (2026-05-18) GAP 3 comparison harness merged to `main`.
- [x] (2026-05-18) Initial reference ingest plan written.
- [x] (2026-05-18) Added `benchmarks/veda_compare/convert_veda_reference.py`.
- [x] (2026-05-18) Added focused ingest tests for missing input, normalized pass-through, explicit column mapping, and missing-column failure.

## Surprises & Discoveries

- Observation: No checked-in original VEDA reference exports are available in this repository.
  Evidence: GAP 3 harness validation reports missing-reference `SKIP`; no reference matrix fixture is present outside synthetic tests.

- Observation: The comparison harness only needs two normalized files for the first validation gate.
  Evidence: `benchmarks/veda_compare/README.md` defines `veda_reference_ped_matrix.csv` and optional `veda_reference_dominant_assignments.csv`.

## Decision Log

- Decision: The first ingest tool accepts only already-normalized CSV files or explicit column mappings.
  Rationale: Inferring unknown native VEDA export formats would risk inventing mode numbers, coordinate labels, or PED percentages. Explicit mappings keep the evidence boundary clear.
  Date/Author: 2026-05-18 / Codex

- Decision: Missing raw reference directories or missing matrix CSV files produce `SKIP`, not `PASS`.
  Rationale: Absence of original VEDA data is not evidence of agreement.
  Date/Author: 2026-05-18 / Codex

- Decision: Invalid columns or nonnumeric mode/PED values produce `FAIL`.
  Rationale: Malformed reference rows must block comparison rather than silently producing misleading fixtures.
  Date/Author: 2026-05-18 / Codex

## Outcomes & Retrospective

Initial outcome: `convert_veda_reference.py` can write normalized comparison fixtures from a checked-in CSV with the required schema or from an explicitly mapped CSV. It writes `veda_reference_ingest_summary.json` with `PASS`, `SKIP`, or `FAIL` status. It does not parse native VEDA text or spreadsheet exports without explicit column mapping.

Validation outcome:

- `.\.venv312\Scripts\python.exe -m py_compile benchmarks\veda_compare\convert_veda_reference.py tests\test_veda_reference_ingest.py` completed successfully on 2026-05-18.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_veda_reference_ingest.py tests\test_veda_reference_compare.py -q` returned `8 passed` on 2026-05-18.
- `.\.venv312\Scripts\python.exe benchmarks\veda_compare\convert_veda_reference.py --raw-reference data\veda_reference_missing --out outputs\veda_reference_ingest_skip_probe` returned `conversion_status=SKIP`, `acceptance_status=SKIP`, and `reason=raw_reference_directory_missing` on 2026-05-18.

## Context and Orientation

The GAP 3 comparison harness lives in `benchmarks/veda_compare/compare_veda_outputs.py`. It compares ORCAVEDA `*__veda_like_ped_matrix.csv` and `*__veda_like_ped_audit.csv` outputs against original VEDA reference CSVs.

The ingest tool lives in `benchmarks/veda_compare/convert_veda_reference.py`. Its output directory can be used directly as the `--reference` directory for `compare_veda_outputs.py`.

## Plan of Work

First, add a conservative converter for the current reference schema. Second, add tests that prove missing references skip, normalized references pass, explicit mappings pass, and malformed references fail. Third, document the command and keep native VEDA parsing out of scope until actual exports are available.

## Concrete Steps

Run syntax checks:

    .\.venv312\Scripts\python.exe -m py_compile benchmarks\veda_compare\convert_veda_reference.py tests\test_veda_reference_ingest.py

Run focused tests:

    .\.venv312\Scripts\python.exe -m pytest tests\test_veda_reference_ingest.py tests\test_veda_reference_compare.py -q

Run a missing-reference CLI probe:

    .\.venv312\Scripts\python.exe benchmarks\veda_compare\convert_veda_reference.py --raw-reference data\veda_reference_missing --out outputs\veda_reference_ingest_skip_probe

Expected result: `conversion_status=SKIP`.

## Validation and Acceptance

Accepted when the converter:

- writes `SKIP` for absent raw reference directories;
- writes normalized `veda_reference_ped_matrix.csv` for a valid source;
- writes optional `veda_reference_dominant_assignments.csv` when supplied;
- supports explicit column mappings;
- writes `FAIL` for malformed required columns or nonnumeric values;
- does not modify ORCAVEDA scientific pipeline code or default output schemas.

## Idempotence and Recovery

The converter writes only to the requested output directory. Re-running it may overwrite generated normalized reference CSVs and the ingest summary. If a conversion fails, remove or fix only the source reference fixture or mapping; no ORCAVEDA outputs need to be regenerated.

## Artifacts and Notes

Added artifacts:

- `benchmarks/veda_compare/convert_veda_reference.py`
- `tests/test_veda_reference_ingest.py`
- `docs/veda_reference_fixture_ingest_execplan.md`

## Interfaces and Dependencies

Input:

- `--raw-reference <dir>`
- optional `--matrix-csv <csv>`
- optional `--dominant-csv <csv>`
- optional explicit column mapping flags

Output:

- `veda_reference_ped_matrix.csv`
- optional `veda_reference_dominant_assignments.csv`
- `veda_reference_ingest_summary.json`

No new runtime dependency is added beyond Pandas, already used by ORCAVEDA tests and benchmark utilities.
