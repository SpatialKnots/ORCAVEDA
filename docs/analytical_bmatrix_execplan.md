# Add Analytical B-Matrix Rows

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

GAP 2 adds analytical B-matrix derivatives where the mathematics is simple and verifiable, while preserving the existing finite-difference pipeline. The first user-visible result is a tested hybrid API that can build analytical rows for distance-like stretches and regular angle/bend coordinates, then fall back to finite differences for unsupported rows.

## Progress

- [x] (2026-05-15) Initial GAP 2 branch created from `main`.
- [x] (2026-05-15) Reviewed current `finite_difference_B`, `InternalCoordinate`, and internal-coordinate builders.
- [x] (2026-05-15) Added additive `analytical_B(...)` hybrid builder for distance and angle rows with finite-difference fallback.
- [x] (2026-05-15) Added focused analytical B-matrix tests.
- [x] (2026-05-15) Focused and PED/Wilson safety tests completed.
- [x] (2026-05-15) Added finite-difference vs hybrid analytical B-matrix comparison harness.
- [x] (2026-05-15) Ran default four-file probe and full `data/hess/*.hess` sweep.

## Surprises & Discoveries

- Observation: The current production path imports and uses `finite_difference_B` directly in Stage 3D and Wilson GF paths.
  Evidence: `src/ORCAVEDA_patched_stage3D_v5_0.py`, `src/wilson_gf.py`, and `tests/test_wilson_gf.py`.

- Observation: `InternalCoordinate` already has enough metadata for safe first-wave analytical rows: `kind`, `atoms0`, `source`, and `fn`.
  Evidence: `src/orcaveda_models.py` and `src/internal_coordinates.py`.

- Observation: Near-linear angle rows can create large finite-difference vs analytical deltas even when rank is preserved.
  Evidence: Initial full sweep showed large row deltas for nearly linear angle coordinates in `CH3CN_freq.hess`, `ethyne.hess`, and `propyne.hess`, with angle sine values around `7e-5` to `2.6e-4`.

- Observation: After adding near-linear angle fallback, the full sweep still reports 3 files with small angle-row deltas above `1e-5`, no redundant rank changes, no selected-rank changes, and 4 files where selected basis indices differ at the same rank.
  Evidence: `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full`.

## Decision Log

- Decision: Keep the first analytical B-matrix implementation additive and do not switch default pipeline behavior.
  Rationale: Analytical derivatives can change numerical rank/condition behavior; default scientific outputs must not shift without broader validation.
  Date/Author: 2026-05-15 / Codex

- Decision: Support only distance-like two-atom rows and regular three-atom angle/bend rows in the first patch.
  Rationale: These formulas are compact and directly comparable to the existing finite-difference baseline. Torsions, composed rows, and linear-bend components need separate treatment.
  Date/Author: 2026-05-15 / Codex

- Decision: Treat angle rows with `sin(theta) <= 1.0e-3` as `singular_or_near_linear_angle` and use finite-difference fallback.
  Rationale: The analytical angle derivative has a `1/sin(theta)` factor and is numerically unsafe for nearly linear geometries relative to the finite-difference baseline used by the current production pipeline.
  Date/Author: 2026-05-15 / Codex

- Decision: Add a standalone benchmark harness instead of wiring `analytical_B` into Stage 3D or Wilson GF.
  Rationale: The sweep provides evidence about row deltas and rank/selection behavior without changing scientific outputs.
  Date/Author: 2026-05-15 / Codex

## Outcomes & Retrospective

Validation outcome:

- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py tests\test_b_matrix_analytical.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py -q` returned `3 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` returned `42 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` failed during collection with `ModuleNotFoundError: No module named 'ORCAVEDA_patched_stage3D_v5_0'` because `PYTHONPATH` was not set.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` returned `2 passed`.

Comparison harness outcome:

- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py benchmarks\bmatrix_compare\compare_bmatrix_methods.py tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` returned `6 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py tests\test_ped.py tests\test_wilson_gf.py -q` returned `48 passed`.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` returned `2 passed`.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal` returned `file_count=4`, `files_with_rows_above_tolerance=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` returned `file_count=55`, `files_with_rows_above_tolerance=3`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=4`.
- Remaining rows above `1e-5` are angle rows in `monoethanolamine_dimer_NH_to_O_DFT.hess`, `monoethanolamine_dimer_OH_to_N_DFT.hess`, and `phenyl_isocyanate.hess`; max reported delta is `8.065514768418325e-05`.

## Context and Orientation

The current B matrix is built by `finite_difference_B(coords_A, internals)` in `src/b_matrix.py`. It evaluates each internal coordinate function at plus/minus Cartesian perturbations. Angles are returned in degrees, torsions in radians, and torsion finite differences wrap through `[-pi, pi]`.

The new `analytical_B(coords_A, internals)` returns a B matrix plus diagnostics. It reports per-row methods and fallback reasons. Near-linear angle rows fall back to finite differences through `angle_sin_tol=1.0e-3`.

## Plan of Work

First, add analytical formulas for distance and regular angle rows. Second, keep unsupported and near-linear rows on finite-difference fallback. Third, test analytical rows against the finite-difference baseline on non-singular geometries. Fourth, use the benchmark harness to compare row deltas and rank behavior across real `.hess` fixtures. Fifth, leave production pipeline wiring unchanged until the remaining deltas and selection-index changes are reviewed.

## Concrete Steps

Run focused tests:

    .\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py -q

Run broader PED/Wilson safety tests before considering any future pipeline integration:

    .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q

Run the B-matrix method comparison harness:

    .\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal
    .\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full

## Validation and Acceptance

Accepted for this milestone when analytical distance and regular angle rows match `finite_difference_B` within tight numerical tolerance on focused tests, unsupported and near-linear rows fall back with explicit diagnostics, the comparison harness reports row deltas/rank/condition/selection diagnostics, and no default Stage 3D or Wilson GF code path changes.

## Idempotence and Recovery

The new API is additive. If validation finds a formula issue, callers can continue using `finite_difference_B` unchanged. The tests can be rerun without generated persistent artifacts.

## Artifacts and Notes

Added artifacts:

- `src/b_matrix.py`
- `tests/test_b_matrix_analytical.py`
- `docs/analytical_bmatrix_execplan.md`
- `benchmarks/bmatrix_compare/compare_bmatrix_methods.py`
- `benchmarks/bmatrix_compare/README.md`
- `tests/test_b_matrix_method_compare.py`

## Interfaces and Dependencies

New API:

- `analytical_B(coords_A, internals, eps=EPS_FD_A, singular_tol=1.0e-12, angle_sin_tol=1.0e-3) -> tuple[np.ndarray, dict]`

Existing default API remains unchanged:

- `finite_difference_B(coords_A, internals, eps=EPS_FD_A) -> np.ndarray`
