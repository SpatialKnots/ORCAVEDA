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

## Surprises & Discoveries

- Observation: The current production path imports and uses `finite_difference_B` directly in Stage 3D and Wilson GF paths.
  Evidence: `src/ORCAVEDA_patched_stage3D_v5_0.py`, `src/wilson_gf.py`, and `tests/test_wilson_gf.py`.

- Observation: `InternalCoordinate` already has enough metadata for safe first-wave analytical rows: `kind`, `atoms0`, `source`, and `fn`.
  Evidence: `src/orcaveda_models.py` and `src/internal_coordinates.py`.

## Decision Log

- Decision: Keep the first analytical B-matrix implementation additive and do not switch default pipeline behavior.
  Rationale: Analytical derivatives can change numerical rank/condition behavior; default scientific outputs must not shift without broader validation.
  Date/Author: 2026-05-15 / Codex

- Decision: Support only distance-like two-atom rows and regular three-atom angle/bend rows in the first patch.
  Rationale: These formulas are compact and directly comparable to the existing finite-difference baseline. Torsions, composed rows, and linear-bend components need separate treatment.
  Date/Author: 2026-05-15 / Codex

## Outcomes & Retrospective

Validation outcome:

- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py tests\test_b_matrix_analytical.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py -q` returned `3 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` returned `42 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` failed during collection with `ModuleNotFoundError: No module named 'ORCAVEDA_patched_stage3D_v5_0'` because `PYTHONPATH` was not set.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` returned `2 passed`.

## Context and Orientation

The current B matrix is built by `finite_difference_B(coords_A, internals)` in `src/b_matrix.py`. It evaluates each internal coordinate function at plus/minus Cartesian perturbations. Angles are returned in degrees, torsions in radians, and torsion finite differences wrap through `[-pi, pi]`.

The new `analytical_B(coords_A, internals)` returns a B matrix plus diagnostics. It reports per-row methods and fallback reasons.

## Plan of Work

First, add analytical formulas for distance and angle rows. Second, keep unsupported rows on finite-difference fallback. Third, test analytical rows against the finite-difference baseline on non-singular geometries. Fourth, leave production pipeline wiring unchanged until a broader full-suite numerical comparison is available.

## Concrete Steps

Run focused tests:

    .\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py -q

Run broader PED/Wilson safety tests before considering any future pipeline integration:

    .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q

## Validation and Acceptance

Accepted when analytical distance and angle rows match `finite_difference_B` within tight numerical tolerance, unsupported rows fall back with explicit diagnostics, and no default Stage 3D or Wilson GF code path changes.

## Idempotence and Recovery

The new API is additive. If validation finds a formula issue, callers can continue using `finite_difference_B` unchanged. The tests can be rerun without generated persistent artifacts.

## Artifacts and Notes

Added artifacts:

- `src/b_matrix.py`
- `tests/test_b_matrix_analytical.py`
- `docs/analytical_bmatrix_execplan.md`

## Interfaces and Dependencies

New API:

- `analytical_B(coords_A, internals, eps=EPS_FD_A, singular_tol=1.0e-12) -> tuple[np.ndarray, dict]`

Existing default API remains unchanged:

- `finite_difference_B(coords_A, internals, eps=EPS_FD_A) -> np.ndarray`
