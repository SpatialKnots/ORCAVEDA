# GAP 2 Analytical B-Matrix Review Note

Branch: `codex-gap2-analytical-bmatrix`

## Scope

GAP 2 adds a hybrid analytical B-matrix API while preserving the finite-difference default pipeline. The hybrid method supports:

- distance-like two-atom rows;
- regular angle/bend three-atom rows;
- regular torsion four-atom rows.

Unsupported coordinates, composed coordinates, singular or near-linear angles, high-angle rows near 180 degrees, and singular or near-linear torsions fall back to finite differences with diagnostics.

## Commits

- `6280ee7` Add hybrid analytical B-matrix rows
- `c45d725` Add B-matrix comparison harness
- `5d844dd` Validate hybrid analytical B-matrix diagnostics
- `9860914` Codify analytical B-matrix acceptance policy
- `f2961fd` Add opt-in hybrid analytical B-matrix wiring
- `aee078f` Add guarded analytical torsion B rows
- `da48d79` Prepare analytical B-matrix merge review

## Production Behavior

Default behavior remains `finite_difference_B`.

Hybrid analytical B is opt-in only:

    --b-matrix-method hybrid_analytical

or:

    analyze_orca_ped_like(..., b_matrix_method="hybrid_analytical")

Opt-in runs emit `b_matrix_diagnostics`.

## Acceptance Evidence

Final validation on 2026-05-16:

- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py src\ORCAVEDA_patched_stage3D_v5_0.py src\orcaveda_cli.py tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py tests\test_wilson_gf.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 11 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` -> 43 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `rows_above_tolerance_count=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=5`, `selected_basis_difference_count=10`, `selected_basis_replacement_rank_loss_count=0`.

## Limitations

- Exact selected-basis index identity is not required for the hybrid API; 10 selected-basis differences remain visible and rank-preserving in the full sweep.
- Composed-coordinate analytical rows remain future work.
- Linear-bend analytical components remain future work.
- No default Stage 3D, Wilson GF, or VEDA-like output path was switched to hybrid analytical B.

## Verdict

PASS. GAP 2 is ready for review/merge as an opt-in analytical B-matrix milestone.
