# ORCAVEDA Work Log

---
Task ID: 1
Agent: Super Z
Task: Clone ORCAVEDA repository and review current state for EPM optimization work

Work Log:
- Cloned/updated ORCAVEDA from https://github.com/SpatialKnots/ORCAVEDA.git
- Pulled latest changes (abb242b..a596dc1)
- Read key source files: wilson_gf.py, b_matrix.py, orcaveda_cli.py, AGENTS.md, PLANS.md
- Discovered 292 lines of uncommitted EPM implementation already in wilson_gf.py (likely from prior Codex session)
- The EPM implementation includes Level 1 (wilson_gf_ped_localization_metrics), Level 2 (optimize_wilson_gf_basis_for_epm), and Level 3 (integration in wilson_gf_diagonalization)
- Identified gap: CLI flags not yet added, epm_optimize not wired through run_orca_ped_like()

Stage Summary:
- EPM core implementation exists but is uncommitted
- CLI integration and pipeline wiring needed
- Collaboration framework (COLLABORATION.md) created for multi-agent paradigm

---
Task ID: 2
Agent: Super Z
Task: Add CLI flags for EPM optimization and create collaboration documentation

Work Log:
- Added --epm-optimize, --epm-max-passes, --epm-improvement-tol to orcaveda_cli.py
- Wired CLI flags through to run_orca_ped_like() call
- Created COLLABORATION.md with multi-agent protocol (User + Codex + Super Z)
- Updated AGENTS.md with full EPM implementation context for Codex
- Documented: what EPM is, current implementation status, remaining work, key design decisions
- Explicitly documented that b_matrix.optimize_independent_coordinates_for_ped() is the OLD geometric optimizer and should NOT be connected to the VEDA-like pipeline

Stage Summary:
- CLI integration complete (--epm-optimize and related flags)
- COLLABORATION.md and AGENTS.md updated for cross-agent context
- Remaining: wire epm_optimize through ORCAVEDA_patched_stage3D_v5_0.py, write EPM tests, generate patches

---
Task ID: 3
Agent: Codex
Task: Complete GAP 1 EPM wiring, tests, and validation on codex-gap1-epm-optimization

Work Log:
- Verified Super Z documentation commit `a056146` is present on branch `codex-gap1-epm-optimization`.
- Verified `epm_optimize`, `epm_max_passes`, and `epm_improvement_tol` are wired through `run_orca_ped_like()` / `analyze_orca_ped_like()` into both Wilson GF validation and VEDA-like diagnostics.
- Verified EPM tests exist for `H2O_freq.hess` and `ethene.hess`.
- Did not add CH4-specific test because no CH4 `.hess` fixture exists in `data/hess`.
- Updated `AGENTS.md` to mark the handoff items as completed and keep the old checklist as historical notes.

Validation:
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q` -> 23 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.

---
Task ID: 4
Agent: Codex
Task: Start GAP 2 analytical B-matrix on codex-gap2-analytical-bmatrix

Work Log:
- Created branch `codex-gap2-analytical-bmatrix` from current `main`.
- Reviewed `src/b_matrix.py`, `src/internal_coordinates.py`, `src/orcaveda_models.py`, and B-matrix/PED tests.
- Added additive `analytical_B(...)` hybrid API in `src/b_matrix.py`.
- Implemented analytical rows for distance-like two-atom coordinates and regular angle/bend three-atom coordinates.
- Kept torsions, composed coordinates, linear-bend components, and singular angle rows on finite-difference fallback with diagnostics.
- Added focused tests in `tests/test_b_matrix_analytical.py`.
- Added `docs/analytical_bmatrix_execplan.md`.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py tests\test_b_matrix_analytical.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py -q` -> 3 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` -> 42 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> collection failed with `ModuleNotFoundError: No module named 'ORCAVEDA_patched_stage3D_v5_0'` because `PYTHONPATH` was not set.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.

---
Task ID: 5
Agent: Codex
Task: Add GAP 2 finite-difference vs hybrid analytical B-matrix comparison harness

Work Log:
- Added `benchmarks/bmatrix_compare/compare_bmatrix_methods.py` and README.
- Added `tests/test_b_matrix_method_compare.py`.
- Extended `analytical_B(...)` with `angle_sin_tol=1.0e-3` so singular or near-linear angle rows fall back to finite differences with `singular_or_near_linear_angle` diagnostics.
- Added focused near-linear angle fallback coverage in `tests/test_b_matrix_analytical.py`.
- Updated `docs/analytical_bmatrix_execplan.md` with comparison results and remaining risks.
- Updated stale GAP status table in `COLLABORATION.md`.
- Production Stage 3D / Wilson GF pipeline remains on `finite_difference_B`; no pipeline switch was made.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py benchmarks\bmatrix_compare\compare_bmatrix_methods.py tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 6 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py tests\test_ped.py tests\test_wilson_gf.py -q` -> 48 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal` -> `file_count=4`, `files_with_rows_above_tolerance=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `files_with_rows_above_tolerance=3`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=4`.

Remaining:
- Full sweep still reports small angle-row deltas above `1e-5` in 3 files, max `8.065514768418325e-05`.
- 4 files select different basis indices at the same selected rank.
- Torsion analytical rows remain unimplemented and use finite-difference fallback.
- Do not switch the production pipeline to hybrid analytical B without reviewing these diagnostics.

---
Task ID: 6
Agent: Codex
Task: Session initialization and Caveman mode activation

Work Log:
- Read `COLLABORATION.md` and `WORKLOG.md` at session start.
- User requested `$caveman full`; activated compressed response style for subsequent work.

Validation:
- No tests run; no code changes requested.

---
Task ID: 7
Agent: Codex
Task: GAP 2 analytical B-matrix diagnostic review

Work Log:
- Added row-level atom-index and angle-geometry diagnostics to `benchmarks/bmatrix_compare/compare_bmatrix_methods.py`.
- Added `bmatrix_method_comparison_selected_basis_differences.csv` output for finite-difference vs hybrid selected-basis index changes.
- Updated `tests/test_b_matrix_method_compare.py`, `benchmarks/bmatrix_compare/README.md`, and `docs/analytical_bmatrix_execplan.md`.
- Re-ran default and full B-matrix comparison harnesses.
- Did not change `analytical_B` formulas, `angle_sin_tol`, Stage 3D, Wilson GF, or production B-matrix wiring.

Findings:
- Full sweep: `file_count=55`, `files_with_rows_above_tolerance=3`, `rows_above_tolerance_count=4`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=4`, `selected_basis_difference_count=8`.
- Remaining above-tolerance rows are high-angle coordinates, not below current `angle_sin_tol=1.0e-3`: two monoethanolamine dimer H-bond angles and duplicated phenyl isocyanate NCO bend rows.
- Selected-basis index changes are regular primitive angle-row swaps in aromatic systems at unchanged selected rank.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile benchmarks\bmatrix_compare\compare_bmatrix_methods.py tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 6 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` -> 42 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal` -> `file_count=4`, no rows/rank/selected-basis changes.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> counts listed above.

Remaining:
- Torsion analytical rows remain unimplemented and use finite-difference fallback.
- Production pipeline remains on `finite_difference_B`.
- Do not switch production to hybrid analytical B until high-angle row deltas and selected-basis swaps are explicitly accepted or guarded.

---
Task ID: 8
Agent: Codex
Task: Resolve GAP 2 high-angle analytical B row deltas

Work Log:
- Ran epsilon sensitivity on the four full-sweep rows above `1e-5`.
- Found the deltas shrink below tolerance when finite-difference `eps` is reduced below the default `1.0e-4`, indicating default finite-difference truncation error near high angles rather than an analytical formula mismatch.
- Changed `analytical_B(...)` default `angle_sin_tol` from `1.0e-3` to `2.0e-1`, so near-linear and high-angle rows fall back to the finite-difference baseline.
- Added a real-fixture regression test using `phenyl_isocyanate.hess` NCO bend.
- Updated `benchmarks/bmatrix_compare/README.md` and `docs/analytical_bmatrix_execplan.md`.
- Did not switch Stage 3D, Wilson GF, or production B-matrix wiring.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py tests\test_b_matrix_analytical.py benchmarks\bmatrix_compare\compare_bmatrix_methods.py tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 7 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` -> 42 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal` -> `file_count=4`, `rows_above_tolerance_count=0`, no rank/selected-basis changes.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `files_with_rows_above_tolerance=0`, `rows_above_tolerance_count=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=4`, `selected_basis_difference_count=8`.

Remaining:
- Selected-basis index swaps remain in 4 aromatic fixtures at unchanged selected rank.
- Torsion analytical rows remain unimplemented and use finite-difference fallback.
- Production pipeline remains on `finite_difference_B`.

---
Task ID: 9
Agent: Codex
Task: Trace GAP 2 selected-basis swaps

Work Log:
- Added selected-basis replacement diagnostics to `benchmarks/bmatrix_compare/compare_bmatrix_methods.py`.
- New selected-basis CSV fields include `replacement_rank_preserved`, replacement rank, replacement condition, and replacement minimum singular value for both finite and hybrid matrices.
- Added focused regression coverage in `tests/test_b_matrix_method_compare.py`.
- Updated `benchmarks/bmatrix_compare/README.md` and `docs/analytical_bmatrix_execplan.md`.
- Did not change production selector logic or production B-matrix wiring.

Findings:
- Full sweep still has `files_with_selected_basis_index_change=4` and `selected_basis_difference_count=8`.
- All 8 selected-basis swaps preserve selected rank when rows are substituted both ways.
- `selected_basis_replacement_rank_loss_count=0`.
- Replacement minimum singular values are around `2.494e-1` to `3.654e-1`, far above the `1.0e-6` rank tolerance.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile benchmarks\bmatrix_compare\compare_bmatrix_methods.py tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_method_compare.py -q` -> 3 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal` -> `selected_basis_replacement_rank_loss_count=0`.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `rows_above_tolerance_count=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=4`, `selected_basis_difference_count=8`, `selected_basis_replacement_rank_loss_count=0`.

Remaining:
- Exact selected-basis index identity still differs in 4 aromatic fixtures.
- Production integration needs an explicit acceptance policy for rank-preserving selected-index differences.
- Torsion analytical rows remain unimplemented and use finite-difference fallback.

---
Task ID: 20
Agent: Codex
Task: Session initialization and Caveman full mode activation

Work Log:
- Read `COLLABORATION.md`, `AGENTS.md`, `WORKLOG.md`, and `orcaveda-core` skill instructions at session start.
- User requested `$caveman full`; activated compressed response style for subsequent work.
- No scientific/code behavior changes requested.

Validation:
- `git status --short` before the WORKLOG append -> `?? incoming/`.
- No tests run; no code changes requested.

Remaining:
- Await next user task.

---
Task ID: 18
Agent: Codex
Task: Reapply GAP 3 VEDA reference validation harness on current main

Work Log:
- Created branch `codex/gap3-veda-validation-current` from current `main`.
- Confirmed old `codex-gap3-veda-validation` was stale relative to `main` and would remove GAP 2 analytical B-matrix files if merged directly.
- Reapplied only GAP 3 artifacts from the old branch: `benchmarks/veda_compare/`, `docs/veda_reference_validation_execplan.md`, and `tests/test_veda_reference_compare.py`.
- Updated `docs/full_veda_implementation_execplan.md` and `docs/veda_reference_validation_execplan.md` with current 2026-05-18 validation evidence.
- Did not modify Stage 3D, Wilson GF, analytical B-matrix, output schemas, thresholds, or scientific assignment logic.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile benchmarks\veda_compare\compare_veda_outputs.py tests\test_veda_reference_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_veda_reference_compare.py tests\test_validate_veda_like_outputs.py -q` -> 9 passed.
- `.\.venv312\Scripts\python.exe benchmarks\veda_compare\compare_veda_outputs.py --orcaveda outputs\pytest_veda_like_epm_opt_in_h2o --reference data\veda_reference_missing --out outputs\veda_reference_compare_skip_probe` -> `comparison_status=SKIP`, `acceptance_status=SKIP`, `reason=veda_reference_directory_missing`.

Remaining:
- Original VEDA reference outputs are still absent; no VEDA reproduction claim is supported.
- Review, commit, and optionally push `codex/gap3-veda-validation-current`.

---
Task ID: 19
Agent: Codex
Task: Merge GAP 3 and start VEDA reference fixture ingest

Work Log:
- Merged `codex/gap3-veda-validation-current` into `main` with merge commit `b5bb0d0`.
- Pushed `main` to origin.
- Created branch `codex/veda-reference-fixture-ingest` from updated `main`.
- Added conservative reference ingest tool `benchmarks/veda_compare/convert_veda_reference.py`.
- Added `tests/test_veda_reference_ingest.py`.
- Added `docs/veda_reference_fixture_ingest_execplan.md`.
- Updated `benchmarks/veda_compare/README.md` with ingest commands.
- The ingest tool accepts already-normalized reference CSVs or explicit column mappings only. It does not infer unknown native VEDA export formats and does not fabricate rows.
- Did not modify Stage 3D, Wilson GF, analytical B-matrix, output schemas, thresholds, or scientific assignment logic.

Validation:
- `git pull` on `main` -> already up to date before merge.
- `git merge --no-ff codex/gap3-veda-validation-current -m "Merge GAP 3 VEDA reference validation harness"` -> completed without conflicts.
- `git push` -> pushed `main` from `f8d8285` to `b5bb0d0`.
- `.\.venv312\Scripts\python.exe -m py_compile benchmarks\veda_compare\convert_veda_reference.py tests\test_veda_reference_ingest.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_veda_reference_ingest.py tests\test_veda_reference_compare.py -q` -> 8 passed.
- `.\.venv312\Scripts\python.exe benchmarks\veda_compare\convert_veda_reference.py --raw-reference data\veda_reference_missing --out outputs\veda_reference_ingest_skip_probe` -> `conversion_status=SKIP`, `acceptance_status=SKIP`, `reason=raw_reference_directory_missing`.

Remaining:
- Original VEDA reference outputs are still absent.
- Commit and push `codex/veda-reference-fixture-ingest` if accepted.

---
Task ID: 11
Agent: Codex
Task: Encode GAP 2 selected-basis acceptance policy

Work Log:
- Pushed commit `5d844dd` to `origin/codex-gap2-analytical-bmatrix`.
- Added a full-sweep acceptance-policy regression test for the hybrid analytical B comparison harness.
- Documented that exact selected-basis index identity is not required when row deltas, redundant rank, selected rank, and replacement rank preservation pass.
- Did not switch Stage 3D, Wilson GF, or production B-matrix wiring.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_method_compare.py -q` -> 4 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `rows_above_tolerance_count=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `selected_basis_difference_count=8`, `selected_basis_replacement_rank_loss_count=0`.

Remaining:
- Production integration still needs a separate explicit opt-in/default decision.
- Torsion analytical rows remain unimplemented and use finite-difference fallback.

---
Task ID: 12
Agent: Codex
Task: Add opt-in hybrid analytical B-matrix production wiring

Work Log:
- Added explicit `b_matrix_method` plumbing to the main ORCAVEDA pipeline and CLI.
- Default remains `finite_difference`; opt-in `hybrid_analytical` uses `analytical_B(...)` for primitive B construction.
- Added opt-in `b_matrix_diagnostics` output rows and manifest text.
- Added tests for the CLI argument and H2O opt-in pipeline diagnostics.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\ORCAVEDA_patched_stage3D_v5_0.py src\orcaveda_cli.py tests\test_wilson_gf.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 9 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q` -> 24 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q` -> 19 passed.

Remaining:
- Torsion analytical rows remain unimplemented and use finite-difference fallback.
- Default production behavior must remain finite-difference unless explicitly changed later.

---
Task ID: 13
Agent: Codex
Task: Add guarded analytical torsion B-matrix rows

Work Log:
- Added analytical torsion rows for regular four-atom torsion coordinates using chain-rule differentiation of the existing projected-plane `atan2` torsion definition.
- Added `singular_or_near_linear_torsion` fallback for torsions with small projected torsion-plane sine.
- Initial full sweep without the near-linear torsion guard failed acceptance with `rows_above_tolerance_count=8`; problematic rows were near-linear nitrile/alkyne-axis torsions in `benzonitrile.hess`, `CH3CN_freq.hess`, and `propyne.hess`.
- After adding the guard, full-sweep row/rank acceptance was restored; selected-basis alternatives increased to 10 and remained rank-preserving.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 11 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `rows_above_tolerance_count=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `selected_basis_difference_count=10`, `selected_basis_replacement_rank_loss_count=0`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` -> 43 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.

Remaining:
- Default production behavior must remain finite-difference unless explicitly changed later.

---
Task ID: 14
Agent: Codex
Task: Prepare GAP 2 merge-readiness review package

Work Log:
- Reviewed current GAP 2 docs and source metadata for stale analytical B-matrix wording.
- Updated `docs/analytical_bmatrix_execplan.md`, `benchmarks/bmatrix_compare/README.md`, and opt-in manifest text to include regular torsion analytical rows and near-linear torsion fallback.
- Left older `WORKLOG.md` entries intact as chronological evidence; their earlier "torsion unimplemented" notes are superseded by Task 13.
- Confirmed default production behavior remains `finite_difference`; hybrid analytical B remains opt-in via `b_matrix_method="hybrid_analytical"` / `--b-matrix-method hybrid_analytical`.

Review Summary:
- Current hybrid analytical B supports distance-like rows, regular angle/bend rows, and regular torsion rows.
- Singular or near-linear angles and torsions fall back to finite differences with explicit diagnostics.
- Full-sweep selected-basis exact index identity is not required; 10 selected-basis differences remain visible and rank-preserving.
- No default Stage 3D, Wilson GF, or VEDA-like output path was switched to hybrid analytical B.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py src\ORCAVEDA_patched_stage3D_v5_0.py src\orcaveda_cli.py tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py tests\test_wilson_gf.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 11 passed.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` -> 43 passed.
- `$env:PYTHONPATH='src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` -> 2 passed.
- `.\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full` -> `file_count=55`, `rows_above_tolerance_count=0`, `files_with_redundant_rank_change=0`, `files_with_selected_rank_change=0`, `files_with_selected_basis_index_change=5`, `selected_basis_difference_count=10`, `selected_basis_replacement_rank_loss_count=0`.

Remaining:
- GAP 2 is ready for review/merge as an opt-in analytical B-matrix milestone.
- Future work, separate from this merge package: composed-coordinate analytical rows, linear-bend components, or broader default-switch policy.

---
Task ID: 15
Agent: Codex
Task: Add GAP 2 review handoff note

Work Log:
- `gh` CLI was not available in the local environment, so no pull request was created from the terminal.
- Added `docs/analytical_bmatrix_review_note.md` with branch scope, commit list, production behavior, validation evidence, limitations, and verdict.
- Did not merge to `main`; project protocol reserves merging for the user.

Validation:
- No code behavior changed in this task.
- Relied on the final validation recorded in Task 14.

Remaining:
- User review/merge to `main` if accepted.
- Future work: composed-coordinate analytical rows, linear-bend components, or broader default-switch policy.

---
Task ID: 16
Agent: Codex
Task: Merge GAP 2 analytical B-matrix branch to main

Work Log:
- User explicitly authorized Codex to perform the merge.
- Checked `codex-gap2-analytical-bmatrix` was clean and contained 8 commits ahead of `main`.
- Switched to `main`, pulled from origin, and merged `codex-gap2-analytical-bmatrix` with a non-fast-forward merge commit.
- Pushed `main` to origin.

Validation:
- `git pull` on `main` -> already up to date before merge.
- `git merge --no-ff codex-gap2-analytical-bmatrix -m "Merge GAP 2 analytical B-matrix"` -> completed without conflicts.
- `git push` -> pushed `main` from `9b25516` to `19d533c`.

Remaining:
- Commit and push this WORKLOG merge note.
- Next project step: start GAP 3 VEDA reference validation or a separate GAP 2 follow-up branch for composed/linear-bend analytical rows.

---
Task ID: 17
Agent: Codex
Task: Session initialization and Caveman full mode activation

Work Log:
- Read `COLLABORATION.md`, `AGENTS.md`, `WORKLOG.md`, and `orcaveda-core` skill instructions at session start.
- User requested `$caveman full`; activated compressed response style for subsequent work.
- No scientific/code behavior changes requested.

Validation:
- `git status --short` before the WORKLOG append -> clean.
- No tests run; no code changes requested.

Remaining:
- Await next user task.

---
Task ID: 10
Agent: Codex
Task: Commit GAP 2 analytical B-matrix validation state

Work Log:
- Re-read `COLLABORATION.md` and current `WORKLOG.md` before committing.
- Checked the dirty tree on `codex-gap2-analytical-bmatrix`.
- Confirmed the patch scope remains limited to hybrid `analytical_B`, comparison diagnostics, tests, docs, and worklog.
- Did not switch Stage 3D, Wilson GF, or production B-matrix wiring.

Validation:
- `.\.venv312\Scripts\python.exe -m py_compile src\b_matrix.py tests\test_b_matrix_analytical.py benchmarks\bmatrix_compare\compare_bmatrix_methods.py tests\test_b_matrix_method_compare.py` -> completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_b_matrix_analytical.py tests\test_b_matrix_method_compare.py -q` -> 8 passed.

Remaining:
- Exact selected-basis index identity still differs in 4 aromatic fixtures.
- Production integration needs an explicit acceptance policy for rank-preserving selected-index differences.
- Torsion analytical rows remain unimplemented and use finite-difference fallback.
