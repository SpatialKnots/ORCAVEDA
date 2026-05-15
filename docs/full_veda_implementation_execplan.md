# Implement Comparable VEDA-Like PED for ORCAVEDA

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA currently provides Stage 3D v5.0, PED v1/v2 diagnostic outputs, Wilson GF-style PED audits, composed-coordinate PED diagnostics, and an opt-in Wilson GF diagonalization validation prototype. These layers are useful evidence for vibrational interpretation, but they do not yet constitute a comparable implementation of the original VEDA workflow.

The purpose of this plan is to define the remaining work needed for an additive, opt-in, comparable VEDA-like PED implementation. After the work in this plan is complete, a user should be able to run ORCAVEDA with a target interface such as `--veda-like-ped` and receive a self-contained set of `veda_like_*` outputs: a PED matrix, per-mode audit rows, basis diagnostics, mode-correspondence diagnostics, metadata describing the method, and a future comparison report when original VEDA reference outputs are available.

This plan does not claim VEDA equivalence. The forbidden wording for source code, output metadata, reports, documentation, and viewer text is `VEDA-equivalent`, `original VEDA reproduced`, and `strict VEDA PED` unless original VEDA reference outputs are added and a comparison suite passes. Until then, the correct wording is `comparable VEDA-like`, `VEDA-inspired`, `Wilson GF validation prototype`, `Wilson GF-style PED audit`, or `geometric and weighted independent-coordinate assignment audit`, depending on the exact layer.

## Progress

- [x] (2026-05-15) Initial plan written as `docs/full_veda_implementation_execplan.md`.
- [x] (2026-05-15) Current source and historical ExecPlans reviewed for Stage 3D, PED v1/v2, Wilson GF validation, and composed-coordinate diagnostics.
- [x] (2026-05-15) Formalized the first opt-in VEDA-like PED mathematics and output conventions before code changes.
- [x] (2026-05-15) Added first opt-in `--veda-like-ped` backend and `veda_like_*` artifacts on top of the closed Wilson GF validation core.
- [x] (2026-05-15) Extended focused VEDA-like validation beyond H2O with NH3/formaldehyde tests and NH3/formaldehyde/ethene probe outputs.
- [x] (2026-05-15) Ran mixed-molecule and full local `.hess` VEDA-like sweeps and recorded PASS/WARN diagnostics.
- [x] (2026-05-15) Hardened diagnostics for linear/near-linear fixed-conversion WARN cases and H-bond-dominated high-frequency X-H modes.
- [x] (2026-05-15) Added a reproducible `veda_like_*` output validator and used it on a fresh full local sweep.
- [ ] Promote the Wilson GF/PED backend from prototype/diagnostic-only status to a production-ready opt-in backend with explicit numerical gates across the full local `.hess` set.
- [ ] Define and implement the stable nonredundant or optimized internal-coordinate basis used by the VEDA-like backend.
- [ ] Extend composed-coordinate optimization with verifiable EPM-like metrics and conservative promotion rules.
- [x] Implement initial mode correspondence between ORCA normal modes, GF eigenvectors, and PED matrix rows for the closed Wilson GF/PED backend.
- [x] Add initial `veda_like_*` output artifacts behind an explicit opt-in interface.
- [ ] Add a future original-VEDA comparison harness that is inactive until VEDA reference outputs are checked in.

## Surprises & Discoveries

- Observation: Stage 3D remains a geometric and weighted independent-coordinate assignment audit, not a strict VEDA PED implementation.
  Evidence: `PLANS.md`, `AGENTS.md`, `src/ORCAVEDA_patched_stage3D_v5_0.py`, and `src/mode_assignment.py`.

- Observation: `src/ped.py` already contains PED v1/v2 and Wilson GF-style PED audit helpers, and the pipeline emits PED audit tables separately from `assignment_audit`.
  Evidence: `src/ped.py`, `src/ORCAVEDA_patched_stage3D_v5_0.py`, `docs/full_ped_execplan.md`, and `tests/test_ped.py`.

- Observation: `src/wilson_gf.py` exists as an opt-in Wilson GF diagonalization validation path, exposed through `--wilson-gf-validation`.
  Evidence: `src/wilson_gf.py`, `src/orcaveda_cli.py`, `src/ORCAVEDA_patched_stage3D_v5_0.py`, `tests/test_wilson_gf.py`, and `docs/wilson_gf_validation_execplan.md`.

- Observation: Composed-coordinate optimization already exists as a diagnostic VEDA-like precursor. It produces composed PED audit and basis diagnostic tables, but remains separate from Stage 3D final semantics.
  Evidence: `src/b_matrix.py`, `src/internal_coordinates.py`, `src/reports.py`, `src/ORCAVEDA_patched_stage3D_v5_0.py`, and `docs/veda_like_coordinate_optimization_execplan.md`.

- Observation: The target `--veda-like-ped` interface and first `veda_like_*` output schemas are now implemented in source after the first milestone.
  Evidence: `src/orcaveda_cli.py`, `src/wilson_gf.py`, `src/ORCAVEDA_patched_stage3D_v5_0.py`, `tests/test_wilson_gf.py`, and generated `outputs\veda_like_h2o_probe\H2O__veda_like_*.csv/json` artifacts.

- Observation: The initial closed Wilson GF/PED backend gives PASS diagnostics for small molecule probes H2O, NH3, formaldehyde, and ethene, with per-mode `contribution_percent` sums of 100.0 in the long-form PED matrix.
  Evidence: generated `outputs\veda_like_h2o_probe`, `outputs\veda_like_nh3_probe`, `outputs\veda_like_formaldehyde_probe`, `outputs\veda_like_ethene_probe`; `tests/test_wilson_gf.py` now covers NH3 and formaldehyde VEDA-like matrix normalization and high-frequency X-H dominant families.

- Observation: A five-file mixed subset completed with PASS status for acetophenone, benzoic acid, aniline, phenol, and pyridine. Aniline retained expected nonpositive-mode warnings while still producing PASS mode correspondence.
  Evidence: `outputs\veda_like_mixed_probe\acetophenone__benzoic_acid__aniline__plus_2_files__multi_file_5__veda_like_*.csv/json`.

- Observation: The full local `.hess` sweep completed for 55 files. Basis-level diagnostics reported 53 PASS and 2 WARN. No long-form PED matrix mode had a `contribution_percent` sum different from 100.0 after rounding to six decimals.
  Evidence: `outputs\veda_like_full_sweep_live\acetaldehyde__acetamide__acetanilide__plus_52_files__multi_file_55__veda_like_*.csv/json`.

- Observation: Full-sweep WARN files are ethyne and propyne. Ethyne reports `linear_bend_coordinate_used; fixed_conversion_failed; empirical_ratio_only`. Propyne reports `near_linear_bend_coordinate; fixed_conversion_failed; empirical_ratio_only`.
  Evidence: `veda_like_basis_diagnostics.csv` in `outputs\veda_like_full_sweep_live`.

- Observation: Most full-sweep high-frequency dominant rows are X-H stretch families, but two monoethanolamine dimer modes above 2800 cm-1 have a hydrogen-bond coordinate as rank 1 and an N-H stretch as rank 2. This is not evidence of original VEDA equivalence and should remain a diagnostic inspection item before production promotion.
  Evidence: `monoethanolamine_dimer_NH_to_O_DFT.hess` mode 63 and `monoethanolamine_dimer_OH_to_N_DFT.hess` mode 60 in full-sweep `veda_like_ped_audit.csv`.

- Observation: Linear and near-linear fixed-conversion failures now get the explicit warning token `linear_or_near_linear_fixed_conversion_review` in addition to the existing `fixed_conversion_failed` and `empirical_ratio_only` tokens.
  Evidence: `tests/test_wilson_gf.py`; regenerated `outputs\veda_like_edge_diagnostics_probe` for ethyne and propyne.

- Observation: H-bond-dominated high-frequency modes with secondary X-H stretch contributors now get the explicit mode-level warning token `high_frequency_hbond_dominates_xh_stretch_secondary` in the VEDA-like audit and matrix warning fields.
  Evidence: `tests/test_wilson_gf.py`; regenerated `outputs\veda_like_edge_diagnostics_probe` for the two monoethanolamine dimer cases.

- Observation: `tools\validate_veda_like_outputs.py` now provides a repeatable gate for `veda_like_*` artifacts by reading basis diagnostics, mode correspondence, PED matrix, and PED audit CSVs. It reports PASS/WARN/FAIL counts, warning-token counts, artifact row counts, and per-mode matrix normalization failures.
  Evidence: `tools\validate_veda_like_outputs.py` and `tests\test_validate_veda_like_outputs.py`.

## Decision Log

- Decision: Start with a planning-only milestone, not immediate algorithm rewrites.
  Rationale: The current repository already has several partially overlapping PED, Wilson, and composed-coordinate layers. A self-contained integration plan is the smallest safe step before changing scientific output semantics or adding new `veda_like_*` schemas.
  Date/Author: 2026-05-15 / Codex

- Decision: The future implementation must be additive and opt-in.
  Rationale: Existing `assignment_audit`, Stage 3D semantics, PED audit schemas, viewer payloads, and regression expectations are already used as baselines. A comparable VEDA-like backend changes scientific meaning and must not silently replace them.
  Date/Author: 2026-05-15 / Codex

- Decision: Do not require original VEDA reference outputs for this first plan.
  Rationale: No checked-in original VEDA outputs were provided with the request. The plan must prepare a comparison harness, but cannot use absent data as acceptance evidence.
  Date/Author: 2026-05-15 / Codex

- Decision: Reserve `VEDA-equivalent`, `original VEDA reproduced`, and `strict VEDA PED` until a future comparison suite passes against original VEDA reference outputs.
  Rationale: Current evidence supports VEDA-like and VEDA-inspired diagnostics, not equivalence to the original VEDA implementation.
  Date/Author: 2026-05-15 / Codex

## Outcomes & Retrospective

Initial outcome: this document establishes the integration roadmap and acceptance gates for a comparable VEDA-like implementation. No source-code behavior, CLI flags, output schemas, or scientific thresholds were changed by this milestone.

The first implementation milestone after this plan should update this section with the exact mathematical specification selected for the target PED matrix and the exact tests run.

First implementation specification outcome: the initial `--veda-like-ped` backend will reuse the existing Wilson GF diagonalization validation machinery as its numerical core, but will emit separate `veda_like_*` outputs with stricter method-boundary metadata. This is a comparable VEDA-like closed Wilson GF/PED artifact, not original VEDA reproduction.

First implementation outcome: `src/wilson_gf.py` now exposes VEDA-like builders for a top-N PED audit, long-form full PED matrix, mode correspondence, basis diagnostics, and metadata. `src/orcaveda_cli.py` accepts `--veda-like-ped`. `src/ORCAVEDA_patched_stage3D_v5_0.py` writes `__veda_like_ped_audit.csv`, `__veda_like_ped_matrix.csv`, `__veda_like_basis_diagnostics.csv`, `__veda_like_mode_correspondence.csv`, and `__veda_like_metadata.json` only when the flag is enabled. Default runs do not include these tables.

Validation outcome for this milestone:

- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py src\orcaveda_cli.py src\ORCAVEDA_patched_stage3D_v5_0.py tests\test_wilson_gf.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q` returned `15 passed` on the latest run after adding NH3/formaldehyde VEDA-like tests.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q` returned `34 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` returned `21 passed` on the latest run after the missing-Hessian matrix diagnostic row patch.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\veda_like_h2o_probe --veda-like-ped` completed in CLI mode.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\NH3.hess --outdir outputs\veda_like_nh3_probe --veda-like-ped` completed in CLI mode.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\formaldehyde.hess --outdir outputs\veda_like_formaldehyde_probe --veda-like-ped` completed in CLI mode.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\ethene.hess --outdir outputs\veda_like_ethene_probe --veda-like-ped` completed in CLI mode.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\acetophenone.hess data\hess\benzoic_acid.hess data\hess\aniline.hess data\hess\phenol.hess data\hess\pyridine.hess --outdir outputs\veda_like_mixed_probe --veda-like-ped` completed in CLI mode.
- Full local `.hess` sweep using all 55 files under `data\hess` completed in CLI mode with output under `outputs\veda_like_full_sweep_live`.
- Generated H2O artifacts under `outputs\veda_like_h2o_probe`: `H2O__veda_like_ped_audit.csv`, `H2O__veda_like_ped_matrix.csv`, `H2O__veda_like_basis_diagnostics.csv`, `H2O__veda_like_mode_correspondence.csv`, and `H2O__veda_like_metadata.json`.
- The H2O long-form PED matrix has 9 rows, sums to 100.0 percent for modes 6, 7, and 8, and mode correspondence status is `PASS`.
- Generated NH3, formaldehyde, and ethene artifacts under their `outputs\veda_like_*_probe` directories. Their mode-correspondence status is `PASS`, warnings are empty, and their long-form PED matrix sums to 100.0 percent for each mapped mode.
- NH3 high-frequency modes have dominant `N-H stretch` rows in the VEDA-like audit. Formaldehyde high-frequency modes have dominant `C-H stretch` rows in the VEDA-like audit. Ethene high-frequency modes have dominant `C-H stretch` rows in the VEDA-like audit.
- Mixed subset result: all five files have PASS mode correspondence. Aniline reports nonpositive-mode warnings already visible in basis diagnostics.
- Full-sweep result: 55 files processed; basis diagnostics report 53 PASS and 2 WARN; mode correspondence rows report 1560 PASS and 21 WARN; full PED matrix normalization has zero detected per-mode deviations from 100.0 after rounding to six decimals.
- Full-sweep limitation: ethyne and propyne remain WARN because the fixed SI conversion failed and only empirical ratio evidence is available. Two monoethanolamine dimer high-frequency modes are dominated by H-bond coordinates with N-H stretch as rank 2, so high-frequency X-H semantics need review before claiming production readiness.

Edge-diagnostic hardening outcome:

- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py tests\test_wilson_gf.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q` returned `19 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q` returned `21 passed`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\ethyne.hess data\hess\propyne.hess data\hess\monoethanolamine_dimer_NH_to_O_DFT.hess data\hess\monoethanolamine_dimer_OH_to_N_DFT.hess --outdir outputs\veda_like_edge_diagnostics_probe --veda-like-ped` completed in CLI mode.
- The edge probe reports `WARN` for ethyne and propyne with `linear_or_near_linear_fixed_conversion_review`.
- The edge probe reports `PASS` for both monoethanolamine dimers and adds `high_frequency_hbond_dominates_xh_stretch_secondary` on the two H-bond-dominated high-frequency modes.
- The edge-probe long-form PED matrix has zero detected per-mode normalization deviations from 100.0 after rounding to six decimals.

Validator outcome:

- `.\.venv312\Scripts\python.exe -m py_compile tools\validate_veda_like_outputs.py tests\test_validate_veda_like_outputs.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_validate_veda_like_outputs.py -q` returned `5 passed` after adding warning-token allowlist gate coverage.
- A fresh full local `.hess` sweep completed under `outputs\veda_like_full_sweep_validated_live`.
- `.\.venv312\Scripts\python.exe tools\validate_veda_like_outputs.py outputs\veda_like_full_sweep_validated_live --json-out outputs\veda_like_full_sweep_validated_live\veda_like_validation_summary.json --csv-out outputs\veda_like_full_sweep_validated_live\veda_like_normalization_failures.csv` returned a `WARN` summary, not `FAIL`.
- Validator summary for the fresh full sweep: 55 files, basis diagnostics 53 PASS and 2 WARN, mode correspondence 1560 PASS and 21 WARN, PED audit 12431 PASS and 156 WARN, 0 normalization failures, warning-token counts `fixed_conversion_failed=440`, `linear_or_near_linear_fixed_conversion_review=440`, and `high_frequency_hbond_dominates_xh_stretch_secondary=136`.
- `tools\validate_veda_like_outputs.py` now reports both `validation_status` and `acceptance_status`. Without an allowlist, the fresh full sweep reports `validation_status=WARN` and `acceptance_status=WARN`.
- With an explicit allowlist for the current full-sweep warning tokens, the fresh full sweep reports `validation_status=WARN` and `acceptance_status=PASS`. The allowed tokens were `empirical_ratio_only`, `fixed_conversion_failed`, `high_frequency_hbond_dominates_xh_stretch_secondary`, `linear_bend_coordinate_used`, `linear_or_near_linear_fixed_conversion_review`, `near_linear_bend_coordinate`, `nonpositive_gf_eigenvalues_within_expected_vibrational_space`, `nonpositive_orca_modes_within_expected_vibrational_space`, `positive_gf_eigenvalue_count_below_expected_vibrational_rank`, and `positive_orca_mode_count_below_expected_vibrational_rank`.

## Context and Orientation

The trusted project baseline is Stage 3D v5.0. The main legacy entrypoint is `src/ORCAVEDA_patched_stage3D_v5_0.py`. Modular source lives under `src/`.

Relevant current files:

- `src/orca_parser.py`: parses ORCA `.hess` data, including atoms, coordinates, frequencies, normal modes, IR intensities, and Cartesian Hessian data when available. The ORCA block parser must continue to support single-column block headers.
- `src/orcaveda_models.py`: defines shared data structures such as `HessData`, `InternalCoordinate`, and composed-coordinate metadata.
- `src/internal_coordinates.py`: builds primitive, functional-group, and first-wave composed internal coordinates.
- `src/b_matrix.py`: builds B matrices, selects independent-coordinate bases, runs EPM-like PED basis optimization, and handles composed-coordinate candidate matrices.
- `src/ped.py`: computes PED v1, force-aware PED v2, Wilson G/F diagnostics, and Wilson GF-style PED audit rows.
- `src/wilson_gf.py`: implements the opt-in Wilson GF diagonalization validation prototype.
- `src/mode_assignment.py`: builds Stage 3D assignment semantics and protects assignment-family wording.
- `src/reports.py`: builds PED/Stage 3D agreement, PED final assignment, composed diagnostic policy, and viewer payloads.
- `tests/test_ped.py`, `tests/test_wilson_gf.py`, `tests/test_stage3d_outputs.py`, and `tests/test_regression_baseline_outputs.py`: focused safety tests for the current layers.

Protected scientific constraints:

- Preserve ORCA normal-mode orientation as `normal_mode_vector = normal_modes[:, mode]`.
- Do not replace it with `normal_modes[mode, :]`.
- Preserve single-column ORCA block matrix header parsing.
- High-frequency C-H, N-H, and O-H modes must not remain unassigned without diagnostics.
- Do not hide singular bases, rank loss, ill-conditioning, failed mode correspondence, or missing reference data behind broad `try/except`.

Working terms:

- Stage 3D: ORCAVEDA's geometric and weighted independent-coordinate assignment audit.
- PED v1/v2: current ORCAVEDA diagnostic PED layers, not original VEDA reproduction.
- Wilson GF-style PED audit: current diagnostic projection using Wilson-style G/F ingredients.
- Wilson GF validation prototype: opt-in closed GF eigenproblem validation currently exposed through `--wilson-gf-validation`.
- Composed-coordinate optimization: VEDA-inspired coordinate mixing and EPM-like localization diagnostics.
- VEDA-like PED target: the future opt-in implementation described by this plan, with explicit method metadata and comparison-ready outputs.

## Plan of Work

First, write the mathematical specification for the target VEDA-like PED layer. This specification must define the coordinate basis, B-matrix convention, mass metric, internal force reconstruction, GF eigenproblem, normal-mode or GF-mode indexing, PED matrix orientation, percent normalization, sign handling, and handling of redundant or near-singular coordinates. It must state whether the target matrix is coordinate-by-mode or mode-by-coordinate and how that maps to CSV rows.

Initial mathematical specification:

- Normal-mode orientation remains protected as `normal_modes[:, mode]`. The first VEDA-like backend uses GF eigenvectors for PED rows and maps them to ORCA modes by sorted positive GF eigenvalues versus sorted positive ORCA frequencies; it does not transpose ORCA normal modes.
- The coordinate basis is the validation basis selected by `wilson_gf_diagonalization(...)`: the Stage 3D independent basis unless the existing Wilson GF conditioning fallback selects a better full-rank mass-metric basis. The selected indices are reported.
- The B matrix is the finite-difference internal-coordinate B matrix in Angstrom geometry units. Bend-like rows are scaled by `pi / 180` through `wilson_coordinate_scales(...)`.
- The mass metric is Wilson `G = B M^-1 B^T`, using atomic masses repeated over Cartesian components.
- The internal force matrix is reconstructed with the Wilson mass-metric back-transform `J = M^-1 B^T G^-1`, then `F_internal = J^T F_cart J`, with ORCA Bohr-based force constants converted to Angstrom-based units.
- The eigenproblem is solved as a symmetric `G^(1/2) F G^(1/2)` problem. Positive eigenvalues are sorted and paired with sorted positive ORCA frequencies.
- The PED contribution for coordinate `i` in GF mode `k` is `q_i * (F_internal q)_i`. The signed fraction is this term divided by the absolute-term normalization denominator. The reported contribution percent is `100 * abs(q_i * (F_internal q)_i) / sum_j(abs(q_j * (F_internal q)_j))`.
- The PED matrix is emitted in long form with one row per mapped mode and selected internal coordinate. Its orientation is explicitly recorded as `mode_rows_by_coordinate_columns_long_form`.
- The audit table is a top-N view of the matrix, with ranks per mode.
- The mode-correspondence table reports ORCA mode index, ORCA frequency, GF eigenvector index, GF eigenvalue, reconstructed frequency, relative error, validation status, and warnings.
- Singular basis, missing Cartesian Hessian, mode-count mismatch, fixed-conversion failure, nonpositive modes, and ill-conditioned G/F matrices remain visible as `PASS`, `WARN`, or `FAIL` diagnostics. They must not be hidden or converted into successful rows.

Second, harden the Wilson GF/PED backend into a production-ready opt-in component. The current `src/wilson_gf.py` is a validation prototype. Promotion requires stable public functions, clear diagnostics for rank and conditioning, deterministic basis selection, explicit linear or near-linear molecule handling, and regression tests showing default outputs remain unchanged when the VEDA-like flag is absent.

Third, define the nonredundant or optimized internal-coordinate basis used by the VEDA-like backend. The implementation must state whether it starts from the Stage 3D basis, the Wilson GF validation basis, the composed-coordinate PED basis, or a new VEDA-like basis selector. It must record rank, condition number, selected indices, rejected coordinates, and fallback decisions.

Fourth, extend composed-coordinate optimization beyond the current first-wave diagnostics only where evidence supports it. The EPM-like metric must be explicitly defined and reported. Any composed-coordinate promotion must preserve rank, preserve explicit X-H semantics, and report when a composed candidate improves localization, worsens semantic agreement, or is rejected.

Fifth, implement mode correspondence. The future `veda_like_mode_correspondence.csv` must map ORCA mode indices, ORCA frequencies, GF eigenvector indices, reconstructed frequencies, PED row indices, and correspondence status. If fixed conversion or mode counts fail, the row must say so explicitly.

Sixth, add output artifacts behind an explicit opt-in interface. The target interface is now `--veda-like-ped`. The first implemented output files are `veda_like_ped_audit.csv`, `veda_like_ped_matrix.csv`, `veda_like_basis_diagnostics.csv`, `veda_like_mode_correspondence.csv`, and `veda_like_metadata.json`. The initial implementation has been validated on H2O; broader molecule validation remains future work.

Seventh, add a future comparison harness for original VEDA outputs. This harness should remain inactive or skip with a clear message until original VEDA reference outputs are checked in. Once references exist, it should compare coordinate labels, mode correspondence, dominant PED contributors, matrix dimensions, percentage normalization, and any EPM-like metrics without requiring bit-for-bit identity unless the method definition justifies it.

## Concrete Steps

1. Confirm current source state and keep unrelated dirty files isolated.

       git status --short
       rg "normal_modes" src tests
       rg "veda_like|--veda-like-ped" src tests docs

   Expected result: note whether the VEDA-like target flag and target output schemas are absent from source, and confirm normal-mode orientation remains column-oriented.

2. Draft the VEDA-like mathematical spec in this file before code changes.

   The spec must include:

       normal-mode orientation
       mass weighting convention
       B-matrix convention
       G and F construction
       eigenproblem form
       PED contribution formula
       matrix orientation
       normalization and rounding rules
       singular-basis diagnostics
       mode-correspondence algorithm

3. Add focused tests for the specification before pipeline integration.

       .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q

   If new test files are added, keep them focused on small synthetic matrices plus H2O/NH3 fixtures before broad molecules.

4. Implement the opt-in backend in source.

   Likely touched files:

       src\wilson_gf.py
       src\ped.py
       src\b_matrix.py
       src\ORCAVEDA_patched_stage3D_v5_0.py
       src\orcaveda_cli.py
       tests\test_wilson_gf.py
       tests\test_ped.py

   The implementation must preserve existing output tables when `--veda-like-ped` is absent.

5. Emit target artifacts for H2O and NH3 first.

       .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\veda_like_h2o_probe --veda-like-ped
       .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\NH3.hess --outdir outputs\veda_like_nh3_probe --veda-like-ped

   Expected artifacts:

       *__veda_like_ped_audit.csv
       *__veda_like_ped_matrix.csv
       *__veda_like_basis_diagnostics.csv
       *__veda_like_mode_correspondence.csv
       *__veda_like_metadata.json

6. Run core regression tests after pipeline integration.

       .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q
       .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q

7. Run a local `.hess` sweep only after H2O/NH3 are stable.

       $hess = Get-ChildItem -LiteralPath 'data\hess' -Filter '*.hess' | Sort-Object Name | ForEach-Object { $_.FullName }
       .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\veda_like_full_sweep_live --veda-like-ped

   Record PASS/WARN/FAIL counts, rank losses, correspondence failures, fixed-conversion failures, high-frequency X-H diagnostics, and output schema validation results.

8. Add future VEDA comparison only when reference outputs exist.

       .\.venv312\Scripts\python.exe benchmarks\veda_compare\compare_veda_outputs.py --orcaveda outputs\veda_like_full_sweep_live --reference <checked-in-veda-reference-dir> --out outputs\veda_like_reference_comparison_live

   This command is a target placeholder. Do not implement a comparison that pretends missing VEDA reference outputs exist.

## Validation and Acceptance

This markdown-only milestone is accepted when `docs/full_veda_implementation_execplan.md` exists, follows the required `PLANS.md` section order, records current implementation boundaries, lists target interfaces as future work, and prohibits unsupported VEDA-equivalence wording.

The future code implementation is accepted only when:

- `--veda-like-ped` or the final chosen opt-in flag exists and is documented in code and metadata.
- Default ORCAVEDA outputs are unchanged when the opt-in flag is absent.
- `veda_like_ped_audit.csv` reports per-mode dominant contributors with normalized percentages and diagnostics.
- `veda_like_ped_matrix.csv` reports the full comparison-ready PED matrix with explicit orientation.
- `veda_like_basis_diagnostics.csv` reports basis size, rank, condition, selected coordinates, rejected coordinates, and fallback status.
- `veda_like_mode_correspondence.csv` maps ORCA modes, GF eigenvectors, and PED rows with PASS/WARN/FAIL status.
- `veda_like_metadata.json` states method name, method boundary, units, thresholds, source `.hess`, ORCAVEDA version label, and forbidden-equivalence disclaimer.
- H2O and NH3 have chemically sensible dominant bend/stretch contributors.
- At least one mixed-mode molecule such as acetophenone, benzoic acid, aniline, phenol, or pyridine has inspected evidence.
- High-frequency C-H, N-H, and O-H modes are assigned or carry explicit diagnostics.
- Focused tests pass:

       .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_wilson_gf.py -q

- Core regression tests pass:

       .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q

- A local `.hess` sweep completes or records explicit, reproducible failure diagnostics.
- Original VEDA comparison is reported only after reference outputs are available.

## Idempotence and Recovery

This document can be edited repeatedly as facts change. It does not alter source behavior.

Future generated outputs must go under `outputs/` with descriptive directories such as `outputs/veda_like_h2o_probe` and `outputs/veda_like_full_sweep_live`. They may be deleted and regenerated.

If `--veda-like-ped` implementation fails for a molecule, the pipeline must keep producing existing Stage 3D, PED, Wilson GF-style, composed PED, and viewer outputs unless the failure is in shared parsing. The VEDA-like outputs should report a diagnostic row where possible instead of silently skipping.

If rank loss, severe ill-conditioning, failed mode correspondence, or missing Cartesian Hessian data prevents VEDA-like computation, record a `WARN` or `FAIL` status with the exact reason. Do not broaden thresholds or change basis semantics silently.

If the comparison harness cannot find original VEDA references, it must report that reference outputs are absent and skip comparison. It must not fabricate reference rows, expected contributors, or PASS statuses.

## Artifacts and Notes

Initial artifact from this milestone:

- `docs/full_veda_implementation_execplan.md`

Existing related historical plans:

- `docs/full_ped_execplan.md`
- `docs/wilson_gf_validation_execplan.md`
- `docs/veda_like_coordinate_optimization_execplan.md`

Target future output artifacts:

- `veda_like_ped_audit.csv`
- `veda_like_ped_matrix.csv`
- `veda_like_basis_diagnostics.csv`
- `veda_like_mode_correspondence.csv`
- `veda_like_metadata.json`

Forbidden wording until future VEDA references and comparison pass:

- `VEDA-equivalent`
- `original VEDA reproduced`
- `strict VEDA PED`

Current allowed wording:

- `comparable VEDA-like PED`
- `VEDA-inspired composed-coordinate optimization`
- `EPM-like metric`
- `Wilson GF validation prototype`
- `Wilson GF-style PED audit`
- `geometric and weighted independent-coordinate assignment audit`

## Interfaces and Dependencies

Planned CLI interface:

- `--veda-like-ped`: opt-in flag for the first comparable VEDA-like closed Wilson GF/PED backend.

Planned output schema names:

- `veda_like_ped_audit`
- `veda_like_ped_matrix`
- `veda_like_basis_diagnostics`
- `veda_like_mode_correspondence`
- `veda_like_metadata`

Existing interfaces that must remain stable by default:

- `assignment_audit.csv`
- `ped_audit.csv`
- `ped_v2_force_audit.csv`
- `wilson_ped_audit.csv`
- `composed_ped_audit.csv`
- `composed_ped_v2_force_audit.csv`
- `composed_wilson_ped_audit.csv`
- `ped_stage3d_agreement.csv`
- `ped_final_assignment.csv`
- `composed_ped_policy_diagnostics.csv`
- `wilson_gf_validation.csv`, when `--wilson-gf-validation` is used
- `wilson_gf_ped_audit.csv`, when `--wilson-gf-validation` is used
- `wilson_gf_basis_diagnostics.csv`, when `--wilson-gf-validation` is used

No new dependency is allowed for the first implementation unless this plan is updated with a specific rationale, source impact, and validation impact. The expected dependency set remains NumPy, Pandas, and existing ORCAVEDA code.
