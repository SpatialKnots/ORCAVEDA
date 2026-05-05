# Implement Full PED for ORCAVEDA

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA currently has Stage 3D, a geometric and weighted independent-coordinate assignment audit. It is useful for assigning vibrational modes, but it is not a strict VEDA PED or full Wilson GF PED implementation. This plan adds a separate PED layer that decomposes each normal mode into internal-coordinate contributions and reports percentage-like mode composition in a scientifically traceable way.

After this work, a user should be able to run ORCAVEDA on an ORCA `.hess` file and receive, alongside the existing Stage 3D assignment audit, a PED audit showing the dominant bond, angle, torsion, and out-of-plane contributions for each vibrational mode. This should improve interpretation of mixed modes such as carbonyl/ring coupling, amide NH2 deformation, carboxylic acid O-H bending context, and aromatic ring modes without pretending that PED fixes frequency errors or method/scale-factor mismatch.

The initial implementation must be additive. Stage 3D remains available and must not be renamed as full PED. PED results should be emitted as separate outputs until they are validated enough to participate in combined assignment decisions.

## Progress

- [x] (2026-05-05) Initial implementation plan written.
- [x] (2026-05-05) Audited existing `.hess` parser, normal-mode orientation, masses, Hessian, and internal-coordinate data flow.
- [x] (2026-05-05) Specified PED v1 as normalized B-matrix internal-coordinate projection, not force-constant Wilson GF PED.
- [x] (2026-05-05) Implemented additive PED module and focused unit tests.
- [x] (2026-05-05) Integrated separate `ped_audit` output table into the ORCAVEDA pipeline.
- [x] (2026-05-05) Validated PED v1 on synthetic modes, water pipeline output, golden RDKit/NIST focused tests, and Stage 3D regression wrappers.
- [x] (2026-05-05) Ran PED benchmark review on water, ammonia, acetophenone, benzoic acid, aniline, phenol, and pyridine.
- [x] (2026-05-05) Decided that PED v1 should remain evidence/diagnostic only for now; no automatic final-label override yet.
- [x] (2026-05-05) Added PED-aware benchmark diagnostics with optional `--ped-audit`, raw/scaled multiscale correspondence, and PED semantic windows.
- [x] (2026-05-05) Ran PED-aware benchmark validation on the current 10-molecule assignment benchmark.
- [x] (2026-05-05) Implemented PED v2 as a force-aware diagnostic using parsed ORCA `$hessian`.
- [x] (2026-05-05) Validated PED v2 on synthetic force-weighting, H2O pipeline output, and the 10-molecule benchmark.
- [x] (2026-05-05) Implemented Wilson GF-style PED audit with G matrix, reconstructed internal F matrix, and mode-projected potential-energy terms.
- [x] (2026-05-05) Validated Wilson PED on H2O, synthetic G/F tests, and the 10-molecule benchmark.
- [x] (2026-05-06) Switched the interactive viewer assignment text to the strongest available PED layer, with priority `wilson_ped_audit`, then `ped_v2_force_audit`, then `ped_audit`, while preserving Stage 3D assignment text as a separate diagnostic field.

## Surprises & Discoveries

- Observation: `HessData` currently stores `filename`, `atoms`, `masses`, `coords_A`, `frequencies_cm1`, `ir_intensities`, `normal_modes`, `temperature_K`, and `frequency_scale_factor`; it does not store a separate Cartesian Hessian or internal force-constant matrix.
  Evidence: `src/orcaveda_models.py` and `src/orca_parser.py`.

- Observation: The parser reads `$atoms`, `$vibrational_frequencies`, `$normal_modes`, and `$ir_spectrum`; normal modes are parsed as a `(3N, 3N)` matrix and existing assignment code uses `hess.normal_modes[:, mode]`.
  Evidence: `src/orca_parser.py`, `src/mode_assignment.py`, and `rg "normal_modes" src tests`.

- Observation: A finite-difference B-matrix implementation already exists and is used by Stage 3D independent-coordinate selection.
  Evidence: `src/b_matrix.py` and `src/ORCAVEDA_patched_stage3D_v5_0.py`.

- Observation: PED v1 has been added as an output-only layer. It writes a separate `ped_audit` table and does not change `assignment_audit` wording.
  Evidence: `src/ped.py` and the `ped_audit` table integration in `src/ORCAVEDA_patched_stage3D_v5_0.py`.

- Observation: PED v1 is useful for exposing chemically meaningful contributors in several benchmark problem cases, but it is not yet reliable enough to override Stage 3D automatically. Acetophenone C=O is correctly exposed by both Stage 3D and PED on the semantically selected mode, but raw frequency matching still needs the multiscale/scale-factor layer. Benzoic acid mixed ring-acid rows show useful O-H torsion/bend and carbonyl context, but some in-plane mixed rows remain aromatic C-H dominated. Aniline reveals strong NH2 scissor evidence in one key row and N-H bend/C-N context in others, but not uniformly. Phenol and pyridine mostly benefit from PED as a more detailed explanation beneath broad labels.
  Evidence: `outputs/ped_benchmark_review/ped_benchmark_summary.csv`, `outputs/ped_benchmark_review/ped_problem_cases_summary.csv`, and `outputs/ped_benchmark_review/ped_key_cases_compact.csv`.

- Observation: The PED-aware comparator preserves legacy Stage 3D status columns and adds PED-only diagnostics. On the 10-molecule benchmark, raw-primary comparison produced 28 PASS, 63 WARN, and 2 FAIL; scaled-primary with explicit scale factor 0.96 produced 27 PASS, 66 WARN, and 0 FAIL. PED semantic status on the raw-primary report produced 33 PASS, 55 WARN, and 5 FAIL. PED supported benchmark semantics in 5 raw-primary rows where the Stage 3D/frequency layer was not fully clean.
  Evidence: `outputs/ped_aware_benchmark_10/ped_aware_comparison_raw.csv`, `outputs/ped_aware_benchmark_10/ped_aware_comparison_scaled.csv`, and `outputs/ped_aware_benchmark_10/ped_aware_summary_raw.csv`.

- Observation: ORCA `.hess` files in `data/hess` contain a `$hessian` block. The parser now stores it as `HessData.cartesian_hessian`, and the pipeline emits `ped_v2_force_audit` when that matrix is available.
  Evidence: `src/orca_parser.py`, `src/orcaveda_models.py`, `src/ped.py`, and generated `outputs/ped_v2_h2o_probe/H2O__ped_v2_force_audit.csv`.

- Observation: PED v2 changes contribution percentages but does not change the legacy benchmark status, because it is still a diagnostic layer. On the 10-molecule benchmark, PED v2 raw-primary legacy status remained 28 PASS, 63 WARN, and 2 FAIL. PED v2 semantic status was 32 PASS, 55 WARN, and 6 FAIL.
  Evidence: `outputs/ped_v2_benchmark_10/ped_v2_aware_comparison_raw.csv` and `outputs/ped_v2_benchmark_10/ped_v2_summary_raw.csv`.

- Observation: Wilson PED produces sharper diagnostic percentages for simple modes, but it is stricter on mixed benchmark rows. On H2O, Wilson PED gives about 100% H-O-H bend for the bend mode and about 50/50 O-H stretch for stretch modes. On the 10-molecule benchmark, Wilson semantic status was 31 PASS, 51 WARN, and 11 FAIL. The additional FAIL rows mostly reflect mixed benchmark expectations where the Wilson top terms do not include every expected semantic class.
  Evidence: `outputs/wilson_ped_h2o_probe/H2O__wilson_ped_audit.csv`, `outputs/wilson_ped_benchmark_10/wilson_ped_comparison_raw.csv`, and `outputs/wilson_ped_benchmark_10/wilson_ped_summary_raw.csv`.

- Observation: The interactive viewer now uses PED-derived assignment text when a PED table is available. The most complete available source is chosen in this order: Wilson PED, PED v2 force-aware, PED v1 B-matrix projection, then Stage 3D fallback. Stage 3D assignment and supporting coordinates remain visible in mode details for comparison.
  Evidence: `src/reports.py`, `src/ORCAVEDA_patched_stage3D_v5_0.py`, `tests/test_interactive_spectrum_viewer.py`, and generated `outputs/ped_frontend_monoethanolamine/monoethanolamine_DFT_therm__spectrum_data.json`.

## Decision Log

- Decision: Implement PED as a separate additive layer rather than replacing Stage 3D.
  Rationale: Stage 3D is already a working assignment-audit layer with existing outputs and tests. Full PED changes scientific meaning and should be introduced with separate diagnostics and validation.
  Date/Author: 2026-05-05 / Codex

- Decision: Start with a B-matrix/internal-coordinate projection PED v1, then evaluate force-constant PED as v2.
  Rationale: A B-matrix projection is the smallest safe route to traceable internal-coordinate percentages. A force-constant Wilson GF treatment is more chemically formal but has higher numerical risk, especially around redundant coordinates and pseudoinverses.
  Date/Author: 2026-05-05 / Codex

- Decision: Treat frequency accuracy and interpretation accuracy as separate validation axes.
  Rationale: PED should improve mode interpretation and mixed-mode decomposition. It does not by itself correct raw calculated frequencies, scaling factors, phase effects, solvent shifts, adsorption effects, or local experimental/calculated mode correspondence.
  Date/Author: 2026-05-05 / Codex

- Decision: PED v1 uses the selected independent internal-coordinate basis by default.
  Rationale: The redundant coordinate pool contains aliases such as primitive and functional-group template rows. Using the selected basis avoids immediate double-counting while still preserving traceability to the Stage 3D basis diagnostics.
  Date/Author: 2026-05-05 / Codex

- Decision: Do not allow PED v1 to automatically override final Stage 3D assignment labels yet.
  Rationale: Benchmark review shows PED v1 is strongest as an explanatory evidence layer. Automatic overrides would be unsafe when mode correspondence is unresolved, when the PED top contributor is diffuse, or when the benchmark label is broader than ORCAVEDA's local-coordinate wording. Future overrides may be allowed only when PED has a strong diagnostic contributor, frequency/mode correspondence is acceptable or multiscale-resolved, and Stage 3D/PED disagreement is explicitly reported.
  Date/Author: 2026-05-05 / Codex

- Decision: Implement PED v2 as `force-aware normalized B-matrix projection`, not as full Wilson GF PED.
  Rationale: The ORCA Cartesian Hessian is available and can strengthen mode-specific internal-coordinate diagnostics. A full Wilson GF/VEDA-equivalent PED still needs a more formal internal-coordinate G/F treatment and stricter validation of redundant-coordinate handling. The v2 label therefore remains diagnostic and explicit about its limits.
  Date/Author: 2026-05-05 / Codex

- Decision: Implement Wilson PED as a separate `wilson_ped_audit` output instead of merging it into `ped_v2_force_audit`.
  Rationale: Wilson PED uses a different normalization convention and reports G/F diagnostics. Keeping it separate preserves reproducibility and allows v1, v2, and Wilson outputs to be compared directly.
  Date/Author: 2026-05-05 / Codex

- Decision: Let the frontend use PED as the primary displayed interpretation when PED output is available, but keep CSV `assignment_audit` as the Stage 3D audit and expose Stage 3D beside PED in viewer details.
  Rationale: The user-facing spectrum table should show the strongest current mode-composition evidence, while reproducibility still requires keeping Stage 3D and PED layers distinguishable. This avoids silently rewriting historical Stage 3D outputs.
  Date/Author: 2026-05-06 / Codex

## Outcomes & Retrospective

Initial implementation outcome: `src/ped.py` now computes PED v1 as normalized B-matrix projections onto each ORCA normal mode using `normal_modes[:, mode]`. `analyze_general_hess_files(...)` now emits a separate `ped_audit` table. This is additive and does not change Stage 3D assignments.

Validation outcomes:

- `.\.venv312\Scripts\python.exe -m py_compile src\ped.py src\ORCAVEDA_patched_stage3D_v5_0.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_ped_tmp` returned `3 passed` before the pipeline regression was added.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_ped_tmp2` returned `4 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_nist_ir_matching.py tests\test_nist_ir_compare.py -q --basetemp outputs\pytest_ped_core_tmp` returned `19 passed`.
- `$env:PYTHONPATH='C:\Users\unive\Documents\Projects\orcaveda\src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_ped_stage3d_tmp` returned `2 passed`.

Water probe outcome: running `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\ped_h2o_probe` generated `H2O__ped_audit.csv`. The positive-frequency rows show the bend mode dominated by `ang(H2-O1-H3)` at about 99.7% and the two stretch modes dominated by the two O-H stretches at about 50/50.

Benchmark review outcome: running ORCAVEDA on `H2O_freq.hess`, `NH3.hess`, `acetophenone.hess`, `benzoic_acid.hess`, `aniline.hess`, `phenol.hess`, and `pyridine.hess` generated `outputs/ped_benchmark_review/*__ped_audit.csv`. The raw-primary benchmark comparison over the requested molecules produced 62 rows: 7 PASS, 53 WARN, and 2 FAIL. The scaled-primary comparison with explicit scale factor 0.96 produced 62 rows: 8 PASS, 54 WARN, and 0 FAIL. This supports keeping frequency/mode correspondence separate from PED interpretation.

PED-aware 10-molecule validation outcome: running ORCAVEDA on water, ammonia, acetaldehyde, acetamide, acetophenone, aniline, benzene, benzoic acid, phenol, and pyridine generated `outputs/ped_aware_benchmark_10/*__ped_audit.csv`. The comparator was run with `--ped-audit`, `--windows-cm1 50,100,200,500`, `--scale-factor 0.96`, and both raw/scaled primary axes. Focused tests for PED and multiscale comparator returned `22 passed`.

PED v2 outcome: `src/ped.py` now includes `compute_ped_v2_force_aware(...)` and `build_ped_v2_force_audit_dataframe(...)`. `src/orca_parser.py` now parses `$hessian` into `HessData.cartesian_hessian`. `analyze_general_hess_files(...)` writes `ped_v2_force_audit` when the Hessian is available. On H2O, v2 preserves the expected bend/stretch identities while shifting force-aware percentages: the bend mode is about 92.7% H-O-H bend, and the two stretch modes remain about 50/50 O-H stretch.

PED v2 validation commands:

- `.\.venv312\Scripts\python.exe -m py_compile src\orca_parser.py src\orcaveda_models.py src\ped.py src\ORCAVEDA_patched_stage3D_v5_0.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_ped_v2_tmp` returned `6 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_golden_rdkit_outputs.py tests\test_nist_ir_matching.py tests\test_nist_ir_compare.py -q --basetemp outputs\pytest_ped_v2_core_tmp` returned `24 passed`.
- `$env:PYTHONPATH='C:\Users\unive\Documents\Projects\orcaveda\src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_ped_v2_stage3d_tmp` returned `2 passed`.

Wilson PED outcome: `src/ped.py` now includes `build_wilson_g_matrix(...)`, `reconstruct_internal_force_matrix(...)`, and `build_wilson_ped_audit_dataframe(...)`. `analyze_general_hess_files(...)` writes `wilson_ped_audit` when `$hessian` is available. The comparator accepts both `ped_rank` and `wilson_rank` schemas for PED-aware validation.

Wilson PED validation commands:

- `.\.venv312\Scripts\python.exe -m py_compile src\ped.py src\ORCAVEDA_patched_stage3D_v5_0.py` completed successfully.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_wilson_ped_tmp` returned `9 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_vibrational_assignment_multiscale.py tests\test_ped.py -q --basetemp outputs\pytest_wilson_ped_core_tmp` returned `12 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_golden_rdkit_outputs.py tests\test_nist_ir_matching.py tests\test_nist_ir_compare.py -q --basetemp outputs\pytest_wilson_ped_focus_tmp` returned `28 passed`.
- `$env:PYTHONPATH='C:\Users\unive\Documents\Projects\orcaveda\src'; .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_wilson_ped_stage3d_tmp` returned `2 passed`.

## Context and Orientation

The trusted baseline is Stage 3D v5.0. The main legacy entrypoint is `src/ORCAVEDA_patched_stage3D_v5_0.py`. Modular code lives mostly in `src/`.

Relevant files and responsibilities:

- `src/orca_parser.py`: parses ORCA `.hess` content, including atoms, coordinates, frequencies, normal modes, and Hessian-like matrices where available.
- `src/internal_coordinates.py`: constructs internal coordinates such as bonds, angles, torsions, and related geometric descriptors.
- `src/mode_assignment.py`: builds Stage 3D assignment audit rows and final vibrational interpretation labels.
- `src/scale_factor_engine.py`: applies explicit frequency scaling and must remain conceptually separate from PED.
- `benchmarks/vibrational_assignments/`: stores literature/NIST benchmark assignments and comparator tooling for interpretation validation.
- `tests/`: contains regression, chemistry, Stage 3D, NIST, and benchmark comparison tests.

Important project boundaries:

- Stage 3D may be called a geometric and weighted independent-coordinate assignment audit or PED-like audit, but not strict VEDA PED or full Wilson GF PED.
- ORCA normal modes must preserve the orientation `normal_modes[:, mode]`.
- High-frequency C-H, N-H, and O-H modes must not remain unassigned without diagnostics.
- Parser support for single-column ORCA block matrix headers must not regress.

Working terms:

- Internal coordinate: a chemically meaningful coordinate such as bond stretch, angle bend, torsion, out-of-plane bend, or ring deformation proxy.
- B-matrix: derivatives of internal coordinates with respect to Cartesian coordinates.
- PED contribution: a normalized contribution of an internal coordinate or internal-coordinate family to a normal mode under the chosen PED v1 definition.
- Mode correspondence: selecting which calculated mode corresponds to an experimental or benchmark peak. This remains separate from the PED calculation itself.

## Plan of Work

First, audit the current parser and internal-coordinate code. Confirm what data are available from `.hess`: atomic numbers, coordinates, masses, frequencies, normal modes, Hessian or mass-weighted Hessian, and any current normalization assumptions. Verify normal-mode orientation from code and tests.

Second, write a short mathematical specification for PED v1 before patching. This must define whether normal modes are treated as Cartesian displacement vectors or mass-weighted displacement vectors, how internal-coordinate derivatives are computed, how contributions are normalized, how redundant coordinates are handled, and what diagnostics are emitted when a mode cannot be decomposed reliably.

PED v1 specification selected for the first implementation:

- Normal-mode vector convention: `mode_vec = hess.normal_modes[:, mode]`.
- Normal-mode normalization: divide each mode vector by its Euclidean norm.
- Internal-coordinate derivative method: use existing finite-difference B-matrix rows from `src/b_matrix.py`.
- Coordinate-row normalization: divide each B row by its Euclidean norm before projection.
- Projection: `projection_i = dot(B_unit_i, mode_unit)`.
- Contribution weight: `weight_i = projection_i ** 2`.
- Percent normalization: `percent_i = 100 * weight_i / sum(weight)`.
- Coordinate basis: use selected independent internal-coordinate indices by default.
- Diagnostics: report zero/invalid normal mode vectors, no internal-coordinate basis, no valid B rows, zero projection weight, and diffuse top contributions.
- Boundary: this is not force-constant Wilson GF PED because no force-constant matrix is used.

Third, implement a new module, likely `src/ped.py`, with minimal public functions:

- `build_b_matrix(...)`
- `project_mode_to_internal_coordinates(...)`
- `compute_ped(...)`
- `summarize_ped_mode(...)`

The implementation should use existing internal-coordinate generation where possible. It should not silently alter Stage 3D outputs.

Fourth, add output generation. The first stable artifact should be a separate `ped_audit.csv` or `ped_summary.csv`, plus optional JSON metadata describing coordinate counts, normalization, rejected coordinates, and numerical warnings. If later integration adds PED columns to `assignment_audit.csv`, that schema change must be explicit and tested.

Fifth, validate on small molecules before complex molecules. Use water and ammonia as sanity checks for stretch/bend behavior, then benzene for ring/mixed/degenerate modes, then acetophenone and benzoic acid for cases where benchmark comparison already showed mode-correspondence and mixed-mode interpretation issues.

Sixth, compare PED assignments to the benchmark comparator. A PED success is not a smaller frequency delta. A PED success is better agreement of functional group, motion family, and mixed-mode composition once the mode correspondence problem is handled.

Seventh, decide how PED should affect final assignment wording. The first production-safe option is to report Stage 3D and PED side by side. A later combined assignment can let PED reinforce, warn against, or override Stage 3D labels only when thresholds and tests are explicit.

## Concrete Steps

1. Inspect current data flow.

       git status --short
       rg "normal_modes" src tests
       rg "internal" src/internal_coordinates.py src/mode_assignment.py tests
       rg "hessian|mass|frequency|freq" src/orca_parser.py src tests

   Expected result: a concise note in this plan listing available arrays, shapes, units, and assumptions.

2. Confirm parser and Stage 3D baseline still pass focused tests before PED edits.

       .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q
       .\.venv312\Scripts\python.exe -m pytest tests\test_golden_rdkit_outputs.py tests\test_nist_ir_matching.py tests\test_nist_ir_compare.py -q

   If any command fails for unrelated existing reasons, record exact failure and use focused substitutes.

3. Draft PED v1 specification in this file under `Artifacts and Notes`.

   The spec must state:

   - normal-mode vector convention;
   - mass-weighting convention;
   - internal-coordinate derivative method;
   - normalization equation;
   - redundant-coordinate handling;
   - failure diagnostics.

4. Implement `src/ped.py` with unit-tested low-level geometry derivatives.

   Initial tests should cover:

   - bond stretch derivative direction;
   - angle bend finite-difference consistency;
   - torsion derivative finite-difference consistency or explicit diagnostic if torsion support is deferred;
   - contribution normalization per mode.

5. Integrate PED output into the ORCAVEDA pipeline as an optional/additive artifact.

   The first CLI or pipeline behavior should be explicit, for example a flag or output file creation that does not change final assignment labels by default.

6. Validate on selected `.hess` files.

       .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\ped_h2o_probe
       .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\NH3.hess --outdir outputs\ped_nh3_probe
       .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\acetophenone.hess --outdir outputs\ped_acetophenone_probe

   Adjust exact commands after confirming the current entrypoint and CLI options.

7. Compare PED outputs against benchmark assignments.

       .\.venv312\Scripts\python.exe benchmarks\vibrational_assignments\compare_orcaveda_assignments.py --audit <assignment_audit.csv> --out <comparison.csv>

   Extend comparator only after PED output schema is stable.

## Validation and Acceptance

The implementation is acceptable only when all of the following are true:

- PED is implemented in source code, not just described in comments or reports.
- Stage 3D remains available and is not mislabeled as full PED.
- PED output includes per-mode top internal-coordinate contributors and normalized percentages or explicitly named normalized weights.
- Output metadata states the PED definition used by ORCAVEDA.
- At least water and ammonia have focused tests showing chemically sensible dominant coordinates.
- At least one mixed-mode benchmark case, preferably acetophenone C=O or benzoic acid carboxylic acid/ring coupling, is used to show whether PED improves interpretation.
- Tests that touch parser, internal coordinates, mode assignment, and PED run and their exact results are recorded.

A strong validation target is:

- no regression in existing Stage 3D/NIST/multiscale benchmark tests;
- PED top contributor agrees with simple benchmark assignments for water and ammonia;
- PED provides explicit mixed-mode decomposition for acetophenone, benzoic acid, aniline, phenol, or pyridine without pretending to solve frequency scaling.

Current acceptance status:

- PED is implemented in `src/ped.py`.
- Stage 3D remains separate and is not relabeled as full PED.
- `ped_audit` reports per-mode ranked contributors and normalized percentages.
- Water has a pipeline regression test and generated probe output.
- Ammonia and complex benchmark molecules still need dedicated PED interpretation review.

## Idempotence and Recovery

The parser audit, test commands, and benchmark comparisons are safe to repeat. Generated outputs should go under `outputs/` with unique descriptive folders such as `outputs/ped_h2o_probe`.

If the PED prototype produces unstable percentages because of redundant internal coordinates, keep the prototype behind a clearly named experimental function or flag and do not integrate it into final assignment labels.

If tests fail after adding PED, first isolate whether the failure is in the new PED module or an unintended Stage 3D regression. Revert only the local PED changes if necessary; do not revert unrelated benchmark or user changes in the working tree.

If force-constant PED proves necessary, record that as a new decision and add a v2 milestone rather than expanding v1 silently.

## Artifacts and Notes

Initial agent/task split:

- Sci-agent owns the PED mathematical definition, scientific boundary language, and interpretation-risk review.
- Backend-agent owns parser inspection, module design, implementation, and tests.
- Sci-agent plus Backend-agent jointly review shapes, units, normal-mode orientation, and benchmark outcomes.
- Frontend-agent is deferred until PED output is stable; later it can expose PED percentages in the HTML viewer.

Initial expected output artifacts:

- `src/ped.py`
- `tests/test_ped.py`
- `ped_audit.csv` or `ped_summary.csv` in generated output folders
- optional `ped_summary.json` metadata
- updated benchmark comparison only after the PED schema is stable

Generated PED benchmark review artifacts:

- `outputs/ped_benchmark_review/stage3d_benchmark_comparison.csv`
- `outputs/ped_benchmark_review/stage3d_benchmark_comparison_scaled_primary.csv`
- `outputs/ped_benchmark_review/ped_benchmark_summary.csv`
- `outputs/ped_benchmark_review/ped_problem_cases_summary.csv`
- `outputs/ped_benchmark_review/ped_key_cases_compact.csv`
- `outputs/ped_benchmark_review/ped_key_cases_compact.md`

Generated PED-aware 10-molecule validation artifacts:

- `outputs/ped_aware_benchmark_10/ped_aware_comparison_raw.csv`
- `outputs/ped_aware_benchmark_10/ped_aware_comparison_scaled.csv`
- `outputs/ped_aware_benchmark_10/ped_aware_summary_raw.csv`
- `outputs/ped_aware_benchmark_10/ped_aware_summary_scaled.csv`
- `outputs/ped_aware_benchmark_10/ped_aware_key_cases_raw.csv`

Generated PED v2 artifacts:

- `outputs/ped_v2_h2o_probe/H2O__ped_v2_force_audit.csv`
- `outputs/ped_v2_benchmark_10/ped_v2_aware_comparison_raw.csv`
- `outputs/ped_v2_benchmark_10/ped_v2_aware_comparison_scaled.csv`
- `outputs/ped_v2_benchmark_10/ped_v2_summary_raw.csv`
- `outputs/ped_v2_benchmark_10/ped_v2_summary_scaled.csv`

Generated Wilson PED artifacts:

- `outputs/wilson_ped_h2o_probe/H2O__wilson_ped_audit.csv`
- `outputs/wilson_ped_benchmark_10/wilson_ped_comparison_raw.csv`
- `outputs/wilson_ped_benchmark_10/wilson_ped_comparison_scaled.csv`
- `outputs/wilson_ped_benchmark_10/wilson_ped_summary_raw.csv`
- `outputs/wilson_ped_benchmark_10/wilson_ped_summary_scaled.csv`
- `outputs/wilson_ped_benchmark_10/wilson_ped_result_table_by_molecule.csv`
- `outputs/wilson_ped_benchmark_10/wilson_ped_key_result_table.csv`

Generated PED frontend artifacts:

- `outputs/ped_frontend_monoethanolamine/monoethanolamine_DFT_therm__interactive_spectrum.html`
- `outputs/ped_frontend_monoethanolamine/monoethanolamine_DFT_therm__spectrum_data.json`
- `outputs/ped_frontend_monoethanolamine/monoethanolamine_DFT_therm__wilson_ped_audit.csv`

Open scientific questions for the first implementation pass:

- Are ORCA normal modes in the current parser Cartesian displacements, mass-weighted displacements, or already normalized in a way requiring conversion?
- Does the `.hess` parser expose enough Hessian/force-constant data for force-constant PED, or should v1 deliberately avoid that claim?
- Which internal-coordinate set gives the most stable and chemically interpretable decomposition without excessive redundancy?
- Should torsions be included in v1 or introduced after bond/angle/out-of-plane validation?
- What threshold should separate primary, secondary, and minor PED contributors in assignment wording?

## Interfaces and Dependencies

The new PED module should depend only on existing ORCAVEDA data structures, NumPy, and standard-library code unless a dependency is explicitly justified and added to `requirements.txt`.

Potential function contracts:

- `compute_ped(atoms, coordinates, normal_modes, internal_coordinates, masses=None, options=None) -> list[PEDModeResult]`
- `PEDModeResult.mode_index`
- `PEDModeResult.frequency_cm1`
- `PEDModeResult.contributions`
- `PEDModeResult.normalization_sum`
- `PEDModeResult.warnings`

Potential CSV columns:

- `mode`
- `frequency_cm-1`
- `ped_rank`
- `internal_coordinate`
- `coordinate_family`
- `atoms`
- `contribution_percent`
- `ped_warning`

The exact schema must be finalized before implementation and covered by tests. Units must remain `cm-1` for frequencies, Angstrom-based coordinates for geometry where applicable, and explicit dimensionless percentages for normalized PED contributions.
