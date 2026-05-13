# VEDA-Like Composed Coordinate PED Optimization

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA already parses ORCA `.hess` files, builds primitive and functional-group internal coordinates, selects a rank-preserving independent basis, computes Stage 3D assignment audits, writes Wilson GF-style PED diagnostics, and now has a conservative EPM-like basis optimizer that swaps existing coordinates to improve PED localization. The next scientific upgrade is to implement the central mathematical idea described by Jamroz 2013 for VEDA: optimize the internal-coordinate description itself by constructing composed local coordinates, such as sums and differences of same-type stretches, bends, torsions, and out-of-plane coordinates.

After this work, a user should be able to run ORCAVEDA on a `.hess` file and receive, alongside the existing Stage 3D and Wilson/PED outputs, a composed-coordinate PED optimization audit. The audit should show which composed coordinates were created, whether the rank was preserved, how the EPM-like score changed, how many diffuse PED modes improved, and which final assignments became better supported by PED. This must improve scientific interpretation without claiming full VEDA equivalence until the implemented mathematics and validation justify that wording.

## Progress

- [x] (2026-05-06 16:31+05:00) Existing project agent roles, ORCAVEDA skills, current PED code, and Jamroz 2013 text extraction were reviewed.
- [x] (2026-05-06 16:31+05:00) Initial ExecPlan written before starting implementation of composed-coordinate optimization.
- [x] (2026-05-07 09:53+05:00) Implement initial composed internal-coordinate data model and B-row composition invariant helper.
- [x] (2026-05-07 09:57+05:00) Implement same-type and X-H/heavy coordinate grouping helper.
- [x] (2026-05-07 10:02+05:00) Implement first narrow X-H pair symmetric/asymmetric stretch candidate generator.
- [x] (2026-05-07 10:06+05:00) Implement composed PED candidate B-matrix builder helper.
- [x] (2026-05-07 10:09+05:00) Implement rank-preserving composed PED basis selection helper.
- [x] (2026-05-07 10:15+05:00) Integrate composed PED basis as an experimental PED-only diagnostics layer without changing `assignment_audit`.
- [x] (2026-05-07 10:20+05:00) Add separate composed PED/PED v2/Wilson audit tables without changing Stage 3D baseline.
- [x] (2026-05-07 10:25+05:00) Run full-golden composed PED audit diagnostics over all local `data/hess/*.hess`.
- [x] (2026-05-07 10:34+05:00) Expose composed PED evidence in the interactive viewer as a separate selectable/visible evidence layer.
- [x] (2026-05-07 11:05+05:00) Add benchmark comparator support for composed PED audit as a separate diagnostic source and run policy review.
- [x] (2026-05-07 11:43+05:00) Harden composed X-H stretch labels and composed benchmark context summaries; rerun focused tests and raw/scaled subset comparator.
- [x] (2026-05-07 11:55+05:00) Repeat full-golden composed benchmark comparison after semantic hardening; raw and scaled comparisons had zero `worsens_semantic_match` rows.
- [x] (2026-05-07 12:18+05:00) Add conservative composed PED diagnostic policy fields to viewer payload and benchmark comparator without changing final assignment policy.
- [x] (2026-05-07 12:39+05:00) Export conservative composed PED diagnostic policy as `composed_ped_policy_diagnostics.csv` and validate full-golden behavior.
- [x] (2026-05-07 13:04+05:00) Add viewer peak-table filter and compact indicator for composed diagnostic hints; validate focused tests and full-golden smoke outputs.
- [x] (2026-05-07 13:32+05:00) Harden primitive torsion labels so terminal C-H/N-H/O-H neighbors drive torsion family semantics; rerun full-golden audit.
- [x] (2026-05-12 14:25+05:00) Add diagnostic triage categories for remaining composed PED hints and rerun full-golden audit.
- [x] (2026-05-12 15:06+05:00) Inspect the 8 `motion_family_coordinate_generation_target` rows and 3 high-frequency ethylene oxide rows; classify ethylene oxide as X-H stretch recovery rather than unresolved high-frequency conflict.
- [x] (2026-05-12 15:31+05:00) Add top-contributor provenance fields to composed PED policy diagnostics and verify the 8 DMF/ethene target rows are primitive-row optimizer substitutions.
- [x] (2026-05-12 15:47+05:00) Split primitive-row optimizer substitutions into their own composed PED triage category.
- [x] (2026-05-12 16:02+05:00) Add stable `composed_ped_evidence_origin` field for future optimizer constraints.
- [x] (2026-05-12 16:30+05:00) Add warning-only diagnostic for primitive-row optimizer substitution evidence without changing final assignments.
- [x] (2026-05-12 16:42+05:00) Run focused tests, subset check, and full-golden acceptance for the warning-only diagnostic.
- [x] (2026-05-12 17:10+05:00) Add an opt-in narrow primitive-substitution constraint experiment and validate the DMF/ethene/ethylene oxide subset.
- [x] (2026-05-13 10:36+05:00) Run full-golden acceptance for the opt-in primitive-substitution constraint; rank was preserved, final labels were unchanged, and primitive-substitution warnings were removed.
- [x] (2026-05-13 10:45+05:00) Decide to keep the primitive-substitution constraint opt-in for now; default promotion needs a separate explicit plan because confirmation candidates decreased.
- [ ] (next) Plan the next scientific expansion: broader composed-coordinate generators or a default-promotion validation plan, without changing final assignment policy silently.

## Surprises & Discoveries

- Observation: Existing project roles are sufficient for this upgrade. Backend-agent owns code paths and schema changes, Sci-agent owns mathematical and chemical validity, and Frontend-agent is only needed after backend evidence exists.
  Evidence: `AGENTS.md` defines Backend-agent, Frontend-agent, and Sci-agent with the necessary ownership boundaries.

- Observation: ORCAVEDA already has the first safe precursor to VEDA-like optimization: `optimize_independent_coordinates_for_ped(...)` in `src/b_matrix.py` swaps existing independent coordinates to improve PED localization while preserving rank.
  Evidence: Local source inspection of `src/b_matrix.py` and previous full golden run output in `outputs/ped_basis_optimizer_full_golden_20260506`.

- Observation: Jamroz 2013 describes a stronger method than the current ORCAVEDA optimizer. VEDA constructs composed local coordinates from sums or differences, usually preserving coordinate type and separating CH-like motions from heavy-atom motions.
  Evidence: Extracted text in `outputs/jamroz_2013_text.txt`, especially the sections "Different ways for optimization of PED analysis" and "Block diagram".

- Observation: On the current curated benchmark rows, composed Wilson evidence did not improve semantic status relative to `ped_final_assignment`, even though it improved localization for one acetophenone carbonyl row. Several rows worsened as a separate evidence source because composed top terms lost explicit C-H or O-H semantic classes.
  Evidence: `outputs/composed_ped_benchmark_20260507/composed_ped_comparison_raw.csv` and `outputs/composed_ped_benchmark_20260507/composed_ped_comparison_scaled.csv`.

- Observation: The initial worsened rows had two separable causes. Composed X-H stretch coordinates were labeled as generic `bond stretch`, losing C-H/N-H/O-H classes, and the composed benchmark summary used only four top contributors, dropping weak but important carboxylic-context contributors such as `C-C-O bend` or `O-H torsion`.
  Evidence: `outputs/composed_ped_benchmark_20260507/composed_ped_comparison_raw.csv`, `outputs/composed_ped_benchmark_20260507/composed_ped_comparison_scaled.csv`, and focused inspection of `composed_ped_top_contributors` for acetaldehyde, acetamide, benzoic acid, and aniline rows.

- Observation: The DMF/ethene `C-H torsion` composed disagreements are not caused by generated composed torsion coordinates. The first-wave composed generator creates only X-H stretch sum/difference candidates. The dominant C-H torsions in the composed Wilson audit are primitive torsion rows selected by the existing EPM-like basis optimizer after composed candidates are added.
  Evidence: `src/internal_coordinates.py` rejects composed torsion components in `make_composed_internal_coordinate(...)` and only `build_xh_pair_composed_stretch_candidates(...)` is integrated for first-wave composed generation; `outputs/composed_xh_recovery_triage_full_golden_20260512\*__composed_wilson_ped_audit.csv` shows primitive `tor(...)` rows as the DMF/ethene top contributors.

## Decision Log

- Decision: Do not create new permanent project agent roles for this upgrade.
  Rationale: The current `AGENTS.md` roles already cover implementation, scientific review, and viewer review. Adding a new role now would add process overhead without increasing correctness. A future "Benchmark-agent" can be considered only if curated experimental benchmark maintenance becomes a distinct recurring task.
  Date/Author: 2026-05-06 / Codex

- Decision: Use existing ORCAVEDA skills rather than creating new skills before implementation.
  Rationale: `$orcaveda-core` covers ORCA `.hess`, B-matrix, internal coordinates, PED, Stage 3D, and regression safety. `$orcaveda-nist-ir` will be used later only if validation touches NIST matching. `$orcaveda-viewer` will be used only when viewer wording or UI changes are made.
  Date/Author: 2026-05-06 / Codex

- Decision: Implement composed-coordinate optimization as an additive PED/Wilson layer, not as a replacement for Stage 3D.
  Rationale: Stage 3D is the trusted baseline geometric assignment audit. The composed-coordinate optimizer changes scientific interpretation and must be introduced with diagnostics, rank checks, and benchmark evidence before influencing final assignment labels more aggressively.
  Date/Author: 2026-05-06 / Codex

- Decision: The first production implementation will use deterministic same-type sum/difference composed coordinates with coefficients limited to `+1` and `-1`.
  Rationale: Jamroz 2013 notes that VEDA can use non-unit coefficients, but unit coefficients are more transparent, easier to test, and safer for a first reproducible implementation. Non-unit coefficient optimization remains a later milestone.
  Date/Author: 2026-05-06 / Codex

- Decision: Optimization metrics will be reported as ORCAVEDA EPM-like diagnostics, not literal VEDA `<EPm>`, until the full PED matrix orientation and coordinate/mode indexing semantics are aligned and tested.
  Rationale: Jamroz describes EPM as the arithmetic average of maximal elements of the PED matrix. ORCAVEDA currently has a mode-centric localization score and Wilson energy terms, but not the full VEDA coordinate-construction workflow.
  Date/Author: 2026-05-06 / Codex

- Decision: Start composed-coordinate implementation with metadata, a manual composed-coordinate factory, and a pure B-row composition helper.
  Rationale: This proves the core invariant that a composed coordinate's B row is the coefficient-weighted sum of primitive B rows without changing Stage 3D assignment behavior, automatic coordinate generation, rank selection, output schemas, or scientific thresholds.
  Date/Author: 2026-05-07 / Codex

- Decision: Name the first hydrogen-containing grouping class `XH_like`, not `CH_like`.
  Rationale: The initial grouping deliberately separates all hydrogen-involving motions, including C-H, O-H, and N-H, from heavy-atom-only motions. Calling this group `CH_like` would obscure O-H and N-H semantics, while `XH_like` preserves the intended safety boundary for future composed-coordinate generation.
  Date/Author: 2026-05-07 / Codex

- Decision: Keep composed-coordinate evidence as viewer/diagnostic evidence only after the first benchmark comparator run.
  Rationale: The benchmark comparator found no semantic improvements over the baseline `ped_final_assignment` source on the current curated benchmark rows. The composed layer can confirm existing PASS rows and expose localization changes, but the evidence does not yet justify warning policy changes or automatic fallback into final assignments.
  Date/Author: 2026-05-07 / Codex

- Decision: Fix composed semantic hardening only in evidence labels and benchmark diagnostics, not final assignment policy.
  Rationale: Mapping `composed_*_XH_stretch(center:H,H)` to C-H/N-H/O-H stretch preserves explicit X-H chemistry without changing the B matrix, basis selection, or assignment policy. Using six composed benchmark contributors matches the existing final-assignment evidence depth closely enough to avoid false semantic regressions from truncated diagnostic summaries.
  Date/Author: 2026-05-07 / Codex

- Decision: Add composed PED policy hints as diagnostic fields only.
  Rationale: Full-golden evidence supports using composed PED as a conservative confirmation/warning candidate, but not as an automatic final-label fallback. The implemented fields (`composed_ped_policy_hint`, `composed_ped_localization_delta_percent`, `composed_ped_semantic_status`, and `composed_ped_semantic_reason`) are surfaced in viewer payloads and benchmark comparison outputs while preserving `assignment_audit`, `ped_stage3d_agreement`, `ped_final_assignment`, B matrices, and rank selection.
  Date/Author: 2026-05-07 / Codex

- Decision: Do not restrict the composed PED optimizer candidate pool globally to composed coordinates only.
  Rationale: A local experiment on DMF/ethene/ethylene oxide removed the original DMF/ethene torsion disagreements, but the full-golden run broadened composed/baseline disagreements from 19 to 181 rows. This failed the minimal safe patch criterion, so the optimizer-scope change was reverted. The safer retained change is diagnostic triage only.
  Date/Author: 2026-05-12 / Codex

- Decision: Keep the primitive-substitution constraint opt-in after full-golden acceptance; do not make it default in this plan.
  Rationale: The opt-in repair passed the planned safety checks: zero rank loss, zero high-frequency unassigned final assignments, no final-label changes, and zero primitive-substitution warnings. However, it reduced `confirmation_candidate` rows from 140 to 126, so default promotion is a policy and validation decision rather than a mechanical acceptance step.
  Date/Author: 2026-05-13 / Codex

## Outcomes & Retrospective

Initial composed-coordinate model outcome: `src/orcaveda_models.py` now adds `ComposedCoordinateTerm` and optional trailing composition metadata on `InternalCoordinate`. `src/internal_coordinates.py` now has `make_composed_internal_coordinate(...)` for manual composed-coordinate construction from existing primitive coordinates. `src/b_matrix.py` now has `compose_b_row(...)`, which builds a coefficient-weighted B row from already-computed primitive B rows. This is additive and does not change Stage 3D assignment behavior, automatic coordinate generation, rank selection, output schemas, or PED final-label policy.

Initial composed-coordinate validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_ped_baseline_tmp` returned `30 passed, 1 skipped` before edits.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_brow_tmp` returned `13 passed` after the focused composed-coordinate patch.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_ped_after_patch_tmp` returned `31 passed, 1 skipped` after the patch.

Current limitation after the first model patch: automatic composed-coordinate generation, grouping rules, rank-preserving composed PED selection, composed-coordinate output tables, and EPM-like before/after reporting remained not implemented in source.

Initial coordinate-grouping outcome: `src/internal_coordinates.py` now has `classify_coordinate_optimization_group(...)`, which classifies internal coordinates by motion family (`stretch`, `bend`, `torsion`, or `other`) and by whether the coordinate contains hydrogen (`XH_like`) or only heavy atoms (`heavy`). This is a pure helper for future composed-coordinate candidate generation and is not connected to the pipeline or output tables.

Coordinate-grouping validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_coordinate_grouping_tmp` returned `14 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_coordinate_grouping_baseline_tmp` returned `32 passed, 1 skipped`.

Current limitation: automatic composed-coordinate generation, rank-preserving composed PED selection, composed-coordinate output tables, and EPM-like before/after reporting remain not implemented in source.

Initial X-H pair generator outcome: `src/internal_coordinates.py` now has `build_xh_pair_composed_stretch_candidates(...)`. It creates symmetric and asymmetric sum/difference candidates only for two X-H stretch coordinates sharing the same heavy atom and having distinct hydrogens. Duplicate primitive X-H stretches for the same atom pair are resolved deterministically by priority, name, and index. The helper returns candidate `InternalCoordinate` objects but is not connected to Stage 3D, PED basis selection, pipeline outputs, or CSV schemas.

X-H pair generator validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_xh_pair_generator_tmp` returned `15 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_xh_pair_generator_baseline_tmp` returned `33 passed, 1 skipped`.

Current limitation: the first generator is available as a pure helper only. Rank-preserving composed PED selection, composed-coordinate output tables, EPM-like before/after reporting, and full-golden behavior remain not implemented in source.

Initial composed candidate B-matrix outcome: `src/b_matrix.py` now has `build_composed_candidate_b_matrix(...)`. It appends composed-coordinate rows to an existing primitive B matrix, returns primitive internals followed by composed internals, and reports diagnostic counts including generation-rule counts. It is side-effect free and does not perform rank selection, modify Stage 3D internals, or write pipeline outputs.

Composed candidate B-matrix validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_candidate_b_matrix_tmp` returned `16 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_candidate_b_matrix_baseline_tmp` returned `34 passed, 1 skipped`.

Current limitation: rank-preserving composed PED selection, composed-coordinate output tables, EPM-like before/after reporting, and full-golden behavior remain not implemented in source.

Initial rank-preserving composed PED basis outcome: `src/b_matrix.py` now has `select_rank_preserving_composed_ped_basis(...)`. It starts from an accepted primitive basis, computes the starting rank, runs the existing EPM-like PED optimizer over the primitive plus composed candidate matrix, and rejects rank loss. The helper reports required rank, starting and optimized condition diagnostics, selected composed candidate indices, and composed selected count. It is still a pure helper and is not connected to Stage 3D, pipeline outputs, or CSV schemas.

Rank-preserving composed PED basis validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_rank_preserving_composed_basis_tmp` returned `17 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_rank_preserving_composed_basis_baseline_tmp` returned `35 passed, 1 skipped`.

Current limitation: this selection helper has only synthetic/H2O-style unit coverage. Pipeline integration, composed-coordinate output tables, EPM-like report fields, full-golden behavior, and any effect on final assignment policy remain not implemented in source.

Initial pipeline diagnostics outcome: `src/ORCAVEDA_patched_stage3D_v5_0.py` now builds X-H pair composed candidates, constructs a primitive plus composed candidate B matrix, and runs rank-preserving composed PED basis selection for each analyzed `.hess` file. The result is written only to new diagnostic tables: `composed_ped_basis_diagnostics` and `composed_ped_basis`. Existing `assignment_audit`, `ped_audit`, `ped_v2_force_audit`, `wilson_ped_audit`, `ped_stage3d_agreement`, and `ped_final_assignment` continue to use the existing basis behavior at this milestone.

Pipeline diagnostics validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_pipeline_tmp` returned `17 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_pipeline_baseline_tmp` returned `35 passed, 1 skipped`.

Current limitation: composed basis selection is visible as diagnostics only. It is not yet used to build PED, force-aware PED, or Wilson PED audit rows, and it does not affect final assignment policy. Full-golden behavior remains not run for this milestone.

Initial composed PED audit outcome: `src/ORCAVEDA_patched_stage3D_v5_0.py` now writes separate experimental composed-coordinate audit tables: `composed_ped_audit`, `composed_ped_v2_force_audit`, and `composed_wilson_ped_audit`. These tables are built from the primitive plus composed candidate B matrix and the rank-preserving composed PED basis selection. The existing `ped_audit`, `ped_v2_force_audit`, `wilson_ped_audit`, `assignment_audit`, `ped_stage3d_agreement`, and `ped_final_assignment` remain on their prior basis behavior.

Composed PED audit validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_ped_audits_tmp` returned `17 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_ped_audits_baseline_tmp` returned `35 passed, 1 skipped`.

Current limitation: composed audit tables are not yet consumed by `ped_stage3d_agreement`, `ped_final_assignment`, viewer payloads, or benchmark comparators. Full-golden behavior remains not run for this milestone.

Full-golden composed diagnostics outcome: running the current pipeline on all local `data/hess/*.hess` completed successfully and wrote outputs under `outputs/composed_ped_audits_full_golden_20260507`. The input set contained 55 `.hess` files, and `source_map.csv` plus `composed_ped_basis_diagnostics.csv` each contained 55 rows, so all input files were represented in the generated tables. The aggregated output prefix was `acetaldehyde__acetamide__acetanilide__plus_52_files__multi_file_55`.

Full-golden command:

- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_ped_audits_full_golden_20260507`, where `@hess` was the PowerShell-expanded sorted list from `Get-ChildItem -LiteralPath 'data\hess' -Filter '*.hess'`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.

Full-golden audit results:

- `composed_candidate_count` total: 366.
- `composed_selected_count` total: 166.
- Rank losses: 0 rows; every row preserved `optimized_rank >= required_rank`.
- Baseline high-frequency rows at or above 2500 cm-1 in `assignment_audit`: 350.
- Baseline high-frequency unassigned rows: 0.
- Rows with improved `optimized_mean_top_percent`: 52 of 55.
- Rows with improved `optimized_mean_top_percent` and at least one selected composed row: 43 of 55.
- Rows improved without selected composed rows, reflecting primitive-only optimization in the same helper path: `benzaldimine.hess`, `benzaldoxime.hess`, `benzene.hess`, `benzoic_acid.hess`, `benzonitrile.hess`, `nitrobenzene.hess`, `phenol.hess`, `pyridine.hess`, and `pyrrole.hess`.
- Rows without improvement: `ethyne.hess`, `H2O2_freq.hess`, and `phenyl_isocyanate.hess`.
- Files with no generated X-H pair candidates: `benzaldimine.hess`, `benzaldoxime.hess`, `benzene.hess`, `benzoic_acid.hess`, `benzonitrile.hess`, `ethyne.hess`, `H2O2_freq.hess`, `nitrobenzene.hess`, `phenol.hess`, `phenyl_isocyanate.hess`, `pyridine.hess`, and `pyrrole.hess`.

Full-golden verdict: PASS for diagnostic integration. The composed-coordinate layer is safe as a separate evidence surface on this local corpus. It is not yet validated as an input to `ped_stage3d_agreement`, `ped_final_assignment`, viewer payloads, or benchmark comparators.

Initial viewer evidence-layer outcome: `src/reports.py` now accepts `composed_wilson_ped_audit`, `composed_ped_v2_force_audit`, and `composed_ped_audit` in `build_spectrum_payload(...)`. The payload carries composed PED interpretation fields separately from baseline PED fields. The interactive HTML viewer adds an `Evidence Layer` selector in the selected-mode card and displays both baseline PED evidence and composed PED evidence without changing `final_assignment`, `assignment`, `ped_stage3d_agreement`, or `ped_final_assignment` policy. `src/ORCAVEDA_patched_stage3D_v5_0.py` passes the composed audit tables into the spectrum payload.

Viewer evidence-layer validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py -q --basetemp outputs\pytest_composed_viewer_tmp` returned `6 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_viewer_ped_tmp` returned `17 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_viewer_baseline_tmp` returned `35 passed, 1 skipped`.

Current limitation: visual browser inspection of the generated HTML was not run at this milestone. The viewer tests verify HTML text and JSON payload fields, but not pixel layout or manual interaction in a browser.

Initial benchmark comparator outcome: `benchmarks/vibrational_assignments/compare_orcaveda_assignments.py` now accepts `--composed-ped-audit` and writes separate composed-coordinate diagnostic columns next to the existing baseline PED and `ped_final_assignment` benchmark columns. The new columns include `composed_ped_top_contributors`, `composed_ped_top_family`, `composed_ped_top_percent`, `composed_ped_classes`, `composed_ped_semantic_status`, `composed_ped_semantic_reason`, `stage3d_composed_ped_overlap_classes`, `stage3d_composed_ped_warning`, `composed_vs_baseline_ped_status`, `composed_ped_localization_delta_percent`, and `composed_ped_policy_hint`. These columns do not change the legacy `status`, `reason`, `orcaveda_assignment`, or `ped_final_assignment` comparison path.

Benchmark comparator validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_vibrational_assignment_multiscale.py -q --basetemp outputs\pytest_composed_comparator_tmp` returned `6 passed`.
- Raw benchmark command with `--ped-final-assignment` and `--composed-ped-audit outputs\composed_ped_audits_full_golden_20260507\acetaldehyde__acetamide__acetanilide__plus_52_files__multi_file_55__composed_wilson_ped_audit.csv` wrote `outputs\composed_ped_benchmark_20260507\composed_ped_comparison_raw.csv` and reported `27 PASS` and `66 WARN` for the primary baseline comparison.
- Scaled benchmark command with the same composed audit source wrote `outputs\composed_ped_benchmark_20260507\composed_ped_comparison_scaled.csv` and reported `26 PASS` and `67 WARN` for the primary baseline comparison.

Benchmark policy review:

- Raw comparison: baseline PED semantic status was `32 PASS` and `61 WARN`; composed PED semantic status was `28 PASS`, `61 WARN`, and `4 FAIL`. `composed_vs_baseline_ped_status` reported `87 same_semantic_match`, `5 worsens_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 improves_semantic_match`.
- Scaled comparison: baseline PED semantic status was `32 PASS` and `61 WARN`; composed PED semantic status was `27 PASS`, `61 WARN`, and `5 FAIL`. `composed_vs_baseline_ped_status` reported `86 same_semantic_match`, `6 worsens_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 improves_semantic_match`.
- The only localization-only gain in both raw and scaled outputs was acetophenone near 1683 cm-1, where composed evidence raised the top contributor from baseline `C-C-C bend` at 45.933% to composed `C=O stretch` at 63.043354%, while both sources remained semantic PASS.
- Worsened composed evidence rows included acetaldehyde C-H bend rows, acetamide C-H stretch rows, benzoic-acid O-H/context rows, and one scaled aniline mixed row. The common risk is that composed-coordinate terms can improve or preserve localization but lose explicit chemical labels needed by the benchmark class extractor.

Current policy: composed evidence remains suitable as a viewer evidence layer and benchmark diagnostic. It is not yet supported as an automatic warning/confirmation policy layer or as fallback into final assignments. A future fallback experiment should require new evidence where baseline PED is diffuse or unavailable and composed evidence improves semantic status without increasing C-H, N-H, O-H, or carboxylic-context failures.

Composed-coordinate semantic hardening outcome: `src/mode_assignment.py` now maps `composed_symmetric_XH_stretch(...)` and `composed_asymmetric_XH_stretch(...)` to explicit C-H, N-H, or O-H stretch families based on the heavy-atom center encoded in the composed coordinate name. `benchmarks/vibrational_assignments/compare_orcaveda_assignments.py` now summarizes composed PED evidence with the top six contributors, so weak but chemically important carboxylic-context contributors are not silently dropped from the diagnostic comparison. This changes composed evidence labeling and benchmark diagnostics only; it does not change `assignment_audit`, `ped_stage3d_agreement`, `ped_final_assignment`, B matrices, rank selection, or final-label policy.

Composed-coordinate semantic hardening validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_semantics_ped_tmp` returned `18 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_vibrational_assignment_multiscale.py -q --basetemp outputs\pytest_composed_semantics_compare_tmp` returned `7 passed`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\acetaldehyde.hess data\hess\acetamide.hess data\hess\aniline.hess data\hess\benzoic_acid.hess --outdir outputs\composed_semantics_subset_20260507` completed in CLI mode.
- Raw comparator for the same four-molecule problem subset wrote `outputs\composed_semantics_subset_20260507\composed_ped_comparison_raw.csv`. On the previous full-golden benchmark filtered to these molecules, `composed_vs_baseline_ped_status` had `5 worsens_semantic_match` and `40 same_semantic_match`; on the regenerated subset it had `0 worsens_semantic_match` and `45 same_semantic_match`, with 48 rows lacking composed comparison because those benchmark rows matched files outside the regenerated composed subset.
- Scaled comparator for the same four-molecule problem subset wrote `outputs\composed_semantics_subset_20260507\composed_ped_comparison_scaled.csv`. On the previous full-golden benchmark filtered to these molecules, `composed_vs_baseline_ped_status` had `6 worsens_semantic_match` and `39 same_semantic_match`; on the regenerated subset it had `0 worsens_semantic_match` and `45 same_semantic_match`, with 48 rows lacking composed comparison for rows outside the regenerated subset.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_semantics_focused_tmp` returned `33 passed, 1 skipped`.

Current limitation after semantic hardening: this is a subset validation, not a refreshed full-golden comparison. The composed layer remains diagnostic/viewer evidence only and still must not drive `ped_final_assignment` or policy warnings until a full-golden rerun shows no new C-H, N-H, O-H, carbonyl, or carboxylic-context regressions.

Full-golden semantic-hardening rerun outcome: running the current pipeline on all local `data/hess/*.hess` after the semantic hardening patch completed successfully and wrote outputs under `outputs/composed_semantics_full_golden_20260507`. The input set contained 55 `.hess` files, and `composed_ped_basis_diagnostics.csv` contained 55 rows. Rank loss remained zero. The total composed candidate count was 366, and the total selected composed coordinate count was 166. `assignment_audit.csv` contained 350 high-frequency rows at or above 2500 cm-1, with zero high-frequency unassigned rows.

Full-golden semantic-hardening commands:

- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_semantics_full_golden_20260507`, where `@hess` was the PowerShell-expanded sorted list from `Get-ChildItem -LiteralPath 'data\hess' -Filter '*.hess'`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Raw comparator with `--ped-final-assignment outputs\composed_semantics_full_golden_20260507\acetaldehyde__acetamide__acetanilide__plus_52_files__multi_file_55__ped_final_assignment.csv` and `--composed-ped-audit outputs\composed_semantics_full_golden_20260507\acetaldehyde__acetamide__acetanilide__plus_52_files__multi_file_55__composed_wilson_ped_audit.csv` wrote `outputs\composed_semantics_full_golden_20260507\composed_ped_comparison_raw.csv` and reported primary baseline counts of `27 PASS` and `66 WARN`.
- Scaled comparator with the same composed audit source wrote `outputs\composed_semantics_full_golden_20260507\composed_ped_comparison_scaled.csv` and reported primary baseline counts of `26 PASS` and `67 WARN`.

Full-golden semantic-hardening composed comparison:

- Raw `composed_vs_baseline_ped_status`: `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`.
- Scaled `composed_vs_baseline_ped_status`: `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`.
- The only localization-gain row in both raw and scaled comparisons was acetophenone near 1683 cm-1, mode 42. Baseline PED semantic status was `PASS`; composed PED semantic status was also `PASS`. The top contribution increased from baseline `45.933%` to composed `63.043354%`, and the composed top contributor remained `C=O stretch [carbonyl_CO_stretch(C12=O17)]`.

Current policy after full-golden semantic-hardening rerun: composed evidence is now cleaner as a diagnostic layer, with no raw or scaled semantic regressions on the current curated benchmark. It still remains a separate viewer/benchmark diagnostic source. Any graduation into warning/confirmation or fallback policy must be planned separately and must preserve `assignment_audit`, `ped_stage3d_agreement`, and `ped_final_assignment` contracts unless explicitly changed.

Conservative composed diagnostic policy outcome: `src/reports.py` now computes viewer-only composed diagnostic fields by comparing baseline PED evidence to composed PED evidence. `composed_confirms_with_better_localization` is emitted only when baseline PED and composed PED are semantically compatible and composed top contribution exceeds baseline top contribution by more than 10 percentage points. When composed evidence is available but differs from baseline or fills a diffuse/unclassified baseline case, the policy emits diagnostic hints rather than final-label fallbacks. The interactive viewer displays the selected composed evidence policy hint and localization delta, plus composed semantic status/reason. `benchmarks/vibrational_assignments/compare_orcaveda_assignments.py` now uses the same conservative naming for benchmark policy hints; benchmark semantic status remains expected-aware, while viewer semantic status is baseline-PED compatibility only.

Conservative composed diagnostic policy validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py -q --basetemp outputs\pytest_composed_policy_viewer_tmp` returned `7 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_vibrational_assignment_multiscale.py -q --basetemp outputs\pytest_composed_policy_compare_tmp` returned `7 passed`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_policy_full_golden_20260507`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Raw comparator wrote `outputs\composed_policy_full_golden_20260507\composed_ped_comparison_raw.csv` and reported primary baseline counts of `27 PASS` and `66 WARN`; `composed_vs_baseline_ped_status` was `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`; `composed_ped_policy_hint` was `92 viewer_evidence_only` and `1 composed_confirms_with_better_localization`.
- Scaled comparator wrote `outputs\composed_policy_full_golden_20260507\composed_ped_comparison_scaled.csv` and reported primary baseline counts of `26 PASS` and `67 WARN`; `composed_vs_baseline_ped_status` was `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`; `composed_ped_policy_hint` was `92 viewer_evidence_only` and `1 composed_confirms_with_better_localization`.
- Full-golden diagnostics remained stable: `composed_ped_basis_diagnostics.csv` had 55 rows, zero rank loss, 366 composed candidates, and 166 selected composed coordinates. `assignment_audit.csv` had 350 high-frequency rows at or above 2500 cm-1 and zero high-frequency unassigned rows.
- `Select-String` checks on the regenerated `spectrum_data.json` and `interactive_spectrum.html` found the new composed policy fields and viewer labels.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_policy_focused_tmp` returned `34 passed, 1 skipped`.

Current limitation after conservative composed diagnostic policy: the policy is intentionally diagnostic. It does not write composed policy fields into `ped_final_assignment.csv`, and it does not let composed PED override or fallback into final labels. A future plan can decide whether these hints should be exported as a dedicated CSV table or used for a stricter warning/confirmation layer.

Composed policy CSV export outcome: `src/reports.py` now has `build_composed_ped_policy_diagnostics_table(...)`, and `src/ORCAVEDA_patched_stage3D_v5_0.py` writes a separate `composed_ped_policy_diagnostics.csv` table for every pipeline run. The table columns are `Filename`, `mode`, `frequency_cm-1`, `ped_assignment`, `ped_top_percent`, `composed_ped_assignment`, `composed_ped_top_percent`, `composed_ped_localization_delta_percent`, `composed_ped_semantic_status`, `composed_ped_semantic_reason`, and `composed_ped_policy_hint`. This table is diagnostic only and is not consumed by `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment`.

Composed policy CSV validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_policy_csv_ped_tmp` returned `18 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py tests\test_vibrational_assignment_multiscale.py -q --basetemp outputs\pytest_composed_policy_csv_viewer_compare_tmp` returned `14 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_policy_csv_full_golden_20260507`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Full-golden `composed_ped_policy_diagnostics.csv` contained 1920 rows and all required columns. Policy hint counts were `1739 viewer_evidence_only`, `145 composed_confirms_with_better_localization`, `20 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Raw comparator wrote `outputs\composed_policy_csv_full_golden_20260507\composed_ped_comparison_raw.csv` and reported `27 PASS`, `66 WARN`, `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`. Benchmark policy hints were `92 viewer_evidence_only` and `1 composed_confirms_with_better_localization`.
- Scaled comparator wrote `outputs\composed_policy_csv_full_golden_20260507\composed_ped_comparison_scaled.csv` and reported `26 PASS`, `67 WARN`, `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`. Benchmark policy hints were `92 viewer_evidence_only` and `1 composed_confirms_with_better_localization`.
- Critical benchmark rows involving X-H, carboxyl/carboxylic acid, or carbonyl had zero `diagnostic_hint_composed_differs_from_baseline` rows in both raw and scaled comparisons.
- Full-golden diagnostics remained stable: `composed_ped_basis_diagnostics.csv` had 55 rows, zero rank loss, 366 composed candidates, and 166 selected composed coordinates. `assignment_audit.csv` had 350 high-frequency rows at or above 2500 cm-1 and zero high-frequency unassigned rows.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_policy_csv_focused_tmp` returned `34 passed, 1 skipped`.

Current limitation after CSV export: the separate policy table now exposes all composed diagnostic hints, including 20 baseline/composed differences outside the current critical benchmark rows. Those rows are diagnostic inspection targets, not final-label changes. The next safe UI step is a viewer filter/indicator for non-`viewer_evidence_only` hints.

Viewer composed-hint filter outcome: `src/reports.py` now adds a peak-table `composedHintFilter` with `All modes`, `Composed hints`, `Better localization`, and `Differs from baseline` options. The peak table has a compact `Composed Hint` column that labels non-`viewer_evidence_only` rows, including `Better localization`, `Differs from baseline`, `Baseline diffuse`, and `Composed only`. The filter uses only viewer JSON fields and does not change assignments, spectrum rendering data, CSV scientific tables, `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment`.

Viewer composed-hint filter validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py -q --basetemp outputs\pytest_composed_hint_filter_viewer_tmp` returned `7 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_hint_filter_ped_tmp` returned `18 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_hint_filter_focused_tmp` returned `34 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_hint_filter_full_golden_20260507`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Full-golden HTML smoke checks found `composedHintFilter`, `Composed Hint`, `Composed hints`, `Better localization`, `Differs from baseline`, and the empty-state text `No modes match the selected composed hint filter.` in the regenerated interactive spectrum HTML.
- Full-golden `composed_ped_policy_diagnostics.csv` remained stable with 1920 rows and hint counts of `1739 viewer_evidence_only`, `145 composed_confirms_with_better_localization`, `20 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Full-golden diagnostics remained stable: 55 composed basis diagnostic rows, zero rank loss, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.
- Raw and scaled benchmark comparators under `outputs\composed_hint_filter_full_golden_20260507` both reported `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`; benchmark policy hints remained `92 viewer_evidence_only` and `1 composed_confirms_with_better_localization`; critical X-H/carboxyl/carbonyl rows had zero `diagnostic_hint_composed_differs_from_baseline` in raw and scaled outputs.

Current limitation after viewer filter: browser pixel/interaction automation was not run for this UI-only patch. The validation covers generated HTML strings, JSON payload fields, pipeline generation, focused tests, and full-golden smoke outputs. The filter is intentionally table-only and does not filter the spectrum plot.

Primitive torsion semantic hardening outcome: `src/mode_assignment.py` now labels primitive torsions containing terminal hydrogens by the heavy atom adjacent to the terminal H in the torsion label. For example, `tor(C1-N3-C4-H12)` and `tor(C1-O2-C6-H5)` now map to `C-H torsion`, while true `tor(C1-C2-O3-H4)` and `tor(C1-C2-N3-H4)` retain `O-H torsion` and `N-H torsion`. This fixes misleading `N-H torsion` and `O-H torsion` labels in molecules without N-H or O-H bonds without changing PED math, coordinate selection, B matrices, final assignment policy, or output schemas.

Primitive torsion semantic hardening validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_torsion_semantics_ped_tmp` returned `19 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py -q --basetemp outputs\pytest_torsion_semantics_viewer_compare_tmp` returned `14 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_torsion_semantics_focused_tmp` returned `35 passed, 1 skipped`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\torsion_semantics_full_golden_20260507`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Full-golden `composed_ped_policy_diagnostics.csv` remained 1920 rows. Policy hint counts changed from the previous `1739 viewer_evidence_only`, `145 composed_confirms_with_better_localization`, `20 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified` to `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- The DMSO 1049 cm-1 S=O row no longer appears among `diagnostic_hint_composed_differs_from_baseline`. DMF, dimethyl carbonate, ethylene oxide, ethanol, and monoethanolamine dimer rows now report terminal-C-H torsions as `C-H torsion` instead of misleading `N-H torsion` or `O-H torsion` where the terminal H is adjacent to carbon.
- Full-golden diagnostics remained stable: 55 composed basis diagnostic rows, zero rank loss, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.
- Raw and scaled benchmark comparators under `outputs\torsion_semantics_full_golden_20260507` both reported `92 same_semantic_match`, `1 localization_gain_without_semantic_improvement`, and `0 worsens_semantic_match`; benchmark policy hints remained `92 viewer_evidence_only` and `1 composed_confirms_with_better_localization`; critical X-H/carboxyl/carbonyl rows had zero `diagnostic_hint_composed_differs_from_baseline` in raw and scaled outputs.

Current limitation after primitive torsion semantic hardening: remaining `diagnostic_hint_composed_differs_from_baseline` rows are mostly genuine bend/torsion/stretch competition or high-frequency stretch localization disagreements, not obvious X-H label bugs. They should be treated as candidate coordinate-generation or benchmark-fixture review targets, not as final-label changes.

Composed hint triage outcome: `src/reports.py` now adds `triage_composed_ped_diagnostic_hint(...)`, and `composed_ped_policy_diagnostics.csv` now includes `composed_ped_triage_category` and `composed_ped_triage_recommendation`. These fields are diagnostic only. They are not written into `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment`, and they do not affect viewer final labels.

Composed hint triage validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py tests\test_ped.py -q --basetemp outputs\pytest_composed_triage_tmp` returned `27 passed` before failing one unrelated headless Chrome smoke test at browser startup with `GPU process isn't usable`. The failure occurred in `test_interactive_spectrum_viewer_headless_chrome_smoke_without_cdn`, before exercising the new composed triage logic.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py tests\test_ped.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_composed_triage_no_chrome_tmp` returned `27 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_vibrational_assignment_multiscale.py -q --basetemp outputs\pytest_composed_triage_compare_tmp` returned `7 passed`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_composed_triage_focused_tmp` returned `36 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_triage_full_golden_20260512`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Full-golden `composed_ped_policy_diagnostics.csv` remained 1920 rows. Policy hint counts remained `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Full-golden triage counts were `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 motion_family_coordinate_generation_target`, and `3 high_frequency_motion_family_review`.
- Full-golden diagnostics remained stable: 55 composed basis diagnostic rows, all 55 rank-preserved, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.

Current limitation after composed hint triage: this patch only makes the remaining diagnostic rows easier to inspect. It does not implement new composed-coordinate generators. The next implementation step should focus on the 8 DMF/ethene-like `motion_family_coordinate_generation_target` rows and 3 ethylene oxide `high_frequency_motion_family_review` rows; the 8 negative-delta `baseline_preferred_composed_lower_localization` rows should not be promoted into composed policy.

DMF/ethene and high-frequency ethylene oxide review outcome: focused inspection of the 8 `motion_family_coordinate_generation_target` rows showed that composed Wilson PED selected primitive C-H torsion rows, not generated composed torsion rows. A trial patch that restricted composed-basis optimization to composed candidates only passed focused PED tests and removed the original DMF/ethene targets on a three-molecule subset, but full-golden validation broadened the diagnostic disagreement surface and was rejected. The retained implementation improves triage only: high-frequency rows where composed PED recovers an explicit X-H stretch from a non-stretch baseline are now categorized as `high_frequency_xh_stretch_recovery`, with recommendation `keep_composed_xh_stretch_as_diagnostic_evidence`.

DMF/ethene and high-frequency ethylene oxide validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q --basetemp outputs\pytest_composed_optimizer_scope_ped_tmp` returned `20 passed` during the rejected optimizer-scope experiment.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\DMF_freq.hess data\hess\ethene.hess data\hess\ethylene_oxide.hess --outdir outputs\composed_optimizer_scope_subset_20260512_b` returned exit code 0 and printed `Terminal mode detected -> using CLI`. The subset had zero original DMF/ethene targets but still one new low-frequency ethylene oxide target; the three high-frequency ethylene oxide rows became `high_frequency_xh_stretch_recovery`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_optimizer_scope_full_golden_20260512`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 but produced `181 diagnostic_hint_composed_differs_from_baseline` rows, including `46 motion_family_coordinate_generation_target` and `13 high_frequency_xh_stretch_recovery`. This was rejected as too broad.
- After reverting the optimizer-scope experiment, `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_interactive_spectrum_viewer.py tests\test_vibrational_assignment_multiscale.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_composed_xh_recovery_triage_tmp` returned `34 passed, 1 deselected`.
- Final full-golden run `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_xh_recovery_triage_full_golden_20260512` returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Final full-golden policy counts remained `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Final full-golden triage counts were `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 motion_family_coordinate_generation_target`, and `3 high_frequency_xh_stretch_recovery`. There were zero `high_frequency_motion_family_review` rows.
- Final full-golden diagnostics remained stable: 55 composed basis diagnostic rows, all 55 rank-preserved, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.

Current limitation after high-frequency X-H recovery triage: the 8 DMF/ethene rows remain genuine coordinate-generation or optimizer-diagnostic targets. They should not be fixed by a global optimizer candidate-pool restriction. A future safe patch should be narrower, likely adding a diagnostic that distinguishes primitive-row optimizer substitutions from composed-coordinate evidence, or a molecule-agnostic generator/constraint that can be validated without increasing full-golden disagreement counts.

Composed PED provenance outcome: `src/ped.py` now writes `generation_rule` in PED v1, PED v2, and Wilson PED audit rows. `src/reports.py` now carries the top-ranked composed PED contributor provenance into `composed_ped_policy_diagnostics.csv` as `composed_ped_top_source`, `composed_ped_top_internal_coordinate`, `composed_ped_top_coord_index`, `composed_ped_top_generation_rule`, and `composed_ped_top_is_composed_coordinate`. These fields are diagnostic only and are intentionally absent from `assignment_audit`, `ped_stage3d_agreement`, and `ped_final_assignment`.

Composed PED provenance validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_interactive_spectrum_viewer.py tests\test_vibrational_assignment_multiscale.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_composed_provenance_focused_tmp` returned `34 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\DMF_freq.hess data\hess\ethene.hess data\hess\ethylene_oxide.hess --outdir outputs\composed_provenance_subset_20260512` returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- In `outputs\composed_provenance_subset_20260512\*__composed_ped_policy_diagnostics.csv`, all 8 DMF/ethene `motion_family_coordinate_generation_target` rows had `composed_ped_top_source == primitive`, `composed_ped_top_is_composed_coordinate == False`, and top coordinates `tor(C1-N3-C4-H12)`, `tor(C5-N3-C4-H12)`, `tor(H3-C1-C2-H6)`, `tor(H3-C1-C2-H5)`, or `tor(H4-C1-C2-H6)`.
- The three high-frequency ethylene oxide `high_frequency_xh_stretch_recovery` rows had explicit top provenance: modes 17 and 19 were topped by composed X-H stretch rows with `generation_rule == xh_pair_sum_difference`, while mode 18 was topped by primitive `r(C2-H4)`.
- Final full-golden run `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_provenance_full_golden_20260512`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Final full-golden policy counts remained unchanged: `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Final full-golden triage counts remained unchanged: `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 motion_family_coordinate_generation_target`, and `3 high_frequency_xh_stretch_recovery`.
- All 8 full-golden `motion_family_coordinate_generation_target` rows had `composed_ped_top_source == primitive` and zero had `composed_ped_top_source == composed_coordinate`.
- Final full-golden diagnostics remained stable: 55 composed basis diagnostic rows, zero rank loss, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.
- Comparing `ped_final_assignment.csv` from `outputs\composed_xh_recovery_triage_full_golden_20260512` and `outputs\composed_provenance_full_golden_20260512` over `Filename`, `mode`, `frequency_cm-1`, `final_assignment`, `final_assignment_source`, `final_assignment_policy`, and `final_assignment_warning` found 1920 rows in both files and no changed columns.

Current limitation after composed PED provenance: provenance confirms the 8 DMF/ethene target rows are primitive-row optimizer substitutions, not composed-coordinate top contributors. The next safe experiment should use this provenance to add an optimizer diagnostic or a narrower constraint that flags or limits primitive substitution artifacts without reducing useful composed X-H stretch recovery and without increasing full-golden disagreement counts.

Primitive-row substitution triage outcome: `src/reports.py` now classifies composed PED motion-family disagreements with positive localization delta and a non-composed top contributor as `primitive_row_optimizer_substitution`, with recommendation `inspect_optimizer_substitution_before_coordinate_generation`. This keeps `diagnostic_hint_composed_differs_from_baseline` policy status unchanged but distinguishes optimizer substitution artifacts from missing composed-coordinate generators. It does not change PED math, B matrices, rank selection, `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment`.

Primitive-row substitution triage validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py tests\test_ped.py tests\test_vibrational_assignment_multiscale.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_primitive_substitution_triage_tmp` returned `34 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\DMF_freq.hess data\hess\ethene.hess data\hess\ethylene_oxide.hess --outdir outputs\primitive_substitution_triage_subset_20260512` returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- In the subset output, all 8 DMF/ethene `diagnostic_hint_composed_differs_from_baseline` rows became `primitive_row_optimizer_substitution`, and all 8 had `composed_ped_top_source == primitive`.
- The three high-frequency ethylene oxide rows remained `high_frequency_xh_stretch_recovery`.
- Final full-golden run `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\primitive_substitution_triage_full_golden_20260512`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Final full-golden policy counts remained unchanged: `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Final full-golden triage counts were `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 primitive_row_optimizer_substitution`, and `3 high_frequency_xh_stretch_recovery`.
- Final full-golden diagnostics remained stable: 55 composed basis diagnostic rows, zero rank loss, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.
- Comparing `ped_final_assignment.csv` from `outputs\composed_provenance_full_golden_20260512` and `outputs\primitive_substitution_triage_full_golden_20260512` over `Filename`, `mode`, `frequency_cm-1`, `final_assignment`, `final_assignment_source`, `final_assignment_policy`, and `final_assignment_warning` found 1920 rows in both files and no changed columns.

Current limitation after primitive-row substitution triage: the next step is still a method change, not a label change. Candidate approaches include adding an explicit warning when composed evidence is topped by newly selected primitive rows, or testing a narrower optimizer constraint that applies only to primitive substitution artifacts and is accepted only if full-golden policy counts do not broaden.

Composed evidence-origin outcome: `src/reports.py` now writes `composed_ped_evidence_origin` in `composed_ped_policy_diagnostics.csv`. Values are `composed_coordinate_top`, `primitive_substitution_top`, or `baseline_or_no_composed_top`. This is a stable machine-readable origin field derived from top-contributor provenance and is intended for future optimizer constraints or warning-only policy. It does not change policy hints, triage categories, PED math, rank selection, `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment`.

Composed evidence-origin validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py tests\test_ped.py tests\test_vibrational_assignment_multiscale.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_evidence_origin_focused_tmp` returned `35 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\DMF_freq.hess data\hess\ethene.hess data\hess\ethylene_oxide.hess --outdir outputs\evidence_origin_subset_20260512` returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- In the subset output, all 8 DMF/ethene primitive-row substitution rows had `composed_ped_evidence_origin == primitive_substitution_top`.
- In the subset output, ethylene oxide high-frequency modes 17 and 19 had `composed_ped_evidence_origin == composed_coordinate_top`, while mode 18 had `primitive_substitution_top` because the top row was primitive `r(C2-H4)`.
- Final full-golden run `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\evidence_origin_full_golden_20260512`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Final full-golden policy counts remained unchanged: `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Final full-golden triage counts remained unchanged: `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 primitive_row_optimizer_substitution`, and `3 high_frequency_xh_stretch_recovery`.
- Final full-golden origin counts were `1440 primitive_substitution_top`, `150 composed_coordinate_top`, and `330 baseline_or_no_composed_top`.
- Final full-golden diagnostics remained stable: 55 composed basis diagnostic rows, zero rank loss, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.
- Comparing `ped_final_assignment.csv` from `outputs\primitive_substitution_triage_full_golden_20260512` and `outputs\evidence_origin_full_golden_20260512` over `Filename`, `mode`, `frequency_cm-1`, `final_assignment`, `final_assignment_source`, `final_assignment_policy`, and `final_assignment_warning` found 1920 rows in both files and no changed columns.

Current limitation after evidence-origin field: many viewer-evidence rows naturally have `primitive_substitution_top` because composed PED basis selection can still be topped by primitive rows. This field should not be interpreted as a failure by itself. It becomes actionable only together with policy/triage fields such as `primitive_row_optimizer_substitution`.

Primitive substitution warning-only outcome: `src/reports.py` now adds `classify_composed_ped_warning(...)` and writes `composed_ped_warning` plus `composed_ped_warning_reason` in `composed_ped_policy_diagnostics.csv` and the interactive viewer payload. The warning is intentionally narrow: it is emitted only when `composed_ped_evidence_origin == primitive_substitution_top` and `composed_ped_triage_category == primitive_row_optimizer_substitution`. Primitive-topped viewer-evidence rows without that triage remain unflagged. This is diagnostic only and does not change policy hints, triage categories, PED math, rank selection, `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment`.

Primitive substitution warning-only validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py tests\test_ped.py tests\test_vibrational_assignment_multiscale.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_primitive_warning_focused_tmp` returned `36 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\DMF_freq.hess data\hess\ethene.hess data\hess\ethylene_oxide.hess --outdir outputs\primitive_warning_subset_20260512` returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- In the subset output, `composed_ped_warning` emitted `8 primitive_row_optimizer_substitution_warning` rows and `67` blank rows. The 8 warnings were the known DMF/ethene primitive-row substitutions, and the 3 ethylene oxide high-frequency rows remained `high_frequency_xh_stretch_recovery`.
- Full-golden run `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\primitive_warning_full_golden_20260512`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Final full-golden policy counts remained unchanged: `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Final full-golden triage counts remained unchanged: `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 primitive_row_optimizer_substitution`, and `3 high_frequency_xh_stretch_recovery`.
- Final full-golden warning counts were `1912` blank rows and `8 primitive_row_optimizer_substitution_warning` rows. The warning rows were limited to `DMF_freq.hess` (`5`) and `ethene.hess` (`3`).
- Final full-golden evidence-origin counts remained `1440 primitive_substitution_top`, `330 baseline_or_no_composed_top`, and `150 composed_coordinate_top`.
- Final full-golden diagnostics remained stable: 55 composed basis diagnostic rows, zero rank loss, 366 composed candidates, 166 selected composed coordinates, 350 high-frequency rows at or above 2500 cm-1, and zero high-frequency unassigned rows.
- Comparing `ped_final_assignment.csv` from `outputs\evidence_origin_full_golden_20260512` and `outputs\primitive_warning_full_golden_20260512` over `Filename`, `mode`, `frequency_cm-1`, `final_assignment`, `final_assignment_source`, `final_assignment_policy`, and `final_assignment_warning` found 1920 rows in both files and no changed columns.

Current limitation after warning-only patch: no optimizer behavior changed. The new warning only identifies primitive-row optimizer substitution evidence for review; it is not a constraint and does not resolve those rows.

Diagnostic baseline acceptance criteria before optimizer-constraint experiments:

- Full-golden policy counts must not broaden beyond `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Full-golden triage counts must not get worse than `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 primitive_row_optimizer_substitution`, and `3 high_frequency_xh_stretch_recovery`.
- `composed_ped_warning` count must be `8` or fewer.
- Composed PED basis rank loss must remain `0`.
- High-frequency rows at or above 2500 cm-1 must have `0` unassigned final assignments.
- `ped_final_assignment.csv` must remain unchanged over `Filename`, `mode`, `frequency_cm-1`, `final_assignment`, `final_assignment_source`, `final_assignment_policy`, and `final_assignment_warning` unless a new explicit plan justifies a final-label policy change.

Primitive-substitution constraint experiment outcome: `src/ORCAVEDA_patched_stage3D_v5_0.py` now has an opt-in experimental repair controlled by `experimental_composed_primitive_substitution_constraint`, exposed on the CLI as `--experimental-composed-primitive-substitution-constraint`. The default behavior is unchanged. When enabled, the pipeline first runs the normal composed PED basis optimization, builds temporary composed policy diagnostics, identifies only rows matching `primitive_row_optimizer_substitution`, `primitive_substitution_top`, positive localization delta, and `motion_family_mismatch`, then attempts to revert the top primitive replacement to the corresponding original Stage 3D basis row if rank is preserved. This does not globally forbid primitive rows, and it does not change `assignment_audit`, `ped_stage3d_agreement`, or `ped_final_assignment` policy.

Primitive-substitution constraint subset validation:

- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_interactive_spectrum_viewer.py -q -k "not headless_chrome_smoke_without_cdn" --basetemp outputs\pytest_experimental_constraint_focused_tmp` returned `29 passed, 1 deselected`.
- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\DMF_freq.hess data\hess\ethene.hess data\hess\ethylene_oxide.hess --outdir outputs\primitive_constraint_subset_20260512 --experimental-composed-primitive-substitution-constraint` returned exit code 0 and printed `Terminal mode detected -> using CLI`.
- Subset `composed_ped_warning` counts changed from the diagnostic baseline `8 primitive_row_optimizer_substitution_warning` rows to `75` blank rows and zero warnings.
- Subset triage counts became `67 viewer_evidence_only`, `5 confirmation_candidate`, and `3 high_frequency_xh_stretch_recovery`; the prior 8 primitive-row substitution rows were removed.
- The three ethylene oxide high-frequency X-H recovery rows were preserved. Modes 17 and 19 remained `composed_coordinate_top`; mode 18 remained `primitive_substitution_top` but stayed classified as `high_frequency_xh_stretch_recovery`, not as a substitution warning.
- Subset diagnostics reported rank preserved for all 3 files, `8` primitive-substitution targets before repair, `0` after repair, and `3` rank-safe reverted primitive indices.
- Subset high-frequency rows at or above 2500 cm-1 had `0` unassigned final assignments.
- Comparing subset `ped_final_assignment.csv` from `outputs\primitive_warning_subset_20260512` and `outputs\primitive_constraint_subset_20260512` over final-label columns found 75 rows in both files and no changed columns.

Primitive-substitution constraint full-golden validation:

- `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\primitive_constraint_full_golden_live --experimental-composed-primitive-substitution-constraint`, where `@hess` was the PowerShell-expanded sorted list from `data\hess`, first timed out at the 120 second tool limit and was then rerun with a longer timeout. The rerun returned exit code 0 in 101.2 seconds and printed `Terminal mode detected -> using CLI`.
- Full-golden `composed_ped_basis_diagnostics.csv` had 55 rows and zero rank loss. Constraint status counts were `53 no_targets` and `2 applied`.
- Constraint diagnostics reported 8 primitive-substitution targets before repair, 0 after repair, and 3 reverted primitive indices.
- Full-golden policy counts changed from the warning baseline `1745 viewer_evidence_only`, `140 composed_confirms_with_better_localization`, `19 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified` to `1767 viewer_evidence_only`, `126 composed_confirms_with_better_localization`, `11 diagnostic_hint_composed_differs_from_baseline`, and `16 diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified`.
- Full-golden triage counts changed from the warning baseline `1745 viewer_evidence_only`, `140 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, `8 primitive_row_optimizer_substitution`, and `3 high_frequency_xh_stretch_recovery` to `1767 viewer_evidence_only`, `126 confirmation_candidate`, `16 baseline_gap_candidate`, `8 baseline_preferred_composed_lower_localization`, and `3 high_frequency_xh_stretch_recovery`.
- `composed_ped_warning` counts changed from 8 `primitive_row_optimizer_substitution_warning` rows in `outputs\primitive_warning_full_golden_20260512` to 1920 blank warning rows in `outputs\primitive_constraint_full_golden_live`.
- Full-golden evidence-origin counts were `1442 primitive_substitution_top`, `148 composed_coordinate_top`, and `330 baseline_or_no_composed_top`.
- Comparing `ped_final_assignment.csv` from `outputs\primitive_warning_full_golden_20260512` and `outputs\primitive_constraint_full_golden_live` over `Filename`, `mode`, `frequency_cm-1`, `final_assignment`, `final_assignment_source`, `final_assignment_policy`, and `final_assignment_warning` found 1920 rows in both files and no changed columns.
- High-frequency rows at or above 2500 cm-1 had 0 unassigned final assignments.

Current limitation after full-golden constraint validation: the opt-in repair passed the planned acceptance checks, but it is still not default behavior. The run reduced diagnostic disagreement and removed primitive-substitution warnings, but it also reduced `confirmation_candidate` rows from 140 to 126. Any promotion from opt-in experiment to default behavior needs a separate explicit decision and validation plan.

## Context and Orientation

The relevant project files are:

- `src/orcaveda_models.py`: defines `HessData`, `InternalCoordinate`, and related data classes.
- `src/internal_coordinates.py`: builds primitive and functional-group internal coordinates from atoms, coordinates, bonds, fragments, hydrogen bonds, and functional groups.
- `src/b_matrix.py`: computes finite-difference B matrices, rank/condition diagnostics, independent-coordinate selection, pivoted Cholesky fallback, and the current EPM-like basis swap optimizer.
- `src/ped.py`: computes PED v1, PED v2 force-aware diagnostics, Wilson G, reconstructed internal F, and Wilson GF-style PED audits.
- `src/mode_assignment.py`: builds Stage 3D assignment audits from selected independent coordinates and contains assignment-family semantics.
- `src/reports.py`: builds PED/Stage3D agreement and PED-driven final assignment tables plus viewer payloads.
- `src/ORCAVEDA_patched_stage3D_v5_0.py`: legacy-compatible main pipeline that wires together parser, chemistry, internals, B matrix, Stage 3D, PED, reports, and CSV outputs.
- `tests/test_ped.py`: focused PED and Wilson tests.
- `tests/test_golden_rdkit_outputs.py`, `tests/test_stage3d_outputs.py`, and `tests/test_regression_baseline_outputs.py`: broader safety checks.
- `docs/full_ped_execplan.md`: current PED history and evidence log.

The key scientific terms are:

- Internal coordinate: a local coordinate such as a bond stretch, angle bend, torsion, out-of-plane motion, or functional-group-local coordinate.
- Primitive coordinate: a single internal coordinate generated directly from geometry or a functional-group template.
- Composed coordinate: a linear combination of primitive coordinates, such as a symmetric CH3 stretch or ring breathing coordinate.
- B matrix: derivative matrix mapping Cartesian displacements to internal-coordinate displacements.
- PED: potential energy distribution, a percentage-like decomposition of normal-mode potential energy into internal-coordinate contributions.
- Wilson GF-style PED audit: ORCAVEDA's current force-aware internal-coordinate diagnostic using `G = B M^-1 B^T`, reconstructed internal F, and mode-projected potential-energy terms.
- EPM-like score: ORCAVEDA diagnostic inspired by Jamroz's EPM parameter. It measures whether modes have dominant PED contributors and whether contributions are diffuse.

Jamroz 2013 describes VEDA as automatically generating internal coordinates, checking linear independence, mixing coordinates by geometry and atom type, preserving motion categories such as stretch/bend/torsion, separating CH-like and not-CH coordinates, optimizing the PED matrix to increase maximal elements, optionally using ADM extraction, reducing coordinate complexity, allowing user/frozen coordinates, and acknowledging that PED interpretation is not unique.

## Plan of Work

First, introduce an explicit representation for composed coordinates. The smallest safe path is to keep `InternalCoordinate` usable as-is and add an optional companion structure or wrapper that stores a composed coordinate's components, coefficients, category, source, and human-readable label. The composed coordinate's function evaluates the weighted sum of the component coordinate values. Its B row can be produced either by finite-differencing the composed function or by linearly combining already-computed primitive B rows. The linear-combination route is preferred for reproducibility and speed because a composed coordinate's derivative is exactly the same weighted sum of component derivatives under the local coordinate definitions.

Second, classify coordinate candidates into optimization groups. Groups must respect Jamroz's first safe rules: combine only same motion type in the initial implementation, and separate H-containing coordinates from heavy-atom-only coordinates. The first categories should be `stretch_CH_like`, `stretch_heavy`, `bend_CH_like`, `bend_heavy`, `torsion_CH_like`, `torsion_heavy`, and `oop` once out-of-plane coordinates exist. OH and NH coordinates may be treated as CH-like for optimizer grouping, but this must be recorded in diagnostics.

Third, add deterministic composed coordinate generators for chemically important motifs. The first wave should target motifs where sums/differences are scientifically clear and easy to validate: H2O symmetric/asymmetric O-H stretches, NH2 symmetric/asymmetric N-H stretches and scissor, CH2 symmetric/asymmetric stretches and bends, CH3 symmetric/asymmetric stretch and deformation combinations, carboxylate/nitro/sulfone symmetric/asymmetric X-O stretches, carbonyl-adjacent C-C/C-O mixed local coordinates only after same-type rules are implemented, and ring breathing/ring stretch sums for small aromatic and aliphatic rings.

Fourth, compute EPM-like metrics before and after adding composed coordinates. Metrics must include mean top PED contribution, median top PED contribution, diffuse mode fraction, rank, condition number, number of composed coordinates selected, number of primitive coordinates replaced, and the exact changed coordinates. Add both mode-centric and coordinate-centric views when possible. Do not call the metric literal VEDA EPM until its definition matches the paper's orientation and denominator.

Fifth, integrate composed candidates into PED basis optimization without changing Stage 3D. The pipeline should build the primitive internal coordinate pool, select the existing Stage 3D independent basis, build composed candidates, combine primitive and composed candidate B rows for PED only, run rank-preserving selection/optimization for PED, and write separate outputs. Stage 3D `assignment_audit` must continue using the original independent basis until a later plan explicitly changes it.

Sixth, update PED and reports. Add `composed_ped_basis.csv`, `composed_coordinate_definitions.csv`, and summary columns showing composed optimization status. Add columns to PED audit rows that identify whether a top contributor is primitive or composed and, for composed rows, list component coordinates and coefficients. Viewer wording must call this "composed-coordinate PED optimization" or "VEDA-inspired EPM-like optimization", not "VEDA-equivalent".

Seventh, validate incrementally. Start with synthetic molecules where expected sums and differences are obvious, then water and ammonia, then benzene and acetophenone, then the full golden `.hess` set. Compare not only EPM-like metrics but also PED-driven final label policy counts, Stage 3D fallback counts, high-frequency unassigned modes, rank loss, and selected NIST matching diagnostics when relevant.

Eighth, harden composed-coordinate semantics before policy integration. Use the benchmark rows where composed evidence worsened as test fixtures and improve only label propagation and diagnostic text first. The target is to preserve explicit X-H labels (`C-H`, `N-H`, `O-H`), carbonyl labels (`C=O`), and carboxylic-acid context when a composed coordinate is made from primitive coordinates that already carry those semantics. Re-run the comparator and require fewer or zero `worsens_semantic_match` rows before considering warning/confirmation or fallback policies.

## Concrete Steps

1. Create data model helpers.

   Edit `src/internal_coordinates.py` or a new focused module such as `src/composed_coordinates.py`. Add a structure that can represent a composed coordinate with:

       name
       kind
       atoms0
       priority
       source = "composed_coordinate"
       components = [(coord_index, coefficient), ...]
       category
       generation_rule

   Add tests that a composed coordinate B row equals the coefficient-weighted sum of component B rows.

2. Add coordinate grouping.

   Implement functions that classify an `InternalCoordinate` into motion type and H/heavy category. Use atom indices and atom symbols, not only label strings. Keep the grouping function pure and unit-test it with CH, OH, NH, heavy-heavy, bend, and torsion examples.

3. Add first composed generators.

   Start with water, ammonia, methyl, methylene, nitro, sulfone, carboxylate-like paired stretches, and aromatic ring stretch sums. Keep each generator deterministic and conservative. Do not generate cross-type mixed coordinates in the first pass.

4. Add composed candidate B matrix builder.

   Build a stacked PED candidate matrix:

       B_ped_candidates = vstack([B_primitive, B_composed])

   Preserve a mapping from candidate row to primitive or composed coordinate definition. Do not mutate the primitive internal-coordinate list used by Stage 3D.

5. Extend optimizer.

   Reuse `select_independent_coordinates` and `optimize_independent_coordinates_for_ped` on the PED candidate matrix. Add constraints so same-type and H/heavy grouping rules are respected during generated candidate construction, and so rank never drops relative to the Stage 3D independent basis.

6. Integrate in pipeline.

   In `src/ORCAVEDA_patched_stage3D_v5_0.py`, keep:

       selected_idx -> Stage 3D assignment_audit

   Add:

       ped_candidate_internals, B_ped_candidates, composed_definitions = build_composed_ped_candidates(...)
       ped_selected_idx, ped_report = optimize_independent_coordinates_for_ped(...)

   Use `ped_candidate_internals` and `B_ped_candidates` only for PED outputs.

7. Add reports.

   Write:

       *__composed_coordinate_definitions.csv
       *__composed_ped_basis.csv

   Add manifest wording that states this is VEDA-inspired composed-coordinate optimization, not a full VEDA implementation.

8. Add viewer and final-label safeguards only after backend evidence exists.

   Use `$orcaveda-viewer` when touching viewer text. Add a source label such as `Composed-coordinate Wilson PED audit`. Keep Stage 3D visible.

9. Benchmark.

   Run full golden and compare before/after:

       EPM-like mean top percent
       diffuse fraction
       PED final vs Stage 3D fallback counts
       high-frequency unassigned modes
       rank losses
       selected NIST matching status where already available

## Validation and Acceptance

Focused unit tests must pass:

    .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q

Core regression tests must pass after pipeline integration:

    .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_golden_rdkit_outputs.py tests\test_vibrational_assignment_multiscale.py tests\test_interactive_spectrum_viewer.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q --basetemp outputs\pytest_composed_ped_optimizer_focused_tmp

The full golden batch must complete:

    $hess = Get-ChildItem -LiteralPath 'data\hess' -Filter '*.hess' | Sort-Object Name | ForEach-Object { $_.FullName }
    .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py @hess --outdir outputs\composed_ped_optimizer_full_golden_live

Acceptance criteria for the first implementation:

- No rank losses in the composed PED basis relative to the Stage 3D independent basis.
- No high-frequency C-H, N-H, or O-H modes become unassigned.
- New composed-coordinate tables are written and machine-readable.
- At least water or ammonia shows expected symmetric/asymmetric X-H composed coordinates in PED contributors.
- On the full golden set, mean top PED contributor or diffuse mode fraction improves relative to the primitive optimized PED basis, or the output clearly reports no improvement without silently degrading assignments.
- Viewer wording, if touched, does not call the method VEDA-equivalent.

Acceptance criteria for a later stronger milestone:

- CH3/CH2, aromatic ring, nitro, sulfone, and carboxyl-like composed coordinates improve PED localization on representative golden molecules.
- PED-only benchmark comparison improves or remains neutral, with failures explained by semantic or experimental-reference limitations rather than missing contributors.
- The report includes enough evidence to decide whether composed PED can expand PED-driven final label policy.

## Idempotence and Recovery

All generated outputs must go under `outputs/` with descriptive directory names and can be deleted and regenerated. The implementation must be additive: if composed candidate generation fails for a molecule, the pipeline should not hide the failure with broad `try/except`; it should either raise a clear error in tests or record a specific diagnostic only where the failure is expected and recoverable.

If composed-coordinate optimization causes rank loss, the pipeline must fall back to the current primitive optimized PED basis and record a warning such as `composed_ped_basis_rejected_rank_loss`. It must not silently use a rank-deficient composed basis.

If a new composed-coordinate generator creates chemically questionable labels, remove or disable only that generator and keep the rest of the composed-coordinate infrastructure. Do not revert unrelated PED, Stage 3D, NIST, or viewer changes.

## Artifacts and Notes

Evidence already available before this plan:

- `outputs/jamroz_2013_text.txt`: extracted text from `Jamroz 2013.pdf`.
- `outputs/ped_basis_optimizer_full_golden_20260506`: previous full golden run for primitive EPM-like basis swap optimization.
- `docs/full_ped_execplan.md`: PED development history and prior validation notes.

Key Jamroz 2013 evidence from extracted text:

- VEDA automatically proposes local mode coordinates and optimizes them to obtain maximal PED matrix elements, called EPM.
- PED analysis requires `3N-6` linearly independent local mode coordinates for a nonlinear molecule.
- VEDA optimization replaces introductory coordinates with sum or difference coordinates while preserving movement type.
- VEDA separates CH-like and not-CH coordinates for clarity; OH and NH can be included in CH-like treatment.
- VEDA allows user editing, freezing satisfactory coordinates, and reoptimizing the rest.
- Jamroz explicitly states PED interpretation is not unique and requires expert judgement.

## Interfaces and Dependencies

The implementation must preserve these existing interfaces unless a later plan explicitly changes them:

- `InternalCoordinate(name, kind, atoms0, priority, fn, source="primitive")`
- `finite_difference_B(coords_A, internals)`
- `select_independent_coordinates(B, internals, target_rank, tol_abs=...)`
- `optimize_independent_coordinates_for_ped(B, internals, selected_idx, normal_modes, mode_indices, target_rank=...)`
- `build_ped_audit_dataframe(...)`
- `build_ped_v2_force_audit_dataframe(...)`
- `build_wilson_ped_audit_dataframe(...)`
- `build_stage3d_assignment_audit(...)`
- output tables `assignment_audit`, `ped_audit`, `ped_v2_force_audit`, `wilson_ped_audit`, `ped_stage3d_agreement`, and `ped_final_assignment`

New additive interfaces proposed by this plan:

- `ComposedCoordinateDefinition` or equivalent structure.
- `classify_coordinate_optimization_group(...)`.
- `build_composed_coordinate_candidates(...)`.
- `build_composed_B_rows(...)`.
- output table `composed_coordinate_definitions`.
- output table `composed_ped_basis`.

No network access is required for the first backend implementation. NIST validation is optional and should use `$orcaveda-nist-ir` only after local PED and golden validations pass.
