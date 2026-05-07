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
- [ ] (next) Harden composed-coordinate semantic labels for X-H, carbonyl, and carboxylic-context evidence before any policy integration.

## Surprises & Discoveries

- Observation: Existing project roles are sufficient for this upgrade. Backend-agent owns code paths and schema changes, Sci-agent owns mathematical and chemical validity, and Frontend-agent is only needed after backend evidence exists.
  Evidence: `AGENTS.md` defines Backend-agent, Frontend-agent, and Sci-agent with the necessary ownership boundaries.

- Observation: ORCAVEDA already has the first safe precursor to VEDA-like optimization: `optimize_independent_coordinates_for_ped(...)` in `src/b_matrix.py` swaps existing independent coordinates to improve PED localization while preserving rank.
  Evidence: Local source inspection of `src/b_matrix.py` and previous full golden run output in `outputs/ped_basis_optimizer_full_golden_20260506`.

- Observation: Jamroz 2013 describes a stronger method than the current ORCAVEDA optimizer. VEDA constructs composed local coordinates from sums or differences, usually preserving coordinate type and separating CH-like motions from heavy-atom motions.
  Evidence: Extracted text in `outputs/jamroz_2013_text.txt`, especially the sections "Different ways for optimization of PED analysis" and "Block diagram".

- Observation: On the current curated benchmark rows, composed Wilson evidence did not improve semantic status relative to `ped_final_assignment`, even though it improved localization for one acetophenone carbonyl row. Several rows worsened as a separate evidence source because composed top terms lost explicit C-H or O-H semantic classes.
  Evidence: `outputs/composed_ped_benchmark_20260507/composed_ped_comparison_raw.csv` and `outputs/composed_ped_benchmark_20260507/composed_ped_comparison_scaled.csv`.

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

Next planned stage: composed-coordinate semantic hardening. The benchmark comparator showed that composed coordinates can improve localization while losing explicit chemical labels such as C-H, O-H, C=O, or carboxylic acid context. The next patch should inspect the worsened rows in `outputs/composed_ped_benchmark_20260507`, improve composed-coordinate naming and `coordinate_family` propagation from primitive components, add focused tests for X-H and carbonyl/carboxylic labels, regenerate a small benchmark subset, and only then repeat the raw/scaled comparator. This stage must not change `ped_final_assignment`, `ped_stage3d_agreement`, or `assignment_audit` policy.

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
