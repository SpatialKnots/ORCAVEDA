# VEDA-Like Composed Coordinate PED Optimization

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA already parses ORCA `.hess` files, builds primitive and functional-group internal coordinates, selects a rank-preserving independent basis, computes Stage 3D assignment audits, writes Wilson GF-style PED diagnostics, and now has a conservative EPM-like basis optimizer that swaps existing coordinates to improve PED localization. The next scientific upgrade is to implement the central mathematical idea described by Jamroz 2013 for VEDA: optimize the internal-coordinate description itself by constructing composed local coordinates, such as sums and differences of same-type stretches, bends, torsions, and out-of-plane coordinates.

After this work, a user should be able to run ORCAVEDA on a `.hess` file and receive, alongside the existing Stage 3D and Wilson/PED outputs, a composed-coordinate PED optimization audit. The audit should show which composed coordinates were created, whether the rank was preserved, how the EPM-like score changed, how many diffuse PED modes improved, and which final assignments became better supported by PED. This must improve scientific interpretation without claiming full VEDA equivalence until the implemented mathematics and validation justify that wording.

## Progress

- [x] (2026-05-06 16:31+05:00) Existing project agent roles, ORCAVEDA skills, current PED code, and Jamroz 2013 text extraction were reviewed.
- [x] (2026-05-06 16:31+05:00) Initial ExecPlan written before starting implementation of composed-coordinate optimization.
- [ ] Implement composed internal-coordinate data model and B-row composition.
- [ ] Implement same-type and CH/not-CH coordinate grouping.
- [ ] Implement first composed coordinate generators for high-value chemical motifs.
- [ ] Implement EPM metrics and optimization gates.
- [ ] Integrate composed PED basis into Wilson/PED outputs without changing Stage 3D baseline.
- [ ] Validate on synthetic cases, H2O/NH3, aromatic/ring cases, and full golden `.hess`.

## Surprises & Discoveries

- Observation: Existing project roles are sufficient for this upgrade. Backend-agent owns code paths and schema changes, Sci-agent owns mathematical and chemical validity, and Frontend-agent is only needed after backend evidence exists.
  Evidence: `AGENTS.md` defines Backend-agent, Frontend-agent, and Sci-agent with the necessary ownership boundaries.

- Observation: ORCAVEDA already has the first safe precursor to VEDA-like optimization: `optimize_independent_coordinates_for_ped(...)` in `src/b_matrix.py` swaps existing independent coordinates to improve PED localization while preserving rank.
  Evidence: Local source inspection of `src/b_matrix.py` and previous full golden run output in `outputs/ped_basis_optimizer_full_golden_20260506`.

- Observation: Jamroz 2013 describes a stronger method than the current ORCAVEDA optimizer. VEDA constructs composed local coordinates from sums or differences, usually preserving coordinate type and separating CH-like motions from heavy-atom motions.
  Evidence: Extracted text in `outputs/jamroz_2013_text.txt`, especially the sections "Different ways for optimization of PED analysis" and "Block diagram".

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

## Outcomes & Retrospective

No implementation outcomes yet. This plan was created before starting the composed-coordinate optimization patch. Update this section after each milestone with the exact commands run, generated files, observed EPM changes, fallback changes, and any scientific limitations discovered.

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
