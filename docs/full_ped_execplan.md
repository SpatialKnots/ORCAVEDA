# Implement Full PED for ORCAVEDA

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA currently has Stage 3D, a geometric and weighted independent-coordinate assignment audit. It is useful for assigning vibrational modes, but it is not a strict VEDA PED or full Wilson GF PED implementation. This plan adds a separate PED layer that decomposes each normal mode into internal-coordinate contributions and reports percentage-like mode composition in a scientifically traceable way.

After this work, a user should be able to run ORCAVEDA on an ORCA `.hess` file and receive, alongside the existing Stage 3D assignment audit, a PED audit showing the dominant bond, angle, torsion, and out-of-plane contributions for each vibrational mode. This should improve interpretation of mixed modes such as carbonyl/ring coupling, amide NH2 deformation, carboxylic acid O-H bending context, and aromatic ring modes without pretending that PED fixes frequency errors or method/scale-factor mismatch.

The initial implementation must be additive. Stage 3D remains available and must not be renamed as full PED. PED results should be emitted as separate outputs until they are validated enough to participate in combined assignment decisions.

## Progress

- [x] (2026-05-05) Initial implementation plan written.
- [ ] Audit existing `.hess` parser, normal-mode orientation, masses, Hessian, and internal-coordinate data flow.
- [ ] Specify the exact PED v1 mathematical definition and output semantics.
- [ ] Implement additive PED module and focused unit tests.
- [ ] Validate PED v1 on small molecules and benchmark molecules.
- [ ] Decide whether and how PED should influence final assignment wording.

## Surprises & Discoveries

- Observation: None yet.
  Evidence: Initial planning only; no PED source implementation has been inspected or changed for this plan.

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

## Outcomes & Retrospective

No outcomes yet. Update this section after each major milestone, especially after the parser audit, first PED prototype, and benchmark validation.

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
