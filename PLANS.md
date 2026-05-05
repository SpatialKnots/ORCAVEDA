# ORCAVEDA ExecPlan Rules

This file defines how to write and maintain execution plans for ORCAVEDA. An execution plan, or ExecPlan, is a living design document that a coding agent or human contributor can follow from a clean checkout to deliver a working, observable change.

Use an ExecPlan for non-trivial work: scientific logic, ORCA `.hess` parsing, mode assignment, NIST IR matching, regression infrastructure, report schemas, interactive viewers, or broad refactors. For tiny edits, direct implementation is fine, but the same evidence and validation rules still apply.

## Non-Negotiable Rules

Every ExecPlan must be self-contained. The reader must not need chat history, prior plans, or unstated repository knowledge. Define ORCAVEDA-specific terms in plain language the first time they appear.

Every ExecPlan must produce demonstrably working behavior. Explain what the user can do after the change, what command to run, and what output, file, plot, report, or UI behavior proves success.

Every ExecPlan is a living document. Update it as work proceeds, when facts change, when tests fail, and when design decisions are made. A future agent must be able to restart from only the current ExecPlan and the repository.

Evidence comes only from source code, uploaded or checked-in files, generated outputs, terminal logs, and test results. User claims and prior chat context are hypotheses until verified. Do not invent data, functions, files, outputs, constants, tests, or successful runs.

Do not claim a test passed unless it actually ran. Syntax checks and imports prove only syntax and importability; they do not validate chemistry, assignments, NIST matching, or viewer behavior.

Prefer the smallest safe patch. Do not silently change units, thresholds, output schemas, scientific logic, file naming, or report semantics. If a plan proposes such a change, state it explicitly and explain why it is necessary.

## Required ExecPlan Format

When an ExecPlan is pasted into chat, wrap the entire plan in one fenced code block labeled `md`. Do not nest triple backtick fences inside it; indent commands, code snippets, transcripts, and diffs instead.

When an ExecPlan is saved as a Markdown file whose whole content is the plan, omit the surrounding triple backticks.

Use plain prose. Prefer short narrative paragraphs over large tables. Checklists are required only in `Progress`.

Every ExecPlan must include these sections, in this order:

1. `Purpose / Big Picture`
2. `Progress`
3. `Surprises & Discoveries`
4. `Decision Log`
5. `Outcomes & Retrospective`
6. `Context and Orientation`
7. `Plan of Work`
8. `Concrete Steps`
9. `Validation and Acceptance`
10. `Idempotence and Recovery`
11. `Artifacts and Notes`
12. `Interfaces and Dependencies`

`Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` are mandatory living sections. Keep them current during implementation.

## ORCAVEDA Scientific Boundaries

The current trusted project baseline is Stage 3D v5.0. The main legacy entrypoint is `src/ORCAVEDA_patched_stage3D_v5_0.py`. The modular implementation lives in `src/`, including `orca_parser.py`, `internal_coordinates.py`, `mode_assignment.py`, `reports.py`, `chemistry.py`, `chemistry_rdkit_backend.py`, `scale_factor_engine.py`, and `src/nist_ir/`.

Stage 3D is a geometric and weighted independent-coordinate assignment audit. It may be called a PED-like assignment-audit layer. Do not call it strict VEDA PED, full Wilson GF PED, or publication-grade universal benchmark validation unless that method is actually implemented and validated.

For ORCA `.hess` normal modes, preserve the orientation rule:

    normal_mode_vector = normal_modes[:, mode]

Never replace it with:

    normal_modes[mode, :]

The ORCA block matrix parser must support single-column headers using this behavior:

    if len(parts) >= 1 and all(re.fullmatch(r"\d+", x) for x in parts):
        current_cols = list(map(int, parts))
        continue

Do not regress this parser behavior.

High-frequency X-H modes must not remain unassigned without diagnostics. X-H means C-H, N-H, and O-H stretching regions. Protected X-H logic must not create false C-H, N-H, or O-H relabeling.

Functional group claims must be traceable to code, RDKit output, parsed molecular data, or generated reports. If a group or assignment is absent from source output, say `Not reported in source.` If a feature is absent from code, say `Not implemented in source.`

NIST IR references must distinguish suitable IR curve references from non-curve records such as absorption index records. Non-suitable references may remain in manifests, but must not be used as normal peak-matching references.

Viewer plans must treat the HTML output as a user-facing scientific artifact. Layout, labels, hover behavior, mode details, and final peak wording must be verified visually and, where possible, with automated browser checks.

## Planning Requirements for ORCAVEDA Work

Begin with the user-visible purpose. Say what someone will be able to do after the change, for example compare a calculated spectrum to a NIST reference, inspect mode assignments in an HTML viewer, or run a regression set and see PASS/FAIL results.

Orient the reader to the relevant repository paths. Name files relative to the repository root. Explain how the touched modules fit together before prescribing edits.

Resolve ambiguity inside the plan. If a design choice affects scientific meaning, output schemas, thresholds, or naming, record the decision and rationale in `Decision Log`.

Milestones must be independently verifiable. Each milestone should end with a command, output file, test result, or UI observation that proves progress.

When prototyping is needed, label the milestone as `prototyping`. Keep prototypes additive, testable, and disposable. State the criterion for promoting or discarding the prototype.

## Validation Defaults

Choose validation based on the touched area. Always run focused tests for the changed subsystem, then run broader tests when the change can affect shared behavior.

For ORCA parsing, internal coordinates, mode assignment, and core reports, prefer:

    .\.venv312\Scripts\python.exe -m pytest tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q

For chemistry semantics and functional groups, prefer:

    .\.venv312\Scripts\python.exe -m pytest tests\test_chemistry_backend.py tests\test_supramolecular_chemistry.py tests\test_golden_rdkit_outputs.py -q

For NIST IR matching and reference handling, prefer:

    .\.venv312\Scripts\python.exe -m pytest tests\test_nist_ir_matching.py tests\test_nist_ir_compare.py tests\test_nist_ir_pipeline.py -q

For interactive spectrum viewer changes, prefer:

    .\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py -q

For scale-factor work, prefer:

    .\.venv312\Scripts\python.exe -m pytest tests\test_scale_factor_engine.py -q

For full regression infrastructure, run:

    .\.venv312\Scripts\python.exe run_regression_tests.py --outdir outputs\regression_live --expectations expectations\regression_expectations_stage3D_v5_0.json

If a listed command cannot run, record why in the ExecPlan and use the closest safe substitute. Never report PASS without the actual command and result.

## Skeleton

Use this skeleton for new ORCAVEDA ExecPlans:

    # <Short, action-oriented title>

    This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

    This plan follows `PLANS.md` at the repository root.

    ## Purpose / Big Picture

    Explain what the user gains and how they can see it working.

    ## Progress

    - [ ] (YYYY-MM-DD HH:MMZ) Initial plan written.

    ## Surprises & Discoveries

    - Observation: None yet.
      Evidence: Initial planning only.

    ## Decision Log

    - Decision: Initial approach selected.
      Rationale: Explain why this approach is the smallest safe path.
      Date/Author: YYYY-MM-DD / Codex

    ## Outcomes & Retrospective

    No outcomes yet. Update this after each major milestone and at completion.

    ## Context and Orientation

    Describe the relevant current files, outputs, tests, and terms as if the reader is new to ORCAVEDA.

    ## Plan of Work

    Describe the sequence of edits in prose. Name files and functions precisely when needed.

    ## Concrete Steps

    State exact commands and working directory. Include concise expected output when useful.

    ## Validation and Acceptance

    State the tests, generated files, viewer behavior, or command outputs that prove the change works.

    ## Idempotence and Recovery

    Explain which steps are safe to repeat and how to recover from partial failure.

    ## Artifacts and Notes

    Add concise transcripts, report snippets, or diffs that prove success.

    ## Interfaces and Dependencies

    Name any APIs, schemas, file formats, dependencies, function signatures, or output contracts that must exist after implementation.

## Reporting While Implementing an ExecPlan

Keep responses concise and evidence-first. After a patch or milestone, report:

    Changed:
    Tests run:
    Limitations:
    Verdict:

Use PASS, FAIL, or WARN when summarizing validation. Do not paste long audit tables unless explicitly requested.

When revising an ExecPlan, update all affected sections and add a short note in `Decision Log` or `Outcomes & Retrospective` explaining what changed and why.
