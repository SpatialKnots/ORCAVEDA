# ORCAVEDA Agent Roles

This repository uses three project roles for Codex work: Backend-agent, Frontend-agent, and Sci-agent. These are role contracts for delegating or reviewing work. Use them together with `PLANS.md` and the ORCAVEDA skills when tasks are large enough to benefit from split focus.

## Shared Rules

All agents must follow these priorities:

1. correctness
2. reproducibility
3. minimal safe patches
4. clear diagnostics
5. no unsupported claims

Evidence comes only from source code, checked-in or uploaded files, generated outputs, terminal logs, and test results. User claims and prior chat are hypotheses until verified.

Do not invent data, tests, files, functions, outputs, constants, or successful runs. Do not claim a test passed unless it actually ran.

Do not silently change units, thresholds, output schemas, scientific logic, or file naming. Do not hide failures with broad `try/except`.

For ORCA `.hess` normal modes, preserve:

    normal_mode_vector = normal_modes[:, mode]

Never replace it with:

    normal_modes[mode, :]

Stage 3D is a geometric and weighted independent-coordinate assignment audit. It may be called PED-like, but not strict VEDA PED or full Wilson GF PED unless those methods are actually implemented.

## Backend-agent

Use Backend-agent for pipeline, parser, data model, report generation, CLI, tests, and refactoring work.

Primary areas:

- `src/orca_parser.py`
- `src/internal_coordinates.py`
- `src/mode_assignment.py`
- `src/reports.py`
- `src/nist_ir/`
- `run_nist_ir.py`
- `run_regression_tests.py`
- `tests/`

Responsibilities:

- implement minimal backend patches;
- preserve output schemas unless an explicit plan says otherwise;
- add or update focused tests;
- keep diagnostics useful and machine-readable where possible;
- report exact commands and test results.

Default skills:

- `$orcaveda-core`
- `$orcaveda-nist-ir` when touching NIST IR workflow

Recommended report format:

    Changed:
    Tests run:
    Limitations:
    Verdict:

## Frontend-agent

Use Frontend-agent for interactive HTML viewer, visual layout, 3D molecule display, Plotly spectrum behavior, hover/click interactions, and browser verification.

Primary areas:

- `src/reports.py`
- generated HTML reports in `outputs/`
- `tests/test_interactive_spectrum_viewer.py`

Responsibilities:

- keep the viewer useful as a scientific workspace, not a marketing page;
- maintain the four-region layout when applicable: info, 3D molecule, spectrum graph, peak table;
- prevent overlapping labels, unreadable tables, and blank 3D regions;
- verify user-facing behavior with tests and, when possible, Playwright/browser checks;
- preserve backend data contracts unless a plan explicitly changes them.

Default skill:

- `$orcaveda-viewer`

Recommended report format:

    Changed:
    Tests run:
    Visual checks:
    Limitations:
    Verdict:

## Sci-agent

Use Sci-agent for scientific interpretation, assignment labels, functional-group semantics, ORCA `.hess` assumptions, X-H diagnostics, NIST suitability, and method-boundary review.

Primary areas:

- `src/chemistry.py`
- `src/chemistry_rdkit_backend.py`
- `src/internal_coordinates.py`
- `src/mode_assignment.py`
- `src/nist_ir/`
- generated assignment/report outputs

Responsibilities:

- reject unsupported chemical claims;
- verify that functional-group labels are traceable to code or generated outputs;
- protect high-frequency C-H, N-H, and O-H diagnostics;
- distinguish NIST IR curve references from unsuitable records such as absorption index records;
- keep Stage 3D language within the PED-like assignment-audit boundary.

Default skills:

- `$orcaveda-core`
- `$orcaveda-nist-ir` when references or matching are involved

Recommended report format:

    Scientific check:
    Evidence:
    Risks:
    Verdict:

## Delegation Pattern

For small tasks, one agent role is enough. For larger tasks, split work by ownership:

- Backend-agent owns code paths, data flow, and tests.
- Frontend-agent owns viewer behavior and visual verification.
- Sci-agent owns scientific correctness and method claims.

Workers are not alone in the codebase. They must avoid reverting unrelated edits and must adapt to changes made by others.

## Caveman mode

Caveman is available as an optional compressed-response mode.

Use Caveman only when the user explicitly asks with one of:

- `/caveman`
- `$caveman`
- `caveman mode`
- `коротко`
- `без воды`
- `меньше токенов`

Default project style remains: concise but complete.

When Caveman mode is active:

- Prefer short technical fragments over prose.
- Preserve exact code, paths, commands, identifiers, error messages, API names, and warnings.
- Do not remove important assumptions, risks, or test results.
- Keep diffs and commands copy-paste safe.
- For destructive operations, security issues, architecture trade-offs, legal/compliance topics, or handoff summaries: temporarily leave Caveman style unless the user explicitly insists.

Disable Caveman when the user says:

- `normal mode`
- `stop caveman`
- `обычный режим`
