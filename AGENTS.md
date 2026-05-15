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

ORCAVEDA's current assignment layer is a geometric and weighted independent-coordinate audit. It may be called PED-like, but not strict VEDA PED or full Wilson GF PED unless those methods are actually implemented.

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
- keep assignment-audit language within the PED-like boundary.

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

## Multi-Agent Collaboration

ORCAVEDA is developed by three agents: User (owner), Codex (primary implementer), and Super Z (architecture/review). See `COLLABORATION.md` for the full protocol, branch strategy, and handoff rules.

Every session must start by reading `COLLABORATION.md` and `WORKLOG.md`. Every session must end by appending to `WORKLOG.md`.

## GAP 1: EPM Optimization — Implementation Context

### What Is EPM?

EPM (Eigenvector Projection Maximization) is VEDA's iterative basis-optimization strategy. It swaps internal coordinates in the nonredundant basis to maximize dominant PED contributions. The real VEDA uses force-constant-weighted PED from the GF eigenproblem, not geometric projections.

### Current Implementation (as of 2026-05-15)

The GF-PED-aware EPM optimizer is implemented in `src/wilson_gf.py` (lines ~1148-1408):

1. **`wilson_gf_ped_localization_metrics()`** — GF-PED-aware score function.
   - Builds G and F_internal for the trial basis.
   - Solves the GF eigenproblem.
   - Computes PED percentages per positive mode using `|q_i * (F*q)_i| / sum|q_j * (F*q)_j|`.
   - Returns localization_score = mean_top + 0.25 * median_top - 10.0 * diffuse_fraction.
   - This is the **correct** EPM metric (force-constant-weighted), unlike the geometric-only `ped_basis_localization_metrics()` in `b_matrix.py`.

2. **`optimize_wilson_gf_basis_for_epm()`** — Iterative swap optimizer.
   - Takes `initial_basis_idx`, tries swapping each position with all candidates.
   - Safety checks: G-matrix rank preservation, condition number cap (1e10).
   - Each swap must improve localization_score by at least `improvement_tol` (default 0.1).
   - Up to `max_passes` (default 3) passes over all positions.
   - Returns optimized basis indices and a report dict with swap log.

3. **Integration in `wilson_gf_diagonalization()`** — Opt-in via `epm_optimize=False`.
   - After `_select_conditioned_wilson_basis()`, if `epm_optimize=True`, runs the swap optimizer.
   - Verifies G-rank is preserved before accepting the optimized basis.
   - Adds `"epm_basis_optimized"` or `"epm_optimization_rejected_rank_loss"` to warnings.

4. **CLI flags** — `--epm-optimize`, `--epm-max-passes`, `--epm-improvement-tol` in `orcaveda_cli.py`.

### GAP 1 Status

Current status after Codex validation:

- Completed: EPM core, CLI flags, pipeline wiring through `run_orca_ped_like()`, and focused tests for `H2O_freq.hess` plus `ethene.hess`.
- Not done: CH4-specific test, because no CH4 `.hess` fixture exists in `data/hess`.
- Optional future work: emit a dedicated EPM swap-log CSV and benchmark large-molecule runtime.

Historical notes from the Super Z handoff follow; completed items may be listed there for traceability.

### Historical GAP 1 Handoff Notes

- **Wire `epm_optimize` through `run_orca_ped_like()`**: The `ORCAVEDA_patched_stage3D_v5_0.py` calls `wilson_gf_diagonalization()` at lines 2085 and 2192 without passing `epm_optimize`. These call sites need to forward the CLI flag.
- **Write EPM-specific tests**: Test `wilson_gf_ped_localization_metrics()` and `optimize_wilson_gf_basis_for_epm()` with small molecules (H2O, CH4).
- **Diagnostic CSV output**: Optionally emit EPM swap log as a CSV alongside existing VEDA-like outputs.
- **Performance benchmark**: The swap optimizer solves a full GF eigenproblem per candidate per position. For large molecules this may be slow. Consider caching or early-exit heuristics.

### Key Design Decisions

1. **Swap-based, not combinatorial**: O(n_basis * n_candidates) per pass, not C(n, k). This is tractable for typical molecules (3N-6 <= ~30 basis coordinates, ~50-100 candidates).

2. **GF-PED-aware scoring, not geometric**: The old `ped_basis_localization_metrics()` in `b_matrix.py` uses `(B_unit @ mode_unit)^2` — geometric projections only, no force constants. The new `wilson_gf_ped_localization_metrics()` uses the actual GF eigenproblem with F_internal, which is what VEDA does.

3. **Safety-first**: G-rank and condition-number checks at every swap. If rank is lost, the optimization is rejected entirely.

4. **Opt-in, never default**: `epm_optimize=False` by default. This cannot break existing workflows.

### Relationship to `b_matrix.py` Optimizer

`b_matrix.optimize_independent_coordinates_for_ped()` is the OLD geometric EPM optimizer. It is NOT connected to the VEDA-like pipeline and uses the WRONG objective function. Do NOT connect it to `wilson_gf_diagonalization()`. The new `optimize_wilson_gf_basis_for_epm()` replaces it for GF-PED-aware optimization.
