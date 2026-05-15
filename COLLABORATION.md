# ORCAVEDA Collaboration Protocol

This document defines the collaboration protocol for the three agents working on ORCAVEDA:
- **User** (project owner, domain expert in quantum chemistry / VEDA)
- **Codex** (OpenAI coding agent, primary implementation worker)
- **Super Z** (GLM-based assistant, architecture / documentation / review)

## Why This File Exists

Super Z and Codex are stateless: they lose all context between sessions. This file, together with `AGENTS.md` and `WORKLOG.md`, acts as shared memory. Every agent must read these files before starting work and update them after finishing.

## Core Rules

1. **Every session starts with reading**: `COLLABORATION.md`, `AGENTS.md`, `WORKLOG.md`.
2. **Every session ends with writing**: append to `WORKLOG.md` with what was done.
3. **No agent modifies another agent's in-progress work without coordination.**
4. **User is the only one who merges to `main`.**
5. **Evidence over claims**: do not report PASS without actual test results.

## Branch Strategy

```
main (stable, user-merged)
  ├── gap1/epm-optimization        ← EPM basis optimizer integration
  ├── gap2/analytical-bmatrix      ← Analytical B-matrix derivatives
  ├── gap3/veda-validation         ← VEDA reference output validation
  └── docs/agents-context          ← Documentation and collaboration files
```

- Each GAP gets its own branch.
- Only one agent works on a branch at a time.
- User creates branches, agents work within assigned branches.
- User merges after review.

## Current GAP Status

| GAP | Description | Branch | Status | Agent | Last Updated |
|-----|-------------|--------|--------|-------|-------------|
| GAP 1 | EPM optimization (GF-PED-aware basis swap) | `codex-gap1-epm-optimization` | Complete and merged to `main` | Codex + Super Z | 2026-05-15 |
| GAP 2 | Analytical B-matrix | `gap2/analytical-bmatrix` | Not started | — | — |
| GAP 3 | VEDA reference validation | `codex-gap3-veda-validation` | In progress | Codex | 2026-05-15 |

## Agent Capabilities

### Codex
- Can read/write the full repository
- Can run tests and commit
- Best for: iterative implementation, refactoring, running test suites
- Limitation: may not have deep chemistry context

### Super Z
- Can clone repo, write code, generate patches
- Cannot push to GitHub or commit directly
- Best for: architecture design, code review, documentation, complex math
- Limitation: no persistence between sessions, works in isolated `/home/z/my-project/`

### User
- Full repository access, merge authority
- Domain expertise in VEDA and quantum chemistry
- Reviews and accepts/rejects all changes
- Creates branches and manages GitHub

## Handoff Protocol

When switching agents (e.g., Codex rate limit hit, Super Z takes over):

1. Current agent appends to `WORKLOG.md` with: what was done, what's next, any blockers.
2. Next agent reads `WORKLOG.md`, `COLLABORATION.md`, `AGENTS.md` before starting.
3. Next agent checks for uncommitted changes (`git diff`).
4. Next agent resumes from the last incomplete item.

## Patch Workflow (Super Z -> User)

Since Super Z cannot push to GitHub directly:

1. Super Z generates `.patch` files from local changes.
2. User downloads patches and applies with `git apply`.
3. User reviews, tests, and commits.

Patch generation command:
```bash
git diff HEAD > /home/z/my-project/download/orcaveda_epm_YYYYMMDD.patch
```

## Communication Conventions

- Use English in code, docstrings, and commit messages.
- Use the user's preferred language in conversation (Russian in current sessions).
- Document all design decisions in `AGENTS.md` or `WORKLOG.md`, not just in chat.
- When uncertain about chemistry, ask the user. Do not guess.

## File Index

| File | Purpose | Who Updates |
|------|---------|-------------|
| `COLLABORATION.md` | This file: collaboration protocol | Any agent (user approves) |
| `AGENTS.md` | Agent roles and project rules | Primarily Codex/User |
| `WORKLOG.md` | Chronological work log | Every agent after each session |
| `PLANS.md` | ExecPlan format and rules | Primarily Codex/User |
