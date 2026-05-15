# ORCAVEDA UI Visual Refactor - Codex Instructions

## Goal

Upgrade ORCAVEDA's frontend from an internal utility look to a serious computational spectroscopy workstation.

This is a visual and interaction-quality refactor only.

Do not redesign the scientific workflow. Do not change scientific logic, schemas, units, thresholds, assignments, file naming, or artifact contracts.

Primary targets:

- `src/reports.py` - self-contained interactive spectrum viewer
- `src/web_app.py` - upload/result wrapper around the viewer

## Non-Negotiable Constraints

Do not change:

- `write_interactive_spectrum_viewer(...)` public API
- `build_spectrum_payload(...)` schema
- generated JSON schema
- output artifact names
- ORCA `.hess` parsing
- normal mode orientation
- scaling logic
- NIST matching logic
- assignment or PED-like diagnostic logic
- units, thresholds, warnings, or policy strings

Keep terminology scientifically bounded:

- Current assignment layer may be called geometric / weighted independent-coordinate assignment audit.
- Do not call it strict VEDA PED or full Wilson GF PED unless actually implemented.
- Preserve diagnostics. Do not hide failures with broad `try/except`.

## Product Identity

ORCAVEDA should feel like:

- computational spectroscopy workstation
- scientific evidence system
- vibrational-analysis platform
- research-grade diagnostic tool

It should not feel like:

- startup landing page
- generic SaaS dashboard
- Bootstrap admin panel
- educational chemistry demo
- decorative AI product page

Target first impression:

> I am operating a serious computational spectroscopy platform.

## Visual Direction

Use a calm, analytical workstation aesthetic.

Preferred style:

- dark scientific interface
- high information density
- restrained contrast
- precise panel hierarchy
- spectral-analysis visual language
- minimal decorative effects

Avoid:

- pure black backgrounds
- bright white text everywhere
- neon cyberpunk colors
- gradient blobs or orbs
- cartoon molecule imagery
- oversized hero sections
- playful animation

## Theme Tokens

Use CSS variables inside the existing self-contained HTML/CSS.

Recommended dark theme:

```css
--bg: #0b1020;
--panel: #121a2b;
--panel-elevated: #182235;
--panel-soft: #0f1726;

--border: rgba(255,255,255,0.08);
--border-strong: rgba(255,255,255,0.14);

--text: #f3f4f6;
--text-secondary: #9ca3af;
--text-muted: #6b7280;

--accent-blue: #60a5fa;
--accent-cyan: #22d3ee;
--accent-violet: #a78bfa;
--accent-green: #34d399;
--accent-red: #f87171;
--accent-amber: #fbbf24;
```

Use muted accents only. No saturated neon.

## Pipeline Color Semantics

Use colors consistently for UI metadata:

| Pipeline Stage | Accent |
|---|---|
| ORCA `.hess` | blue |
| Normal modes | cyan |
| Coordinate projection | violet |
| Assignments | green |
| Reports / warnings | red or amber |

Allowed uses:

- badges
- tabs
- artifact pills
- diagnostics severity
- selected mode highlights
- chart overlays
- section accents

Do not imply scientific confidence from color unless the underlying data explicitly supports it.

## Typography

Use system-safe fonts unless external assets already exist.

Recommended UI font:

```css
font-family: Inter, "Segoe UI", system-ui, sans-serif;
```

Recommended diagnostic / numeric font:

```css
font-family: "JetBrains Mono", "Cascadia Mono", Consolas, monospace;
```

Apply monospace to:

- frequencies
- intensities
- mode ids
- diagnostics
- logs
- artifact metadata
- numerical tables

Do not import web fonts unless explicitly approved.

## Layout Requirements

Keep the summary-first workstation layout:

1. Compact top toolbar / file selector
2. Summary panel with selected-mode quick facts
3. Main workspace:
   - IR spectrum
   - 3D molecule viewer
4. Peak table
5. Advanced diagnostics collapsed by default

The viewer must remain practical:

- spectrum and 3D areas must be visually prominent
- selected mode must be immediately readable
- advanced PED/NIST/raw diagnostics must remain accessible
- no text overlap
- no clipped controls
- no blank 3D region

## Main Viewer Requirements

The interactive spectrum viewer is the core scientific workspace.

Keep:

- spectrum canvas
- scale factor control
- HWHM control
- y-axis mode control
- show sticks toggle
- invert x-axis toggle
- click/hover mode selection
- peak table selection
- 3Dmol rendering and native fallback

Move or keep collapsed by default:

- NIST reference selector
- scale engine selector
- matching layer selector
- engine tables
- raw PED/composed PED diagnostic details

## Selected Mode Panel

Default visible fields only:

- Mode
- Scaled Frequency
- Original Frequency
- IR Intensity
- Final Assignment
- Final Assignment Source
- Final Assignment Policy
- Final Assignment Warning
- Warnings

Advanced tabs:

- `Summary`
- `Evidence`
- `NIST / Scaling`
- `Raw diagnostics`

Evidence tab may show:

- Stage 3D assignment
- baseline PED-like diagnostic fields
- composed PED-like diagnostic fields
- top contributors
- agreement status
- policy warning
- supporting coordinates

Keep the method language conservative.

## Peak Table

Default columns:

- Mode
- Frequency
- Intensity
- Final Assignment
- Warning

Behavior:

- sort by scaled frequency
- hover updates selected mode
- click pins/updates selected mode
- active row is visually clear
- table remains scrollable
- long labels wrap cleanly

Do not show dense diagnostic columns by default.

## Diagnostics

Diagnostics are important, but not first-screen clutter.

Use collapsed inspectors for:

- advanced diagnostics
- raw logs
- NIST/scaling tables
- evidence details

Diagnostics style:

- dark terminal-like panel
- monospace text
- severity accents

Severity mapping:

- info: muted blue/cyan
- success: green
- warning: amber
- critical/error: muted red

Do not remove diagnostics to make UI cleaner.

## Panels

Use workstation panels, not marketing cards.

Recommended:

```css
background: var(--panel);
border: 1px solid var(--border);
border-radius: 8px;
box-shadow:
  0 8px 24px rgba(0,0,0,0.22),
  inset 0 1px 0 rgba(255,255,255,0.03);
```

Use 6-10px radius by default. Use 12px only for large major panels if it improves hierarchy.

Avoid:

- nested card stacks
- huge rounded cards
- decorative shadows
- glossy gradients

## Buttons And Controls

Controls should feel precise and instrument-grade.

Use:

- compact sizing
- clear hover states
- muted accent borders
- consistent spacing
- visible focus states

Avoid:

- oversized SaaS buttons
- playful animations
- heavy gradients
- icon-only controls without labels/tooltips

## Artifacts

Artifact links should look like scientific resources/modules.

Style as compact action pills:

- interactive viewer
- JSON data
- XLSX report
- run manifest
- NIST reference sets

Do not rename artifact keys or URLs.

## Motion

Use minimal motion only.

Allowed:

- hover transitions
- focus transitions
- collapse/expand transitions if simple
- subtle selected-row highlight

Avoid:

- pulsing glows
- floating panels
- large animated backgrounds
- startup-style motion

Motion must not distract from analysis.

## Technical Implementation Rules

Use the existing architecture:

- self-contained HTML/CSS/JS generated in `src/reports.py`
- small CSS/JS edits
- no React
- no Tailwind
- no shadcn/ui
- no new build pipeline
- no external font/CDN dependency unless explicitly approved

Keep generated HTML standalone.

Use CSS variables and local classes.

Prefer minimal safe patches.

## Accessibility And Responsiveness

Required:

- readable contrast in dark theme
- keyboard-visible focus states
- no text overlap
- no clipped labels
- responsive stacking below tablet width
- stable dimensions for spectrum, 3D viewer, and table
- long scientific labels wrap safely

Check at least:

- desktop around `1400x1100`
- narrower viewport around `900px`
- mobile-ish width if feasible

## Verification

Run:

```powershell
.\.venv312\Scripts\python.exe -m pytest tests\test_interactive_spectrum_viewer.py -q
.\.venv312\Scripts\python.exe -m pytest tests\test_web_app.py tests\test_web_import.py -q
git diff --check
```

If Chrome/browser visual verification is available:

- generate viewer HTML
- open viewer
- verify:
  - summary-first layout exists
  - spectrum visible and nonblank
  - 3D viewer/fallback visible and nonblank
  - peak table readable
  - advanced diagnostics collapsed by default
  - tabs work
  - no obvious overlap/clipping

If visual verification is blocked, report why.

## Acceptance Criteria

The refactor is successful when:

- UI looks like a serious computational spectroscopy workstation
- first screen gives clear summary + spectrum + 3D + peak table
- advanced diagnostics remain accessible but not visually dominant
- no scientific data contract changes
- tests pass
- generated viewer remains self-contained
- terminology remains scientifically conservative

## Required Report Format

After implementation, report:

```text
Changed:
Tests run:
Visual checks:
Limitations:
Verdict:
```
