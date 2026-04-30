from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Dict, Sequence

import pandas as pd


def safe_output_stem(name: str) -> str:
    stem = Path(str(name)).name
    for suffix in ("_freq.hess", ".hess", "_freq.out", ".out"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    stem = re.sub(r"[^A-Za-z0-9._+-]+", "_", stem).strip("_")
    return stem or "ORCAVEDA_output"


def output_prefix_for_hess_paths(paths: Sequence[str | Path]) -> str:
    stems = [safe_output_stem(str(path)) for path in paths]
    if not stems:
        return "ORCAVEDA_output"
    if len(stems) == 1:
        return stems[0]
    joined = "__".join(stems[:3])
    if len(stems) > 3:
        joined += f"__plus_{len(stems) - 3}_files"
    return f"{joined}__multi_file_{len(stems)}"


def normalize_sheet_name(name: str) -> str:
    bad = set("[]:*?/\\")
    clean = "".join("_" if c in bad else c for c in name)
    return clean[:31] if len(clean) > 31 else clean


def load_nist_reference_set(manifest_path: str | Path) -> Dict[str, object]:
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_dir = manifest_path.parent

    spectra = []
    for item in manifest.get("reference_spectra", []):
        csv_path = manifest_dir / Path(item["csv"]).name if not Path(item["csv"]).is_absolute() else Path(item["csv"])
        meta_path = manifest_dir / Path(item["meta_json"]).name if not Path(item["meta_json"]).is_absolute() else Path(item["meta_json"])
        spectrum_df = pd.read_csv(csv_path, encoding="utf-8")
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        spectra.append(
            {
                "index": str(item.get("index", "")),
                "phase_tag": str(item.get("phase_tag", "")),
                "phase_label": str(item.get("phase_label", "")),
                "selection_priority": int(item.get("selection_priority", 0)),
                "description": str(item.get("description", "")),
                "y_units": str(meta.get("jcamp_metadata", {}).get("YUNITS", "")),
                "state": str(meta.get("jcamp_metadata", {}).get("STATE", "")),
                "csv": str(csv_path),
                "points": [
                    {
                        "x": float(row["wavenumber_cm-1"]),
                        "y": float(row["intensity"]),
                    }
                    for _, row in spectrum_df.iterrows()
                ],
            }
        )

    return {
        "inchikey": str(manifest.get("inchikey", "")),
        "canonical_smiles": str(manifest.get("canonical_smiles", "")),
        "reference_spectra": spectra,
        "preferred_reference": manifest.get("preferred_reference", {}),
    }


def attach_nist_reference_set(
    payload: Dict[str, object],
    manifest_path: str | Path,
    *,
    file_title: str | None = None,
) -> Dict[str, object]:
    updated = dict(payload)
    refs = dict(updated.get("nist_reference_sets", {}))
    target_title = file_title
    if target_title is None:
        files = updated.get("files", [])
        if len(files) == 1:
            target_title = str(files[0].get("title", ""))
    if not target_title:
        raise ValueError("file_title is required when payload has multiple files")

    refs[str(target_title)] = load_nist_reference_set(manifest_path)
    updated["nist_reference_sets"] = refs
    return updated


def write_xlsx_report(report_tables: Dict[str, pd.DataFrame], xlsx_path: str | Path) -> Path:
    xlsx_path = Path(xlsx_path)
    xlsx_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        import xlsxwriter  # noqa: F401
        engine = "xlsxwriter"
    except ModuleNotFoundError:
        try:
            import openpyxl  # noqa: F401
            engine = "openpyxl"
        except ModuleNotFoundError:
            print("WARNING: neither xlsxwriter nor openpyxl is installed; XLSX report skipped.")
            return xlsx_path

    with pd.ExcelWriter(xlsx_path, engine=engine) as writer:
        if engine == "xlsxwriter":
            workbook = writer.book
            header_fmt = workbook.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1})
            warning_fmt = workbook.add_format({"bg_color": "#FFF2CC"})
            critical_fmt = workbook.add_format({"bg_color": "#F4CCCC"})
            number_fmt = workbook.add_format({"num_format": "0.000"})
            sci_fmt = workbook.add_format({"num_format": "0.00E+00"})
        else:
            header_fmt = warning_fmt = critical_fmt = number_fmt = sci_fmt = None

        for raw_name, df in report_tables.items():
            if df is None:
                continue
            if not isinstance(df, pd.DataFrame):
                df = pd.DataFrame(df)

            sheet = normalize_sheet_name(str(raw_name).replace(".csv", ""))
            df.to_excel(writer, index=False, sheet_name=sheet)
            if engine != "xlsxwriter":
                continue

            ws = writer.sheets[sheet]
            for col_idx, col_name in enumerate(df.columns):
                ws.write(0, col_idx, col_name, header_fmt)
            ws.freeze_panes(1, 0)
            ws.autofilter(0, 0, max(len(df), 1), max(len(df.columns) - 1, 0))

            for col_idx, col_name in enumerate(df.columns):
                series = df[col_name].astype(str) if len(df) else pd.Series([str(col_name)])
                width_quantile = series.str.len().quantile(0.90) if len(series) else 10
                width_basis = int(width_quantile) if pd.notna(width_quantile) else len(str(col_name))
                width = min(max(len(str(col_name)), width_basis) + 2, 60)
                ws.set_column(col_idx, col_idx, width)

            for col_idx, col_name in enumerate(df.columns):
                low = str(col_name).lower()
                if "condition" in low:
                    ws.set_column(col_idx, col_idx, 14, sci_fmt)
                elif any(key in low for key in ["score", "freq", "intensity", "angle", "rha", "rda"]):
                    ws.set_column(col_idx, col_idx, 14, number_fmt)

            if "warnings" in df.columns and len(df):
                wcol = df.columns.get_loc("warnings")
                ws.conditional_format(1, wcol, len(df), wcol, {"type": "text", "criteria": "containing", "value": "negative", "format": critical_fmt})
                ws.conditional_format(1, wcol, len(df), wcol, {"type": "text", "criteria": "containing", "value": "near_degenerate", "format": warning_fmt})

            if "system_flags" in df.columns and len(df):
                fcol = df.columns.get_loc("system_flags")
                ws.conditional_format(1, fcol, len(df), fcol, {"type": "text", "criteria": "not containing", "value": "", "format": warning_fmt})

    return xlsx_path


def build_spectrum_payload(
    hess_list,
    assignment_audit: pd.DataFrame | None = None,
    *,
    nist_reference_sets: Dict[str, object] | None = None,
) -> Dict[str, object]:
    from chemistry import (
        build_connectivity as chemistry_build_connectivity,
        classify_system as chemistry_classify_system,
        formula_string as chemistry_formula_string,
        split_fragments as chemistry_split_fragments,
    )

    audit_df = assignment_audit if isinstance(assignment_audit, pd.DataFrame) else pd.DataFrame()
    files = []

    for hess in hess_list:
        rows = []
        file_audit = audit_df[audit_df.get("Filename", pd.Series(dtype=str)).astype(str) == str(hess.filename)].copy() if not audit_df.empty and "Filename" in audit_df.columns else pd.DataFrame()
        audit_by_mode = {}
        if not file_audit.empty and "mode" in file_audit.columns:
            for _, row in file_audit.iterrows():
                try:
                    audit_by_mode[int(row["mode"])] = row
                except Exception:
                    continue

        for mode, (freq, intensity) in enumerate(zip(hess.frequencies_cm1, hess.ir_intensities)):
            if not pd.notna(freq) or float(freq) <= 0.0:
                continue
            audit_row = audit_by_mode.get(mode)
            rows.append(
                {
                    "mode": int(mode),
                    "frequency_cm1": float(freq),
                    "intensity": float(intensity),
                    "assignment": str(audit_row.get("functional_group_assignment", "")) if audit_row is not None else "",
                    "top_internal_coordinates": str(audit_row.get("top_internal_coordinates", "")) if audit_row is not None else "",
                    "warnings": str(audit_row.get("warnings", "")) if audit_row is not None else "",
                }
            )

        rows.sort(key=lambda row: row["frequency_cm1"])
        bonds = chemistry_build_connectivity(hess.atoms, hess.coords_A)
        geometry_atoms = [
            {
                "index": int(idx),
                "element": str(atom),
                "x": float(coord[0]),
                "y": float(coord[1]),
                "z": float(coord[2]),
            }
            for idx, (atom, coord) in enumerate(zip(hess.atoms, hess.coords_A))
        ]
        geometry_bonds = [
            {
                "i": int(i),
                "j": int(j),
                "distance_A": float(distance),
            }
            for i, j, distance in bonds
        ]
        fragments = chemistry_split_fragments(len(hess.atoms), bonds)
        positive_freqs = [row["frequency_cm1"] for row in rows]
        files.append(
            {
                "filename": str(hess.filename),
                "title": safe_output_stem(str(hess.filename)),
                "summary": {
                    "formula": chemistry_formula_string(hess.atoms),
                    "natoms": int(len(hess.atoms)),
                    "total_modes": int(len(hess.frequencies_cm1)),
                    "positive_mode_count": int(len(rows)),
                    "system_type": chemistry_classify_system(fragments),
                    "fragment_count": int(len(fragments)),
                    "frequency_min_cm1": float(min(positive_freqs)) if positive_freqs else None,
                    "frequency_max_cm1": float(max(positive_freqs)) if positive_freqs else None,
                },
                "modes": rows,
                "geometry": {
                    "atoms": geometry_atoms,
                    "bonds": geometry_bonds,
                },
            }
        )

    return {
        "viewer_title": "ORCAVEDA Interactive IR Spectrum",
        "default_scale_factor": 1.0,
        "default_lorentz_hwhm": 12.0,
        "files": files,
        "nist_reference_sets": nist_reference_sets or {},
    }


def write_interactive_spectrum_viewer(
    payload: Dict[str, object],
    html_path: str | Path,
    json_path: str | Path | None = None,
) -> Path:
    html_path = Path(html_path)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    if json_path is not None:
        json_path = Path(json_path)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    payload_json = json.dumps(payload, ensure_ascii=False)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ORCAVEDA Interactive IR Spectrum</title>
  <script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
  <style>
    :root {{
      --bg: #f5f2eb;
      --panel: #fffdf8;
      --ink: #1f2933;
      --muted: #677483;
      --accent: #0f766e;
      --accent-2: #9a3412;
      --line: #0b5670;
      --grid: #d8d2c5;
      --sticks: rgba(15, 118, 110, 0.34);
      --border: #d5cec1;
      --soft: #f0ece3;
      --shadow: 0 18px 44px rgba(31, 41, 51, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(15, 118, 110, 0.08), transparent 24%),
        radial-gradient(circle at bottom right, rgba(154, 52, 18, 0.08), transparent 28%),
        linear-gradient(180deg, #f8f6f1 0%, var(--bg) 100%);
      color: var(--ink);
      font-family: "Segoe UI", "Trebuchet MS", sans-serif;
    }}
    .wrap {{
      max-width: 1600px;
      margin: 0 auto;
      padding: 14px 18px 28px;
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      gap: 12px;
      margin-bottom: 14px;
      padding: 10px 14px;
      background: rgba(255,255,255,0.82);
      border: 1px solid var(--border);
      border-radius: 18px;
      box-shadow: var(--shadow);
      flex-wrap: wrap;
    }}
    .toolbar h1 {{
      margin: 0;
      font-size: 20px;
      letter-spacing: 0.01em;
    }}
    .toolbar p {{
      margin: 0;
      color: var(--muted);
      font-size: 12px;
      flex: 1 1 320px;
    }}
    .toolbar select {{
      min-width: 280px;
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 8px 12px;
      background: white;
      color: var(--ink);
      font-size: 13px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      grid-template-rows: auto auto;
      gap: 16px;
      align-items: start;
    }}
    .panel-spectrum {{
      grid-column: 1 / -1;
    }}
    .info3d-stack {{
      display: grid;
      grid-template-rows: auto auto;
      gap: 16px;
      min-height: 0;
      align-content: start;
    }}
    .panel {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 20px;
      box-shadow: var(--shadow);
      overflow: hidden;
      display: flex;
      flex-direction: column;
      min-height: 0;
    }}
    .panel-head {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 12px 16px;
      border-bottom: 1px solid var(--border);
      background: linear-gradient(180deg, rgba(255,255,255,0.96), rgba(247,244,237,0.96));
    }}
    .panel-head h2 {{
      margin: 0;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--muted);
    }}
    .panel-body {{
      padding: 14px;
      min-height: 0;
      flex: 1 1 auto;
    }}
    .panel-body.info-body {{
      overflow: visible;
    }}
    .panel-body.viewer-body {{
      display: flex;
      padding: 16px;
    }}
    .panel-body.spectrum-body {{
      display: flex;
      flex-direction: column;
      gap: 8px;
      overflow: hidden;
    }}
    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 8px 12px;
      margin-bottom: 10px;
    }}
    .kv {{
      padding: 8px 10px;
      border: 1px solid #ebe5d9;
      border-radius: 14px;
      background: rgba(255,255,255,0.72);
      font-size: 12px;
      line-height: 1.28;
      overflow-wrap: anywhere;
      word-break: break-word;
    }}
    .kv strong {{
      display: block;
      margin-bottom: 2px;
      color: var(--muted);
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }}
    .mode-card {{
      border: 1px solid var(--border);
      border-radius: 16px;
      background: linear-gradient(180deg, rgba(255,255,255,0.88), rgba(246,243,236,0.92));
      padding: 10px;
      overflow-wrap: anywhere;
      word-break: break-word;
    }}
    .mode-card h3 {{
      margin: 0 0 10px 0;
      font-size: 12px;
      color: var(--accent);
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }}
    .mode-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 8px 10px;
      font-size: 12px;
    }}
    .mode-grid .wide {{
      grid-column: 1 / -1;
    }}
    .controls {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 6px 12px;
      margin-bottom: 4px;
    }}
    .control {{
      display: grid;
      gap: 4px;
    }}
    .control label {{
      color: var(--muted);
      font-size: 10px;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      font-weight: 700;
    }}
    .control input[type="range"], .control select {{
      width: 100%;
      transform: scaleY(0.88);
      transform-origin: center;
    }}
    .control select {{
      padding: 6px 10px;
      border: 1px solid var(--border);
      border-radius: 999px;
      background: white;
      color: var(--accent);
      font-size: 11px;
    }}
    .value {{
      color: var(--accent);
      font-weight: 700;
    }}
    .checkrow {{
      display: flex;
      gap: 14px;
      align-items: center;
      color: var(--muted);
      font-size: 11px;
      margin-bottom: 2px;
      flex-wrap: wrap;
    }}
        .ghost-btn {{
          border: 1px solid var(--border);
          border-radius: 999px;
          padding: 5px 10px;
          background: white;
      color: var(--accent);
      font-size: 11px;
          cursor: pointer;
        }}
        .button-row {{
          display: flex;
          gap: 8px;
          align-items: center;
          flex-wrap: wrap;
        }}
        .fit-summary {{
          color: var(--muted);
          font-size: 11px;
          line-height: 1.3;
          margin-bottom: 2px;
        }}
    #chart {{
      width: 100%;
      height: 100%;
      display: block;
      border-radius: 16px;
      border: 1px solid var(--border);
      background: linear-gradient(180deg, rgba(255,255,255,0.96), rgba(248,245,239,0.98));
    }}
    .chart-wrap {{
      position: relative;
      flex: 1 1 auto;
      min-height: 0;
      height: 470px;
      border-radius: 16px;
      overflow: hidden;
    }}
    .chart-tooltip {{
      position: absolute;
      z-index: 5;
      min-width: 180px;
      max-width: 320px;
      padding: 8px 10px;
      border: 1px solid rgba(15, 118, 110, 0.22);
      border-radius: 12px;
      background: rgba(255, 253, 248, 0.96);
      box-shadow: 0 10px 26px rgba(31, 41, 51, 0.16);
      color: var(--ink);
      font-size: 11px;
      line-height: 1.3;
      pointer-events: none;
      opacity: 0;
      transform: translateY(6px);
      transition: opacity 120ms ease, transform 120ms ease;
    }}
    .chart-tooltip.visible {{
      opacity: 1;
      transform: translateY(0);
    }}
    .chart-tooltip strong {{
      display: block;
      margin-bottom: 3px;
      color: var(--accent);
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}
    #moleculeViewer {{
      position: relative;
      flex: 1 1 auto;
      width: 100%;
      min-height: 500px;
      border: 1px solid var(--border);
      border-radius: 16px;
      overflow: hidden;
      background:
        radial-gradient(circle at 30% 20%, rgba(15, 118, 110, 0.08), transparent 28%),
        radial-gradient(circle at 70% 80%, rgba(217, 119, 6, 0.08), transparent 32%),
        linear-gradient(180deg, rgba(255,255,255,0.96), rgba(246,244,239,0.98));
    }}
    #moleculeViewer canvas {{
      width: 100% !important;
      height: 100% !important;
      display: block;
      border: 0;
      background: transparent;
    }}
    .viewer-actions {{
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }}
    .viewer-actions select,
    .viewer-actions button {{
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 6px 10px;
      background: white;
      color: var(--accent);
      font-size: 11px;
    }}
    .table-wrap {{
      overflow: auto;
      height: 520px;
      max-height: 520px;
      border: 1px solid var(--border);
      border-radius: 16px;
      background: rgba(255,255,255,0.72);
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
      table-layout: fixed;
    }}
    thead th {{
      position: sticky;
      top: 0;
      background: #faf7f2;
      z-index: 1;
    }}
    th, td {{
      padding: 8px 7px;
      border-bottom: 1px solid #ebe5d9;
      text-align: left;
      vertical-align: top;
      overflow-wrap: anywhere;
      word-break: break-word;
    }}
    th {{
      color: var(--muted);
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }}
    tbody tr {{
      cursor: pointer;
    }}
    tbody tr:hover,
    tbody tr.active {{
      background: rgba(15, 118, 110, 0.08);
    }}
    .chart-wrap.panning {{
      cursor: grabbing;
    }}
    .hint {{
      color: var(--muted);
      font-size: 11px;
      margin-top: 2px;
      line-height: 1.25;
    }}
    @media (max-width: 1080px) {{
      .grid {{
        grid-template-columns: 1fr;
        grid-template-rows: auto;
      }}
      .panel-spectrum {{
        grid-column: auto;
      }}
      .info3d-stack {{
        grid-template-rows: auto;
      }}
      .summary-grid, .mode-grid, .controls {{
        grid-template-columns: 1fr;
      }}
      .chart-wrap {{
        height: 340px;
      }}
      #moleculeViewer {{
        min-height: 420px;
      }}
      .table-wrap {{
        height: 420px;
        max-height: 420px;
      }}
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <div class="toolbar">
      <h1>ORCAVEDA Viewer</h1>
      <p>Interactive spectrum viewer with file summary, 3D molecular view, clean peak table, and mode-specific detail panel.</p>
      <select id="fileSelect"></select>
    </div>

    <div class="grid">
      <section class="panel panel-spectrum">
        <div class="panel-head"><h2>Interactive Spectrum</h2></div>
        <div class="panel-body spectrum-body">
          <div class="controls">
            <div class="control">
              <label for="scaleFactor">Scale Factor <span id="scaleValue" class="value"></span></label>
              <input id="scaleFactor" type="range" min="0.900" max="1.050" step="0.001" value="1.000">
            </div>
            <div class="control">
              <label for="hwhm">Lorentz HWHM (cm-1) <span id="hwhmValue" class="value"></span></label>
              <input id="hwhm" type="range" min="2" max="40" step="0.5" value="12">
            </div>
            <div class="control">
              <label for="yMode">Y-axis Mode</label>
              <select id="yMode">
                <option value="transmittance" selected>Transmittance</option>
                <option value="absorbance">Absorbance</option>
              </select>
            </div>
            <div class="control">
              <label for="nistReference">NIST Reference</label>
              <select id="nistReference"></select>
            </div>
            <div class="control">
              <label>NIST Fit</label>
              <div class="button-row">
                <button id="autoFitScale" type="button" class="ghost-btn">Auto-fit scale</button>
                <button id="resetZoom" type="button" class="ghost-btn">Reset Zoom</button>
              </div>
            </div>
          </div>
          <div id="fitSummary" class="fit-summary">Choose a NIST reference and press Auto-fit scale to estimate the best frequency scaling.</div>
          <div class="checkrow">
            <label><input id="showSticks" type="checkbox" checked> Show sticks</label>
            <label><input id="invertAxis" type="checkbox" checked> Invert x-axis</label>
          </div>
          <div class="chart-wrap">
            <canvas id="chart" width="1200" height="360"></canvas>
            <div id="chartTooltip" class="chart-tooltip"></div>
          </div>
          <div class="hint">Click a stick or row in the peak table to update mode-specific details. Hovering rows also previews the mode.</div>
        </div>
      </section>

      <div class="info3d-stack">
        <section class="panel">
          <div class="panel-head"><h2>File & Molecule Info</h2></div>
          <div class="panel-body info-body">
            <div id="summaryGrid" class="summary-grid"></div>
            <div class="mode-card">
              <h3>Selected Mode</h3>
              <div id="modeDetails" class="mode-grid"></div>
            </div>
          </div>
        </section>

        <section class="panel">
          <div class="panel-head">
            <h2>3D Molecule Viewer</h2>
            <div class="viewer-actions">
              <select id="molStyle">
                <option value="ballstick" selected>Ball &amp; Stick</option>
                <option value="stick">Stick</option>
                <option value="line">Line</option>
                <option value="sphere">Sphere</option>
              </select>
              <button id="reset3d" type="button">Reset View</button>
            </div>
          </div>
          <div class="panel-body viewer-body">
            <div id="moleculeViewer"></div>
          </div>
        </section>
      </div>

      <section class="panel">
        <div class="panel-head"><h2>Peak Table</h2></div>
        <div class="panel-body">
          <div class="table-wrap">
            <table>
              <thead>
                <tr>
                  <th>Mode</th>
                  <th>Scaled Frequency</th>
                  <th>IR Intensity</th>
                  <th>Final Assignment</th>
                </tr>
              </thead>
              <tbody id="peakTable"></tbody>
            </table>
          </div>
        </div>
      </section>
    </div>
  </div>

  <script>
    const payload = {payload_json};
    const fileSelect = document.getElementById("fileSelect");
    const scaleFactor = document.getElementById("scaleFactor");
    const hwhm = document.getElementById("hwhm");
    const yMode = document.getElementById("yMode");
    const nistReference = document.getElementById("nistReference");
    const autoFitScale = document.getElementById("autoFitScale");
    const resetZoom = document.getElementById("resetZoom");
    const fitSummary = document.getElementById("fitSummary");
    const showSticks = document.getElementById("showSticks");
    const invertAxis = document.getElementById("invertAxis");
    const scaleValue = document.getElementById("scaleValue");
    const hwhmValue = document.getElementById("hwhmValue");
    const summaryGrid = document.getElementById("summaryGrid");
    const modeDetails = document.getElementById("modeDetails");
    const peakTable = document.getElementById("peakTable");
    const canvas = document.getElementById("chart");
    const ctx = canvas.getContext("2d");
    const chartTooltip = document.getElementById("chartTooltip");
    const chartWrap = canvas.parentElement;
    const moleculeViewerHost = document.getElementById("moleculeViewer");
    const molStyle = document.getElementById("molStyle");
    const reset3d = document.getElementById("reset3d");

    let currentIndex = 0;
    let currentRender = null;
    let selectedMode = null;
    let moleculeViewer = null;
    let currentView = null;
    let pinnedTooltipMode = null;
    let isPanning = false;
    let panStart = null;

    function resizeChart() {{
      const wrap = canvas.parentElement;
      if (!wrap) return;
      const dpr = window.devicePixelRatio || 1;
      const targetWidth = Math.max(640, Math.round(wrap.clientWidth * dpr));
      const targetHeight = Math.max(420, Math.round(wrap.clientHeight * dpr));
      if (canvas.width !== targetWidth || canvas.height !== targetHeight) {{
        canvas.width = targetWidth;
        canvas.height = targetHeight;
      }}
    }}

    payload.files.forEach((file, idx) => {{
      const opt = document.createElement("option");
      opt.value = String(idx);
      opt.textContent = file.filename;
      fileSelect.appendChild(opt);
    }});

    function getCurrentFile() {{
      return payload.files[currentIndex] || {{ filename: "", summary: {{}}, modes: [], geometry: {{ atoms: [] }} }};
    }}

    function getCurrentReferenceSet() {{
      const file = getCurrentFile();
      const sets = payload.nist_reference_sets || {{}};
      return sets[file.title] || null;
    }}

    function populateReferenceOptions() {{
      const refSet = getCurrentReferenceSet();
      nistReference.innerHTML = "";
      const noneOpt = document.createElement("option");
      noneOpt.value = "";
      noneOpt.textContent = "None";
      nistReference.appendChild(noneOpt);

      if (!refSet || !Array.isArray(refSet.reference_spectra)) {{
        nistReference.disabled = true;
        return;
      }}
      nistReference.disabled = false;
      for (const item of refSet.reference_spectra) {{
        const opt = document.createElement("option");
        opt.value = String(item.index);
        opt.textContent = `${{item.phase_label || item.phase_tag || item.index}}`;
        nistReference.appendChild(opt);
      }}
      const preferred = refSet.preferred_reference || null;
      nistReference.value = preferred ? String(preferred.index || "") : "";
    }}

    function lorentz(x, x0, intensity, gamma) {{
      const g2 = gamma * gamma;
      const dx = x - x0;
      return intensity * (g2 / (dx * dx + g2));
    }}

    function ensureMoleculeViewer() {{
      if (moleculeViewer || typeof $3Dmol === "undefined") return moleculeViewer;
      moleculeViewer = $3Dmol.createViewer(moleculeViewerHost, {{ backgroundColor: "rgba(0,0,0,0)" }});
      return moleculeViewer;
    }}

    function geometryToXyz(file) {{
      const geometry = file.geometry || {{ atoms: [] }};
      const lines = [String(geometry.atoms.length), file.filename || "ORCAVEDA molecule"];
      for (const atom of geometry.atoms) {{
        lines.push(`${{atom.element}} ${{atom.x.toFixed(6)}} ${{atom.y.toFixed(6)}} ${{atom.z.toFixed(6)}}`);
      }}
      return lines.join("\\n");
    }}

    function currentMolStyle() {{
      const value = molStyle.value;
      if (value === "line") return {{ line: {{ linewidth: 2.4, colorscheme: "Jmol" }} }};
      if (value === "sphere") return {{ sphere: {{ scale: 0.34, colorscheme: "Jmol" }} }};
      if (value === "stick") return {{ stick: {{ radius: 0.18, colorscheme: "Jmol" }} }};
      return {{
        stick: {{ radius: 0.18, colorscheme: "Jmol" }},
        sphere: {{ scale: 0.28, colorscheme: "Jmol" }}
      }};
    }}

    function renderMolecule() {{
      const file = getCurrentFile();
      const viewer = ensureMoleculeViewer();
      if (!viewer) {{
        moleculeViewerHost.innerHTML = '<div style="padding:18px;color:#677483;font:14px Segoe UI;">3Dmol.js failed to load. Check internet access or vendor a local 3Dmol-min.js file.</div>';
        return;
      }}
      viewer.clear();
      viewer.addModel(geometryToXyz(file), "xyz");
      viewer.setStyle({{}}, currentMolStyle());
      viewer.zoomTo();
      viewer.resize();
      viewer.render();
    }}

    function defaultRange(file, scale) {{
      if (!file.modes.length) return [0, 4000];
      const scaled = file.modes.map(mode => mode.frequency_cm1 * scale);
      const minScaled = Math.min(...scaled);
      const maxScaled = Math.max(...scaled);
      const left = Math.max(0, Math.floor(minScaled / 25) * 25 - 50);
      const right = Math.min(4200, Math.ceil(maxScaled / 25) * 25 + 50);
      return [left, Math.max(left + 100, right)];
    }}

    function clampView(file, scale, x1, x2) {{
      const [d1, d2] = defaultRange(file, scale);
      const minSpan = 120;
      const maxSpan = Math.max(minSpan, d2 - d1);
      let left = Number.isFinite(x1) ? x1 : d1;
      let right = Number.isFinite(x2) ? x2 : d2;
      if (right - left < minSpan) {{
        const center = (left + right) / 2 || (d1 + d2) / 2;
        left = center - minSpan / 2;
        right = center + minSpan / 2;
      }}
      if (right - left > maxSpan) {{
        left = d1;
        right = d2;
      }}
      if (left < d1) {{
        right += d1 - left;
        left = d1;
      }}
      if (right > d2) {{
        left -= right - d2;
        right = d2;
      }}
      left = Math.max(d1, left);
      right = Math.min(d2, right);
      if (right - left < minSpan) {{
        right = Math.min(d2, left + minSpan);
        left = Math.max(d1, right - minSpan);
      }}
      return [left, right];
    }}

    function transformYValue(value, mode, maxIntensity) {{
      const clipped = Math.max(0, Number(value) || 0);
      const norm = maxIntensity > 0 ? Math.min(1, clipped / maxIntensity) : 0;
      if (mode === "absorbance") return Math.log10(1 + 9 * norm);
      return Math.max(0.02, 1 - 0.92 * norm);
    }}

    function getSelectedReferenceSpectrum() {{
      const refSet = getCurrentReferenceSet();
      if (!refSet || !Array.isArray(refSet.reference_spectra)) return null;
      const selected = String(nistReference.value || "");
      if (!selected) return null;
      return refSet.reference_spectra.find(item => String(item.index) === selected) || null;
    }}

    function convertReferenceY(value, yUnits, axisMode) {{
      const unit = String(yUnits || "").toUpperCase();
      const clipped = Math.max(0, Number(value) || 0);
      if (axisMode === "absorbance") {{
        if (unit.includes("ABSORB")) return clipped;
        let t = clipped;
        if (t > 1.5) t /= 100.0;
        t = Math.max(1e-6, Math.min(1.0, t));
        return -Math.log10(t);
      }}
      if (unit.includes("TRANS")) {{
        let t = clipped;
        if (t > 1.5) t /= 100.0;
        return Math.max(0.0, Math.min(1.0, t));
      }}
      return Math.pow(10, -clipped);
    }}

    function buildReferenceOverlay(referenceSpectrum, x1, x2, axisMode) {{
      if (!referenceSpectrum || !Array.isArray(referenceSpectrum.points)) return null;
      const points = referenceSpectrum.points
        .filter(pt => Number(pt.x) >= x1 && Number(pt.x) <= x2)
        .sort((a, b) => Number(a.x) - Number(b.x));
      if (!points.length) return null;
      const rawY = points.map(pt => convertReferenceY(pt.y, referenceSpectrum.y_units, axisMode));
      const minY = Math.min(...rawY);
      const maxY = Math.max(...rawY);
      const span = Math.max(1e-9, maxY - minY);
      let ys;
      let yMin;
      let yMax;
      if (axisMode === "absorbance") {{
        ys = rawY.map(y => (y - minY) / span);
        yMin = 0.0;
        yMax = 1.0;
      }} else {{
        ys = rawY.map(y => 0.02 + 0.98 * ((y - minY) / span));
        yMin = 0.02;
        yMax = 1.0;
      }}
      return {{
        xs: points.map(pt => Number(pt.x)),
        ys,
        yMin,
        yMax,
        label: referenceSpectrum.phase_label || referenceSpectrum.phase_tag || `Index ${{referenceSpectrum.index}}`,
      }};
    }}

    function inferReferencePeakDirection(referenceSpectrum) {{
      const unit = String(referenceSpectrum?.y_units || "").toUpperCase();
      return unit.includes("ABSORB") ? "max" : "min";
    }}

    function smoothSeries(values, radius = 2) {{
      if (!Array.isArray(values) || values.length < 3) return values.slice();
      const out = [];
      for (let i = 0; i < values.length; i += 1) {{
        let total = 0;
        let count = 0;
        for (let j = Math.max(0, i - radius); j <= Math.min(values.length - 1, i + radius); j += 1) {{
          total += Number(values[j]) || 0;
          count += 1;
        }}
        out.push(total / Math.max(1, count));
      }}
      return out;
    }}

    function pickReferencePeaks(referenceSpectrum, limit = 16, minSpacingCm1 = 18) {{
      if (!referenceSpectrum || !Array.isArray(referenceSpectrum.points) || referenceSpectrum.points.length < 5) return [];
      const points = referenceSpectrum.points
        .map(pt => ({{ x: Number(pt.x), y: Number(pt.y) }}))
        .filter(pt => Number.isFinite(pt.x) && Number.isFinite(pt.y))
        .sort((a, b) => a.x - b.x);
      if (points.length < 5) return [];

      const ys = smoothSeries(points.map(pt => pt.y), 2);
      const direction = inferReferencePeakDirection(referenceSpectrum);
      const candidates = [];
      for (let i = 2; i < points.length - 2; i += 1) {{
        const y0 = ys[i - 1];
        const y1 = ys[i];
        const y2 = ys[i + 1];
        const isExtremum = direction === "min"
          ? (y1 <= y0 && y1 < y2)
          : (y1 >= y0 && y1 > y2);
        if (!isExtremum) continue;

        const leftWindow = ys.slice(Math.max(0, i - 8), i);
        const rightWindow = ys.slice(i + 1, Math.min(ys.length, i + 9));
        if (!leftWindow.length || !rightWindow.length) continue;
        const shoulder = (Math.max(...leftWindow) + Math.max(...rightWindow)) / 2;
        const valley = (Math.min(...leftWindow) + Math.min(...rightWindow)) / 2;
        const prominence = direction === "min" ? shoulder - y1 : y1 - valley;
        if (!(prominence > 0)) continue;
        candidates.push({{
          x: points[i].x,
          y: points[i].y,
          prominence,
        }});
      }}

      candidates.sort((a, b) => b.prominence - a.prominence);
      const selected = [];
      for (const candidate of candidates) {{
        if (selected.some(existing => Math.abs(existing.x - candidate.x) < minSpacingCm1)) continue;
        selected.push(candidate);
        if (selected.length >= limit) break;
      }}
      return selected.sort((a, b) => a.x - b.x);
    }}

    function scoreScaleAgainstReference(file, referenceSpectrum, scale) {{
      const referencePeaks = pickReferencePeaks(referenceSpectrum, 16, 18);
      if (!referencePeaks.length) return null;

      const modes = (file.modes || [])
        .map(mode => ({{
          ...mode,
          scaled: Number(mode.frequency_cm1) * scale,
          weight: Math.max(1e-6, Number(mode.intensity) || 0),
        }}))
        .sort((a, b) => b.weight - a.weight);
      if (!modes.length) return null;

      const rankedPeaks = referencePeaks
        .map(peak => ({{
          ...peak,
          weight: Math.max(1e-6, peak.prominence),
        }}))
        .sort((a, b) => b.weight - a.weight);

      const tolerance = 35.0;
      const unmatchedPenalty = 28.0;
      const usedModes = new Set();
      const matches = [];
      let weightedDelta = 0;
      let totalWeight = 0;
      let unmatchedCount = 0;

      for (const peak of rankedPeaks) {{
        let best = null;
        for (const mode of modes) {{
          if (usedModes.has(mode.mode)) continue;
          const delta = Math.abs(mode.scaled - peak.x);
          if (delta > tolerance) continue;
          const score = delta - 0.015 * Math.log10(1 + mode.weight);
          if (!best || score < best.score) {{
            best = {{ mode, delta, score }};
          }}
        }}

        const w = peak.weight;
        totalWeight += w;
        if (!best) {{
          unmatchedCount += 1;
          weightedDelta += w * unmatchedPenalty;
          continue;
        }}

        usedModes.add(best.mode.mode);
        weightedDelta += w * best.delta;
        matches.push({{
          peak_x: peak.x,
          mode: best.mode.mode,
          scaled: best.mode.scaled,
          delta: best.delta,
          weight: w,
          assignment: best.mode.assignment || "unassigned",
        }});
      }}

      const score = weightedDelta / Math.max(1e-6, totalWeight) + unmatchedCount * 2.5;
      const meanAbsDelta = matches.length
        ? matches.reduce((acc, row) => acc + row.delta, 0) / matches.length
        : Infinity;
      return {{
        score,
        scale,
        matchedCount: matches.length,
        totalPeaks: rankedPeaks.length,
        unmatchedCount,
        meanAbsDelta,
        matches,
      }};
    }}

    function autoFitScaleAgainstReference() {{
      const file = getCurrentFile();
      const referenceSpectrum = getSelectedReferenceSpectrum();
      if (!file || !referenceSpectrum) {{
        fitSummary.textContent = "Select a NIST reference first, then press Auto-fit scale.";
        return;
      }}

      let best = null;
      for (let scale = 0.9; scale <= 1.0500001; scale += 0.0005) {{
        const rounded = Number(scale.toFixed(4));
        const result = scoreScaleAgainstReference(file, referenceSpectrum, rounded);
        if (!result) continue;
        if (!best || result.score < best.score) best = result;
      }}

      if (!best || !Number.isFinite(best.scale)) {{
        fitSummary.textContent = "Auto-fit could not estimate a stable scale for the selected NIST reference.";
        return;
      }}

      scaleFactor.value = best.scale.toFixed(3);
      fitSummary.textContent = `Best scale ${{best.scale.toFixed(3)}} | mean Δ ${{best.meanAbsDelta.toFixed(1)}} cm-1 | matched ${{best.matchedCount}}/${{best.totalPeaks}}`;
      pinnedTooltipMode = null;
      chartTooltip.classList.remove("visible");
      currentView = null;
      drawSpectrum(false);
    }}

    function buildSpectrum(file, scale, gamma, x1, x2, axisMode) {{
      const n = 1500;
      const xs = [];
      const rawYs = [];
      let rawMax = 0;
      for (let i = 0; i < n; i += 1) {{
        const x = x1 + (x2 - x1) * i / (n - 1);
        let y = 0;
        for (const mode of file.modes) {{
          y += lorentz(x, mode.frequency_cm1 * scale, Math.max(0, mode.intensity), gamma);
        }}
        xs.push(x);
        rawYs.push(y);
        rawMax = Math.max(rawMax, y);
      }}
      const ys = rawYs.map(y => transformYValue(y, axisMode, rawMax));
      const yMin = axisMode === "transmittance" ? Math.min(...ys) : 0;
      const yMax = axisMode === "transmittance" ? 1.0 : Math.max(...ys, 1e-6);
      return {{ xs, ys, yMin, yMax, rawMax }};
    }}

    function tx(x, x1, x2, left, width, inverted) {{
      const t = (x - x1) / (x2 - x1);
      return inverted ? left + width * (1 - t) : left + width * t;
    }}

    function ty(y, yMin, yMax, top, height) {{
      const span = Math.max(1e-9, yMax - yMin);
      return top + height * (1 - (y - yMin) / span);
    }}

    function formatMaybe(value, digits = 1) {{
      if (value == null || Number.isNaN(Number(value))) return "n/a";
      return Number(value).toFixed(digits);
    }}

    function updateSummary(file) {{
      const s = file.summary || {{}};
      const items = [
        ["Filename", file.filename || "n/a"],
        ["Formula", s.formula || "n/a"],
        ["System", s.system_type || "n/a"],
        ["Atoms", s.natoms ?? "n/a"],
        ["Total Modes", s.total_modes ?? "n/a"],
        ["Positive Modes", s.positive_mode_count ?? "n/a"],
        ["Fragments", s.fragment_count ?? "n/a"],
        ["Range", s.frequency_min_cm1 != null && s.frequency_max_cm1 != null ? `${{Number(s.frequency_min_cm1).toFixed(1)}} – ${{Number(s.frequency_max_cm1).toFixed(1)}} cm-1` : "n/a"]
      ];
      summaryGrid.innerHTML = items.map(([k, v]) => `<div class="kv"><strong>${{k}}</strong>${{v}}</div>`).join("");
    }}

    function updateModeDetails(file, scale) {{
      const mode = file.modes.find(row => row.mode === selectedMode) || file.modes[0];
      if (!mode) {{
        modeDetails.innerHTML = '<div class="wide">No positive-frequency modes available.</div>';
        return;
      }}
      modeDetails.innerHTML = `
        <div class="kv"><strong>Mode</strong>${{mode.mode}}</div>
        <div class="kv"><strong>Scaled Frequency</strong>${{(mode.frequency_cm1 * scale).toFixed(2)}} cm-1</div>
        <div class="kv"><strong>Original Frequency</strong>${{mode.frequency_cm1.toFixed(2)}} cm-1</div>
        <div class="kv"><strong>IR Intensity</strong>${{Number(mode.intensity).toFixed(4)}}</div>
        <div class="kv wide"><strong>Final Assignment</strong>${{mode.assignment || "unassigned"}}</div>
        <div class="kv wide"><strong>Warnings</strong>${{mode.warnings || "none"}}</div>
        <div class="kv wide"><strong>Supporting Coordinates</strong>${{mode.top_internal_coordinates || "n/a"}}</div>
      `;
    }}

    function updatePeakTable(file, scale, x1, x2) {{
      const rows = file.modes
        .map(mode => ({{ ...mode, scaled: mode.frequency_cm1 * scale }}))
        .filter(mode => mode.scaled >= x1 && mode.scaled <= x2)
        .sort((a, b) => a.scaled - b.scaled);
      peakTable.innerHTML = "";
      for (const mode of rows) {{
        const tr = document.createElement("tr");
        if (mode.mode === selectedMode) tr.classList.add("active");
        tr.innerHTML = `
          <td>${{mode.mode}}</td>
          <td>${{mode.scaled.toFixed(1)}}</td>
          <td>${{Number(mode.intensity).toFixed(3)}}</td>
          <td>${{mode.assignment || "unassigned"}}</td>
        `;
        tr.addEventListener("mouseenter", () => {{
          selectedMode = mode.mode;
          updateModeDetails(file, scale);
          highlightSelectedRow();
        }});
        tr.addEventListener("click", () => {{
          selectedMode = mode.mode;
          updateModeDetails(file, scale);
          highlightSelectedRow();
          drawSpectrum(false);
        }});
        peakTable.appendChild(tr);
      }}
    }}

    function highlightSelectedRow() {{
      for (const row of peakTable.querySelectorAll("tr")) row.classList.remove("active");
      const target = Array.from(peakTable.querySelectorAll("tr")).find(row => String(row.children[0]?.textContent || "") === String(selectedMode));
      if (target) target.classList.add("active");
    }}

    function modeScreenX(mode, render) {{
      const center = mode.frequency_cm1 * render.scale;
      return tx(center, render.x1, render.x2, render.margin.left, render.plotW, invertAxis.checked);
    }}

    function nearestModeAtCanvasX(px, threshold = 16) {{
      if (!currentRender) return null;
      let best = null;
      let bestDx = Infinity;
      for (const mode of currentRender.file.modes) {{
        const center = mode.frequency_cm1 * currentRender.scale;
        if (center < currentRender.x1 || center > currentRender.x2) continue;
        const sx = modeScreenX(mode, currentRender);
        const dx = Math.abs(px - sx);
        if (dx < bestDx) {{
          bestDx = dx;
          best = mode;
        }}
      }}
      return bestDx <= threshold ? best : null;
    }}

    function showTooltip(mode, clientX, clientY) {{
      if (!mode || !currentRender) {{
        chartTooltip.classList.remove("visible");
        return;
      }}
      chartTooltip.innerHTML = `
        <strong>Mode ${{mode.mode}}</strong>
        <div>${{(mode.frequency_cm1 * currentRender.scale).toFixed(1)}} cm-1</div>
        <div>IR: ${{Number(mode.intensity).toFixed(3)}}</div>
        <div>${{mode.assignment || "unassigned"}}</div>
      `;
      const rect = canvas.parentElement.getBoundingClientRect();
      const tooltipWidth = 240;
      const offsetX = 14;
      const offsetY = 12;
      let left = clientX - rect.left + offsetX;
      let top = clientY - rect.top + offsetY;
      if (left + tooltipWidth > rect.width - 8) left = rect.width - tooltipWidth - 8;
      if (left < 8) left = 8;
      const maxTop = rect.height - 92;
      if (top > maxTop) top = maxTop;
      if (top < 8) top = 8;
      chartTooltip.style.left = `${{left}}px`;
      chartTooltip.style.top = `${{top}}px`;
      chartTooltip.style.width = `${{tooltipWidth}}px`;
      chartTooltip.classList.add("visible");
    }}

    function showTooltipNearMode(mode) {{
      if (!mode || !currentRender) {{
        chartTooltip.classList.remove("visible");
        return;
      }}
      const px = modeScreenX(mode, currentRender);
      const py = currentRender.margin.top + Math.min(32, currentRender.plotH * 0.12);
      const rect = chartWrap.getBoundingClientRect();
      showTooltip(mode, rect.left + (px / canvas.width) * rect.width, rect.top + (py / canvas.height) * rect.height);
    }}

    function drawSpectrum(allowDefaultSelection = true) {{
      resizeChart();
      const file = getCurrentFile();
      const scale = Number(scaleFactor.value);
      const gamma = Number(hwhm.value);
      const axisMode = yMode.value || "intensity";
      const [x1, x2] = clampView(file, scale, currentView?.x1, currentView?.x2);
      currentView = {{ x1, x2 }};
      scaleValue.textContent = scale.toFixed(3);
      hwhmValue.textContent = gamma.toFixed(1);

      const width = canvas.width;
      const height = canvas.height;
      const margin = {{ left: 76, right: 24, top: 20, bottom: 64 }};
      const plotW = width - margin.left - margin.right;
      const plotH = height - margin.top - margin.bottom;

      const render = buildSpectrum(file, scale, gamma, x1, x2, axisMode);
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const referenceOverlay = buildReferenceOverlay(referenceSpectrum, x1, x2, axisMode);
      currentRender = {{ ...render, x1, x2, scale, gamma, axisMode, margin, plotW, plotH, file, referenceSpectrum, referenceOverlay }};
      ctx.clearRect(0, 0, width, height);
      ctx.fillStyle = "#fffdf8";
      ctx.fillRect(0, 0, width, height);

      ctx.strokeStyle = "#ddd6c8";
      ctx.lineWidth = 1;
      for (let i = 0; i <= 6; i += 1) {{
        const y = margin.top + plotH * i / 6;
        ctx.beginPath();
        ctx.moveTo(margin.left, y);
        ctx.lineTo(margin.left + plotW, y);
        ctx.stroke();
      }}

      const tickCount = 8;
      ctx.fillStyle = "#677483";
      ctx.font = "12px Segoe UI";
      ctx.textAlign = "center";
      for (let i = 0; i <= tickCount; i += 1) {{
        const value = x1 + (x2 - x1) * i / tickCount;
        const px = margin.left + plotW * i / tickCount;
        const shown = invertAxis.checked ? x2 - (x2 - x1) * i / tickCount : value;
        ctx.beginPath();
        ctx.strokeStyle = "#e8e1d4";
        ctx.moveTo(px, margin.top);
        ctx.lineTo(px, margin.top + plotH);
        ctx.stroke();
        ctx.fillText(Math.round(shown).toString(), px, margin.top + plotH + 22);
      }}

      ctx.strokeStyle = "#0b5670";
      ctx.lineWidth = 2.2;
      ctx.beginPath();
      render.xs.forEach((x, idx) => {{
        const px = tx(x, x1, x2, margin.left, plotW, invertAxis.checked);
        const py = ty(render.ys[idx], render.yMin, render.yMax, margin.top, plotH);
        if (idx === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }});
      ctx.stroke();

      if (referenceOverlay) {{
        ctx.strokeStyle = "rgba(154, 52, 18, 0.88)";
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        referenceOverlay.xs.forEach((x, idx) => {{
          const px = tx(x, x1, x2, margin.left, plotW, invertAxis.checked);
          const py = ty(referenceOverlay.ys[idx], referenceOverlay.yMin, referenceOverlay.yMax, margin.top, plotH);
          if (idx === 0) ctx.moveTo(px, py);
          else ctx.lineTo(px, py);
        }});
        ctx.stroke();

        ctx.fillStyle = "rgba(154, 52, 18, 0.88)";
        ctx.font = "12px Segoe UI";
        ctx.textAlign = "left";
        ctx.fillText(`NIST: ${{referenceOverlay.label}}`, margin.left + 6, margin.top + 16);
      }}

      if (showSticks.checked) {{
        for (const mode of file.modes) {{
          const center = mode.frequency_cm1 * scale;
          if (center < x1 || center > x2) continue;
          const px = tx(center, x1, x2, margin.left, plotW, invertAxis.checked);
          const stickValue = transformYValue(Math.max(0, mode.intensity), axisMode, Math.max(render.rawMax, 1e-6));
          const baseline = axisMode === "transmittance" ? render.yMax : render.yMin;
          const stickTop = ty(stickValue, render.yMin, render.yMax, margin.top, plotH);
          const stickBase = ty(baseline, render.yMin, render.yMax, margin.top, plotH);
          ctx.strokeStyle = mode.mode === selectedMode ? "rgba(154,52,18,0.85)" : "rgba(15,118,110,0.34)";
          ctx.lineWidth = mode.mode === selectedMode ? 2.4 : 1.2;
          ctx.beginPath();
          ctx.moveTo(px, stickBase);
          ctx.lineTo(px, stickTop);
          ctx.stroke();
        }}
      }}

      ctx.strokeStyle = "#1f2933";
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.moveTo(margin.left, margin.top);
      ctx.lineTo(margin.left, margin.top + plotH);
      ctx.lineTo(margin.left + plotW, margin.top + plotH);
      ctx.stroke();

      ctx.fillStyle = "#677483";
      ctx.textAlign = "center";
      ctx.fillText("Wavenumber (cm-1)", margin.left + plotW / 2, height - 16);
      ctx.save();
      ctx.translate(24, margin.top + plotH / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(axisMode === "absorbance" ? "Absorbance (a.u.)" : "Transmittance", 0, 0);
      ctx.restore();

      if (allowDefaultSelection && (selectedMode == null) && file.modes.length) {{
        selectedMode = file.modes[Math.min(file.modes.length - 1, Math.floor(file.modes.length * 0.7))].mode;
      }}
      updateSummary(file);
      updateModeDetails(file, scale);
      updatePeakTable(file, scale, x1, x2);
      highlightSelectedRow();
    }}

    canvas.addEventListener("click", (event) => {{
      if (!currentRender) return;
      const rect = canvas.getBoundingClientRect();
      const px = (event.clientX - rect.left) * (canvas.width / rect.width);
      const best = nearestModeAtCanvasX(px, 22);
      if (best) {{
        selectedMode = best.mode;
        pinnedTooltipMode = best.mode;
        updateModeDetails(currentRender.file, currentRender.scale);
        highlightSelectedRow();
        drawSpectrum(false);
        showTooltipNearMode(best);
      }} else {{
        pinnedTooltipMode = null;
        chartTooltip.classList.remove("visible");
      }}
    }});

    canvas.addEventListener("wheel", (event) => {{
      if (!currentRender) return;
      event.preventDefault();
      const rect = canvas.getBoundingClientRect();
      const px = (event.clientX - rect.left) * (canvas.width / rect.width);
      const plotLeft = currentRender.margin.left;
      const plotRight = currentRender.margin.left + currentRender.plotW;
      const clampedPx = Math.max(plotLeft, Math.min(plotRight, px));
      const tPlot = (clampedPx - plotLeft) / currentRender.plotW;
      const t = invertAxis.checked ? 1 - tPlot : tPlot;
      const center = currentRender.x1 + (currentRender.x2 - currentRender.x1) * t;
      const span = currentRender.x2 - currentRender.x1;
      const zoomFactor = event.deltaY < 0 ? 0.86 : 1.16;
      const nextSpan = span * zoomFactor;
      const left = center - nextSpan * t;
      const right = center + nextSpan * (1 - t);
      const [x1, x2] = clampView(currentRender.file, currentRender.scale, left, right);
      currentView = {{ x1, x2 }};
      drawSpectrum(false);
      if (pinnedTooltipMode != null) {{
        const mode = currentRender.file.modes.find(row => row.mode === pinnedTooltipMode);
        showTooltipNearMode(mode);
      }}
    }}, {{ passive: false }});

    canvas.addEventListener("mousedown", (event) => {{
      if (!currentRender || event.button !== 0) return;
      const rect = canvas.getBoundingClientRect();
      const px = (event.clientX - rect.left) * (canvas.width / rect.width);
      const plotLeft = currentRender.margin.left;
      const plotRight = currentRender.margin.left + currentRender.plotW;
      if (px < plotLeft || px > plotRight) return;
      isPanning = true;
      chartWrap.classList.add("panning");
      panStart = {{
        clientX: event.clientX,
        x1: currentRender.x1,
        x2: currentRender.x2,
      }};
    }});

    window.addEventListener("mousemove", (event) => {{
      if (!isPanning || !currentRender || !panStart) return;
      const rect = canvas.getBoundingClientRect();
      const dxPx = (event.clientX - panStart.clientX) * (canvas.width / rect.width);
      const span = panStart.x2 - panStart.x1;
      const delta = (dxPx / currentRender.plotW) * span * (invertAxis.checked ? 1 : -1);
      const [x1, x2] = clampView(currentRender.file, currentRender.scale, panStart.x1 + delta, panStart.x2 + delta);
      currentView = {{ x1, x2 }};
      drawSpectrum(false);
      if (pinnedTooltipMode != null) {{
        const mode = currentRender.file.modes.find(row => row.mode === pinnedTooltipMode);
        showTooltipNearMode(mode);
      }}
    }});

    window.addEventListener("mouseup", () => {{
      isPanning = false;
      panStart = null;
      chartWrap.classList.remove("panning");
    }});

    reset3d.addEventListener("click", () => {{
      const viewer = ensureMoleculeViewer();
      if (!viewer) return;
      viewer.zoomTo();
      viewer.render();
    }});

    resetZoom.addEventListener("click", () => {{
      currentView = null;
      pinnedTooltipMode = null;
      chartTooltip.classList.remove("visible");
      drawSpectrum(false);
    }});
    autoFitScale.addEventListener("click", autoFitScaleAgainstReference);

    molStyle.addEventListener("input", renderMolecule);
    fileSelect.addEventListener("input", () => {{
      currentIndex = Number(fileSelect.value || 0);
      selectedMode = null;
      currentView = null;
      pinnedTooltipMode = null;
      chartTooltip.classList.remove("visible");
      populateReferenceOptions();
      drawSpectrum(true);
      renderMolecule();
    }});
    [scaleFactor, hwhm, yMode, showSticks, invertAxis, nistReference].forEach(el => {{
      el.addEventListener("input", () => {{
        if (el === yMode || el === nistReference) {{
          pinnedTooltipMode = null;
          chartTooltip.classList.remove("visible");
          if (el === nistReference) {{
            fitSummary.textContent = "Reference changed. Press Auto-fit scale to estimate the best frequency scaling for this spectrum.";
          }}
        }}
        drawSpectrum(false);
      }});
    }});
    window.addEventListener("resize", () => {{
      resizeChart();
      drawSpectrum(false);
      if (moleculeViewer) {{
        moleculeViewer.resize();
        moleculeViewer.render();
      }}
    }});

    scaleFactor.value = String(payload.default_scale_factor || 1.0);
    hwhm.value = String(payload.default_lorentz_hwhm || 12.0);
    fileSelect.value = "0";
    populateReferenceOptions();
    drawSpectrum(true);
    renderMolecule();
  </script>
</body>
</html>
"""
    html_path.write_text(html, encoding="utf-8")
    return html_path
