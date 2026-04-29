from __future__ import annotations

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
