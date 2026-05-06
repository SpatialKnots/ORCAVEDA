from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Dict, Sequence

import pandas as pd


def classify_reference_suitability(
    *,
    y_units: str = "",
    phase_tag: str = "",
    description: str = "",
) -> tuple[bool, str]:
    y_text = str(y_units or "").strip().lower()
    phase_text = str(phase_tag or "").strip().lower()
    desc_text = str(description or "").strip().lower()
    combined = " ".join(part for part in (y_text, phase_text, desc_text) if part)

    if "absorption index" in combined:
        return False, "absorption_index_reference"
    if not y_text:
        return False, "missing_y_units"
    if "transmittance" in y_text or "absorbance" in y_text:
        return True, "ir_curve"
    return False, f"unsupported_y_units:{y_text}"


def build_nist_spectrum_url(
    *,
    nist_id: str = "",
    index: str = "",
    fallback_page_url: str = "",
    jcamp_url: str = "",
) -> str:
    nist_id = str(nist_id or "").strip()
    index = str(index or "").strip()
    if nist_id and index:
        return f"https://webbook.nist.gov/cgi/inchi?ID={nist_id}&Type=IR-SPEC&Index={index}#IR-SPEC"
    if fallback_page_url:
        return str(fallback_page_url)
    return str(jcamp_url or "")


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
    manifest_page_url = str(manifest.get("nist_page_url", ""))

    spectra = []
    for item in manifest.get("reference_spectra", []):
        csv_path = manifest_dir / Path(item["csv"]).name if not Path(item["csv"]).is_absolute() else Path(item["csv"])
        meta_path = manifest_dir / Path(item["meta_json"]).name if not Path(item["meta_json"]).is_absolute() else Path(item["meta_json"])
        spectrum_df = pd.read_csv(csv_path, encoding="utf-8")
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        y_units = str(meta.get("jcamp_metadata", {}).get("YUNITS", ""))
        phase_tag = str(item.get("phase_tag", ""))
        description = str(item.get("description", ""))
        suitable_for_matching, suitability_reason = classify_reference_suitability(
            y_units=y_units,
            phase_tag=phase_tag,
            description=description,
        )
        nist_spectrum_url = build_nist_spectrum_url(
            nist_id=str(item.get("nist_id", "")),
            index=str(item.get("index", "")),
            fallback_page_url=manifest_page_url,
            jcamp_url=str(item.get("jcamp_url", "")),
        )
        spectra.append(
            {
                "index": str(item.get("index", "")),
                "phase_tag": phase_tag,
                "phase_label": str(item.get("phase_label", "")),
                "selection_priority": int(item.get("selection_priority", 0)),
                "description": description,
                "y_units": y_units,
                "state": str(meta.get("jcamp_metadata", {}).get("STATE", "")),
                "nist_spectrum_url": nist_spectrum_url,
                "suitable_for_matching": bool(suitable_for_matching),
                "suitability_reason": str(suitability_reason),
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
        "nist_page_url": manifest_page_url,
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

    ref_set = load_nist_reference_set(manifest_path)
    target_file = next(
        (item for item in updated.get("files", []) if str(item.get("title", "")) == str(target_title)),
        None,
    )
    if target_file is not None:
        from nist_ir.compare import (
            assignment_modes_to_dataframe,
            build_scale_engine_payload,
            reference_points_to_peaks,
        )

        assignment_audit = assignment_modes_to_dataframe(target_file.get("modes", []), scale_factor=1.0)
        if not assignment_audit.empty:
            for spectrum in ref_set.get("reference_spectra", []):
                reference_context = {
                    "phase_tag": spectrum.get("phase_tag", ""),
                    "phase_label": spectrum.get("phase_label", ""),
                    "state": spectrum.get("state", ""),
                    "description": spectrum.get("description", ""),
                    "y_units": spectrum.get("y_units", ""),
                }
                if not bool(spectrum.get("suitable_for_matching", True)):
                    spectrum["scale_engine_payload"] = None
                    spectrum["matching_status"] = "skipped_unsuitable_reference"
                    spectrum["matching_status_reason"] = str(spectrum.get("suitability_reason", "unsuitable_reference"))
                    continue
                reference_peaks = reference_points_to_peaks(
                    spectrum.get("points", []),
                    top_n=16,
                    min_separation_cm1=18.0,
                    reference_context=reference_context,
                )
                spectrum["scale_engine_payload"] = build_scale_engine_payload(
                    reference_peaks,
                    assignment_audit,
                    reference_context=reference_context,
                )

    refs[str(target_title)] = ref_set
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


def _ped_rank_column(ped_df: pd.DataFrame) -> str:
    if "wilson_rank" in ped_df.columns:
        return "wilson_rank"
    if "ped_rank" in ped_df.columns:
        return "ped_rank"
    return ""


def _ped_method_column(ped_df: pd.DataFrame) -> str:
    if "wilson_ped_method" in ped_df.columns:
        return "wilson_ped_method"
    if "ped_method" in ped_df.columns:
        return "ped_method"
    return ""


def _select_ped_diagnostic_layer(
    *,
    wilson_ped_audit: pd.DataFrame | None = None,
    ped_v2_force_audit: pd.DataFrame | None = None,
    ped_audit: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, str]:
    if isinstance(wilson_ped_audit, pd.DataFrame) and not wilson_ped_audit.empty:
        return wilson_ped_audit, "Wilson GF-style PED audit"
    if isinstance(ped_v2_force_audit, pd.DataFrame) and not ped_v2_force_audit.empty:
        return ped_v2_force_audit, "PED v2 force-aware diagnostic"
    if isinstance(ped_audit, pd.DataFrame) and not ped_audit.empty:
        return ped_audit, "PED v1 B-matrix projection diagnostic"
    return pd.DataFrame(), ""


def _assignment_semantic_classes(text: object) -> set[str]:
    value = str(text or "").lower()
    replacements = {
        "c=o": "carbonyl co stretch",
        "c=s": "thiocarbonyl cs stretch",
        "c#c": "alkyne cc stretch",
        "c=c": "alkene cc stretch",
        "c#n": "cn stretch",
        "c=n": "imine cn stretch",
        "s-h": "sh",
        "o-h": "oh",
        "n-h": "nh",
        "c-h": "ch",
        "c-o": "co",
        "c-n": "cn",
        "c-c": "cc",
        "h-o-h": "hoh",
        "nh2": "nh2",
        "ch2": "ch2",
    }
    for needle, replacement in replacements.items():
        value = value.replace(needle, replacement)
    classes: set[str] = set()
    checks = {
        "carbonyl": ["carbonyl", "c=o"],
        "thiocarbonyl": ["thiocarbonyl", "c=s"],
        "nitro": ["nitro", "no2", "n-o"],
        "thiol": ["thiol", "sh", "s-h"],
        "thioamide": ["thioamide"],
        "isocyanate": ["isocyanate", "nco", "n=c=o"],
        "alkyne": ["alkyne", "c#c"],
        "alkene": ["alkene", "vinylic", "c=c"],
        "sulfone": ["sulfone", "o=s=o"],
        "thioether": ["thioether"],
        "imine": ["imine", "c=n"],
        "oxime": ["oxime"],
        "acyl_chloride": ["acyl chloride", "c-cl"],
        "lactone": ["lactone"],
        "carbonate": ["carbonate"],
        "anhydride": ["anhydride"],
        "epoxide": ["epoxide"],
        "acetal": ["acetal"],
        "oh": ["oh", "hydroxyl", "alcohol", "phenolic", "carboxylic"],
        "sh": ["sh", "thiol"],
        "nh": ["nh", "amine"],
        "nh2": ["nh2"],
        "ch": ["ch", "methyl", "methylene", "methine", "aromatic"],
        "ch2": ["ch2", "methylene"],
        "ring": ["ring", "aromatic", "phenyl", "pyridine"],
        "co": ["co", "carbonyl", "carboxylic", "alcohol"],
        "cn": ["cn", "amine", "n-c", "c-n"],
        "cc": ["cc", "ring", "aromatic"],
        "stretch": ["stretch"],
        "bend": ["bend", "scissor", "deformation", "umbrella"],
        "torsion": ["torsion", "out-of-plane", "oop"],
        "mixed": ["mixed"],
    }
    for cls, terms in checks.items():
        if any(term in value for term in terms):
            classes.add(cls)
    return classes


def classify_ped_stage3d_agreement(
    stage3d_assignment: object,
    ped_assignment: object,
    ped_top_percent: object,
    ped_source: str,
) -> tuple[str, str]:
    stage3d_text = str(stage3d_assignment or "").strip()
    ped_text = str(ped_assignment or "").strip()
    if not ped_text:
        return "not_available", "ped_not_available"
    if not stage3d_text or stage3d_text == "unassigned":
        return "adds_context", "stage3d_missing_or_unassigned"

    try:
        top_percent = float(ped_top_percent)
    except (TypeError, ValueError):
        top_percent = 0.0

    stage3d_classes = _assignment_semantic_classes(stage3d_text)
    ped_classes = _assignment_semantic_classes(ped_text)
    overlap = stage3d_classes & ped_classes
    motion_classes = {"stretch", "bend", "torsion"}
    stage3d_motion = stage3d_classes & motion_classes
    ped_motion = ped_classes & motion_classes
    if stage3d_motion and ped_motion and not (stage3d_motion & ped_motion):
        overlap = set()
    warnings = []
    if top_percent < 25.0:
        warnings.append("diffuse_ped_contributions")
    if ped_source:
        warnings.append("ped_diagnostic_basis")

    if overlap:
        if top_percent >= 50.0 and ped_classes.issubset(stage3d_classes | {"mixed"}):
            return "confirms", "; ".join(warnings)
        return "adds_context", "; ".join(warnings)
    if top_percent < 25.0:
        return "diffuse", "; ".join(warnings)
    warnings.append("ped_stage3d_semantic_disagreement")
    return "disagrees", "; ".join(warnings)


def decide_ped_driven_final_assignment(
    stage3d_assignment: object,
    ped_assignment: object,
    ped_source: object,
    ped_agreement_status: object,
    ped_policy_warning: object = "",
    ped_top_percent: object = 0.0,
) -> tuple[str, str, str, str]:
    stage3d_text = str(stage3d_assignment or "").strip()
    ped_text = str(ped_assignment or "").strip()
    source = str(ped_source or "").strip()
    status = str(ped_agreement_status or "").strip()
    warning = str(ped_policy_warning or "").strip()
    try:
        top_percent = float(ped_top_percent)
    except (TypeError, ValueError):
        top_percent = 0.0

    if ped_text and status in {"confirms", "adds_context"}:
        policy = "ped_confirms_stage3d" if status == "confirms" else "ped_adds_context"
        return ped_text, source or "PED diagnostic", policy, warning
    if ped_text and (not stage3d_text or stage3d_text == "unassigned"):
        return ped_text, source or "PED diagnostic", "ped_used_when_stage3d_unassigned", warning
    if (
        ped_text
        and status == "disagrees"
        and top_percent >= 60.0
        and "diffuse_ped_contributions" not in warning
        and "torsion" in _assignment_semantic_classes(stage3d_text)
        and "bend" in _assignment_semantic_classes(ped_text)
    ):
        calibrated_warning = "stage3d_torsion_reclassified_by_high_confidence_ped_bend"
        final_warning = f"{warning}; {calibrated_warning}" if warning else calibrated_warning
        return ped_text, source or "PED diagnostic", "ped_reclassifies_stage3d_torsion", final_warning

    final = stage3d_text or ped_text or "unassigned"
    if status == "disagrees":
        policy = "stage3d_fallback_due_to_ped_disagreement"
    elif status == "diffuse":
        policy = "stage3d_fallback_due_to_diffuse_ped"
    elif status == "not_available":
        policy = "stage3d_fallback_ped_not_available"
    else:
        policy = "stage3d_fallback"
    final_warning = warning
    if ped_text and status in {"disagrees", "diffuse"}:
        suffix = "final_label_kept_stage3d"
        final_warning = f"{final_warning}; {suffix}" if final_warning else suffix
    return final, "Stage 3D assignment audit", policy, final_warning


def _build_ped_viewer_mode_summary(ped_df: pd.DataFrame | None, *, top_n: int = 6) -> Dict[tuple[str, int], Dict[str, object]]:
    if not isinstance(ped_df, pd.DataFrame) or ped_df.empty:
        return {}
    rank_col = _ped_rank_column(ped_df)
    if not rank_col or "Filename" not in ped_df.columns or "mode" not in ped_df.columns:
        return {}

    df = ped_df.copy()
    df["mode"] = pd.to_numeric(df["mode"], errors="coerce")
    df[rank_col] = pd.to_numeric(df[rank_col], errors="coerce")
    if "contribution_percent" in df.columns:
        df["contribution_percent"] = pd.to_numeric(df["contribution_percent"], errors="coerce").fillna(0.0)
    else:
        df["contribution_percent"] = 0.0
    df = df[df["mode"].notna()].copy()
    df = df[(df[rank_col] >= 1) & (df[rank_col] <= int(top_n))].copy()
    if df.empty:
        return {}

    method_col = _ped_method_column(df)
    summaries: Dict[tuple[str, int], Dict[str, object]] = {}
    for (filename, mode), group in df.sort_values(["Filename", "mode", rank_col]).groupby(["Filename", "mode"], dropna=False):
        terms = []
        family_totals: Dict[str, float] = {}
        for _, row in group.iterrows():
            family = str(row.get("coordinate_family", "") or "").strip() or "internal coordinate"
            coord = str(row.get("internal_coordinate", "") or "").strip()
            pct = float(row.get("contribution_percent", 0.0) or 0.0)
            terms.append(f"{family} {pct:.1f}% [{coord}]")
            family_totals[family] = family_totals.get(family, 0.0) + pct

        ordered = sorted(family_totals.items(), key=lambda item: item[1], reverse=True)
        if not ordered:
            assignment = ""
            top_family = ""
            top_percent = 0.0
        else:
            top_family, top_percent = ordered[0]
            assignment = top_family
            mixed = [(family, pct) for family, pct in ordered[1:] if pct >= 12.0 and family != top_family]
            if mixed:
                assignment = f"{top_family} mixed with {mixed[0][0]}"
                if len(mixed) > 1 and mixed[1][1] >= 10.0:
                    assignment = f"{assignment} and {mixed[1][0]}"

        method = str(group.iloc[0].get(method_col, "") or "") if method_col else ""
        summaries[(str(filename), int(mode))] = {
            "ped_assignment": assignment,
            "ped_top_family": top_family,
            "ped_top_percent": round(float(top_percent), 3),
            "ped_top_contributors": "; ".join(terms),
            "ped_method": method,
        }
    return summaries


def build_ped_stage3d_agreement_table(
    assignment_audit: pd.DataFrame | None,
    *,
    wilson_ped_audit: pd.DataFrame | None = None,
    ped_v2_force_audit: pd.DataFrame | None = None,
    ped_audit: pd.DataFrame | None = None,
) -> pd.DataFrame:
    audit_df = assignment_audit if isinstance(assignment_audit, pd.DataFrame) else pd.DataFrame()
    columns = [
        "Source",
        "Filename",
        "mode",
        "frequency_cm-1",
        "stage3d_assignment",
        "ped_assignment",
        "ped_source",
        "ped_top_family",
        "ped_top_percent",
        "ped_agreement_status",
        "ped_policy_warning",
        "ped_top_contributors",
    ]
    if audit_df.empty or "Filename" not in audit_df.columns or "mode" not in audit_df.columns:
        return pd.DataFrame(columns=columns)

    ped_df, ped_source_label = _select_ped_diagnostic_layer(
        wilson_ped_audit=wilson_ped_audit,
        ped_v2_force_audit=ped_v2_force_audit,
        ped_audit=ped_audit,
    )
    ped_by_mode = _build_ped_viewer_mode_summary(ped_df)
    rows = []
    for _, audit_row in audit_df.iterrows():
        try:
            mode = int(audit_row.get("mode"))
        except (TypeError, ValueError):
            continue
        filename = str(audit_row.get("Filename", "") or "")
        ped_row = ped_by_mode.get((filename, mode), {})
        stage3d_assignment = str(audit_row.get("functional_group_assignment", "") or "")
        ped_assignment = str(ped_row.get("ped_assignment", "") or "")
        ped_top_percent = float(ped_row.get("ped_top_percent", 0.0) or 0.0)
        ped_source = ped_source_label if ped_assignment else ""
        status, warning = classify_ped_stage3d_agreement(
            stage3d_assignment,
            ped_assignment,
            ped_top_percent,
            ped_source,
        )
        rows.append(
            {
                "Source": audit_row.get("Source", ""),
                "Filename": filename,
                "mode": mode,
                "frequency_cm-1": audit_row.get("frequency_cm-1", ""),
                "stage3d_assignment": stage3d_assignment,
                "ped_assignment": ped_assignment,
                "ped_source": ped_source,
                "ped_top_family": str(ped_row.get("ped_top_family", "") or ""),
                "ped_top_percent": ped_top_percent,
                "ped_agreement_status": status,
                "ped_policy_warning": warning,
                "ped_top_contributors": str(ped_row.get("ped_top_contributors", "") or ""),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def build_ped_driven_final_assignment_table(agreement_table: pd.DataFrame | None) -> pd.DataFrame:
    agreement_df = agreement_table if isinstance(agreement_table, pd.DataFrame) else pd.DataFrame()
    columns = [
        "Source",
        "Filename",
        "mode",
        "frequency_cm-1",
        "final_assignment",
        "final_assignment_source",
        "final_assignment_policy",
        "final_assignment_warning",
        "stage3d_assignment",
        "ped_assignment",
        "ped_source",
        "ped_agreement_status",
        "ped_policy_warning",
        "ped_top_family",
        "ped_top_percent",
        "ped_top_contributors",
    ]
    if agreement_df.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for _, row in agreement_df.iterrows():
        final_assignment, final_source, final_policy, final_warning = decide_ped_driven_final_assignment(
            row.get("stage3d_assignment", ""),
            row.get("ped_assignment", ""),
            row.get("ped_source", ""),
            row.get("ped_agreement_status", ""),
            row.get("ped_policy_warning", ""),
            row.get("ped_top_percent", 0.0),
        )
        rows.append(
            {
                "Source": row.get("Source", ""),
                "Filename": row.get("Filename", ""),
                "mode": row.get("mode", ""),
                "frequency_cm-1": row.get("frequency_cm-1", ""),
                "final_assignment": final_assignment,
                "final_assignment_source": final_source,
                "final_assignment_policy": final_policy,
                "final_assignment_warning": final_warning,
                "stage3d_assignment": row.get("stage3d_assignment", ""),
                "ped_assignment": row.get("ped_assignment", ""),
                "ped_source": row.get("ped_source", ""),
                "ped_agreement_status": row.get("ped_agreement_status", ""),
                "ped_policy_warning": row.get("ped_policy_warning", ""),
                "ped_top_family": row.get("ped_top_family", ""),
                "ped_top_percent": row.get("ped_top_percent", ""),
                "ped_top_contributors": row.get("ped_top_contributors", ""),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def build_spectrum_payload(
    hess_list,
    assignment_audit: pd.DataFrame | None = None,
    *,
    wilson_ped_audit: pd.DataFrame | None = None,
    ped_v2_force_audit: pd.DataFrame | None = None,
    ped_audit: pd.DataFrame | None = None,
    nist_reference_sets: Dict[str, object] | None = None,
) -> Dict[str, object]:
    from chemistry import (
        build_connectivity as chemistry_build_connectivity,
        classify_system as chemistry_classify_system,
        formula_string as chemistry_formula_string,
        split_fragments as chemistry_split_fragments,
    )

    audit_df = assignment_audit if isinstance(assignment_audit, pd.DataFrame) else pd.DataFrame()
    ped_df, ped_source_label = _select_ped_diagnostic_layer(
        wilson_ped_audit=wilson_ped_audit,
        ped_v2_force_audit=ped_v2_force_audit,
        ped_audit=ped_audit,
    )
    ped_by_mode = _build_ped_viewer_mode_summary(ped_df)
    files = []

    for hess in hess_list:
        rows = []
        file_audit = audit_df[audit_df.get("Filename", pd.Series(dtype=str)).astype(str) == str(hess.filename)].copy() if not audit_df.empty and "Filename" in audit_df.columns else pd.DataFrame()
        audit_by_mode = {}
        if not file_audit.empty and "mode" in file_audit.columns:
            for _, row in file_audit.iterrows():
                try:
                    audit_by_mode[int(row["mode"])] = row
                except (TypeError, ValueError):
                    continue

        for mode, (freq, intensity) in enumerate(zip(hess.frequencies_cm1, hess.ir_intensities)):
            if not pd.notna(freq) or float(freq) <= 0.0:
                continue
            audit_row = audit_by_mode.get(mode)
            ped_row = ped_by_mode.get((str(hess.filename), int(mode)), {})
            stage3d_assignment = str(audit_row.get("functional_group_assignment", "")) if audit_row is not None else ""
            ped_assignment = str(ped_row.get("ped_assignment", "") or "")
            ped_source = ped_source_label if ped_assignment else ""
            ped_agreement_status, ped_policy_warning = classify_ped_stage3d_agreement(
                stage3d_assignment,
                ped_assignment,
                float(ped_row.get("ped_top_percent", 0.0) or 0.0),
                ped_source,
            )
            final_assignment, final_source, final_policy, final_warning = decide_ped_driven_final_assignment(
                stage3d_assignment,
                ped_assignment,
                ped_source,
                ped_agreement_status,
                ped_policy_warning,
                float(ped_row.get("ped_top_percent", 0.0) or 0.0),
            )
            rows.append(
                {
                    "mode": int(mode),
                    "frequency_cm1": float(freq),
                    "intensity": float(intensity),
                    "assignment": final_assignment,
                    "final_assignment": final_assignment,
                    "final_assignment_source": final_source,
                    "final_assignment_policy": final_policy,
                    "final_assignment_warning": final_warning,
                    "stage3d_assignment": stage3d_assignment,
                    "ped_assignment": ped_assignment,
                    "ped_source": ped_source,
                    "ped_top_family": str(ped_row.get("ped_top_family", "") or ""),
                    "ped_top_percent": float(ped_row.get("ped_top_percent", 0.0) or 0.0),
                    "ped_agreement_status": ped_agreement_status,
                    "ped_policy_warning": ped_policy_warning,
                    "ped_top_contributors": str(ped_row.get("ped_top_contributors", "") or ""),
                    "ped_method": str(ped_row.get("ped_method", "") or ""),
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
        "viewer_assignment_source": ped_source_label or "Stage 3D",
        "files": files,
        "nist_reference_sets": nist_reference_sets or {},
    }


def _three_dmol_script_tag(three_dmol_js_path: str | Path | None = None) -> str:
    if three_dmol_js_path is None:
        return '<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>'

    js_path = Path(three_dmol_js_path)
    if not js_path.exists() or not js_path.is_file():
        raise FileNotFoundError(f"3Dmol.js asset not found: {js_path}")
    script_text = js_path.read_text(encoding="utf-8")
    script_text = script_text.replace("</script", "<\\/script")
    return f"<script>\n/* Inlined local 3Dmol.js asset: {js_path.name} */\n{script_text}\n</script>"


def _json_for_script(value: object) -> str:
    text = json.dumps(value, ensure_ascii=False)
    return (
        text.replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def write_interactive_spectrum_viewer(
    payload: Dict[str, object],
    html_path: str | Path,
    json_path: str | Path | None = None,
    three_dmol_js_path: str | Path | None = None,
) -> Path:
    html_path = Path(html_path)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    if json_path is not None:
        json_path = Path(json_path)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    payload_json = _json_for_script(payload)
    three_dmol_script_tag = _three_dmol_script_tag(three_dmol_js_path)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ORCAVEDA Interactive IR Spectrum</title>
  {three_dmol_script_tag}
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
    .engine-table-wrap {{
      border: 1px solid var(--border);
      border-radius: 14px;
      background: rgba(255,255,255,0.78);
      overflow: hidden;
      margin-bottom: 4px;
    }}
    .engine-table-scroll {{
      max-height: 160px;
      overflow: auto;
    }}
    .engine-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 11px;
      line-height: 1.25;
    }}
    .engine-table thead th {{
      position: sticky;
      top: 0;
      z-index: 1;
      background: #f8f4eb;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.04em;
      font-size: 10px;
      text-align: left;
      padding: 7px 8px;
      border-bottom: 1px solid var(--border);
      white-space: nowrap;
    }}
    .engine-table tbody td {{
      padding: 6px 8px;
      border-bottom: 1px solid #efe8da;
      vertical-align: top;
    }}
    .engine-table tbody tr.active {{
      background: rgba(15, 118, 110, 0.10);
    }}
    .engine-table tbody tr:last-child td {{
      border-bottom: 0;
    }}
    .engine-table .num {{
      text-align: right;
      white-space: nowrap;
      font-variant-numeric: tabular-nums;
    }}
    .inline-link-row {{
      margin-top: 5px;
      min-height: 16px;
    }}
    .inline-link-row a {{
      color: var(--accent);
      font-size: 11px;
      text-decoration: none;
      border-bottom: 1px solid transparent;
    }}
    .inline-link-row a:hover {{
      border-bottom-color: currentColor;
    }}
    .inline-link-row a.disabled {{
      color: #93a1af;
      pointer-events: none;
      border-bottom-color: transparent;
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
    .molecule-fallback-canvas {{
      display: block;
      width: 100%;
      height: 100%;
      min-height: 500px;
    }}
    .molecule-fallback-note {{
      position: absolute;
      left: 12px;
      bottom: 10px;
      max-width: calc(100% - 24px);
      color: #677483;
      font-size: 12px;
      background: rgba(255,253,248,0.86);
      border: 1px solid #e4ddd0;
      border-radius: 10px;
      padding: 6px 8px;
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
              <div class="inline-link-row"><a id="nistReferenceLink" href="#" target="_blank" rel="noopener noreferrer">Open on NIST</a></div>
            </div>
            <div class="control">
              <label for="scaleEngine">Scale Engine</label>
              <select id="scaleEngine">
                <option value="manual_static" selected>Manual Static</option>
                <option value="global_ls">Global LS</option>
                <option value="global_weighted_ls">Weighted LS</option>
                <option value="global_huber">Huber</option>
                <option value="piecewise_region">Piecewise Region</option>
                <option value="power_law">Power Law</option>
              </select>
            </div>
            <div class="control">
              <label for="matchingLayer">Matching Layer</label>
              <select id="matchingLayer">
                <option value="nearest">Nearest</option>
                <option value="extended" selected>Extended</option>
                <option value="high_confidence">High-Confidence</option>
              </select>
            </div>
            <div class="control">
              <label>NIST Fit</label>
              <div class="button-row">
                <button id="autoFitScale" type="button" class="ghost-btn">Auto-fit scale</button>
                <button id="resetZoom" type="button" class="ghost-btn">Reset Zoom</button>
              </div>
            </div>
          </div>
          <div id="fitSummary" class="fit-summary">Choose a NIST reference and engine. Manual Static keeps the slider and Auto-fit; the other engines use precomputed matched-peak fits.</div>
          <div id="matchingLayerSummary" class="fit-summary">Matching layers will appear here after you choose a NIST reference.</div>
          <div class="engine-table-wrap">
            <div class="engine-table-scroll">
              <table class="engine-table">
                <thead>
                  <tr>
                    <th>Layer</th>
                    <th class="num">Matched</th>
                    <th class="num">Coverage</th>
                    <th class="num">Mean %Δ</th>
                  </tr>
                </thead>
                <tbody id="matchingLayerTableBody"></tbody>
              </table>
            </div>
          </div>
          <div class="engine-table-wrap">
            <div class="engine-table-scroll">
              <table class="engine-table">
                <thead>
                  <tr>
                    <th>Engine</th>
                    <th class="num">Nearest %Δ</th>
                    <th class="num">High-Conf %Δ</th>
                    <th class="num">Extended %Δ</th>
                  </tr>
                </thead>
                <tbody id="engineLayerMatrixBody"></tbody>
              </table>
            </div>
          </div>
          <div class="engine-table-wrap">
            <div class="engine-table-scroll">
              <table class="engine-table">
                <thead>
                  <tr>
                    <th>Engine</th>
                    <th class="num">Mean %Δ</th>
                    <th class="num">RMS %Δ</th>
                    <th class="num">Max %Δ</th>
                    <th class="num">Matched</th>
                  </tr>
                </thead>
                <tbody id="engineTableBody"></tbody>
              </table>
            </div>
          </div>
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
                  <th>Mode Interpretation</th>
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
    const nistReferenceLink = document.getElementById("nistReferenceLink");
    const scaleEngine = document.getElementById("scaleEngine");
    const matchingLayer = document.getElementById("matchingLayer");
    const autoFitScale = document.getElementById("autoFitScale");
    const resetZoom = document.getElementById("resetZoom");
    const fitSummary = document.getElementById("fitSummary");
    const matchingLayerSummary = document.getElementById("matchingLayerSummary");
    const matchingLayerTableBody = document.getElementById("matchingLayerTableBody");
    const engineLayerMatrixBody = document.getElementById("engineLayerMatrixBody");
    const engineTableBody = document.getElementById("engineTableBody");
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

    function clearElement(el) {{
      while (el.firstChild) el.removeChild(el.firstChild);
    }}

    function appendCell(row, text, className = "") {{
      const td = document.createElement("td");
      if (className) td.className = className;
      td.textContent = String(text ?? "");
      row.appendChild(td);
      return td;
    }}

    function appendEmptyRow(tbody, colspan, text) {{
      const tr = document.createElement("tr");
      const td = appendCell(tr, text);
      td.colSpan = colspan;
      td.style.color = "#677483";
      td.style.padding = "8px";
      tbody.appendChild(tr);
    }}

    function appendKv(container, label, value, wide = false) {{
      const div = document.createElement("div");
      div.className = wide ? "kv wide" : "kv";
      const strong = document.createElement("strong");
      strong.textContent = String(label ?? "");
      div.appendChild(strong);
      div.appendChild(document.createTextNode(String(value ?? "")));
      container.appendChild(div);
      return div;
    }}

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

    function updateNistReferenceLink() {{
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const href = String(referenceSpectrum?.nist_spectrum_url || "");
      if (!href) {{
        nistReferenceLink.href = "#";
        nistReferenceLink.textContent = "NIST link unavailable";
        nistReferenceLink.classList.add("disabled");
        return;
      }}
      nistReferenceLink.href = href;
      nistReferenceLink.textContent = "Open on NIST";
      nistReferenceLink.classList.remove("disabled");
    }}

    function populateReferenceOptions() {{
      const refSet = getCurrentReferenceSet();
      clearElement(nistReference);
      const noneOpt = document.createElement("option");
      noneOpt.value = "";
      noneOpt.textContent = "None";
      nistReference.appendChild(noneOpt);

      if (!refSet || !Array.isArray(refSet.reference_spectra)) {{
        nistReference.disabled = true;
        updateNistReferenceLink();
        return;
      }}
      nistReference.disabled = false;
      for (const item of refSet.reference_spectra) {{
        const opt = document.createElement("option");
        opt.value = String(item.index);
        const phase = item.phase_label || item.phase_tag || "unknown phase";
        const units = item.y_units ? `; ${{item.y_units}}` : "";
        const description = item.description ? `; ${{item.description}}` : "";
        const suitability = item.suitable_for_matching === false ? `; skipped: ${{item.suitability_reason || "unsuitable_reference"}}` : "";
        opt.textContent = `Index ${{item.index}} - ${{phase}}${{units}}${{description}}${{suitability}}`;
        nistReference.appendChild(opt);
      }}
      const preferred = refSet.preferred_reference || null;
      nistReference.value = preferred ? String(preferred.index || "") : "";
      updateNistReferenceLink();
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

    function elementColor(element) {{
      return {{
        H: "#f8fafc",
        C: "#374151",
        N: "#2563eb",
        O: "#dc2626",
        S: "#ca8a04",
        P: "#f97316",
        F: "#16a34a",
        Cl: "#22c55e",
        Br: "#92400e",
        I: "#7c3aed",
      }}[String(element || "")] || "#64748b";
    }}

    function renderMoleculeFallback(file) {{
      clearElement(moleculeViewerHost);
      moleculeViewerHost.dataset.renderer = "native-fallback";
      const canvasEl = document.createElement("canvas");
      canvasEl.className = "molecule-fallback-canvas";
      const rect = moleculeViewerHost.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      const width = Math.max(420, Math.round((rect.width || 640) * dpr));
      const height = Math.max(360, Math.round((rect.height || 500) * dpr));
      canvasEl.width = width;
      canvasEl.height = height;
      moleculeViewerHost.appendChild(canvasEl);

      const note = document.createElement("div");
      note.className = "molecule-fallback-note";
      note.textContent = "Native 2D projection used because 3Dmol.js is unavailable.";
      moleculeViewerHost.appendChild(note);

      const g = canvasEl.getContext("2d");
      g.fillStyle = "#fffdf8";
      g.fillRect(0, 0, width, height);
      const atoms = Array.isArray(file?.geometry?.atoms) ? file.geometry.atoms : [];
      const bonds = Array.isArray(file?.geometry?.bonds) ? file.geometry.bonds : [];
      if (!atoms.length) {{
        g.fillStyle = "#677483";
        g.font = `${{14 * dpr}}px Segoe UI`;
        g.fillText("No geometry atoms reported.", 24 * dpr, 36 * dpr);
        return;
      }}

      const xs = atoms.map(atom => Number(atom.x || 0));
      const ys = atoms.map(atom => Number(atom.y || 0));
      const minX = Math.min(...xs);
      const maxX = Math.max(...xs);
      const minY = Math.min(...ys);
      const maxY = Math.max(...ys);
      const spanX = Math.max(1e-6, maxX - minX);
      const spanY = Math.max(1e-6, maxY - minY);
      const pad = 52 * dpr;
      const scale = Math.min((width - pad * 2) / spanX, (height - pad * 2) / spanY);
      const projected = atoms.map(atom => ({{
        x: pad + (Number(atom.x || 0) - minX) * scale,
        y: height - pad - (Number(atom.y || 0) - minY) * scale,
        z: Number(atom.z || 0),
        element: String(atom.element || "X"),
      }}));

      g.lineCap = "round";
      g.lineWidth = 4 * dpr;
      for (const bond of bonds) {{
        const a = projected[Number(bond.i)];
        const b = projected[Number(bond.j)];
        if (!a || !b) continue;
        g.strokeStyle = "rgba(55, 65, 81, 0.45)";
        g.beginPath();
        g.moveTo(a.x, a.y);
        g.lineTo(b.x, b.y);
        g.stroke();
      }}

      const orderedAtoms = projected
        .map((atom, idx) => ({{ ...atom, idx }}))
        .sort((a, b) => a.z - b.z);
      for (const atom of orderedAtoms) {{
        const radius = (atom.element === "H" ? 10 : 15) * dpr;
        g.beginPath();
        g.arc(atom.x, atom.y, radius, 0, Math.PI * 2);
        g.fillStyle = elementColor(atom.element);
        g.fill();
        g.strokeStyle = "rgba(31, 41, 51, 0.42)";
        g.lineWidth = 1.5 * dpr;
        g.stroke();
        g.fillStyle = atom.element === "H" ? "#1f2933" : "#ffffff";
        g.font = `${{10 * dpr}}px Segoe UI`;
        g.textAlign = "center";
        g.textBaseline = "middle";
        g.fillText(atom.element, atom.x, atom.y);
      }}
    }}

    function renderMolecule() {{
      const file = getCurrentFile();
      const viewer = ensureMoleculeViewer();
      if (!viewer) {{
        renderMoleculeFallback(file);
        return;
      }}
      moleculeViewerHost.dataset.renderer = "3dmol";
      viewer.clear();
      viewer.addModel(geometryToXyz(file), "xyz");
      viewer.setStyle({{}}, currentMolStyle());
      viewer.zoomTo();
      viewer.resize();
      viewer.render();
    }}

    function getScaledModes(file, scale, engineName = getSelectedScaleEngine(), engineFit = getSelectedScaleEngineFit()) {{
      return (file.modes || [])
        .map(mode => ({{
          ...mode,
          scaled: transformFrequencyByEngine(mode.frequency_cm1, scale, engineName, engineFit),
        }}))
        .filter(mode => Number.isFinite(mode.scaled) && mode.scaled > 0);
    }}

    function defaultRange(file, scale) {{
      if (!file.modes.length) return [0, 4000];
      const scaled = getScaledModes(file, scale).map(mode => mode.scaled);
      if (!scaled.length) return [0, 4000];
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

    function getSelectedScaleEngine() {{
      return String(scaleEngine.value || "manual_static");
    }}

    function getSelectedMatchingLayer() {{
      return String(matchingLayer?.value || "extended");
    }}

    function getSelectedMatchingLayerPayload(basePayload) {{
      if (!basePayload) return null;
      const layerName = getSelectedMatchingLayer();
      const layer = basePayload?.matching_layers?.[layerName] || null;
      if (!layer) return basePayload;
      return {{
        ...basePayload,
        engine_table: Array.isArray(layer.engine_table) ? layer.engine_table : [],
        engine_fits: layer.engine_fits || {{}},
        matched_pairs: Array.isArray(layer.matched_pairs) ? layer.matched_pairs : [],
        matched_count: Number(layer.matched_count ?? 0),
        total_reference_peaks: Number(layer.total_reference_peaks ?? 0),
        matching_layer_name: layerName,
      }};
    }}

    function getSelectedScaleEnginePayload() {{
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const payload = referenceSpectrum?.scale_engine_payload || null;
      return getSelectedMatchingLayerPayload(payload);
    }}

    function getSelectedScaleEngineFit() {{
      const engineName = getSelectedScaleEngine();
      const payload = getSelectedScaleEnginePayload();
      if (!payload || !payload.engine_fits) return null;
      return payload.engine_fits[engineName] || null;
    }}

    function transformFrequencyByEngine(freq, scale, engineName, engineFit) {{
      const omega = Number(freq);
      if (!Number.isFinite(omega) || omega <= 0) return Number.NaN;
      if (!engineName || engineName === "manual_static" || !engineFit || !engineFit.parameters) {{
        return omega * scale;
      }}

      const params = engineFit.parameters || {{}};
      if (engineName === "global_ls" || engineName === "global_weighted_ls" || engineName === "global_huber") {{
        const k = Number(params.k);
        return Number.isFinite(k) ? omega * k : omega * scale;
      }}
      if (engineName === "power_law") {{
        const a = Number(params.a);
        const b = Number(params.b);
        return (Number.isFinite(a) && Number.isFinite(b) && omega > 0)
          ? omega * (a * Math.pow(omega, b))
          : omega * scale;
      }}
      if (engineName === "piecewise_region") {{
        const regions = Array.isArray(params.regions) ? params.regions : [];
        for (const region of regions) {{
          const range = region.range || [];
          const lo = Number(range[0]);
          const hi = Number(range[1]);
          const k = Number(region.k);
          if (!Number.isFinite(lo) || !Number.isFinite(hi) || !Number.isFinite(k)) continue;
          if (omega >= lo && omega < hi) return omega * k;
        }}
      }}
      return omega * scale;
    }}

    function updateScaleControlsState() {{
      const engineName = getSelectedScaleEngine();
      const manual = engineName === "manual_static";
      scaleFactor.disabled = !manual;
      autoFitScale.disabled = !manual;
    }}

    function updateMatchingLayerSummary() {{
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const payload = referenceSpectrum?.scale_engine_payload || null;
      const layers = payload?.matching_layers || {{}};
      const hc = layers.high_confidence || {{}};
      const ext = layers.extended || {{}};
      const total = Number(ext.total_reference_peaks ?? hc.total_reference_peaks ?? 0);
      const hcMatched = Number(hc.matched_count ?? 0);
      const extMatched = Number(ext.matched_count ?? 0);
      const active = getSelectedMatchingLayer() === "high_confidence" ? "High-Confidence" : "Extended";
      if (!total) {{
        matchingLayerSummary.textContent = "Matching layers will appear here after you choose a NIST reference.";
        return;
      }}
      matchingLayerSummary.textContent = `Active layer: ${{active}} | High-confidence matches: ${{hcMatched}}/${{total}} | Extended matches: ${{extMatched}}/${{total}}`;
    }}

    function prettyMatchingLayerName(layer) {{
      const mapping = {{
        nearest: "Nearest",
        high_confidence: "High-Confidence",
        extended: "Extended",
      }};
      return mapping[String(layer || "")] || String(layer || "n/a");
    }}

    function updateMatchingLayerTable() {{
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const basePayload = referenceSpectrum?.scale_engine_payload || null;
      clearElement(matchingLayerTableBody);
      const rows = Array.isArray(basePayload?.matching_layer_overview) ? basePayload.matching_layer_overview : [];
      if (!rows.length) {{
        appendEmptyRow(matchingLayerTableBody, 4, "No matching-layer summary is available for the selected NIST reference.");
        return;
      }}
      const activeLayer = getSelectedMatchingLayer();
      for (const row of rows) {{
        const tr = document.createElement("tr");
        if (String(row.layer) === String(activeLayer)) tr.classList.add("active");
        const total = Number(row.total_reference_peaks ?? 0);
        const matched = Number(row.matched_count ?? 0);
        const coverage = Number(row.coverage ?? NaN);
        const mean = Number(row.mean_percent_deviation ?? NaN);
        appendCell(tr, prettyMatchingLayerName(row.layer));
        appendCell(tr, `${{matched}}/${{total}}`, "num");
        appendCell(tr, `${{(coverage * 100).toFixed(1)}}%`, "num");
        appendCell(tr, `${{mean.toFixed(2)}}%`, "num");
        tr.addEventListener("click", () => {{
          if (row.layer && matchingLayer.value !== String(row.layer)) {{
            matchingLayer.value = String(row.layer);
            pinnedTooltipMode = null;
            chartTooltip.classList.remove("visible");
            drawSpectrum(false);
          }}
        }});
        matchingLayerTableBody.appendChild(tr);
      }}
    }}

    function updateEngineLayerMatrix() {{
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const basePayload = referenceSpectrum?.scale_engine_payload || null;
      clearElement(engineLayerMatrixBody);
      const rows = Array.isArray(basePayload?.engine_layer_matrix) ? basePayload.engine_layer_matrix : [];
      if (!rows.length) {{
        appendEmptyRow(engineLayerMatrixBody, 4, "No engine-by-layer comparison is available for the selected NIST reference.");
        return;
      }}
      const activeEngine = getSelectedScaleEngine();
      for (const row of rows) {{
        const tr = document.createElement("tr");
        if (String(row.engine) === String(activeEngine)) tr.classList.add("active");
        const nearestMean = Number(row.nearest_mean_percent_deviation ?? NaN);
        const highMean = Number(row.high_confidence_mean_percent_deviation ?? NaN);
        const extMean = Number(row.extended_mean_percent_deviation ?? NaN);
        appendCell(tr, prettyEngineName(row.engine));
        appendCell(tr, `${{nearestMean.toFixed(2)}}%`, "num");
        appendCell(tr, `${{highMean.toFixed(2)}}%`, "num");
        appendCell(tr, `${{extMean.toFixed(2)}}%`, "num");
        tr.addEventListener("click", () => {{
          if (row.engine && scaleEngine.value !== String(row.engine)) {{
            scaleEngine.value = String(row.engine);
            pinnedTooltipMode = null;
            chartTooltip.classList.remove("visible");
            drawSpectrum(false);
          }}
        }});
        engineLayerMatrixBody.appendChild(tr);
      }}
    }}

    function updateFitSummaryForCurrentContext() {{
      const engineName = getSelectedScaleEngine();
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const payload = getSelectedScaleEnginePayload();
      if (!referenceSpectrum) {{
        fitSummary.textContent = engineName === "manual_static"
          ? "Choose a NIST reference and press Auto-fit scale to estimate the best frequency scaling."
          : "Choose a NIST reference to use precomputed matched-peak scale-engine fits.";
        updateMatchingLayerSummary();
        updateMatchingLayerTable();
        return;
      }}

      if (engineName === "manual_static") {{
        const bestScale = Number(payload?.default_manual_scale);
        fitSummary.textContent = Number.isFinite(bestScale)
          ? `Reference loaded. Manual Static is active; Auto-fit can refine the slider around ${{bestScale.toFixed(3)}}.`
          : "Reference loaded. Manual Static is active; press Auto-fit scale to estimate the best frequency scaling.";
        updateMatchingLayerSummary();
        updateMatchingLayerTable();
        return;
      }}

      const fit = getSelectedScaleEngineFit();
      if (!fit) {{
        fitSummary.textContent = `No precomputed fit is available for ${{engineName}} on the selected reference spectrum.`;
        updateMatchingLayerSummary();
        updateMatchingLayerTable();
        return;
      }}
      const met = fit.metrics || {{}};
      fitSummary.textContent = `${{fit.engine}} | mean %Δ ${{Number(met.mean_percent_deviation ?? NaN).toFixed(2)}}% | RMS %Δ ${{Number(met.rmse_percent_deviation ?? NaN).toFixed(2)}}% | matched ${{Number(fit.matched_count ?? 0)}}`;
      updateMatchingLayerSummary();
      updateMatchingLayerTable();
    }}

    function prettyEngineName(engine) {{
      const mapping = {{
        manual_static: "Manual Static",
        global_ls: "Global LS",
        global_weighted_ls: "Weighted LS",
        global_huber: "Huber",
        piecewise_region: "Piecewise Region",
        power_law: "Power Law",
      }};
      return mapping[String(engine || "")] || String(engine || "n/a");
    }}

    function updateEngineTable() {{
      const payload = getSelectedScaleEnginePayload();
      const activeEngine = getSelectedScaleEngine();
      clearElement(engineTableBody);

      const rows = Array.isArray(payload?.engine_table) ? payload.engine_table : [];
      if (!rows.length) {{
        appendEmptyRow(engineTableBody, 5, "No engine comparison is available for the selected NIST reference.");
        return;
      }}

      for (const row of rows) {{
        const tr = document.createElement("tr");
        if (String(row.engine) === String(activeEngine)) tr.classList.add("active");
        appendCell(tr, prettyEngineName(row.engine));
        appendCell(tr, `${{Number(row.mean_percent_deviation ?? NaN).toFixed(2)}}%`, "num");
        appendCell(tr, `${{Number(row.rmse_percent_deviation ?? NaN).toFixed(2)}}%`, "num");
        appendCell(tr, `${{Number(row.max_percent_deviation ?? NaN).toFixed(2)}}%`, "num");
        appendCell(tr, Number(row.matched_count ?? 0), "num");
        tr.addEventListener("click", () => {{
          if (row.engine && scaleEngine.value !== String(row.engine)) {{
            scaleEngine.value = String(row.engine);
            pinnedTooltipMode = null;
            chartTooltip.classList.remove("visible");
            drawSpectrum(false);
          }}
        }});
        engineTableBody.appendChild(tr);
      }}
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
      if (getSelectedScaleEngine() !== "manual_static") {{
        updateFitSummaryForCurrentContext();
        return;
      }}
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

      const fitMessage = `Best scale ${{best.scale.toFixed(3)}} | mean Δ ${{best.meanAbsDelta.toFixed(1)}} cm-1 | matched ${{best.matchedCount}}/${{best.totalPeaks}}`;
      scaleFactor.value = best.scale.toFixed(3);
      pinnedTooltipMode = null;
      chartTooltip.classList.remove("visible");
      currentView = null;
      drawSpectrum(false);
      fitSummary.textContent = fitMessage;
    }}

    function buildSpectrum(file, scale, gamma, x1, x2, axisMode, engineName, engineFit) {{
      const n = 1500;
      const xs = [];
      const rawYs = [];
      const scaledModes = getScaledModes(file, scale, engineName, engineFit);
      let rawMax = 0;
      for (let i = 0; i < n; i += 1) {{
        const x = x1 + (x2 - x1) * i / (n - 1);
        let y = 0;
        for (const mode of scaledModes) {{
          y += lorentz(x, mode.scaled, Math.max(0, mode.intensity), gamma);
        }}
        xs.push(x);
        rawYs.push(y);
        rawMax = Math.max(rawMax, y);
      }}
      const ys = rawYs.map(y => transformYValue(y, axisMode, rawMax));
      const yMin = axisMode === "transmittance" ? Math.min(...ys) : 0;
      const yMax = axisMode === "transmittance" ? 1.0 : Math.max(...ys, 1e-6);
      return {{ xs, ys, yMin, yMax, rawMax, scaledModes }};
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
      clearElement(summaryGrid);
      for (const [key, value] of items) appendKv(summaryGrid, key, value);
    }}

    function updateModeDetails(file, scale, engineName = getSelectedScaleEngine(), engineFit = getSelectedScaleEngineFit()) {{
      const scaledModes = getScaledModes(file, scale, engineName, engineFit);
      const mode = scaledModes.find(row => row.mode === selectedMode) || scaledModes[0];
      clearElement(modeDetails);
      if (!mode) {{
        const empty = document.createElement("div");
        empty.className = "wide";
        empty.textContent = "No positive-frequency modes available.";
        modeDetails.appendChild(empty);
        return;
      }}
      appendKv(modeDetails, "Mode", mode.mode);
      appendKv(modeDetails, "Scaled Frequency", `${{mode.scaled.toFixed(2)}} cm-1`);
      appendKv(modeDetails, "Original Frequency", `${{mode.frequency_cm1.toFixed(2)}} cm-1`);
      appendKv(modeDetails, "IR Intensity", Number(mode.intensity).toFixed(4));
      appendKv(modeDetails, "Final Assignment", mode.final_assignment || mode.assignment || "unassigned", true);
      appendKv(modeDetails, "Final Assignment Source", mode.final_assignment_source || "Stage 3D assignment audit");
      appendKv(modeDetails, "Final Assignment Policy", mode.final_assignment_policy || "n/a", true);
      appendKv(modeDetails, "Final Assignment Warning", mode.final_assignment_warning || "n/a", true);
      appendKv(modeDetails, "Stage 3D Assignment", mode.stage3d_assignment || "n/a", true);
      appendKv(modeDetails, "PED Diagnostic Interpretation", mode.ped_assignment || "n/a", true);
      appendKv(modeDetails, "PED Agreement Status", mode.ped_agreement_status || "n/a", true);
      appendKv(modeDetails, "PED Policy Warning", mode.ped_policy_warning || "n/a", true);
      appendKv(modeDetails, "PED Top Contributor", mode.ped_top_family ? `${{mode.ped_top_family}} (${{Number(mode.ped_top_percent || 0).toFixed(1)}}%)` : "n/a", true);
      appendKv(modeDetails, "PED Contributors", mode.ped_top_contributors || "n/a", true);
      appendKv(modeDetails, "PED Method", mode.ped_method || "n/a", true);
      appendKv(modeDetails, "Warnings", mode.warnings || "none", true);
      appendKv(modeDetails, "Stage 3D Supporting Coordinates", mode.top_internal_coordinates || "n/a", true);
    }}

    function updatePeakTable(file, scale, x1, x2, engineName = getSelectedScaleEngine(), engineFit = getSelectedScaleEngineFit()) {{
      const rows = getScaledModes(file, scale, engineName, engineFit)
        .filter(mode => mode.scaled >= x1 && mode.scaled <= x2)
        .sort((a, b) => a.scaled - b.scaled);
      clearElement(peakTable);
      for (const mode of rows) {{
        const tr = document.createElement("tr");
        if (mode.mode === selectedMode) tr.classList.add("active");
        appendCell(tr, mode.mode);
        appendCell(tr, mode.scaled.toFixed(1));
        appendCell(tr, Number(mode.intensity).toFixed(3));
        appendCell(tr, mode.assignment || "unassigned");
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
      const center = Number.isFinite(mode.scaled)
        ? mode.scaled
        : transformFrequencyByEngine(mode.frequency_cm1, render.scale, render.engineName, render.engineFit);
      return tx(center, render.x1, render.x2, render.margin.left, render.plotW, invertAxis.checked);
    }}

    function nearestModeAtCanvasX(px, threshold = 16) {{
      if (!currentRender) return null;
      let best = null;
      let bestDx = Infinity;
      for (const mode of currentRender.scaledModes || []) {{
        const center = mode.scaled;
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
      clearElement(chartTooltip);
      const strong = document.createElement("strong");
      strong.textContent = `Mode ${{mode.mode}}`;
      chartTooltip.appendChild(strong);
      for (const text of [
        `${{Number(mode.scaled ?? transformFrequencyByEngine(mode.frequency_cm1, currentRender.scale, currentRender.engineName, currentRender.engineFit)).toFixed(1)}} cm-1`,
        `IR: ${{Number(mode.intensity).toFixed(3)}}`,
        mode.assignment || "unassigned",
      ]) {{
        const div = document.createElement("div");
        div.textContent = text;
        chartTooltip.appendChild(div);
      }}
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
      const engineName = getSelectedScaleEngine();
      const engineFit = getSelectedScaleEngineFit();
      const [x1, x2] = clampView(file, scale, currentView?.x1, currentView?.x2);
      currentView = {{ x1, x2 }};
      updateScaleControlsState();
      updateFitSummaryForCurrentContext();
      updateEngineLayerMatrix();
      updateEngineTable();
      scaleValue.textContent = engineName === "manual_static" ? scale.toFixed(3) : "engine";
      hwhmValue.textContent = gamma.toFixed(1);

      const width = canvas.width;
      const height = canvas.height;
      const margin = {{ left: 76, right: 24, top: 20, bottom: 64 }};
      const plotW = width - margin.left - margin.right;
      const plotH = height - margin.top - margin.bottom;

      const render = buildSpectrum(file, scale, gamma, x1, x2, axisMode, engineName, engineFit);
      const referenceSpectrum = getSelectedReferenceSpectrum();
      const referenceOverlay = buildReferenceOverlay(referenceSpectrum, x1, x2, axisMode);
      currentRender = {{ ...render, x1, x2, scale, gamma, axisMode, margin, plotW, plotH, file, referenceSpectrum, referenceOverlay, engineName, engineFit }};
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
        for (const mode of render.scaledModes) {{
          const center = mode.scaled;
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
      updateModeDetails(file, scale, engineName, engineFit);
      updatePeakTable(file, scale, x1, x2, engineName, engineFit);
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
    [scaleFactor, hwhm, yMode, showSticks, invertAxis, nistReference, scaleEngine, matchingLayer].forEach(el => {{
      el.addEventListener("input", () => {{
        if (el === yMode || el === nistReference || el === scaleEngine || el === matchingLayer) {{
          pinnedTooltipMode = null;
          chartTooltip.classList.remove("visible");
          if (el === nistReference) {{
            updateNistReferenceLink();
            fitSummary.textContent = "Reference changed. Manual Static can use Auto-fit; the other engines use precomputed matched-peak fits.";
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
