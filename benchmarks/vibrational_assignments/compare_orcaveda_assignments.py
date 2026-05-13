#!/usr/bin/env python3
"""Compare ORCAVEDA assignment-audit labels with curated benchmark rows.

This is a benchmark diagnostic, not a chemistry oracle. It matches each
benchmark row to the nearest ORCAVEDA mode by frequency and checks whether the
human-readable assignment has overlapping chemistry terms.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Sequence

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scale_factor_engine import apply_scale  # noqa: E402


NUMBER_RE = re.compile(r"-?\d+(?:\.\d+)?")


def first_number(value: object) -> float | None:
    text = "" if pd.isna(value) else str(value)
    match = NUMBER_RE.search(text)
    return float(match.group(0)) if match else None


def target_frequency(row: pd.Series) -> tuple[float | None, str]:
    for column, label in (
        ("calculated_cm1", "benchmark_calculated"),
        ("observed_ir_cm1", "observed_ir"),
        ("observed_raman_cm1", "observed_raman"),
    ):
        freq = first_number(row.get(column, ""))
        if freq is not None:
            return freq, label
    return None, "missing"


def classes_from_text(*parts: object) -> set[str]:
    text = " ".join("" if pd.isna(part) else str(part) for part in parts).lower()

    classes: set[str] = set()
    if re.search(r"\bo[- ]?h\b", text) or "phenolic" in text:
        classes.add("oh")
    if "carboxylic acid" in text or "carboxylic" in text or "cooh" in text:
        classes.add("carboxylic_acid")
        classes.add("acid_context")
    if "carboxylic o-h" in text or "carboxylic oh" in text or "o-h torsion" in text or "o-h bend" in text:
        classes.add("carboxylic_oh_context")
        classes.add("acid_context")
    if re.search(r"\bn[- ]?h2?\b", text) or "amido" in text or "amine" in text:
        classes.add("nh")
    if re.search(r"\bc[- ]?c[- ]?h\b|\bc[- ]?h\b|\bch3\b", text) or "methyl" in text:
        classes.add("ch")
    if "methyl" in text or "ch3" in text:
        classes.add("methyl")
    if "c=o" in text or "carbonyl" in text or "ketone" in text:
        classes.add("carbonyl")
    if "carboxylic c=o" in text or "carboxylic_co" in text:
        classes.add("carboxyl_carbonyl")
        classes.add("acid_context")
    if re.search(r"\bc[- ]?o\b", text) or "phenolic c-o" in text or "carboxylic" in text:
        classes.add("co")
    if "c-c-o bend" in text or "o-c-o bend" in text or "c-o-o" in text or "carboxyl" in text:
        classes.add("carboxyl_deformation")
        classes.add("acid_context")
    if re.search(r"\bc[- ]?n\b", text) or "pyridine" in text or "heteroaromatic" in text:
        classes.add("cn")
    if re.search(r"\bc[- ]?c\b", text):
        classes.add("cc")
    if any(needle in text for needle in ("aromatic", "ring", "benzene", "phenol", "pyridine")):
        classes.add("aromatic_ring")
    if any(needle in text for needle in ("stretch", " str", "nu", "ν")):
        classes.add("stretch")
    if any(needle in text for needle in ("bend", "deform", "sciss", "umbrella", "delta", "δ")):
        classes.add("bend")
    if "out-of-plane" in text or "oop" in text:
        classes.add("out_of_plane")
    if "in-plane" in text or " ip" in text or "_ip" in text:
        classes.add("in_plane")
    if "torsion" in text or "torsional" in text:
        classes.add("torsion")
    if "wag" in text or "inversion" in text:
        classes.add("wag")
    if "mixed" in text or "+" in text or "fermi" in text:
        classes.add("mixed")
    if "aromatic_ring" in classes and "cn" in classes:
        classes.add("heteroaromatic_ring")
    return classes


def semantic_status(expected: set[str], actual: set[str]) -> tuple[str, str]:
    if not expected:
        return "WARN", "no_expected_classes"
    overlap = expected & actual
    if {"stretch", "bend", "torsion", "wag"} & expected and not ({"stretch", "bend", "torsion", "wag"} & actual):
        return "FAIL", "motion_family_mismatch"
    if "carbonyl" in expected and "carbonyl" not in actual:
        return "FAIL", "missing_carbonyl"
    acid_oh_expected = "oh" in expected and "carboxylic_acid" in expected
    if acid_oh_expected and "oh" not in actual and "acid_context" in actual:
        if "mixed" in expected or "carbonyl" in expected or "aromatic_ring" in expected:
            return "WARN", "acid_context_without_explicit_oh"
    if "oh" in expected and "oh" not in actual:
        return "FAIL", "missing_oh"
    if "nh" in expected and "nh" not in actual:
        return "FAIL", "missing_nh"
    if "ch" in expected and "ch" not in actual:
        return "FAIL", "missing_ch"
    if "aromatic_ring" in expected and ("aromatic_ring" not in actual and "heteroaromatic_ring" not in actual):
        return "WARN", "missing_explicit_ring_label"
    if overlap:
        return "PASS", "class_overlap:" + "|".join(sorted(overlap))
    return "WARN", "no_class_overlap"


def pick_benchmark_match(
    molecule_audit: pd.DataFrame,
    freq: float,
    expected: set[str],
    *,
    frequency_column: str,
    max_delta_cm1: float,
) -> tuple[pd.Series, float, str, str, str]:
    nearest_idx = (molecule_audit[frequency_column] - freq).abs().idxmin()
    nearest = molecule_audit.loc[nearest_idx]
    nearest_delta = float(nearest[frequency_column] - freq)

    semantic_window = max(float(max_delta_cm1) * 1.5, 120.0)
    candidates: list[tuple[int, float, float, str, str, pd.Series]] = []
    for _, candidate in molecule_audit.iterrows():
        delta = float(candidate[frequency_column] - freq)
        if abs(delta) > semantic_window:
            continue
        actual = classes_from_text(
            candidate.get("functional_group_assignment", ""),
            candidate.get("top_internal_coordinates", ""),
            candidate.get("top1_coord", ""),
        )
        status, reason = semantic_status(expected, actual)
        if status == "FAIL":
            rank = 2
        elif status == "WARN":
            rank = 1
        else:
            rank = 0
        candidates.append((rank, abs(delta), delta, status, reason, candidate))

    if not candidates:
        actual = classes_from_text(
            nearest.get("functional_group_assignment", ""),
            nearest.get("top_internal_coordinates", ""),
            nearest.get("top1_coord", ""),
        )
        status, reason = semantic_status(expected, actual)
        return nearest, nearest_delta, "nearest_frequency", status, reason

    candidates.sort(key=lambda item: (item[0], item[1]))
    rank, _abs_delta, delta, status, reason, chosen = candidates[0]
    if int(chosen.get("mode", -1)) == int(nearest.get("mode", -2)):
        strategy = "nearest_frequency"
    elif rank < 2:
        strategy = "semantic_within_window"
    else:
        strategy = "nearest_frequency"
        chosen = nearest
        delta = nearest_delta
        actual = classes_from_text(
            nearest.get("functional_group_assignment", ""),
            nearest.get("top_internal_coordinates", ""),
            nearest.get("top1_coord", ""),
        )
        status, reason = semantic_status(expected, actual)
    return chosen, delta, strategy, status, reason


def multiscale_semantic_summary(
    molecule_audit: pd.DataFrame,
    freq: float,
    expected: set[str],
    *,
    frequency_column: str,
    windows_cm1: Sequence[float],
) -> dict[str, object]:
    summary: dict[str, object] = {}
    for window in windows_cm1:
        window_value = float(window)
        half_width = window_value / 2.0
        key = str(int(window_value)) if window_value.is_integer() else str(window_value).replace(".", "p")
        in_window = molecule_audit[(molecule_audit[frequency_column] - freq).abs() <= half_width]
        found_classes: set[str] = set()
        best: tuple[int, float, str, str, pd.Series] | None = None
        for _, candidate in in_window.iterrows():
            actual = classes_from_text(
                candidate.get("functional_group_assignment", ""),
                candidate.get("top_internal_coordinates", ""),
                candidate.get("top1_coord", ""),
            )
            found_classes |= actual
            status, reason = semantic_status(expected, actual)
            rank = 0 if status == "PASS" else (1 if status == "WARN" else 2)
            delta_abs = abs(float(candidate[frequency_column] - freq))
            item = (rank, delta_abs, status, reason, candidate)
            if best is None or item[:2] < best[:2]:
                best = item

        prefix = f"{frequency_column.replace('_frequency_cm-1', '')}_window_{key}"
        summary[f"{prefix}_classes_found"] = "|".join(sorted(found_classes))
        if best is None:
            summary[f"{prefix}_status"] = "FAIL"
            summary[f"{prefix}_reason"] = "no_modes_in_window"
            summary[f"{prefix}_mode"] = ""
            summary[f"{prefix}_frequency_cm-1"] = ""
            summary[f"{prefix}_delta_cm-1"] = ""
            summary[f"{prefix}_assignment"] = ""
            continue

        _rank, _delta_abs, status, reason, candidate = best
        summary[f"{prefix}_status"] = status
        summary[f"{prefix}_reason"] = reason
        summary[f"{prefix}_mode"] = candidate.get("mode", "")
        summary[f"{prefix}_frequency_cm-1"] = float(candidate[frequency_column])
        summary[f"{prefix}_delta_cm-1"] = float(candidate[frequency_column] - freq)
        summary[f"{prefix}_assignment"] = str(candidate.get("functional_group_assignment", ""))
    return summary


def ped_mode_summaries(ped_audit: pd.DataFrame, *, top_n: int = 4) -> pd.DataFrame:
    if ped_audit.empty:
        return pd.DataFrame(
            columns=[
                "Filename",
                "mode",
                "ped_top_contributors",
                "ped_top_family",
                "ped_top_percent",
                "ped_classes",
                "ped_warnings",
            ]
        )

    ped = ped_audit.copy()
    ped["mode"] = pd.to_numeric(ped["mode"], errors="coerce")
    if "ped_rank" not in ped.columns and "wilson_rank" in ped.columns:
        ped["ped_rank"] = ped["wilson_rank"]
    ped["ped_rank"] = pd.to_numeric(ped.get("ped_rank", 0), errors="coerce")
    ped["contribution_percent"] = pd.to_numeric(ped.get("contribution_percent", 0.0), errors="coerce")
    ped = ped[(ped["ped_rank"] >= 1) & (ped["ped_rank"] <= int(top_n))].copy()
    if ped.empty:
        return pd.DataFrame(columns=["Filename", "mode", "ped_top_contributors", "ped_top_family", "ped_top_percent", "ped_classes", "ped_warnings"])

    warning_column = "ped_warnings" if "ped_warnings" in ped.columns else "wilson_ped_warnings"
    ped["ped_term"] = ped.apply(
        lambda row: (
            f"{row.get('coordinate_family', '')} [{row.get('internal_coordinate', '')}] "
            f"{float(row.get('contribution_percent', 0.0) or 0.0):.1f}%"
        ),
        axis=1,
    )

    rows: list[dict[str, object]] = []
    for (filename, mode), group in ped.sort_values(["Filename", "mode", "ped_rank"]).groupby(["Filename", "mode"], dropna=False):
        top = group.iloc[0]
        text_parts = []
        for _, row in group.iterrows():
            text_parts.extend([row.get("coordinate_family", ""), row.get("internal_coordinate", ""), row.get("coordinate_class", "")])
        rows.append(
            {
                "Filename": filename,
                "mode": mode,
                "ped_top_contributors": "; ".join(group["ped_term"].astype(str)),
                "ped_top_family": top.get("coordinate_family", ""),
                "ped_top_percent": float(top.get("contribution_percent", 0.0) or 0.0),
                "ped_classes": "|".join(sorted(classes_from_text(*text_parts))),
                "ped_warnings": "; ".join(
                    sorted({str(value) for value in group.get(warning_column, pd.Series(dtype=str)).dropna() if str(value)})
                ),
            }
        )
    return pd.DataFrame(rows)


def rename_ped_summary(summary: pd.DataFrame, *, prefix: str) -> pd.DataFrame:
    """Rename PED summary columns so alternate evidence stays diagnostic-only."""
    renamed = summary.copy()
    renamed = renamed.rename(
        columns={
            "ped_top_contributors": f"{prefix}_top_contributors",
            "ped_top_family": f"{prefix}_top_family",
            "ped_top_percent": f"{prefix}_top_percent",
            "ped_classes": f"{prefix}_classes",
            "ped_warnings": f"{prefix}_warnings",
        }
    )
    return renamed


def audit_from_ped_final_assignment(final_assignment: pd.DataFrame) -> pd.DataFrame:
    """Adapt ped_final_assignment.csv to the legacy comparator audit schema."""
    final = final_assignment.copy()
    if "frequency_cm-1" not in final.columns and "frequency_cm1" in final.columns:
        final["frequency_cm-1"] = final["frequency_cm1"]
    adapted = pd.DataFrame(
        {
            "Source": final.get("Source", ""),
            "Filename": final.get("Filename", ""),
            "mode": final.get("mode", ""),
            "frequency_cm-1": final.get("frequency_cm-1", ""),
            "functional_group_assignment": final.get("final_assignment", ""),
            "assignment_confidence": final.get("final_assignment_policy", ""),
            "top_internal_coordinates": final.get("ped_top_contributors", ""),
            "top1_coord": final.get("ped_top_family", ""),
            "final_assignment_source": final.get("final_assignment_source", ""),
            "final_assignment_policy": final.get("final_assignment_policy", ""),
            "final_assignment_warning": final.get("final_assignment_warning", ""),
            "stage3d_assignment": final.get("stage3d_assignment", ""),
            "ped_assignment": final.get("ped_assignment", ""),
            "ped_source": final.get("ped_source", ""),
            "ped_agreement_status": final.get("ped_agreement_status", ""),
            "ped_policy_warning": final.get("ped_policy_warning", ""),
            "ped_top_family": final.get("ped_top_family", ""),
            "ped_top_percent": final.get("ped_top_percent", 0.0),
            "ped_top_contributors": final.get("ped_top_contributors", ""),
        }
    )
    return adapted


def ped_coverage_audit(final_assignment: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    final = final_assignment.copy()
    if final.empty:
        columns = [
            "molecule",
            "Filename",
            "modes",
            "ped_final",
            "stage3d_fallback",
            "disagrees",
            "diffuse",
            "not_available",
            "ped_final_fraction",
            "stage3d_fallback_fraction",
        ]
        return pd.DataFrame(columns=columns), pd.DataFrame()

    final["mode"] = pd.to_numeric(final.get("mode", 0), errors="coerce")
    final["frequency_cm-1"] = pd.to_numeric(final.get("frequency_cm-1", 0.0), errors="coerce")
    final = final[final["frequency_cm-1"] > 0.0].copy()
    final["final_assignment_policy"] = final.get("final_assignment_policy", "").fillna("").astype(str)
    final["ped_agreement_status"] = final.get("ped_agreement_status", "").fillna("").astype(str)
    final["Filename"] = final.get("Filename", "").fillna("").astype(str)
    final["molecule"] = final["Filename"].str.replace(r"_freq\.hess$|\.hess$", "", regex=True)

    final["is_ped_final"] = final["final_assignment_policy"].str.startswith("ped_")
    final["is_stage3d_fallback"] = final["final_assignment_policy"].str.startswith("stage3d_fallback")
    rows: list[dict[str, object]] = []
    for (molecule, filename), group in final.groupby(["molecule", "Filename"], dropna=False):
        modes = int(len(group))
        ped_final = int(group["is_ped_final"].sum())
        stage3d_fallback = int(group["is_stage3d_fallback"].sum())
        rows.append(
            {
                "molecule": molecule,
                "Filename": filename,
                "modes": modes,
                "ped_final": ped_final,
                "stage3d_fallback": stage3d_fallback,
                "disagrees": int((group["ped_agreement_status"] == "disagrees").sum()),
                "diffuse": int((group["ped_agreement_status"] == "diffuse").sum()),
                "not_available": int((group["ped_agreement_status"] == "not_available").sum()),
                "ped_final_fraction": ped_final / modes if modes else 0.0,
                "stage3d_fallback_fraction": stage3d_fallback / modes if modes else 0.0,
            }
        )
    return pd.DataFrame(rows), final


def ped_diagnostic(
    expected: set[str],
    stage3d_classes: set[str],
    ped_classes: set[str],
    *,
    stage3d_status: str,
    ped_top_percent: float,
    match_strategy: str,
    delta_cm1: float,
    max_delta_cm1: float,
) -> tuple[str, str, str, str]:
    if not ped_classes:
        return "WARN", "no_ped_classes", "", "ped_not_reported"

    ped_status, ped_reason = semantic_status(expected, ped_classes)
    overlap = "|".join(sorted(stage3d_classes & ped_classes))
    warnings: list[str] = []
    if float(ped_top_percent) < 25.0:
        warnings.append("ped_diffuse_contributions")
    if str(match_strategy) == "semantic_within_window":
        warnings.append("nearest_frequency_differs_from_semantic_mode")
    if abs(float(delta_cm1)) > float(max_delta_cm1):
        warnings.append("selected_mode_outside_primary_tolerance")

    if ped_status == "PASS" and stage3d_status == "PASS":
        agreement = "ped_confirms_stage3d"
    elif ped_status == "PASS" and stage3d_status in {"WARN", "FAIL"}:
        agreement = "ped_supports_benchmark_semantics"
    elif ped_status == "FAIL" and stage3d_status == "PASS":
        agreement = "ped_disagrees_with_stage3d"
    elif overlap:
        agreement = "ped_adds_compatible_context"
    else:
        agreement = "ped_context_uncertain"

    return ped_status, ped_reason, overlap, "; ".join([agreement, *warnings])


def semantic_rank(status: object) -> int:
    text = "" if pd.isna(status) else str(status)
    return {"FAIL": 0, "WARN": 1, "PASS": 2}.get(text, 1)


def composed_policy_hint(
    *,
    ped_status: str,
    ped_reason: str,
    composed_status: str,
    ped_top_percent: float,
    composed_top_percent: float,
) -> str:
    """Summarize possible use of composed evidence without changing assignments."""
    ped_rank = semantic_rank(ped_status)
    composed_rank = semantic_rank(composed_status)
    if composed_rank == ped_rank and composed_status == "PASS" and float(composed_top_percent) > float(ped_top_percent) + 10.0:
        return "composed_confirms_with_better_localization"
    if composed_rank > ped_rank:
        if float(ped_top_percent) < 25.0 or str(ped_reason) == "no_ped_classes":
            return "diagnostic_hint_composed_available_when_baseline_diffuse_or_unclassified"
        return "diagnostic_hint_composed_semantic_improvement"
    if composed_rank == ped_rank and composed_status == "PASS":
        return "viewer_evidence_only"
    if float(composed_top_percent) > float(ped_top_percent) + 10.0:
        return "localization_gain_without_semantic_improvement"
    return "viewer_evidence_only"


def multiscale_ped_summary(
    molecule_audit: pd.DataFrame,
    freq: float,
    expected: set[str],
    *,
    frequency_column: str,
    windows_cm1: Sequence[float],
) -> dict[str, object]:
    summary: dict[str, object] = {}
    for window in windows_cm1:
        window_value = float(window)
        half_width = window_value / 2.0
        key = str(int(window_value)) if window_value.is_integer() else str(window_value).replace(".", "p")
        in_window = molecule_audit[(molecule_audit[frequency_column] - freq).abs() <= half_width]
        found_classes: set[str] = set()
        best: tuple[int, float, str, str, pd.Series] | None = None
        for _, candidate in in_window.iterrows():
            actual = classes_from_text(
                candidate.get("ped_top_contributors", ""),
                candidate.get("ped_top_family", ""),
                candidate.get("ped_classes", ""),
            )
            found_classes |= actual
            status, reason = semantic_status(expected, actual)
            rank = 0 if status == "PASS" else (1 if status == "WARN" else 2)
            delta_abs = abs(float(candidate[frequency_column] - freq))
            item = (rank, delta_abs, status, reason, candidate)
            if best is None or item[:2] < best[:2]:
                best = item

        prefix = f"{frequency_column.replace('_frequency_cm-1', '')}_ped_window_{key}"
        summary[f"{prefix}_classes_found"] = "|".join(sorted(found_classes))
        if best is None:
            summary[f"{prefix}_status"] = "FAIL"
            summary[f"{prefix}_reason"] = "no_modes_in_window"
            summary[f"{prefix}_mode"] = ""
            summary[f"{prefix}_frequency_cm-1"] = ""
            summary[f"{prefix}_delta_cm-1"] = ""
            summary[f"{prefix}_ped_top"] = ""
            continue

        _rank, _delta_abs, status, reason, candidate = best
        summary[f"{prefix}_status"] = status
        summary[f"{prefix}_reason"] = reason
        summary[f"{prefix}_mode"] = candidate.get("mode", "")
        summary[f"{prefix}_frequency_cm-1"] = float(candidate[frequency_column])
        summary[f"{prefix}_delta_cm-1"] = float(candidate[frequency_column] - freq)
        summary[f"{prefix}_ped_top"] = str(candidate.get("ped_top_contributors", ""))
    return summary


def compare(
    benchmark_csv: Path,
    audit_csv: Path,
    out_csv: Path,
    *,
    max_delta_cm1: float,
    windows_cm1: Sequence[float],
    scale_factor: float,
    primary_frequency: str,
    ped_audit_csv: Path | None = None,
    ped_final_assignment_csv: Path | None = None,
    composed_ped_audit_csv: Path | None = None,
) -> pd.DataFrame:
    benchmark = pd.read_csv(benchmark_csv)
    if ped_final_assignment_csv is not None:
        audit = audit_from_ped_final_assignment(pd.read_csv(ped_final_assignment_csv))
    else:
        audit = pd.read_csv(audit_csv)
    audit = audit[pd.to_numeric(audit["frequency_cm-1"], errors="coerce") > 0.0].copy()
    audit["frequency_cm-1"] = pd.to_numeric(audit["frequency_cm-1"], errors="coerce")
    audit["raw_frequency_cm-1"] = audit["frequency_cm-1"]
    audit["scaled_frequency_cm-1"] = apply_scale(audit["raw_frequency_cm-1"].to_numpy(), float(scale_factor))
    if ped_audit_csv is not None:
        ped_summary = ped_mode_summaries(pd.read_csv(ped_audit_csv))
        audit = audit.merge(ped_summary, on=["Filename", "mode"], how="left")
    if composed_ped_audit_csv is not None:
        composed_summary = rename_ped_summary(
            ped_mode_summaries(pd.read_csv(composed_ped_audit_csv), top_n=6),
            prefix="composed_ped",
        )
        audit = audit.merge(composed_summary, on=["Filename", "mode"], how="left")
    for column in ("ped_top_contributors", "ped_top_family", "ped_classes", "ped_warnings"):
        if column not in audit.columns:
            audit[column] = ""
        audit[column] = audit[column].fillna("")
    if "ped_top_percent" not in audit.columns:
        audit["ped_top_percent"] = 0.0
    audit["ped_top_percent"] = pd.to_numeric(audit["ped_top_percent"], errors="coerce").fillna(0.0)
    for column in ("composed_ped_top_contributors", "composed_ped_top_family", "composed_ped_classes", "composed_ped_warnings"):
        if column not in audit.columns:
            audit[column] = ""
        audit[column] = audit[column].fillna("")
    if "composed_ped_top_percent" not in audit.columns:
        audit["composed_ped_top_percent"] = 0.0
    audit["composed_ped_top_percent"] = pd.to_numeric(audit["composed_ped_top_percent"], errors="coerce").fillna(0.0)
    primary_column = "scaled_frequency_cm-1" if primary_frequency == "scaled" else "raw_frequency_cm-1"

    rows: list[dict[str, object]] = []
    for idx, bench in benchmark.iterrows():
        freq, freq_source = target_frequency(bench)
        molecule_audit = audit[audit["Filename"].astype(str) == str(bench["hess_file"])]
        if freq is None or molecule_audit.empty:
            rows.append(
                {
                    "benchmark_row": idx,
                    "molecule": bench.get("molecule", ""),
                    "hess_file": bench.get("hess_file", ""),
                    "status": "WARN",
                    "reason": "missing_target_frequency_or_audit",
                }
            )
            continue

        expected = classes_from_text(
            bench.get("assignment_normalized", ""),
            bench.get("mode_family", ""),
            bench.get("functional_group", ""),
        )
        nearest_idx = (molecule_audit[primary_column] - freq).abs().idxmin()
        nearest = molecule_audit.loc[nearest_idx]
        nearest_delta = float(nearest[primary_column] - freq)
        chosen, delta, match_strategy, status, reason = pick_benchmark_match(
            molecule_audit,
            freq,
            expected,
            frequency_column=primary_column,
            max_delta_cm1=max_delta_cm1,
        )
        assignment = str(chosen.get("functional_group_assignment", ""))
        actual = classes_from_text(
            assignment,
            chosen.get("top_internal_coordinates", ""),
            chosen.get("top1_coord", ""),
        )
        ped_classes = classes_from_text(
            chosen.get("ped_top_contributors", ""),
            chosen.get("ped_top_family", ""),
            chosen.get("ped_classes", ""),
        )
        if abs(delta) > max_delta_cm1 and status == "PASS":
            status = "WARN"
            reason = f"semantic_pass_but_frequency_delta_gt_{max_delta_cm1:g}"
        ped_status, ped_reason, stage3d_ped_overlap, stage3d_ped_warning = ped_diagnostic(
            expected,
            actual,
            ped_classes,
            stage3d_status=status,
            ped_top_percent=float(chosen.get("ped_top_percent", 0.0) or 0.0),
            match_strategy=match_strategy,
            delta_cm1=delta,
            max_delta_cm1=max_delta_cm1,
        )
        composed_ped_classes = classes_from_text(
            chosen.get("composed_ped_top_contributors", ""),
            chosen.get("composed_ped_top_family", ""),
            chosen.get("composed_ped_classes", ""),
        )
        (
            composed_ped_status,
            composed_ped_reason,
            stage3d_composed_ped_overlap,
            stage3d_composed_ped_warning,
        ) = ped_diagnostic(
            expected,
            actual,
            composed_ped_classes,
            stage3d_status=status,
            ped_top_percent=float(chosen.get("composed_ped_top_percent", 0.0) or 0.0),
            match_strategy=match_strategy,
            delta_cm1=delta,
            max_delta_cm1=max_delta_cm1,
        )
        semantic_delta = semantic_rank(composed_ped_status) - semantic_rank(ped_status)
        if semantic_delta > 0:
            composed_vs_baseline = "improves_semantic_match"
        elif semantic_delta < 0:
            composed_vs_baseline = "worsens_semantic_match"
        elif float(chosen.get("composed_ped_top_percent", 0.0) or 0.0) > float(chosen.get("ped_top_percent", 0.0) or 0.0) + 10.0:
            composed_vs_baseline = "localization_gain_without_semantic_improvement"
        else:
            composed_vs_baseline = "same_semantic_match"

        row = {
            "benchmark_row": idx,
            "molecule": bench.get("molecule", ""),
            "hess_file": bench.get("hess_file", ""),
            "confidence": bench.get("confidence", ""),
            "target_frequency_cm-1": freq,
            "target_frequency_source": freq_source,
            "primary_frequency": primary_frequency,
            "scale_factor": float(scale_factor),
            "match_strategy": match_strategy,
            "nearest_mode": nearest.get("mode", ""),
            "nearest_frequency_cm-1": float(nearest[primary_column]),
            "nearest_delta_cm-1": nearest_delta,
            "orcaveda_mode": chosen.get("mode", ""),
            "orcaveda_frequency_cm-1": float(chosen[primary_column]),
            "orcaveda_raw_frequency_cm-1": float(chosen["raw_frequency_cm-1"]),
            "orcaveda_scaled_frequency_cm-1": float(chosen["scaled_frequency_cm-1"]),
            "delta_cm-1": delta,
            "benchmark_assignment": bench.get("assignment_normalized", ""),
            "benchmark_mode_family": bench.get("mode_family", ""),
            "benchmark_functional_group": bench.get("functional_group", ""),
            "orcaveda_assignment": assignment,
            "orcaveda_confidence": chosen.get("assignment_confidence", ""),
            "final_assignment_source": chosen.get("final_assignment_source", ""),
            "final_assignment_policy": chosen.get("final_assignment_policy", ""),
            "final_assignment_warning": chosen.get("final_assignment_warning", ""),
            "stage3d_assignment": chosen.get("stage3d_assignment", ""),
            "ped_assignment": chosen.get("ped_assignment", ""),
            "ped_agreement_status": chosen.get("ped_agreement_status", ""),
            "ped_policy_warning": chosen.get("ped_policy_warning", ""),
            "expected_classes": "|".join(sorted(expected)),
            "actual_classes": "|".join(sorted(actual)),
            "ped_top_contributors": chosen.get("ped_top_contributors", ""),
            "ped_top_family": chosen.get("ped_top_family", ""),
            "ped_top_percent": float(chosen.get("ped_top_percent", 0.0) or 0.0),
            "ped_classes": "|".join(sorted(ped_classes)),
            "ped_semantic_status": ped_status,
            "ped_semantic_reason": ped_reason,
            "stage3d_ped_overlap_classes": stage3d_ped_overlap,
            "stage3d_ped_warning": stage3d_ped_warning,
            "composed_ped_top_contributors": chosen.get("composed_ped_top_contributors", ""),
            "composed_ped_top_family": chosen.get("composed_ped_top_family", ""),
            "composed_ped_top_percent": float(chosen.get("composed_ped_top_percent", 0.0) or 0.0),
            "composed_ped_classes": "|".join(sorted(composed_ped_classes)),
            "composed_ped_semantic_status": composed_ped_status,
            "composed_ped_semantic_reason": composed_ped_reason,
            "stage3d_composed_ped_overlap_classes": stage3d_composed_ped_overlap,
            "stage3d_composed_ped_warning": stage3d_composed_ped_warning,
            "composed_vs_baseline_ped_status": composed_vs_baseline,
            "composed_ped_localization_delta_percent": float(chosen.get("composed_ped_top_percent", 0.0) or 0.0)
            - float(chosen.get("ped_top_percent", 0.0) or 0.0),
            "composed_ped_policy_hint": composed_policy_hint(
                ped_status=ped_status,
                ped_reason=ped_reason,
                composed_status=composed_ped_status,
                ped_top_percent=float(chosen.get("ped_top_percent", 0.0) or 0.0),
                composed_top_percent=float(chosen.get("composed_ped_top_percent", 0.0) or 0.0),
            ),
            "status": status,
            "reason": reason,
        }
        row.update(
            multiscale_semantic_summary(
                molecule_audit,
                freq,
                expected,
                frequency_column="raw_frequency_cm-1",
                windows_cm1=windows_cm1,
            )
        )
        row.update(
            multiscale_semantic_summary(
                molecule_audit,
                freq,
                expected,
                frequency_column="scaled_frequency_cm-1",
                windows_cm1=windows_cm1,
            )
        )
        if ped_audit_csv is not None:
            row.update(
                multiscale_ped_summary(
                    molecule_audit,
                    freq,
                    expected,
                    frequency_column="raw_frequency_cm-1",
                    windows_cm1=windows_cm1,
                )
            )
            row.update(
                multiscale_ped_summary(
                    molecule_audit,
                    freq,
                    expected,
                    frequency_column="scaled_frequency_cm-1",
                    windows_cm1=windows_cm1,
                )
            )
        rows.append(row)
    result = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_csv, index=False)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=Path("benchmarks/vibrational_assignments/assignments.csv"))
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-delta-cm1", type=float, default=80.0)
    parser.add_argument("--windows-cm1", default="50,100,200,500", help="Comma-separated full window widths in cm-1")
    parser.add_argument("--scale-factor", type=float, default=1.0, help="Constant scale factor applied to ORCAVEDA frequencies for scaled diagnostics")
    parser.add_argument("--primary-frequency", choices=("raw", "scaled"), default="raw", help="Frequency axis used for legacy status columns")
    parser.add_argument("--ped-audit", type=Path, default=None, help="Optional ORCAVEDA ped_audit CSV to add PED-aware diagnostics")
    parser.add_argument("--ped-final-assignment", type=Path, default=None, help="Optional ORCAVEDA ped_final_assignment CSV; when set, benchmark labels are compared against PED-driven final labels")
    parser.add_argument("--composed-ped-audit", type=Path, default=None, help="Optional composed-coordinate PED audit CSV; adds separate diagnostic evidence columns only")
    parser.add_argument("--coverage-out", type=Path, default=None, help="Optional CSV path for PED final-label coverage by molecule")
    parser.add_argument("--coverage-detail-out", type=Path, default=None, help="Optional CSV path for per-mode PED final-label coverage detail")
    args = parser.parse_args()

    windows = [float(part.strip()) for part in str(args.windows_cm1).split(",") if part.strip()]
    result = compare(
        args.benchmark,
        args.audit,
        args.out,
        max_delta_cm1=args.max_delta_cm1,
        windows_cm1=windows,
        scale_factor=args.scale_factor,
        primary_frequency=args.primary_frequency,
        ped_audit_csv=args.ped_audit,
        ped_final_assignment_csv=args.ped_final_assignment,
        composed_ped_audit_csv=args.composed_ped_audit,
    )
    if args.ped_final_assignment is not None and args.coverage_out is not None:
        coverage, detail = ped_coverage_audit(pd.read_csv(args.ped_final_assignment))
        args.coverage_out.parent.mkdir(parents=True, exist_ok=True)
        coverage.to_csv(args.coverage_out, index=False)
        if args.coverage_detail_out is not None:
            args.coverage_detail_out.parent.mkdir(parents=True, exist_ok=True)
            detail.to_csv(args.coverage_detail_out, index=False)
    print(f"wrote {args.out}")
    print(result.groupby("status").size().to_string())
    print()
    print(result.groupby(["molecule", "status"]).size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
