from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd

from scale_factor_engine import (
    compare_scaling_models,
    fit_constant_scale,
)
from nist_ir.matching import match_reference_to_orcaveda_v2


def _phase_peak_picking_profile(
    reference_context: dict[str, object] | None = None,
) -> dict[str, float | int | None]:
    phase_text = " ".join(
        str((reference_context or {}).get(key, "") or "").lower()
        for key in ("phase_tag", "phase_label", "state", "description")
    )
    if any(token in phase_text for token in ("liquid", "solution", "solid")):
        return {
            "top_n_factor": 0.75,
            "min_separation_factor": 1.5,
            "relative_min_intensity": 0.12,
        }
    return {
        "top_n_factor": 1.0,
        "min_separation_factor": 1.0,
        "relative_min_intensity": None,
    }


def _prefer_minima(reference_context: dict[str, object] | None = None) -> bool:
    y_units = str((reference_context or {}).get("y_units", "") or "").lower()
    return "transmittance" in y_units


def pick_reference_peaks(
    spectrum: pd.DataFrame,
    *,
    top_n: int = 12,
    min_intensity: float | None = None,
    min_separation_cm1: float = 20.0,
    reference_context: dict[str, object] | None = None,
) -> pd.DataFrame:
    if spectrum.empty:
        return pd.DataFrame(columns=["wavenumber_cm-1", "intensity"])

    profile = _phase_peak_picking_profile(reference_context)
    adjusted_top_n = max(1, int(round(float(top_n) * float(profile["top_n_factor"]))))
    adjusted_min_separation = float(min_separation_cm1) * float(profile["min_separation_factor"])

    df = spectrum.copy()
    df = df.sort_values("wavenumber_cm-1").reset_index(drop=True)
    prefer_minima = _prefer_minima(reference_context)
    if min_intensity is None and profile["relative_min_intensity"] is not None:
        max_intensity = float(df["intensity"].max()) if not df.empty else 0.0
        min_intensity = max_intensity * float(profile["relative_min_intensity"])
    peaks: List[Dict[str, float]] = []
    values = df["intensity"].tolist()
    wns = df["wavenumber_cm-1"].tolist()

    for i in range(1, len(df) - 1):
        y = values[i]
        if prefer_minima:
            if y > values[i - 1] or y > values[i + 1]:
                continue
            if min_intensity is not None:
                max_intensity = float(df["intensity"].max()) if not df.empty else 0.0
                if (max_intensity - y) < min_intensity:
                    continue
            peaks.append({"wavenumber_cm-1": wns[i], "intensity": float((float(df["intensity"].max()) if not df.empty else 0.0) - y)})
            continue

        if y < values[i - 1] or y < values[i + 1]:
            continue
        if min_intensity is not None and y < min_intensity:
            continue
        peaks.append({"wavenumber_cm-1": wns[i], "intensity": y})

    peak_df = pd.DataFrame(peaks)
    if peak_df.empty:
        return peak_df

    peak_df = peak_df.sort_values("intensity", ascending=False).reset_index(drop=True)
    selected = []
    for _, row in peak_df.iterrows():
        wn = float(row["wavenumber_cm-1"])
        if all(abs(wn - item["wavenumber_cm-1"]) > adjusted_min_separation for item in selected):
            selected.append({"wavenumber_cm-1": wn, "intensity": float(row["intensity"])})
        if len(selected) >= adjusted_top_n:
            break

    return pd.DataFrame(selected).sort_values("wavenumber_cm-1", ascending=False).reset_index(drop=True)


def load_orcaveda_assignments(assignment_audit_csv: str | Path, *, scale_factor: float = 0.96) -> pd.DataFrame:
    audit = pd.read_csv(assignment_audit_csv, encoding="utf-8-sig")
    required = {"mode", "frequency_cm-1", "IR_intensity", "functional_group_assignment"}
    missing = required - set(audit.columns)
    if missing:
        raise ValueError(f"Assignment audit missing required columns: {sorted(missing)}")

    audit = audit.copy()
    audit["mode"] = audit["mode"].astype(int)
    audit["frequency_cm-1"] = audit["frequency_cm-1"].astype(float)
    audit["IR_intensity"] = audit["IR_intensity"].astype(float)
    audit = audit[audit["frequency_cm-1"] > 0.0].reset_index(drop=True)
    audit["scaled_frequency_cm-1"] = audit["frequency_cm-1"] * float(scale_factor)
    return audit


def assignment_modes_to_dataframe(modes: Sequence[dict], *, scale_factor: float = 1.0) -> pd.DataFrame:
    if not modes:
        return pd.DataFrame(
            columns=[
                "mode",
                "frequency_cm-1",
                "IR_intensity",
                "functional_group_assignment",
                "scaled_frequency_cm-1",
            ]
        )

    rows = []
    for row in modes:
        freq = float(row.get("frequency_cm1", row.get("frequency_cm-1", 0.0)))
        if freq <= 0.0:
            continue
        rows.append(
            {
                "mode": int(row.get("mode", 0)),
                "frequency_cm-1": freq,
                "IR_intensity": float(row.get("intensity", row.get("IR_intensity", 0.0))),
                "functional_group_assignment": str(row.get("assignment", row.get("functional_group_assignment", ""))),
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "mode",
                "frequency_cm-1",
                "IR_intensity",
                "functional_group_assignment",
                "scaled_frequency_cm-1",
            ]
        )
    df = pd.DataFrame(rows)
    df["scaled_frequency_cm-1"] = df["frequency_cm-1"] * float(scale_factor)
    return df


def reference_points_to_peaks(
    points: Sequence[dict],
    *,
    top_n: int = 12,
    min_intensity: float | None = None,
    min_separation_cm1: float = 20.0,
    reference_context: dict[str, object] | None = None,
) -> pd.DataFrame:
    if not points:
        return pd.DataFrame(columns=["wavenumber_cm-1", "intensity"])
    spectrum = pd.DataFrame(
        [
            {
                "wavenumber_cm-1": float(pt["x"]),
                "intensity": float(pt["y"]),
            }
            for pt in points
            if "x" in pt and "y" in pt
        ]
    )
    return pick_reference_peaks(
        spectrum,
        top_n=top_n,
        min_intensity=min_intensity,
        min_separation_cm1=min_separation_cm1,
        reference_context=reference_context,
    )


def match_reference_to_orcaveda(
    reference_peaks: pd.DataFrame,
    assignment_audit: pd.DataFrame,
) -> pd.DataFrame:
    if reference_peaks.empty or assignment_audit.empty:
        return pd.DataFrame(
            columns=[
                "reference_peak_cm-1",
                "reference_intensity",
                "orcaveda_mode",
                "orcaveda_scaled_frequency_cm-1",
                "delta_cm-1",
                "orcaveda_assignment",
                "orcaveda_intensity",
            ]
        )

    rows = []
    for _, peak in reference_peaks.iterrows():
        wn = float(peak["wavenumber_cm-1"])
        best_idx = (assignment_audit["scaled_frequency_cm-1"] - wn).abs().idxmin()
        best = assignment_audit.loc[best_idx]
        rows.append(
            {
                "reference_peak_cm-1": wn,
                "reference_intensity": float(peak["intensity"]),
                "orcaveda_mode": int(best["mode"]),
                "orcaveda_scaled_frequency_cm-1": float(best["scaled_frequency_cm-1"]),
                "delta_cm-1": float(best["scaled_frequency_cm-1"] - wn),
                "orcaveda_assignment": str(best.get("functional_group_assignment", "")),
                "orcaveda_intensity": float(best["IR_intensity"]),
            }
        )

    return pd.DataFrame(rows)


def build_matched_peak_pairs(
    reference_peaks: pd.DataFrame,
    assignment_audit: pd.DataFrame,
    *,
    method: str = "nearest",
    reference_context: dict[str, object] | None = None,
) -> pd.DataFrame:
    if method in {"scored", "scored_extended", "scored_high_confidence"}:
        include_stages = None
        if method == "scored_high_confidence":
            include_stages = {"primary"}
        matched = match_reference_to_orcaveda_v2(
            reference_peaks,
            assignment_audit,
            include_stages=include_stages,
            reference_context=reference_context,
        )
        if matched.empty:
            return pd.DataFrame(
                columns=[
                    "reference_peak_cm-1",
                    "reference_intensity",
                    "calc_frequency_cm-1",
                    "calc_intensity",
                    "mode",
                    "assignment",
                ]
            )
        return pd.DataFrame(
            {
                "reference_peak_cm-1": matched["reference_peak_cm-1"].astype(float),
                "reference_intensity": matched["reference_intensity"].astype(float),
                "calc_frequency_cm-1": matched["orcaveda_scaled_frequency_cm-1"].astype(float),
                "calc_intensity": matched["orcaveda_intensity"].astype(float),
                "mode": matched["orcaveda_mode"].astype(int),
                "assignment": matched["orcaveda_assignment"].astype(str),
            }
        )

    matched = match_reference_to_orcaveda(reference_peaks, assignment_audit)
    if matched.empty:
        return pd.DataFrame(
            columns=[
                "reference_peak_cm-1",
                "reference_intensity",
                "calc_frequency_cm-1",
                "calc_intensity",
                "mode",
                "assignment",
            ]
        )
    return pd.DataFrame(
        {
            "reference_peak_cm-1": matched["reference_peak_cm-1"].astype(float),
            "reference_intensity": matched["reference_intensity"].astype(float),
            "calc_frequency_cm-1": matched["orcaveda_scaled_frequency_cm-1"].astype(float),
            "calc_intensity": matched["orcaveda_intensity"].astype(float),
            "mode": matched["orcaveda_mode"].astype(int),
            "assignment": matched["orcaveda_assignment"].astype(str),
        }
    )


def compare_scale_engines_on_matched_peaks(
    matched_pairs: pd.DataFrame,
    *,
    piecewise_regions: Sequence[tuple[float, float]] | None = None,
) -> pd.DataFrame:
    if matched_pairs.empty:
        return pd.DataFrame(
            columns=[
                "engine",
                "parameters_json",
                "mean_percent_deviation",
                "rmse_percent_deviation",
                "max_percent_deviation",
                "matched_count",
            ]
        )

    omega = matched_pairs["calc_frequency_cm-1"].astype(float).to_numpy()
    nu = matched_pairs["reference_peak_cm-1"].astype(float).to_numpy()
    weights = matched_pairs["reference_intensity"].astype(float).to_numpy()
    if np.max(weights) > 0:
        weights = weights / np.max(weights)
    weights = np.clip(weights, 1e-6, None)

    models = compare_scaling_models(
        omega,
        nu,
        weights=weights,
        piecewise_regions=piecewise_regions,
    )
    rows = []
    for engine_name, payload in models.items():
        met = payload.get("metrics")
        if not met:
            continue
        params = payload.get("params", {})
        if engine_name in {"global_ls", "global_weighted_ls", "global_huber"}:
            pred = omega * float(params.get("k", 1.0))
        elif engine_name == "power_law":
            a = float(params.get("a", 1.0))
            b = float(params.get("b", 0.0))
            pred = omega * (a * omega**b)
        elif engine_name == "piecewise_region":
            pred = np.full_like(omega, np.nan, dtype=float)
            for region in params.get("regions", []):
                try:
                    lo, hi = region["range"]
                    k = float(region["k"])
                except (KeyError, TypeError, ValueError):
                    continue
                if not np.isfinite(k):
                    continue
                mask = (omega >= float(lo)) & (omega < float(hi))
                pred[mask] = omega[mask] * k
        else:
            pred = omega.copy()

        valid_mask = np.isfinite(pred)
        abs_percent = np.abs((pred[valid_mask] - nu[valid_mask]) / np.clip(np.abs(nu[valid_mask]), 1e-9, None)) * 100.0
        rows.append(
            {
                "engine": engine_name,
                "parameters_json": json.dumps(params, ensure_ascii=False),
                "mean_percent_deviation": float(np.mean(abs_percent)) if abs_percent.size else float("nan"),
                "rmse_percent_deviation": float(np.sqrt(np.mean(abs_percent**2))) if abs_percent.size else float("nan"),
                "max_percent_deviation": float(np.max(abs_percent)) if abs_percent.size else float("nan"),
                "matched_count": int(np.sum(valid_mask)),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["mean_percent_deviation", "rmse_percent_deviation", "max_percent_deviation"],
        ascending=[True, True, True],
    ).reset_index(drop=True)


def build_scale_engine_payload(
    reference_peaks: pd.DataFrame,
    assignment_audit: pd.DataFrame,
    *,
    piecewise_regions: Sequence[tuple[float, float]] | None = None,
    matching_method: str = "nearest",
    reference_context: dict[str, object] | None = None,
) -> dict:
    matched_pairs = build_matched_peak_pairs(reference_peaks, assignment_audit, method=matching_method, reference_context=reference_context)
    nearest_pairs = build_matched_peak_pairs(reference_peaks, assignment_audit, method="nearest", reference_context=reference_context)
    high_confidence_pairs = build_matched_peak_pairs(reference_peaks, assignment_audit, method="scored_high_confidence", reference_context=reference_context)
    extended_pairs = build_matched_peak_pairs(reference_peaks, assignment_audit, method="scored_extended", reference_context=reference_context)
    engine_table = compare_scale_engines_on_matched_peaks(
        matched_pairs,
        piecewise_regions=piecewise_regions,
    )
    nearest_engine_table = compare_scale_engines_on_matched_peaks(
        nearest_pairs,
        piecewise_regions=piecewise_regions,
    )
    high_conf_engine_table = compare_scale_engines_on_matched_peaks(
        high_confidence_pairs,
        piecewise_regions=piecewise_regions,
    )
    extended_engine_table = compare_scale_engines_on_matched_peaks(
        extended_pairs,
        piecewise_regions=piecewise_regions,
    )
    engine_fits: dict[str, dict] = {}
    for _, row in engine_table.iterrows():
        engine_fits[str(row["engine"])] = {
            "engine": str(row["engine"]),
            "parameters": json.loads(str(row["parameters_json"])),
            "metrics": {
                "mean_percent_deviation": float(row["mean_percent_deviation"]),
                "rmse_percent_deviation": float(row["rmse_percent_deviation"]),
                "max_percent_deviation": float(row["max_percent_deviation"]),
            },
            "matched_count": int(row["matched_count"]),
        }

    def _table_to_fits(table: pd.DataFrame) -> dict[str, dict]:
        fits: dict[str, dict] = {}
        for _, row in table.iterrows():
            fits[str(row["engine"])] = {
                "engine": str(row["engine"]),
                "parameters": json.loads(str(row["parameters_json"])),
                "metrics": {
                    "mean_percent_deviation": float(row["mean_percent_deviation"]),
                    "rmse_percent_deviation": float(row["rmse_percent_deviation"]),
                    "max_percent_deviation": float(row["max_percent_deviation"]),
                },
                "matched_count": int(row["matched_count"]),
            }
        return fits

    def _layer_overview(label: str, pairs: pd.DataFrame, total_reference_peaks: int) -> dict:
        matched_count = int(len(pairs))
        coverage = float(matched_count / total_reference_peaks) if total_reference_peaks > 0 else 0.0
        if pairs.empty:
            mean_percent = float("nan")
            rmse_percent = float("nan")
        else:
            abs_percent = (
                (
                    pairs["calc_frequency_cm-1"].astype(float)
                    - pairs["reference_peak_cm-1"].astype(float)
                ).abs()
                / pairs["reference_peak_cm-1"].astype(float).abs().clip(lower=1e-9)
            ) * 100.0
            mean_percent = float(abs_percent.mean()) if not abs_percent.empty else float("nan")
            rmse_percent = float(np.sqrt(np.mean(abs_percent**2))) if not abs_percent.empty else float("nan")
        return {
            "layer": label,
            "matched_count": matched_count,
            "total_reference_peaks": int(total_reference_peaks),
            "coverage": coverage,
            "mean_percent_deviation": mean_percent,
            "rmse_percent_deviation": rmse_percent,
        }

    def _engine_layer_matrix(nearest_table: pd.DataFrame, high_table: pd.DataFrame, ext_table: pd.DataFrame) -> list[dict]:
        by_layer = {
            "nearest": nearest_table,
            "high_confidence": high_table,
            "extended": ext_table,
        }
        engine_names: set[str] = set()
        for table in by_layer.values():
            if not table.empty:
                engine_names.update(str(name) for name in table["engine"].tolist())
        rows: list[dict] = []
        for engine in sorted(engine_names):
            row: dict[str, object] = {"engine": engine}
            for layer_name, table in by_layer.items():
                sub = table[table["engine"].astype(str) == engine]
                if sub.empty:
                    row[f"{layer_name}_mean_percent_deviation"] = float("nan")
                    row[f"{layer_name}_matched_count"] = 0
                else:
                    hit = sub.iloc[0]
                    row[f"{layer_name}_mean_percent_deviation"] = float(hit["mean_percent_deviation"])
                    row[f"{layer_name}_matched_count"] = int(hit["matched_count"])
            rows.append(row)
        return rows

    default_manual_scale = fit_constant_scale(
        matched_pairs["calc_frequency_cm-1"].astype(float).to_numpy(),
        matched_pairs["reference_peak_cm-1"].astype(float).to_numpy(),
        weights=matched_pairs["reference_intensity"].astype(float).to_numpy(),
    ) if not matched_pairs.empty else 1.0

    total_reference_peaks = int(len(reference_peaks))
    return {
        "matched_pairs": matched_pairs.to_dict(orient="records"),
        "nearest_matched_pairs": nearest_pairs.to_dict(orient="records"),
        "high_confidence_matched_pairs": high_confidence_pairs.to_dict(orient="records"),
        "extended_matched_pairs": extended_pairs.to_dict(orient="records"),
        "engine_table": engine_table.to_dict(orient="records"),
        "engine_fits": engine_fits,
        "matching_layer_overview": [
            _layer_overview("nearest", nearest_pairs, total_reference_peaks),
            _layer_overview("high_confidence", high_confidence_pairs, total_reference_peaks),
            _layer_overview("extended", extended_pairs, total_reference_peaks),
        ],
        "engine_layer_matrix": _engine_layer_matrix(
            nearest_engine_table,
            high_conf_engine_table,
            extended_engine_table,
        ),
        "matching_layers": {
            "nearest": {
                "matched_pairs": nearest_pairs.to_dict(orient="records"),
                "engine_table": nearest_engine_table.to_dict(orient="records"),
                "engine_fits": _table_to_fits(nearest_engine_table),
                "matched_count": int(len(nearest_pairs)),
                "total_reference_peaks": total_reference_peaks,
            },
            "high_confidence": {
                "matched_pairs": high_confidence_pairs.to_dict(orient="records"),
                "engine_table": high_conf_engine_table.to_dict(orient="records"),
                "engine_fits": _table_to_fits(high_conf_engine_table),
                "matched_count": int(len(high_confidence_pairs)),
                "total_reference_peaks": total_reference_peaks,
            },
            "extended": {
                "matched_pairs": extended_pairs.to_dict(orient="records"),
                "engine_table": extended_engine_table.to_dict(orient="records"),
                "engine_fits": _table_to_fits(extended_engine_table),
                "matched_count": int(len(extended_pairs)),
                "total_reference_peaks": total_reference_peaks,
            },
        },
        "default_manual_scale": float(default_manual_scale),
    }
