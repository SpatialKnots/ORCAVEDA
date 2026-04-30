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


def pick_reference_peaks(
    spectrum: pd.DataFrame,
    *,
    top_n: int = 12,
    min_intensity: float | None = None,
    min_separation_cm1: float = 20.0,
) -> pd.DataFrame:
    if spectrum.empty:
        return pd.DataFrame(columns=["wavenumber_cm-1", "intensity"])

    df = spectrum.copy()
    df = df.sort_values("wavenumber_cm-1").reset_index(drop=True)
    peaks: List[Dict[str, float]] = []
    values = df["intensity"].tolist()
    wns = df["wavenumber_cm-1"].tolist()

    for i in range(1, len(df) - 1):
        y = values[i]
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
        if all(abs(wn - item["wavenumber_cm-1"]) > min_separation_cm1 for item in selected):
            selected.append({"wavenumber_cm-1": wn, "intensity": float(row["intensity"])})
        if len(selected) >= top_n:
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
) -> pd.DataFrame:
    if method == "scored":
        matched = match_reference_to_orcaveda_v2(reference_peaks, assignment_audit)
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
                except Exception:
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
) -> dict:
    matched_pairs = build_matched_peak_pairs(reference_peaks, assignment_audit, method=matching_method)
    engine_table = compare_scale_engines_on_matched_peaks(
        matched_pairs,
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

    default_manual_scale = fit_constant_scale(
        matched_pairs["calc_frequency_cm-1"].astype(float).to_numpy(),
        matched_pairs["reference_peak_cm-1"].astype(float).to_numpy(),
        weights=matched_pairs["reference_intensity"].astype(float).to_numpy(),
    ) if not matched_pairs.empty else 1.0

    return {
        "matched_pairs": matched_pairs.to_dict(orient="records"),
        "engine_table": engine_table.to_dict(orient="records"),
        "engine_fits": engine_fits,
        "default_manual_scale": float(default_manual_scale),
    }
