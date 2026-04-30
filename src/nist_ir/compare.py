from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd


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
