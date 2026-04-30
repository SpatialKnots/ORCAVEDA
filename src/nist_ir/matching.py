from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Iterable, Sequence

import pandas as pd


@dataclass(frozen=True)
class ExperimentalPeak:
    peak_id: str
    frequency_cm1: float
    intensity: float
    prominence: float | None = None
    region: str | None = None


@dataclass(frozen=True)
class CalculatedMode:
    mode: int
    frequency_cm1: float
    intensity: float
    assignment: str
    assignment_class: str
    warnings: str = ""
    region: str | None = None


@dataclass(frozen=True)
class MatchCandidate:
    peak_id: str
    mode: int
    delta_cm1: float
    freq_penalty: float
    class_penalty: float
    intensity_penalty: float
    warning_penalty: float
    total_cost: float


@dataclass(frozen=True)
class MatchedPair:
    peak_id: str
    mode: int
    experimental_frequency_cm1: float
    calculated_frequency_cm1: float
    delta_cm1: float
    experimental_intensity: float
    calculated_intensity: float
    assignment: str
    assignment_class: str
    region: str
    total_cost: float
    confidence: str


DEFAULT_TOLERANCE_BY_REGION = {
    "low": 35.0,
    "fingerprint": 20.0,
    "double_bond": 18.0,
    "xh_stretch": 30.0,
}

DEFAULT_WEIGHTS = {
    "freq": 1.0,
    "class": 0.6,
    "intensity": 0.2,
    "warning": 0.2,
}


def infer_region(frequency_cm1: float) -> str:
    freq = float(frequency_cm1)
    if freq < 800.0:
        return "low"
    if freq < 1500.0:
        return "fingerprint"
    if freq < 1900.0:
        return "double_bond"
    if freq >= 2800.0:
        return "xh_stretch"
    return "fingerprint"


def infer_assignment_class(assignment: str) -> str:
    text = str(assignment or "").strip().lower()
    if not text:
        return "generic"
    if "intermolecular" in text or "h-bond" in text or "cluster" in text:
        return "intermolecular"
    if "c=o" in text or "carbonyl" in text:
        return "carbonyl"
    if "o-h" in text and "stretch" in text:
        return "oh_stretch"
    if "n-h" in text and "stretch" in text:
        return "nh_stretch"
    if "c-h" in text and "stretch" in text:
        return "ch_stretch"
    if "c-o" in text and "stretch" in text:
        return "co_stretch"
    if "c-n" in text and "stretch" in text:
        return "cn_stretch"
    if "aromatic c-h bend" in text:
        return "aromatic_ch_bend"
    if "aromatic c-h" in text:
        return "aromatic_ch"
    if "torsion" in text:
        return "torsion"
    if "bend" in text:
        return "bend"
    if "stretch" in text:
        return "stretch"
    return "generic"


def build_experimental_peaks(reference_peaks: pd.DataFrame) -> list[ExperimentalPeak]:
    peaks: list[ExperimentalPeak] = []
    if reference_peaks.empty:
        return peaks
    for idx, row in reference_peaks.reset_index(drop=True).iterrows():
        freq = float(row["wavenumber_cm-1"])
        peaks.append(
            ExperimentalPeak(
                peak_id=f"exp_{idx+1}",
                frequency_cm1=freq,
                intensity=float(row.get("intensity", 0.0)),
                prominence=float(row["intensity"]) if pd.notna(row.get("intensity")) else None,
                region=infer_region(freq),
            )
        )
    return peaks


def build_calculated_modes(assignment_audit: pd.DataFrame) -> list[CalculatedMode]:
    modes: list[CalculatedMode] = []
    if assignment_audit.empty:
        return modes
    for _, row in assignment_audit.iterrows():
        freq = float(row["scaled_frequency_cm-1"])
        assignment = str(row.get("functional_group_assignment", ""))
        modes.append(
            CalculatedMode(
                mode=int(row["mode"]),
                frequency_cm1=freq,
                intensity=float(row.get("IR_intensity", 0.0)),
                assignment=assignment,
                assignment_class=infer_assignment_class(assignment),
                warnings=str(row.get("warnings", "")),
                region=infer_region(freq),
            )
        )
    return modes


def region_tolerance(region: str, tolerance_by_region: dict[str, float] | None = None) -> float:
    mapping = dict(DEFAULT_TOLERANCE_BY_REGION)
    if tolerance_by_region:
        mapping.update({str(k): float(v) for k, v in tolerance_by_region.items()})
    return float(mapping.get(str(region), 20.0))


def _experimental_expected_classes(peak: ExperimentalPeak) -> set[str]:
    if peak.region == "double_bond":
        return {"carbonyl", "co_stretch", "cn_stretch", "generic"}
    if peak.region == "xh_stretch":
        return {"oh_stretch", "nh_stretch", "ch_stretch", "stretch", "generic"}
    if peak.region == "low":
        return {"torsion", "intermolecular", "bend", "generic"}
    return {"co_stretch", "cn_stretch", "aromatic_ch_bend", "bend", "stretch", "generic"}


def score_candidate(
    exp_peak: ExperimentalPeak,
    calc_mode: CalculatedMode,
    *,
    tolerance_by_region: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> MatchCandidate | None:
    tol = region_tolerance(exp_peak.region or "fingerprint", tolerance_by_region)
    delta = float(calc_mode.frequency_cm1 - exp_peak.frequency_cm1)
    abs_delta = abs(delta)
    if abs_delta > tol:
        return None

    w = dict(DEFAULT_WEIGHTS)
    if weights:
        w.update({str(k): float(v) for k, v in weights.items()})

    freq_penalty = abs_delta / max(tol, 1e-9)
    expected = _experimental_expected_classes(exp_peak)
    if calc_mode.assignment_class in expected:
        class_penalty = 0.0
    elif calc_mode.assignment_class in {"generic", "stretch", "bend"}:
        class_penalty = 0.25
    else:
        class_penalty = 0.75

    exp_int = max(float(exp_peak.intensity), 1e-9)
    calc_int = max(float(calc_mode.intensity), 1e-9)
    ratio = max(exp_int, calc_int) / min(exp_int, calc_int)
    intensity_penalty = min(1.0, abs(math.log10(ratio)) / 3.0)

    warn_text = calc_mode.warnings.lower()
    if "mixed" in calc_mode.assignment.lower():
        warning_penalty = 0.15
    elif warn_text and warn_text != "none":
        warning_penalty = 0.25
    else:
        warning_penalty = 0.0

    total_cost = (
        w["freq"] * freq_penalty
        + w["class"] * class_penalty
        + w["intensity"] * intensity_penalty
        + w["warning"] * warning_penalty
    )
    return MatchCandidate(
        peak_id=exp_peak.peak_id,
        mode=calc_mode.mode,
        delta_cm1=delta,
        freq_penalty=freq_penalty,
        class_penalty=class_penalty,
        intensity_penalty=intensity_penalty,
        warning_penalty=warning_penalty,
        total_cost=total_cost,
    )


def generate_match_candidates(
    exp_peaks: Sequence[ExperimentalPeak],
    calc_modes: Sequence[CalculatedMode],
    *,
    tolerance_by_region: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> list[MatchCandidate]:
    candidates: list[MatchCandidate] = []
    for peak in exp_peaks:
        for mode in calc_modes:
            candidate = score_candidate(
                peak,
                mode,
                tolerance_by_region=tolerance_by_region,
                weights=weights,
            )
            if candidate is not None:
                candidates.append(candidate)
    return candidates


def assign_match_confidence(total_cost: float, abs_delta_cm1: float, class_penalty: float) -> str:
    if abs_delta_cm1 <= 10.0 and total_cost <= 0.6 and class_penalty <= 0.25:
        return "high"
    if abs_delta_cm1 <= 20.0 and total_cost <= 1.1:
        return "medium"
    return "low"


def solve_peak_matching(
    candidates: Sequence[MatchCandidate],
    exp_peaks: Sequence[ExperimentalPeak],
    calc_modes: Sequence[CalculatedMode],
) -> dict[str, object]:
    peak_by_id = {peak.peak_id: peak for peak in exp_peaks}
    mode_by_id = {mode.mode: mode for mode in calc_modes}
    used_peaks: set[str] = set()
    used_modes: set[int] = set()
    matched_pairs: list[MatchedPair] = []

    for candidate in sorted(candidates, key=lambda item: (item.total_cost, abs(item.delta_cm1), item.mode, item.peak_id)):
        if candidate.peak_id in used_peaks or candidate.mode in used_modes:
            continue
        peak = peak_by_id[candidate.peak_id]
        mode = mode_by_id[candidate.mode]
        used_peaks.add(candidate.peak_id)
        used_modes.add(candidate.mode)
        matched_pairs.append(
            MatchedPair(
                peak_id=peak.peak_id,
                mode=mode.mode,
                experimental_frequency_cm1=peak.frequency_cm1,
                calculated_frequency_cm1=mode.frequency_cm1,
                delta_cm1=candidate.delta_cm1,
                experimental_intensity=peak.intensity,
                calculated_intensity=mode.intensity,
                assignment=mode.assignment,
                assignment_class=mode.assignment_class,
                region=peak.region or mode.region or "fingerprint",
                total_cost=candidate.total_cost,
                confidence=assign_match_confidence(candidate.total_cost, abs(candidate.delta_cm1), candidate.class_penalty),
            )
        )

    unmatched_exp = [asdict(peak) for peak in exp_peaks if peak.peak_id not in used_peaks]
    unmatched_calc = [asdict(mode) for mode in calc_modes if mode.mode not in used_modes]
    return {
        "matched_pairs": matched_pairs,
        "unmatched_experimental": unmatched_exp,
        "unmatched_calculated": unmatched_calc,
    }


def match_reference_to_orcaveda_v2(
    reference_peaks: pd.DataFrame,
    assignment_audit: pd.DataFrame,
    *,
    tolerance_by_region: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    exp_peaks = build_experimental_peaks(reference_peaks)
    calc_modes = build_calculated_modes(assignment_audit)
    candidates = generate_match_candidates(
        exp_peaks,
        calc_modes,
        tolerance_by_region=tolerance_by_region,
        weights=weights,
    )
    solved = solve_peak_matching(candidates, exp_peaks, calc_modes)
    matched_pairs: Iterable[MatchedPair] = solved["matched_pairs"]
    rows = [
        {
            "peak_id": pair.peak_id,
            "reference_peak_cm-1": pair.experimental_frequency_cm1,
            "reference_intensity": pair.experimental_intensity,
            "orcaveda_mode": pair.mode,
            "orcaveda_scaled_frequency_cm-1": pair.calculated_frequency_cm1,
            "delta_cm-1": pair.delta_cm1,
            "orcaveda_assignment": pair.assignment,
            "orcaveda_assignment_class": pair.assignment_class,
            "orcaveda_intensity": pair.calculated_intensity,
            "region": pair.region,
            "total_cost": pair.total_cost,
            "match_confidence": pair.confidence,
        }
        for pair in matched_pairs
    ]
    return pd.DataFrame(rows).sort_values("reference_peak_cm-1", ascending=False).reset_index(drop=True) if rows else pd.DataFrame(
        columns=[
            "peak_id",
            "reference_peak_cm-1",
            "reference_intensity",
            "orcaveda_mode",
            "orcaveda_scaled_frequency_cm-1",
            "delta_cm-1",
            "orcaveda_assignment",
            "orcaveda_assignment_class",
            "orcaveda_intensity",
            "region",
            "total_cost",
            "match_confidence",
        ]
    )
