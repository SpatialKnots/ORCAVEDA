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
    phase_tag: str | None = None
    state: str | None = None


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
    stage: str


DEFAULT_TOLERANCE_BY_REGION = {
    "low": 45.0,
    "fingerprint": 28.0,
    "double_bond": 24.0,
    "xh_stretch": 40.0,
}

SECONDARY_TOLERANCE_BY_REGION = {
    "low": 55.0,
    "fingerprint": 35.0,
    "double_bond": 28.0,
    "xh_stretch": 48.0,
}

BACKFILL_TOLERANCE_BY_REGION = {
    "low": 65.0,
    "fingerprint": 40.0,
    "double_bond": 32.0,
    "xh_stretch": 55.0,
}

AROMATIC_FINGERPRINT_CLASSES = {
    "aromatic_ch_bend",
    "aromatic_ch",
    "co_stretch",
    "cn_stretch",
    "bend",
    "stretch",
    "generic",
}

DOUBLE_BOND_CLASSES = {
    "carbonyl",
    "co_stretch",
    "cn_stretch",
    "stretch",
    "generic",
}

XH_STRETCH_CLASSES = {
    "oh_stretch",
    "nh_stretch",
    "ch_stretch",
    "aromatic_ch",
    "stretch",
    "generic",
}

DEFAULT_WEIGHTS = {
    "freq": 1.0,
    "class": 0.35,
    "intensity": 0.10,
    "warning": 0.10,
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
    if "ring" in text and "stretch" in text:
        return "aromatic_ring"
    if "ring" in text:
        return "aromatic_ring"
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


def build_experimental_peaks(reference_peaks: pd.DataFrame, *, reference_context: dict[str, object] | None = None) -> list[ExperimentalPeak]:
    peaks: list[ExperimentalPeak] = []
    if reference_peaks.empty:
        return peaks
    phase_tag = str((reference_context or {}).get("phase_tag", "") or "")
    state = str((reference_context or {}).get("state", "") or "")
    for idx, row in reference_peaks.reset_index(drop=True).iterrows():
        freq = float(row["wavenumber_cm-1"])
        peaks.append(
            ExperimentalPeak(
                peak_id=f"exp_{idx+1}",
                frequency_cm1=freq,
                intensity=float(row.get("intensity", 0.0)),
                prominence=float(row["intensity"]) if pd.notna(row.get("intensity")) else None,
                region=infer_region(freq),
                phase_tag=phase_tag,
                state=state,
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


def is_condensed_phase(peak: ExperimentalPeak) -> bool:
    phase = str(peak.phase_tag or peak.state or "").lower()
    return any(token in phase for token in ("liquid", "solution", "solid"))


def phase_scaled_tolerance(peak: ExperimentalPeak, base: float) -> float:
    if not is_condensed_phase(peak):
        return base
    region = peak.region or "fingerprint"
    factor = {
        "low": 1.35,
        "fingerprint": 1.28,
        "double_bond": 1.22,
        "xh_stretch": 1.45,
    }.get(region, 1.2)
    return base * factor


def _experimental_expected_classes(peak: ExperimentalPeak) -> set[str]:
    if peak.region == "double_bond":
        expected = {"carbonyl", "co_stretch", "cn_stretch", "stretch", "generic"}
        if is_condensed_phase(peak):
            expected |= {"bend"}
        return expected
    if peak.region == "xh_stretch":
        expected = {"oh_stretch", "nh_stretch", "ch_stretch", "aromatic_ch", "stretch", "generic"}
        if is_condensed_phase(peak):
            expected |= {"intermolecular", "bend"}
        return expected
    if peak.region == "low":
        expected = {"torsion", "intermolecular", "bend", "stretch", "generic"}
        if is_condensed_phase(peak):
            expected |= {"oh_stretch", "nh_stretch"}
        return expected
    expected = {"co_stretch", "cn_stretch", "aromatic_ch_bend", "aromatic_ch", "aromatic_ring", "bend", "stretch", "generic"}
    if is_condensed_phase(peak):
        expected |= {"intermolecular"}
    return expected


def score_candidate(
    exp_peak: ExperimentalPeak,
    calc_mode: CalculatedMode,
    *,
    tolerance_by_region: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> MatchCandidate | None:
    tol = phase_scaled_tolerance(exp_peak, region_tolerance(exp_peak.region or "fingerprint", tolerance_by_region))
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
    elif is_condensed_phase(exp_peak) and calc_mode.assignment_class in {"generic", "stretch", "bend", "intermolecular", "aromatic_ring"}:
        class_penalty = 0.08
    elif calc_mode.assignment_class in {"generic", "stretch", "bend"}:
        class_penalty = 0.12
    else:
        class_penalty = 0.35

    exp_int = max(float(exp_peak.intensity), 1e-9)
    calc_int = max(float(calc_mode.intensity), 1e-9)
    ratio = max(exp_int, calc_int) / min(exp_int, calc_int)
    intensity_penalty = min(1.0, abs(math.log10(ratio)) / 6.0)

    warn_text = calc_mode.warnings.lower()
    if "mixed" in calc_mode.assignment.lower():
        warning_penalty = 0.08
    elif warn_text and warn_text != "none":
        warning_penalty = 0.15
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


def _secondary_region_policy(peak: ExperimentalPeak, mode: CalculatedMode) -> tuple[float, float, set[str] | None]:
    region = peak.region or "fingerprint"
    if region == "fingerprint":
        # Aromatic fingerprint modes need a wider fallback window than generic bends,
        # but only for chemically plausible aromatic / heteroatom-coupled classes.
        tol = 38.0 if mode.assignment_class in AROMATIC_FINGERPRINT_CLASSES else 32.0
        if is_condensed_phase(peak):
            tol *= 1.2
        return (
            tol,
            0.50 if mode.assignment_class in AROMATIC_FINGERPRINT_CLASSES else (0.35 if is_condensed_phase(peak) else 0.28),
            AROMATIC_FINGERPRINT_CLASSES,
        )
    if region == "double_bond":
        tol = 26.0 * (1.2 if is_condensed_phase(peak) else 1.0)
        return (tol, 0.24 if is_condensed_phase(peak) else 0.18, DOUBLE_BOND_CLASSES)
    if region == "xh_stretch":
        tol = 48.0 * (1.25 if is_condensed_phase(peak) else 1.0)
        return (tol, 0.38 if is_condensed_phase(peak) else 0.28, XH_STRETCH_CLASSES | {"intermolecular"})
    return (phase_scaled_tolerance(peak, region_tolerance(region, SECONDARY_TOLERANCE_BY_REGION)), 0.35, None)


def _backfill_region_policy(peak: ExperimentalPeak, mode: CalculatedMode) -> tuple[float, set[str] | None]:
    region = peak.region or "fingerprint"
    if region == "fingerprint":
        tol = 42.0 if mode.assignment_class in AROMATIC_FINGERPRINT_CLASSES else 36.0
        if is_condensed_phase(peak):
            tol *= 1.2
        return (tol, AROMATIC_FINGERPRINT_CLASSES)
    if region == "double_bond":
        return (30.0 * (1.2 if is_condensed_phase(peak) else 1.0), DOUBLE_BOND_CLASSES)
    if region == "xh_stretch":
        return (52.0 * (1.25 if is_condensed_phase(peak) else 1.0), XH_STRETCH_CLASSES | {"intermolecular"})
    return (phase_scaled_tolerance(peak, region_tolerance(region, BACKFILL_TOLERANCE_BY_REGION)), None)


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

    primary_candidates = [
        item
        for item in candidates
        if item.total_cost <= 0.85 and item.class_penalty <= 0.12 and item.freq_penalty <= 0.85
    ]
    secondary_candidates = list(candidates)

    def accept(candidate: MatchCandidate, *, stage: str) -> bool:
        if candidate.peak_id in used_peaks or candidate.mode in used_modes:
            return False
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
                stage=stage,
            )
        )
        return True

    for candidate in sorted(primary_candidates, key=lambda item: (item.total_cost, abs(item.delta_cm1), item.mode, item.peak_id)):
        accept(candidate, stage="primary")

    def secondary_sort_key(candidate: MatchCandidate) -> tuple[float, float, int, str]:
        peak = peak_by_id[candidate.peak_id]
        mode = mode_by_id[candidate.mode]
        fallback_tol, _, _ = _secondary_region_policy(peak, mode)
        secondary_score = (
            abs(candidate.delta_cm1) / max(fallback_tol, 1e-9)
            + 0.5 * candidate.class_penalty
            + 0.08 * candidate.intensity_penalty
            + 0.05 * candidate.warning_penalty
        )
        return (secondary_score, abs(candidate.delta_cm1), mode.mode, peak.peak_id)

    for candidate in sorted(secondary_candidates, key=secondary_sort_key):
        if candidate.peak_id in used_peaks or candidate.mode in used_modes:
            continue
        peak = peak_by_id[candidate.peak_id]
        mode = mode_by_id[candidate.mode]
        fallback_tol, class_limit, allowed_classes = _secondary_region_policy(peak, mode)
        if abs(candidate.delta_cm1) > fallback_tol:
            continue
        if allowed_classes is not None and mode.assignment_class not in allowed_classes:
            continue
        if candidate.class_penalty > class_limit:
            continue
        accept(candidate, stage="secondary")

    # Final completion pass: region-aware nearest backfill for still-uncovered peaks.
    for peak in sorted((item for item in exp_peaks if item.peak_id not in used_peaks), key=lambda item: item.intensity, reverse=True):
        candidate_pool = [
            item
            for item in candidates
            if item.peak_id == peak.peak_id and item.mode not in used_modes
        ]
        filtered_pool: list[MatchCandidate] = []
        for item in candidate_pool:
            mode = mode_by_id[item.mode]
            backfill_tol, allowed_classes = _backfill_region_policy(peak, mode)
            if abs(item.delta_cm1) > backfill_tol:
                continue
            if allowed_classes is not None and mode.assignment_class not in allowed_classes:
                continue
            filtered_pool.append(item)
        candidate_pool = filtered_pool
        if not candidate_pool:
            continue
        best = min(
            candidate_pool,
            key=lambda item: (
                abs(item.delta_cm1),
                item.class_penalty,
                item.intensity_penalty,
                item.warning_penalty,
                item.mode,
            ),
        )
        accept(best, stage="backfill")

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
    include_stages: set[str] | None = None,
    reference_context: dict[str, object] | None = None,
) -> pd.DataFrame:
    exp_peaks = build_experimental_peaks(reference_peaks, reference_context=reference_context)
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
            "match_stage": pair.stage,
        }
        for pair in matched_pairs
        if include_stages is None or pair.stage in include_stages
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
            "match_stage",
        ]
    )
