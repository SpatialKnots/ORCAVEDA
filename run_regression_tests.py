#!/usr/bin/env python3
"""
ORCAVEDA Stage 3D v5.0 formal regression runner.

This runner validates already-generated ORCAVEDA CSV outputs against a JSON
expectation file. It is intentionally output-based and does not modify chemistry
assignment logic.

Required output files in --outdir:
  *__general_summary.csv
  *__functional_groups.csv
  *__assignment_audit.csv
  *__sanity_check.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Any

import pandas as pd


def _find_one(outdir: Path, suffix: str) -> Path:
    matches = sorted(outdir.glob(f"*__{suffix}.csv"))
    if not matches:
        raise FileNotFoundError(f"Missing ORCAVEDA output: *{suffix}.csv in {outdir}")
    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous ORCAVEDA outputs for {suffix}: {matches}")
    return matches[0]


def _status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def _safe_str(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x)


def _bool_series_false_or_false_text(s: pd.Series) -> pd.Series:
    # Handles boolean False, string "False", 0, NaN.
    return s.apply(lambda v: (v is False) or str(v).strip().lower() in {"false", "0", "", "nan"})


def run_checks(outdir: Path, expectations: Dict[str, Any]) -> pd.DataFrame:
    summary = pd.read_csv(_find_one(outdir, "general_summary"))
    fgroups = pd.read_csv(_find_one(outdir, "functional_groups"))
    audit = pd.read_csv(_find_one(outdir, "assignment_audit"))
    sanity = pd.read_csv(_find_one(outdir, "sanity_check"))

    rows: List[Dict[str, Any]] = []
    molecules = expectations.get("molecules", {})

    def add(filename: str, check: str, ok: bool, detail: Any = ""):
        rows.append({
            "Filename": filename,
            "check": check,
            "status": _status(bool(ok)),
            "detail": detail,
        })

    # Global required molecules present.
    observed = set(summary["Filename"].astype(str))
    expected = set(molecules)
    add("GLOBAL", "expected_molecules_present", expected.issubset(observed),
        f"missing={sorted(expected - observed)}; observed={sorted(observed)}")

    # Rank and parser-level summary checks.
    for _, srow in summary.iterrows():
        fn = str(srow["Filename"])
        if fn not in molecules:
            continue
        add(fn, "rank_B_independent_equals_expected_rank_3N_minus_6",
            int(srow["rank_B_independent"]) == int(srow["expected_rank_3N_minus_6"]),
            f"{srow['rank_B_independent']} vs {srow['expected_rank_3N_minus_6']}")
        add(fn, "normal_mode_orientation_rule_reported",
            "normal_modes[:, mode]" in _safe_str(srow.get("normal_mode_orientation_rule", "")),
            _safe_str(srow.get("normal_mode_orientation_rule", "")))
        add(fn, "no_negative_freq_after_first_6_unless_expected",
            int(srow.get("negative_freq_count_after_first_6", 0)) == 0
            or "negative_vibrational_frequency_after_first_6" in molecules[fn].get("expected_warning_flags", []),
            srow.get("negative_freq_count_after_first_6", ""))

    # Molecule-specific checks.
    for fn, exp in molecules.items():
        mol_fg = fgroups[fgroups["Filename"].astype(str) == fn]
        mol_audit = audit[audit["Filename"].astype(str) == fn]
        mol_sanity = sanity[sanity["Filename"].astype(str) == fn]

        fg_set = set(mol_fg["group"].astype(str)) if not mol_fg.empty else set()
        for g in exp.get("must_detect_functional_groups", []):
            add(fn, f"must_detect_functional_group:{g}", g in fg_set, sorted(fg_set))

        assignments = mol_audit["functional_group_assignment"].fillna("").astype(str).tolist()
        assignment_text = " | ".join(assignments)

        for target in exp.get("must_assign_mid_frequency", []):
            ok = any(target in a for a in assignments)
            hits = mol_audit[mol_audit["functional_group_assignment"].fillna("").astype(str).str.contains(target, regex=False)]
            detail = ""
            if not hits.empty:
                best = hits.sort_values("top1_percent", ascending=False).iloc[0]
                detail = f"mode={best['mode']}; freq={best['frequency_cm-1']}; assignment={best['functional_group_assignment']}; top1={best['top1_coord']}; top1_percent={best['top1_percent']}"
            add(fn, f"must_assign_mid_frequency:{target}", ok, detail or "not found")

        for forbidden in exp.get("forbidden_assignments", []):
            ok = forbidden not in assignment_text
            add(fn, f"forbidden_assignment_absent:{forbidden}", ok,
                "present" if not ok else "absent")

        # High-frequency unassigned X-H guard.
        high = mol_audit[pd.to_numeric(mol_audit["frequency_cm-1"], errors="coerce") >= 2800.0]
        high_unassigned = high[high["functional_group_assignment"].fillna("").astype(str).str.lower().eq("unassigned")]
        add(fn, "high_frequency_unassigned_count_is_zero", len(high_unassigned) == 0, len(high_unassigned))

        # Monoethanolamine-specific sanity checks should be skipped for this external set.
        if not mol_sanity.empty and "verdict" in mol_sanity.columns:
            verdicts = set(mol_sanity["verdict"].fillna("").astype(str))
            add(fn, "monoethanolamine_sanity_skipped", verdicts == {"SKIPPED"}, sorted(verdicts))

        # v4.9/v5.0 protected cleanup.
        cols = set(mol_audit.columns)
        if {"protected_xh_used", "protected_xh_normal_mode_norm", "protected_xh_total_power"}.issubset(cols):
            unused = _bool_series_false_or_false_text(mol_audit["protected_xh_used"])
            bad = mol_audit.loc[unused, ["protected_xh_normal_mode_norm", "protected_xh_total_power"]].notna().any(axis=1).sum()
            add(fn, "protected_xh_unused_norm_power_blank", int(bad) == 0, int(bad))

    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="Directory containing ORCAVEDA CSV outputs")
    ap.add_argument("--expectations", required=True, help="Regression expectations JSON")
    ap.add_argument("--output-csv", default="regression_harness_results_stage3D_v5_0.csv")
    ap.add_argument("--output-json", default="regression_harness_summary_stage3D_v5_0.json")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    expectations = json.loads(Path(args.expectations).read_text(encoding="utf-8"))
    results = run_checks(outdir, expectations)

    output_csv = Path(args.output_csv)
    if not output_csv.is_absolute():
        output_csv = outdir / output_csv
    results.to_csv(output_csv, index=False)

    summary = {
        "version": expectations.get("version", "unknown"),
        "rows": int(len(results)),
        "pass": int((results["status"] == "PASS").sum()),
        "fail": int((results["status"] == "FAIL").sum()),
        "status": "PASS" if (results["status"] != "FAIL").all() else "FAIL",
    }
    output_json = Path(args.output_json)
    if not output_json.is_absolute():
        output_json = outdir / output_json
    output_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(results.to_string(index=False))
    print(f"\nWrote: {output_csv}")
    print(f"Wrote: {output_json}")
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
