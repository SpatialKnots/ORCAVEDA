from __future__ import annotations

import argparse
import json
from pathlib import Path

from .compare import load_orcaveda_assignments, match_reference_to_orcaveda, pick_reference_peaks
from .pipeline import nist_ir_from_hess, nist_ir_from_smiles


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch NIST IR JCAMP-DX data from a SMILES string or directly from a .hess structure.")
    parser.add_argument("smiles", nargs="?", help="Input SMILES")
    parser.add_argument("--hess", help="Optional ORCA .hess path; RDKit will infer canonical SMILES/InChI/InChIKey from geometry")
    parser.add_argument("--outdir", default="outputs/nist_ir", help="Output directory")
    parser.add_argument("--nist-id", help="Optional explicit NIST species identifier, e.g. C98862")
    parser.add_argument("--charge", type=int, default=0, help="Formal charge to assume when inferring identifiers from .hess")
    parser.add_argument("--scale-factor", type=float, default=0.96, help="Scale factor for ORCAVEDA comparison")
    parser.add_argument("--compare-audit", help="Optional ORCAVEDA assignment audit CSV to compare against NIST peaks")
    parser.add_argument("--top-peaks", type=int, default=12, help="How many NIST peaks to keep in comparison")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if bool(args.smiles) == bool(args.hess):
        parser.error("Provide exactly one of: positional SMILES or --hess <path>")

    if args.hess:
        results = nist_ir_from_hess(args.hess, outdir, charge=args.charge, nist_id=args.nist_id)
    else:
        results = nist_ir_from_smiles(args.smiles, outdir, nist_id=args.nist_id)
    print(json.dumps(results, ensure_ascii=False, indent=2))

    if not args.compare_audit:
        return

    first_csv = Path(results[0]["csv"])
    reference = pick_reference_peaks(
        __import__("pandas").read_csv(first_csv),
        top_n=args.top_peaks,
    )
    audit = load_orcaveda_assignments(args.compare_audit, scale_factor=args.scale_factor)
    comparison = match_reference_to_orcaveda(reference, audit)
    comparison_path = outdir / "nist_vs_orcaveda_comparison.csv"
    comparison.to_csv(comparison_path, index=False, encoding="utf-8")
    print(f"Comparison saved to: {comparison_path}")


if __name__ == "__main__":
    main()
