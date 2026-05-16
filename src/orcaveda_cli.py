from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Callable, Optional, Sequence


def cli_main(
    run_orca_ped_like: Callable[[Sequence[str], str], object],
    *,
    set_chem_backend: Optional[Callable[[str], object]] = None,
    list_chem_backends: Optional[Callable[[], Sequence[str]]] = None,
    default_chem_backend: Optional[Callable[[], str]] = None,
    chem_backend_env_var: str = "ORCAVEDA_CHEM_BACKEND",
) -> None:
    parser = argparse.ArgumentParser(description="ORCAVEDA: ORCA-native vibrational/PED-like analyzer")
    parser.add_argument("hess", nargs="*", help="Input ORCA .hess files")
    parser.add_argument("--outdir", default="orca_ped_results", help="Output directory")
    parser.add_argument("--chem-backend", help=f"Chemistry backend name. Overrides {chem_backend_env_var} if provided.")
    parser.add_argument("--list-chem-backends", action="store_true", help="Print available chemistry backends and exit.")
    parser.add_argument(
        "--experimental-composed-primitive-substitution-constraint",
        action="store_true",
        help=(
            "Opt-in experiment: repair composed PED primitive-row optimizer substitution warnings "
            "without changing assignment_audit, ped_stage3d_agreement, or ped_final_assignment policy."
        ),
    )
    parser.add_argument(
        "--wilson-gf-validation",
        action="store_true",
        help=(
            "Opt-in diagnostic prototype: emit Wilson GF diagonalization validation CSVs. "
            "Does not change default Stage 3D or PED outputs."
        ),
    )
    parser.add_argument(
        "--veda-like-ped",
        action="store_true",
        help=(
            "Opt-in comparable VEDA-like closed Wilson GF/PED CSVs. "
            "Diagnostic only; does not claim original VEDA reproduction."
        ),
    )
    parser.add_argument(
        "--epm-optimize",
        action="store_true",
        help=(
            "Opt-in EPM-like Wilson GF/PED basis optimization for Wilson GF validation and VEDA-like diagnostics. "
            "Does not change default Stage 3D assignment audit labels."
        ),
    )
    parser.add_argument("--epm-max-passes", type=int, default=2, help="Maximum greedy EPM optimization passes.")
    parser.add_argument(
        "--epm-improvement-tol",
        type=float,
        default=1.0e-6,
        help="Minimum Wilson GF/PED localization-score improvement required for an EPM swap.",
    )
    parser.add_argument(
        "--b-matrix-method",
        choices=("finite_difference", "hybrid_analytical"),
        default="finite_difference",
        help=(
            "B-matrix construction method. Default finite_difference preserves existing outputs; "
            "hybrid_analytical is opt-in and uses analytical distance/regular angle rows with finite-difference fallback."
        ),
    )

    args, _unknown = parser.parse_known_args()
    if args.list_chem_backends:
        if list_chem_backends:
            print(f"Available chemistry backends: {', '.join(list_chem_backends())}")
            if default_chem_backend:
                print(f"Active default chemistry backend: {default_chem_backend()}")
        else:
            print("Chemistry backend registry is not available in this entry point.")
        return

    selected_backend = args.chem_backend or os.environ.get(chem_backend_env_var, "").strip()
    if selected_backend:
        if not set_chem_backend:
            parser.error("Chemistry backend selection is not supported by this entry point.")
        try:
            set_chem_backend(selected_backend)
        except Exception as exc:
            parser.error(str(exc))

    hess_paths = [str(Path(x)) for x in args.hess if str(x).lower().endswith(".hess")]

    if not hess_paths:
        print("No .hess files provided via command line.")
        print("Terminal example: python ORCAVEDA.py file1.hess file2.hess --outdir results")
        print('Colab example: run_orca_ped_like(["/content/file1.hess"], "orca_ped_results")')
        if list_chem_backends and default_chem_backend:
            print(f"Available chemistry backends: {', '.join(list_chem_backends())}")
            print(f"Active default chemistry backend: {default_chem_backend()}")
        return

    run_orca_ped_like(
        hess_paths,
        args.outdir,
        experimental_composed_primitive_substitution_constraint=args.experimental_composed_primitive_substitution_constraint,
        wilson_gf_validation=args.wilson_gf_validation,
        veda_like_ped=args.veda_like_ped,
        epm_optimize=args.epm_optimize,
        epm_max_passes=args.epm_max_passes,
        epm_improvement_tol=args.epm_improvement_tol,
        b_matrix_method=args.b_matrix_method,
    )


def colab_upload_and_run(run_orca_ped_like: Callable[[Sequence[str], str], object]) -> None:
    from google.colab import files

    print("Upload one or more ORCA .hess files")
    uploaded = files.upload()

    hess_paths = []
    for filename in uploaded.keys():
        if filename.lower().endswith(".hess"):
            hess_paths.append(str(Path(filename).resolve()))

    if not hess_paths:
        print("No .hess files uploaded.")
        return

    print("\nDetected .hess files:")
    for path in hess_paths:
        print(" -", path)

    outdir = "orca_ped_results"
    run_orca_ped_like(hess_paths, outdir)
    print(f"\nDone. Results saved to: {outdir}")


def is_google_colab() -> bool:
    try:
        import google.colab  # noqa: F401
        return True
    except ModuleNotFoundError:
        return False
