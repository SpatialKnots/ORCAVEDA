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

    run_orca_ped_like(hess_paths, args.outdir)


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
