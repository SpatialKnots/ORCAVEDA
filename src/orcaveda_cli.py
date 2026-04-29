from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Sequence


def cli_main(run_orca_ped_like: Callable[[Sequence[str], str], object]) -> None:
    parser = argparse.ArgumentParser(description="ORCAVEDA: ORCA-native vibrational/PED-like analyzer")
    parser.add_argument("hess", nargs="*", help="Input ORCA .hess files")
    parser.add_argument("--outdir", default="orca_ped_results", help="Output directory")

    args, _unknown = parser.parse_known_args()
    hess_paths = [str(Path(x)) for x in args.hess if str(x).lower().endswith(".hess")]

    if not hess_paths:
        print("No .hess files provided via command line.")
        print("Terminal example: python ORCAVEDA.py file1.hess file2.hess --outdir results")
        print('Colab example: run_orca_ped_like(["/content/file1.hess"], "orca_ped_results")')
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
