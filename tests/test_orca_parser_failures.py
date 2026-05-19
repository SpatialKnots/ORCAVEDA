from __future__ import annotations

import sys
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from orca_parser import read_orca_hess  # noqa: E402


def test_read_orca_hess_reports_missing_required_sections():
    outdir = ROOT / "outputs" / "pytest_orca_parser_failures"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    broken_hess = outdir / "broken.hess"
    broken_hess.write_text("$atoms\n1\nH 1.0 0.0 0.0 0.0\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"Missing required \.hess sections"):
        read_orca_hess(broken_hess)


def test_read_orca_hess_reports_decode_failures():
    outdir = ROOT / "outputs" / "pytest_orca_parser_failures"
    outdir.mkdir(parents=True, exist_ok=True)
    broken_hess = outdir / "bad_encoding.hess"
    broken_hess.write_bytes(b"$atoms\n1\nH 1.0 0.0 0.0 0.0\n\xff")

    with pytest.raises(ValueError, match=r"Failed to decode \.hess file bad_encoding\.hess as UTF-8"):
        read_orca_hess(broken_hess)
