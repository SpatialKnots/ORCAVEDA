from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir.jcamp_parser import parse_jcamp_text  # noqa: E402


SAMPLE_JDX = """##TITLE=Example
##JCAMP-DX=4.24
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##YFACTOR=0.1
##DELTAX=4.0
##FIRSTX=450.0
##LASTX=486.0
##XYDATA=(X++(Y..Y))
450.0 10 20 30 40 50
470.0 60 70 80 90 100
##END=
"""


def test_parse_jcamp_text_nist_xydata():
    meta, spectrum = parse_jcamp_text(SAMPLE_JDX)
    assert meta["TITLE"] == "Example"
    assert meta["YUNITS"] == "ABSORBANCE"
    assert len(spectrum) == 10
    assert list(spectrum.columns) == ["wavenumber_cm-1", "intensity"]
    assert spectrum.iloc[0]["wavenumber_cm-1"] == 450.0
    assert spectrum.iloc[0]["intensity"] == 1.0
    assert spectrum.iloc[-1]["wavenumber_cm-1"] == 486.0
    assert spectrum.iloc[-1]["intensity"] == 10.0
