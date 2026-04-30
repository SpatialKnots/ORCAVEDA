from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir.pipeline import nist_ir_from_smiles  # noqa: E402


SAMPLE_HTML = """
<html><body>
<a href="/cgi/inchi?ID=C98862&Type=IR-SPEC&Index=1#IR-SPEC">VAPOR (12 MICROLITER AT 150 C)</a>
<a href="/cgi/inchi?ID=C98862&Type=IR-SPEC&Index=0#IR-SPEC">gas</a>
</body></html>
"""

SAMPLE_JDX = """##TITLE=Acetophenone
##JCAMP-DX=4.24
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##YFACTOR=0.01
##DELTAX=2.0
##FIRSTX=1000.0
##LASTX=1008.0
##XYDATA=(X++(Y..Y))
1000.0 1 2 3 4 5
##END=
"""


def test_nist_ir_from_smiles_local_smoke():
    tmp_path = ROOT / "outputs" / "pytest_nist_ir_pipeline"
    if tmp_path.exists():
        shutil.rmtree(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)

    def fake_page_fetch(url: str) -> str:
        assert "InChI=" in url
        return SAMPLE_HTML

    def fake_jcamp_fetch(url: str) -> str:
        assert "JCAMP=C98862" in url
        assert "Type=IR" in url
        return SAMPLE_JDX

    results = nist_ir_from_smiles(
        "CC(=O)c1ccccc1",
        tmp_path,
        fetch_page_text=fake_page_fetch,
        fetch_jcamp_text=fake_jcamp_fetch,
    )

    assert len(results) == 2
    item = results[0]
    assert Path(item["jdx"]).exists()
    assert Path(item["csv"]).exists()
    assert Path(item["meta_json"]).exists()
    assert item["phase_tag"] == "gas"
    assert item["selection_priority"] == 100

    df = pd.read_csv(item["csv"])
    assert list(df.columns) == ["wavenumber_cm-1", "intensity"]
    assert len(df) == 5

    meta = json.loads(Path(item["meta_json"]).read_text(encoding="utf-8"))
    assert meta["nist_id"] == "C98862"
    assert meta["index"] == "0"
    assert meta["jcamp_metadata"]["TITLE"] == "Acetophenone"

    manifest_path = tmp_path / "KWOLFJPFCHCOCG-UHFFFAOYSA-N_reference_set.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["preferred_reference"]["phase_tag"] == "gas"
    assert len(manifest["reference_spectra"]) == 2
