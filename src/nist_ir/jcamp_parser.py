from __future__ import annotations

from io import StringIO
from typing import Dict, List, Tuple

import pandas as pd


def _parse_metadata(lines: List[str]) -> Dict[str, str]:
    meta: Dict[str, str] = {}
    for raw in lines:
        line = raw.strip()
        if not line.startswith("##") or "=" not in line:
            continue
        key, value = line[2:].split("=", 1)
        meta[key.strip()] = value.strip()
    return meta


def _parse_xydata(lines: List[str], yfactor: float, deltax: float) -> Tuple[List[float], List[float]]:
    xs: List[float] = []
    ys: List[float] = []
    capture = False

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("##XYDATA="):
            capture = True
            continue
        if not capture:
            continue
        if line.startswith("##END"):
            break
        if line.startswith("##"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            start_x = float(parts[0])
            y_values = [float(part) * yfactor for part in parts[1:]]
        except ValueError:
            continue

        x = start_x
        for y in y_values:
            xs.append(x)
            ys.append(y)
            x += deltax

    return xs, ys


def parse_jcamp_text(jdx_text: str) -> Tuple[Dict[str, str], pd.DataFrame]:
    text = str(jdx_text or "")
    lines = text.splitlines()
    metadata = _parse_metadata(lines)

    yfactor = float(metadata.get("YFACTOR", "1.0"))
    deltax = float(metadata.get("DELTAX", "1.0"))
    xs, ys = _parse_xydata(lines, yfactor=yfactor, deltax=deltax)

    spectrum = pd.read_csv(
        StringIO(
            "\n".join(
                f"{x},{y}"
                for x, y in zip(xs, ys)
            )
        ),
        names=["wavenumber_cm-1", "intensity"],
    )
    spectrum = spectrum.sort_values("wavenumber_cm-1").reset_index(drop=True)
    return metadata, spectrum
