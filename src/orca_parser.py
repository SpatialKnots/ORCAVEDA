from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import re

from orcaveda_models import HessData


BOHR_TO_ANGSTROM = 0.529177210903


def split_orca_hess_sections(text: str) -> Dict[str, List[str]]:
    lines = text.splitlines()
    starts = [(line.strip(), i) for i, line in enumerate(lines) if line.strip().startswith("$")]
    return {
        name: lines[i + 1 : (starts[k + 1][1] if k + 1 < len(starts) else len(lines))]
        for k, (name, i) in enumerate(starts)
    }


def parse_atoms(lines: List[str]) -> Tuple[List[str], np.ndarray, np.ndarray]:
    natoms = int(lines[0].strip())
    atoms, masses, coords = [], [], []
    for line in lines[1 : 1 + natoms]:
        parts = line.split()
        atoms.append(parts[0])
        masses.append(float(parts[1]))
        coords.append([float(x) * BOHR_TO_ANGSTROM for x in parts[2:5]])
    return atoms, np.array(masses), np.array(coords)


def parse_frequencies(lines: List[str]) -> np.ndarray:
    n = int(lines[0].strip())
    values = np.zeros(n)
    for line in lines[1:]:
        if line.strip():
            i, v = line.split()[:2]
            values[int(i)] = float(v)
    return values


def parse_scalar_section(lines: List[str]) -> Optional[float]:
    for line in lines:
        parts = line.split()
        for part in parts:
            try:
                return float(part.replace("D", "E"))
            except ValueError:
                continue
    return None


def parse_ir_spectrum(lines: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    n = int(lines[0].strip())
    freqs, intensities = [], []
    for line in lines[1 : 1 + n]:
        parts = line.split()
        if len(parts) >= 3:
            freqs.append(float(parts[0]))
            intensities.append(float(parts[2]))
    return np.array(freqs), np.array(intensities)


def parse_block_matrix(lines: List[str]) -> np.ndarray:
    shape = list(map(int, lines[0].split()[:2]))
    if len(shape) == 1:
        nrow = ncol = shape[0]
    else:
        nrow, ncol = shape[:2]
    matrix = np.zeros((nrow, ncol))
    current_cols: List[int] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 1 and all(re.fullmatch(r"\d+", part) for part in parts):
            current_cols = list(map(int, parts))
            continue
        if current_cols and re.fullmatch(r"\d+", parts[0]):
            row = int(parts[0])
            values = [float(value.replace("D", "E")) for value in parts[1:]]
            for col, value in zip(current_cols, values):
                matrix[row, col] = value
    return matrix


def read_orca_hess(path: str | Path) -> HessData:
    path = Path(path)
    sections = split_orca_hess_sections(path.read_text(errors="ignore"))
    required = ["$atoms", "$vibrational_frequencies", "$normal_modes", "$ir_spectrum"]
    missing = [section for section in required if section not in sections]
    if missing:
        raise ValueError(f"Missing required .hess sections in {path.name}: {missing}")

    atoms, masses, coords_A = parse_atoms(sections["$atoms"])
    frequencies = parse_frequencies(sections["$vibrational_frequencies"])
    _, intensities = parse_ir_spectrum(sections["$ir_spectrum"])
    normal_modes = parse_block_matrix(sections["$normal_modes"])
    cartesian_hessian = parse_block_matrix(sections["$hessian"]) if "$hessian" in sections else None
    temp = parse_scalar_section(sections.get("$actual_temperature", []))
    scale = parse_scalar_section(sections.get("$frequency_scale_factor", []))

    n3 = 3 * len(atoms)
    if normal_modes.shape != (n3, n3):
        raise ValueError(f"normal_modes shape mismatch: {normal_modes.shape}, expected {(n3, n3)}")
    if cartesian_hessian is not None and cartesian_hessian.shape != (n3, n3):
        raise ValueError(f"hessian shape mismatch: {cartesian_hessian.shape}, expected {(n3, n3)}")
    if len(frequencies) != n3 or len(intensities) != n3:
        raise ValueError("frequency/intensity length mismatch with 3N")

    return HessData(path.name, atoms, masses, coords_A, frequencies, intensities, normal_modes, temp, scale, cartesian_hessian)
