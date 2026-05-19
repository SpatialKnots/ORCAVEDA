from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from mode_tracking import normalize_tracking_vector  # noqa: E402


def test_normalize_tracking_vector_returns_zero_for_zero_norm_input():
    vec = np.zeros((2, 3))

    normalized = normalize_tracking_vector(vec)

    assert np.array_equal(normalized, np.zeros_like(vec))
    assert normalized.shape == vec.shape


def test_normalize_tracking_vector_mass_weighted_zero_input_is_zero():
    vec = np.zeros((2, 3))

    normalized = normalize_tracking_vector(vec, masses=np.array([1.0, 16.0]), mass_weighted=True)

    assert np.array_equal(normalized, np.zeros_like(vec))
    assert normalized.shape == vec.shape
