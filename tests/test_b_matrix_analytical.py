from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from b_matrix import analytical_B, finite_difference_B  # noqa: E402
from internal_coordinates import angle_fn, distance_fn, make_composed_internal_coordinate, torsion_fn  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402
from orcaveda_models import InternalCoordinate  # noqa: E402


def test_analytical_b_matrix_matches_finite_difference_for_distance_and_angle_rows():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.2399872, 0.927297, 0.0],
        ],
        dtype=float,
    )
    internals = [
        InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1)),
        InternalCoordinate("r(O1-H3)", "stretch", (0, 2), 10, distance_fn(0, 2)),
        InternalCoordinate("ang(H2-O1-H3)", "bend", (1, 0, 2), 30, angle_fn(1, 0, 2)),
    ]

    analytical, diagnostics = analytical_B(coords, internals)
    finite = finite_difference_B(coords, internals)

    assert diagnostics["method_counts"] == {"analytical_distance": 2, "analytical_angle": 1}
    assert diagnostics["fallback_reasons"] == {}
    assert np.allclose(analytical, finite, atol=1.0e-6, rtol=1.0e-6)


def test_analytical_b_matrix_falls_back_for_torsion_and_composed_rows():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    primitive = [
        InternalCoordinate("r(A1-A2)", "stretch", (0, 1), 10, distance_fn(0, 1)),
        InternalCoordinate("r(A2-A3)", "stretch", (1, 2), 10, distance_fn(1, 2)),
    ]
    composed = make_composed_internal_coordinate(
        "composed_stretch_sum",
        "composed_stretch",
        primitive,
        ((0, 1.0), (1, 1.0)),
    )
    internals = [
        *primitive,
        InternalCoordinate("tor(A1-A2-A3-A4)", "torsion", (0, 1, 2, 3), 55, torsion_fn(0, 1, 2, 3)),
        composed,
    ]

    analytical, diagnostics = analytical_B(coords, internals)
    finite = finite_difference_B(coords, internals)

    assert diagnostics["method_counts"] == {
        "analytical_distance": 2,
        "finite_difference_fallback": 2,
    }
    assert diagnostics["fallback_reasons"] == {
        "composed_coordinate": 1,
        "unsupported_coordinate_kind": 1,
    }
    assert np.allclose(analytical, finite, atol=1.0e-6, rtol=1.0e-6)


def test_analytical_b_matrix_falls_back_for_linear_angle_singularity():
    coords = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    internals = [
        InternalCoordinate("ang(A1-A2-A3)", "bend", (0, 1, 2), 30, angle_fn(0, 1, 2)),
    ]

    analytical, diagnostics = analytical_B(coords, internals)
    finite = finite_difference_B(coords, internals)

    assert diagnostics["method_counts"] == {"finite_difference_fallback": 1}
    assert diagnostics["fallback_reasons"] == {"singular_or_near_linear_angle": 1}
    assert np.allclose(analytical, finite, atol=1.0e-6, rtol=1.0e-6)


def test_analytical_b_matrix_falls_back_for_near_linear_angle():
    coords = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0e-4, 0.0],
        ],
        dtype=float,
    )
    internals = [
        InternalCoordinate("ang(A1-A2-A3)", "bend", (0, 1, 2), 30, angle_fn(0, 1, 2)),
    ]

    analytical, diagnostics = analytical_B(coords, internals)
    finite = finite_difference_B(coords, internals)

    assert diagnostics["method_counts"] == {"finite_difference_fallback": 1}
    assert diagnostics["fallback_reasons"] == {"singular_or_near_linear_angle": 1}
    assert np.allclose(analytical, finite, atol=1.0e-6, rtol=1.0e-6)


def test_analytical_b_matrix_falls_back_for_high_angle_baseline_parity():
    hess = read_orca_hess(ROOT / "data" / "hess" / "phenyl_isocyanate.hess")
    internals = [
        InternalCoordinate(
            "ang(N2-C13-O14)",
            "bend",
            (1, 12, 13),
            30,
            angle_fn(1, 12, 13),
        )
    ]

    analytical, diagnostics = analytical_B(hess.coords_A, internals)
    finite = finite_difference_B(hess.coords_A, internals)

    assert diagnostics["method_counts"] == {"finite_difference_fallback": 1}
    assert diagnostics["fallback_reasons"] == {"singular_or_near_linear_angle": 1}
    assert np.allclose(analytical, finite, atol=0.0, rtol=0.0)
