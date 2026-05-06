from __future__ import annotations

import sys
import shutil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from b_matrix import finite_difference_B, optimize_independent_coordinates_for_ped, ped_basis_localization_metrics, svd_rank_condition  # noqa: E402
from ORCAVEDA_patched_stage3D_v5_0 import analyze_orca_ped_like  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402
from orcaveda_models import HessData, InternalCoordinate  # noqa: E402
from ped import (  # noqa: E402
    PED_V1_METHOD,
    PED_V2_METHOD,
    WILSON_PED_METHOD,
    BOHR_TO_ANGSTROM,
    build_ped_audit_dataframe,
    build_ped_v2_force_audit_dataframe,
    build_wilson_g_matrix,
    build_wilson_ped_audit_dataframe,
    compute_ped,
    reconstruct_internal_force_matrix,
)
from internal_coordinates import angle_fn, distance_fn  # noqa: E402
from reports import decide_ped_driven_final_assignment  # noqa: E402


def _hess_with_mode(atoms, coords, mode_vec, freq=3650.0):
    n3 = 3 * len(atoms)
    normal_modes = np.zeros((n3, n3), dtype=float)
    normal_modes[:, 0] = np.asarray(mode_vec, dtype=float).reshape(-1)
    return HessData(
        filename="synthetic.hess",
        atoms=list(atoms),
        masses=np.ones(len(atoms), dtype=float),
        coords_A=np.asarray(coords, dtype=float),
        frequencies_cm1=np.array([freq] + [0.0] * (n3 - 1), dtype=float),
        ir_intensities=np.zeros(n3, dtype=float),
        normal_modes=normal_modes,
    )


def test_ped_v1_uses_column_normal_mode_orientation_for_oh_stretch():
    atoms = ["O", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float)
    internals = [InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1))]
    B = finite_difference_B(coords, internals)
    hess = _hess_with_mode(atoms, coords, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    results = compute_ped(hess, internals, B, [0], source_label="[1]")

    mode0 = results[0]
    assert mode0.contributions[0].internal_coordinate == "r(O1-H2)"
    assert mode0.contributions[0].coordinate_family == "O-H stretch"
    assert mode0.contributions[0].percent == 100.0
    assert mode0.normalization_sum_percent == 100.0
    assert "Wilson GF" in PED_V1_METHOD
    assert "not force-constant" in PED_V1_METHOD


def test_ped_v1_reports_mixed_water_like_symmetric_stretch():
    atoms = ["O", "H", "H"]
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    internals = [
        InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1)),
        InternalCoordinate("r(O1-H3)", "stretch", (0, 2), 10, distance_fn(0, 2)),
        InternalCoordinate("ang(H2-O1-H3)", "bend", (1, 0, 2), 30, angle_fn(1, 0, 2)),
    ]
    B = finite_difference_B(coords, internals)
    hess = _hess_with_mode(
        atoms,
        coords,
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
    )

    df = build_ped_audit_dataframe(hess, internals, B, [0, 1, 2], source_label="[1]", top_n=3)
    mode0 = df[df["mode"] == 0].sort_values("ped_rank")

    assert list(mode0["internal_coordinate"].head(2)) == ["r(O1-H3)", "r(O1-H2)"]
    assert np.allclose(sorted(mode0["contribution_percent"].head(2)), [50.0, 50.0])
    assert mode0["normalization_sum_percent"].iloc[0] == 100.0
    assert set(mode0["coordinate_family"].head(2)) == {"O-H stretch"}


def test_ped_v1_emits_diagnostics_for_zero_projection():
    atoms = ["O", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float)
    internals = [InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1))]
    B = finite_difference_B(coords, internals)
    hess = _hess_with_mode(atoms, coords, [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    df = build_ped_audit_dataframe(hess, internals, B, [0], source_label="[1]", top_n=1)
    mode0 = df[df["mode"] == 0].iloc[0]

    assert mode0["ped_rank"] == 0
    assert "zero_ped_projection_weight" in mode0["ped_warnings"]


def test_pipeline_writes_separate_ped_audit_for_water():
    outdir = ROOT / "outputs" / "pytest_ped_pipeline_h2o"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tables = analyze_orca_ped_like([ROOT / "data" / "hess" / "H2O_freq.hess"], outdir)

    assert "ped_audit" in tables
    assert "ped_v2_force_audit" in tables
    assert "wilson_ped_audit" in tables
    assert "optimized_ped_basis" in tables
    assert "ped_stage3d_agreement" in tables
    assert "ped_final_assignment" in tables
    ped_audit = tables["ped_audit"]
    ped_v2 = tables["ped_v2_force_audit"]
    wilson_ped = tables["wilson_ped_audit"]
    optimized_ped_basis = tables["optimized_ped_basis"]
    agreement = tables["ped_stage3d_agreement"]
    final_assignment = tables["ped_final_assignment"]
    assert not ped_audit.empty
    assert not ped_v2.empty
    assert not wilson_ped.empty
    assert not optimized_ped_basis.empty
    assert not agreement.empty
    assert not final_assignment.empty
    assert next(outdir.glob("*__ped_audit.csv")).is_file()
    assert next(outdir.glob("*__ped_v2_force_audit.csv")).is_file()
    assert next(outdir.glob("*__wilson_ped_audit.csv")).is_file()
    assert next(outdir.glob("*__optimized_ped_basis.csv")).is_file()
    assert next(outdir.glob("*__ped_stage3d_agreement.csv")).is_file()
    assert next(outdir.glob("*__ped_final_assignment.csv")).is_file()
    assert {
        "stage3d_assignment",
        "ped_assignment",
        "ped_source",
        "ped_agreement_status",
        "ped_policy_warning",
    }.issubset(set(agreement.columns))
    assert set(agreement["ped_agreement_status"]).issubset({"confirms", "adds_context", "disagrees", "diffuse", "not_available"})
    assert {
        "final_assignment",
        "final_assignment_source",
        "final_assignment_policy",
        "final_assignment_warning",
    }.issubset(set(final_assignment.columns))
    assert final_assignment["final_assignment"].astype(str).str.len().gt(0).all()

    positive = ped_audit[ped_audit["frequency_cm-1"] > 0.0].copy()
    assert positive["coordinate_family"].astype(str).str.contains("O-H stretch", regex=False).any()
    assert positive["coordinate_family"].astype(str).str.contains("H-O-H bend", regex=False).any()


def test_ped_basis_optimizer_improves_localization_without_losing_rank():
    def dummy_coord(_xyz):
        return 0.0

    internals = [
        InternalCoordinate("mixed_a", "stretch", (0, 1), 10, dummy_coord),
        InternalCoordinate("mixed_b", "stretch", (0, 2), 10, dummy_coord),
        InternalCoordinate("localized_x", "stretch", (0, 3), 20, dummy_coord),
    ]
    B = np.array(
        [
            [1.0, 1.0],
            [1.0, -1.0],
            [1.0, 0.0],
        ],
        dtype=float,
    )
    normal_modes = np.eye(2, dtype=float)
    selected = [0, 1]

    initial = ped_basis_localization_metrics(B, normal_modes, selected, [0])
    optimized, report = optimize_independent_coordinates_for_ped(
        B,
        internals,
        selected,
        normal_modes,
        [0],
        target_rank=2,
    )
    final = ped_basis_localization_metrics(B, normal_modes, optimized, [0])
    rank, _, _ = svd_rank_condition(B[optimized, :])

    assert report["changed"] is True
    assert rank == 2
    assert final["localization_score"] > initial["localization_score"]
    assert 2 in optimized


def test_ped_final_policy_reclassifies_high_confidence_stage3d_torsion_bend_conflict():
    final, source, policy, warning = decide_ped_driven_final_assignment(
        "C-H torsion",
        "C-C-H bend",
        "Wilson GF-style PED audit",
        "disagrees",
        "ped_diagnostic_basis; ped_stage3d_semantic_disagreement",
        91.2,
    )

    assert final == "C-C-H bend"
    assert source == "Wilson GF-style PED audit"
    assert policy == "ped_reclassifies_stage3d_torsion"
    assert "stage3d_torsion_reclassified_by_high_confidence_ped_bend" in warning


def test_ped_final_policy_keeps_stage3d_for_diffuse_or_low_confidence_torsion_bend_conflict():
    low_confidence = decide_ped_driven_final_assignment(
        "C-C-C torsion",
        "aromatic C-H bend mixed with aromatic ring deformation",
        "Wilson GF-style PED audit",
        "disagrees",
        "ped_diagnostic_basis; ped_stage3d_semantic_disagreement",
        31.6,
    )
    diffuse = decide_ped_driven_final_assignment(
        "C-C-C torsion",
        "N-C-C bend mixed with aromatic ring deformation",
        "Wilson GF-style PED audit",
        "disagrees",
        "diffuse_ped_contributions; ped_diagnostic_basis; ped_stage3d_semantic_disagreement",
        75.0,
    )

    assert low_confidence[0] == "C-C-C torsion"
    assert low_confidence[2] == "stage3d_fallback_due_to_ped_disagreement"
    assert diffuse[0] == "C-C-C torsion"
    assert diffuse[2] == "stage3d_fallback_due_to_ped_disagreement"


def test_orca_parser_reads_cartesian_hessian_block():
    hess = read_orca_hess(ROOT / "data" / "hess" / "H2O_freq.hess")

    assert hess.cartesian_hessian is not None
    assert hess.cartesian_hessian.shape == (9, 9)
    assert np.allclose(hess.cartesian_hessian, hess.cartesian_hessian.T, atol=1.0e-3)


def test_ped_v2_force_aware_weights_use_cartesian_hessian():
    def x_coord(xyz):
        return float(xyz[0, 0])

    def y_coord(xyz):
        return float(xyz[0, 1])

    atoms = ["C"]
    coords = np.array([[0.0, 0.0, 0.0]], dtype=float)
    internals = [
        InternalCoordinate("cart_x(C1)", "cartesian_probe", (0,), 10, x_coord),
        InternalCoordinate("cart_y(C1)", "cartesian_probe", (0,), 10, y_coord),
    ]
    B = finite_difference_B(coords, internals)
    hess = _hess_with_mode(atoms, coords, [[1.0, 1.0, 0.0]], freq=1000.0)
    hess.cartesian_hessian = np.diag([9.0, 1.0, 0.0])

    v1 = build_ped_audit_dataframe(hess, internals, B, [0, 1], top_n=2)
    v2 = build_ped_v2_force_audit_dataframe(hess, internals, B, [0, 1], top_n=2)

    assert np.allclose(sorted(v1[v1["mode"] == 0]["contribution_percent"]), [50.0, 50.0])
    mode0_v2 = v2[v2["mode"] == 0].sort_values("ped_rank")
    assert mode0_v2.iloc[0]["internal_coordinate"] == "cart_x(C1)"
    assert np.isclose(float(mode0_v2.iloc[0]["contribution_percent"]), 90.0)
    assert np.isclose(float(mode0_v2.iloc[1]["contribution_percent"]), 10.0)
    assert "force-aware" in PED_V2_METHOD
    assert "not full Wilson GF PED" in PED_V2_METHOD


def test_wilson_g_matrix_for_oh_stretch_matches_mass_inverse_sum():
    atoms = ["O", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float)
    internals = [InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1))]
    B = finite_difference_B(coords, internals)

    G = build_wilson_g_matrix(B, np.array([16.0, 1.0]))

    assert G.shape == (1, 1)
    assert np.isclose(G[0, 0], 1.0 / 16.0 + 1.0, atol=1.0e-6)


def test_wilson_internal_force_reconstruction_for_one_coordinate():
    atoms = ["O", "H"]
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float)
    internals = [InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1))]
    B = finite_difference_B(coords, internals)
    k_internal = 2.5
    f_cart_A = B.T @ np.array([[k_internal]]) @ B
    f_cart_bohr = f_cart_A * (BOHR_TO_ANGSTROM ** 2)

    F_internal = reconstruct_internal_force_matrix(B, f_cart_bohr)

    assert F_internal.shape == (1, 1)
    assert np.isclose(F_internal[0, 0], k_internal, atol=1.0e-6)


def test_wilson_ped_audit_reports_h2o_bend_and_stretches():
    hess = read_orca_hess(ROOT / "data" / "hess" / "H2O_freq.hess")
    internals = [
        InternalCoordinate("r(O1-H2)", "stretch", (0, 1), 10, distance_fn(0, 1)),
        InternalCoordinate("r(O1-H3)", "stretch", (0, 2), 10, distance_fn(0, 2)),
        InternalCoordinate("ang(H2-O1-H3)", "bend", (1, 0, 2), 30, angle_fn(1, 0, 2)),
    ]
    B = finite_difference_B(hess.coords_A, internals)

    df = build_wilson_ped_audit_dataframe(hess, internals, B, [0, 1, 2], top_n=3)
    positive = df[df["frequency_cm-1"] > 0.0].copy()

    assert positive["wilson_ped_method"].astype(str).str.contains("Wilson GF PED", regex=False).all()
    assert positive["coordinate_family"].astype(str).str.contains("H-O-H bend", regex=False).any()
    assert positive["coordinate_family"].astype(str).str.contains("O-H stretch", regex=False).any()
    assert "G = B M^-1 B^T" in WILSON_PED_METHOD
