from __future__ import annotations

import sys
import shutil
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from b_matrix import finite_difference_B, select_independent_coordinates  # noqa: E402
from chemistry import annotate_chemical_system  # noqa: E402
from internal_coordinates import build_internal_coordinates  # noqa: E402
from ORCAVEDA_patched_stage3D_v5_0 import analyze_orca_ped_like  # noqa: E402
from orca_parser import read_orca_hess  # noqa: E402
from wilson_gf import (  # noqa: E402
    VEDA_LIKE_PED_METHOD,
    WILSON_GF_VALIDATION_METHOD,
    build_veda_like_mode_correspondence_dataframe,
    build_veda_like_ped_audit_dataframe,
    build_veda_like_ped_matrix_dataframe,
    build_wilson_gf_basis_diagnostics_dataframe,
    build_wilson_gf_validation_dataframe,
    solve_symmetric_gf_eigenproblem,
    symmetric_sqrt_decomp,
    wilson_gf_closed_ped,
    wilson_gf_diagonalization,
)


def _pipeline_basis(hess_name: str):
    hess = read_orca_hess(ROOT / "data" / "hess" / hess_name)
    annotation = annotate_chemical_system(hess.atoms, hess.coords_A)
    internals = build_internal_coordinates(
        hess.atoms,
        hess.coords_A,
        list(annotation.bonds),
        [list(fragment) for fragment in annotation.fragments],
        list(annotation.interfragment_hbonds),
        list(annotation.functional_groups),
    )
    B = finite_difference_B(hess.coords_A, internals)
    selected_idx, rank, cond, _ = select_independent_coordinates(B, internals, 3 * len(hess.atoms) - 6)
    return hess, internals, B, selected_idx, rank, cond


def _h2o_pipeline_basis():
    return _pipeline_basis("H2O_freq.hess")


def test_symmetric_sqrt_decomp_identity_and_spd_matrix():
    ident_sqrt, ident_eigs = symmetric_sqrt_decomp(np.eye(3))

    assert np.allclose(ident_sqrt, np.eye(3))
    assert np.allclose(ident_eigs, np.ones(3))

    matrix = np.array([[4.0, 0.0], [0.0, 9.0]])
    sqrt_matrix, eigs = symmetric_sqrt_decomp(matrix)

    assert np.allclose(sqrt_matrix @ sqrt_matrix, matrix)
    assert np.allclose(eigs, [4.0, 9.0])


def test_symmetric_sqrt_decomp_rejects_invalid_shapes_and_symmetry():
    with pytest.raises(ValueError, match="square"):
        symmetric_sqrt_decomp(np.ones((2, 3)))

    with pytest.raises(ValueError, match="symmetric"):
        symmetric_sqrt_decomp(np.array([[1.0, 2.0], [0.0, 1.0]]))


def test_solve_symmetric_gf_eigenproblem_known_diagonal_case():
    G = np.diag([4.0, 9.0])
    F = np.diag([25.0, 16.0])

    eigenvalues, eigenvectors, symmetric_gf = solve_symmetric_gf_eigenproblem(G, F)

    assert np.allclose(eigenvalues, [100.0, 144.0])
    assert eigenvectors.shape == (2, 2)
    assert np.allclose(symmetric_gf, np.diag([100.0, 144.0]))


def test_h2o_wilson_gf_validation_records_basis_and_conversion_diagnostics():
    hess, internals, B, selected_idx, rank, cond = _h2o_pipeline_basis()

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    validation = build_wilson_gf_validation_dataframe(result)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)

    assert hess.cartesian_hessian is not None
    assert B.shape == (3, 9)
    assert selected_idx == [0, 1, 2]
    assert rank == 3
    assert np.isfinite(cond)
    assert result.expected_vibrational_rank == 3
    assert result.internal_basis_size == 3
    assert result.g_rank == 3
    assert result.f_rank == 3
    assert len(result.orca_frequencies_cm1) == 3
    assert result.validation_status == "PASS"
    assert result.max_relative_error < 1.0e-4
    assert "fixed_conversion_failed" not in result.warnings
    assert "not VEDA-equivalent" in WILSON_GF_VALIDATION_METHOD

    required = {
        "method",
        "basis_size",
        "expected_vibrational_rank",
        "conversion_method",
        "mapping_method",
        "max_relative_error",
        "warnings",
        "validation_status",
    }
    assert required.issubset(set(validation.columns))
    assert set(validation["method"]) == {WILSON_GF_VALIDATION_METHOD}
    assert validation["conversion_method"].astype(str).str.contains("fixed_SI", regex=False).all()
    assert int(basis.iloc[0]["basis_size"]) == 3
    assert int(basis.iloc[0]["positive_orca_mode_count"]) == 3


def test_h2o_closed_ped_rows_normalize_to_100_percent_per_mode():
    hess, internals, B, selected_idx, _, _ = _h2o_pipeline_basis()
    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)

    ped = wilson_gf_closed_ped(result, hess, internals, B, selected_idx, top_n=3)

    assert not ped.empty
    assert ped["method"].astype(str).str.contains("Wilson GF diagonalization validation prototype", regex=False).all()
    sums = ped.groupby("mode")["contribution_percent"].sum()
    assert np.allclose(sums.to_numpy(dtype=float), np.full(len(sums), 100.0), atol=1.0e-6)
    normalization = ped.groupby("mode")["normalization_sum_percent"].first()
    assert np.allclose(normalization.to_numpy(dtype=float), np.full(len(normalization), 100.0), atol=1.0e-6)


def test_h2o_veda_like_outputs_are_separate_and_normalized():
    hess, internals, B, selected_idx, _, _ = _h2o_pipeline_basis()
    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)

    audit = build_veda_like_ped_audit_dataframe(result, hess, internals, B, selected_idx, top_n=3)
    matrix = build_veda_like_ped_matrix_dataframe(result, hess, internals, B, selected_idx)
    correspondence = build_veda_like_mode_correspondence_dataframe(result)

    assert not audit.empty
    assert not matrix.empty
    assert not correspondence.empty
    assert set(audit["method"]) == {VEDA_LIKE_PED_METHOD}
    assert set(matrix["method"]) == {VEDA_LIKE_PED_METHOD}
    assert "does not reproduce original VEDA" in VEDA_LIKE_PED_METHOD
    assert set(matrix["matrix_orientation"]) == {"mode_rows_by_coordinate_columns_long_form"}
    sums = matrix.groupby("mode")["contribution_percent"].sum()
    assert np.allclose(sums.to_numpy(dtype=float), np.full(len(sums), 100.0), atol=1.0e-6)
    assert set(correspondence["validation_status"]) == {"PASS"}


@pytest.mark.parametrize(
    ("hess_name", "expected_rank", "high_frequency_family"),
    [
        ("NH3.hess", 6, "N-H stretch"),
        ("formaldehyde.hess", 6, "C-H stretch"),
    ],
)
def test_small_molecule_veda_like_outputs_pass_and_keep_xh_stretches(
    hess_name: str,
    expected_rank: int,
    high_frequency_family: str,
):
    hess, internals, B, selected_idx, _, _ = _pipeline_basis(hess_name)
    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)

    audit = build_veda_like_ped_audit_dataframe(
        result,
        hess,
        result.validation_internals or internals,
        result.validation_B if result.validation_B is not None else B,
        result.basis_indices,
        top_n=8,
    )
    matrix = build_veda_like_ped_matrix_dataframe(
        result,
        hess,
        result.validation_internals or internals,
        result.validation_B if result.validation_B is not None else B,
        result.basis_indices,
    )
    correspondence = build_veda_like_mode_correspondence_dataframe(result)

    assert result.expected_vibrational_rank == expected_rank
    assert result.validation_status == "PASS"
    assert result.warnings == ()
    assert set(correspondence["validation_status"]) == {"PASS"}
    assert matrix.shape[0] == expected_rank * expected_rank
    assert np.allclose(
        matrix.groupby("mode")["contribution_percent"].sum().to_numpy(dtype=float),
        np.full(matrix["mode"].nunique(), 100.0),
        atol=1.0e-6,
    )

    dominant = audit[audit["veda_like_rank"] == 1]
    high_frequency = dominant[dominant["frequency_cm-1"] > 2800.0]
    assert not high_frequency.empty
    assert set(high_frequency["coordinate_family"]) == {high_frequency_family}
    assert set(high_frequency["validation_status"]) == {"PASS"}


def test_ethene_wilson_gf_validation_selects_conditioned_basis():
    hess, internals, B, selected_idx, rank, cond = _pipeline_basis("ethene.hess")

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)
    ped = wilson_gf_closed_ped(result, hess, internals, B, result.basis_indices, top_n=8)

    assert rank == 12
    assert cond > 1.0e6
    assert tuple(selected_idx) != result.basis_indices
    assert result.validation_status == "PASS"
    assert result.g_rank == 12
    assert result.f_rank == 12
    assert result.g_condition < 100.0
    assert result.f_condition < 100.0
    assert result.max_relative_error < 1.0e-6
    assert result.warnings == ()
    assert str(basis.iloc[0]["selected_indices"]) == "0;4;5;6;7;9;12;13;14;15;17;18"
    assert not ped.empty
    assert set(ped["validation_status"]) == {"PASS"}


def test_ch3cn_wilson_gf_uses_linear_bend_components_for_near_linear_bend():
    hess, internals, B, selected_idx, rank, _ = _pipeline_basis("CH3CN_freq.hess")

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)

    assert rank == 12
    assert result.validation_status == "PASS"
    assert result.max_relative_error < 1.0e-6
    assert "linear_bend_coordinate_used" in result.warnings
    assert "near_linear_bend_coordinate" not in result.warnings
    assert "fixed_conversion_failed" not in result.warnings
    assert "empirical_ratio_only" not in result.warnings
    assert str(basis.iloc[0]["selected_indices"]) == "1;2;4;5;6;10;11;12;13;14;19;20"
    assert "linear_bend_coordinate_used" in str(basis.iloc[0]["warnings"])


@pytest.mark.parametrize(
    ("hess_name", "geometry_warning"),
    [
        ("ethyne.hess", "linear_bend_coordinate_used"),
        ("propyne.hess", "near_linear_bend_coordinate"),
    ],
)
def test_linear_edge_cases_keep_fixed_conversion_review_warning(hess_name: str, geometry_warning: str):
    hess, internals, B, selected_idx, _, _ = _pipeline_basis(hess_name)

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)
    correspondence = build_veda_like_mode_correspondence_dataframe(result)

    assert result.validation_status == "WARN"
    assert geometry_warning in result.warnings
    assert "fixed_conversion_failed" in result.warnings
    assert "linear_or_near_linear_fixed_conversion_review" in result.warnings
    assert "empirical_ratio_only" in result.warnings
    assert set(correspondence["validation_status"]) == {"WARN"}
    assert "linear_or_near_linear_fixed_conversion_review" in str(basis.iloc[0]["warnings"])


@pytest.mark.parametrize(
    ("hess_name", "mode"),
    [
        ("monoethanolamine_dimer_NH_to_O_DFT.hess", 63),
        ("monoethanolamine_dimer_OH_to_N_DFT.hess", 60),
    ],
)
def test_hbonded_high_frequency_xh_modes_report_secondary_stretch_warning(hess_name: str, mode: int):
    hess, internals, B, selected_idx, _, _ = _pipeline_basis(hess_name)
    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)

    audit = build_veda_like_ped_audit_dataframe(
        result,
        hess,
        result.validation_internals or internals,
        result.validation_B if result.validation_B is not None else B,
        result.basis_indices,
        top_n=8,
    )
    mode_rows = audit[audit["mode"] == mode].sort_values("veda_like_rank")

    assert result.validation_status == "PASS"
    assert not mode_rows.empty
    assert mode_rows.iloc[0]["coordinate_family"] == "H-bond / intermolecular"
    assert mode_rows.iloc[0]["frequency_cm-1"] > 2800.0
    assert mode_rows["coordinate_family"].astype(str).str.contains("N-H stretch", regex=False).any()
    assert mode_rows["warnings"].astype(str).str.contains(
        "high_frequency_hbond_dominates_xh_stretch_secondary",
        regex=False,
    ).all()


def test_aniline_wilson_gf_warns_when_positive_modes_are_below_expected_rank():
    hess, internals, B, selected_idx, _, _ = _pipeline_basis("aniline.hess")

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)

    assert result.expected_vibrational_rank == 36
    assert len(result.orca_frequencies_cm1) == 35
    assert len(result.gf_eigenvalues) == 35
    assert result.validation_status == "PASS"
    assert "positive_orca_mode_count_below_expected_vibrational_rank" in result.warnings
    assert "positive_gf_eigenvalue_count_below_expected_vibrational_rank" in result.warnings
    assert "nonpositive_orca_modes_within_expected_vibrational_space" in result.warnings
    assert "nonpositive_gf_eigenvalues_within_expected_vibrational_space" in result.warnings
    assert result.orca_nonpositive_mode_count == 7
    assert result.orca_min_nonpositive_frequency_cm1 == pytest.approx(-367.3568575130425)
    assert result.gf_nonpositive_eigenvalue_count == 1
    assert "mode_count_mismatch" not in result.warnings
    assert int(basis.iloc[0]["positive_orca_mode_count"]) == 35
    assert int(basis.iloc[0]["positive_gf_eigenvalue_count"]) == 35
    assert int(basis.iloc[0]["orca_nonpositive_mode_count"]) == 7


def test_acetanilide_wilson_gf_uses_large_system_conditioned_basis():
    hess, internals, B, selected_idx, rank, _ = _pipeline_basis("acetanilide.hess")

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)

    assert rank == 51
    assert result.expected_vibrational_rank == 51
    assert result.internal_basis_size == 51
    assert tuple(selected_idx) != result.basis_indices
    assert result.validation_status == "PASS"
    assert result.g_rank == 51
    assert result.f_rank == 51
    assert result.g_condition < 1.0e6
    assert result.f_condition < 1.0e6
    assert result.max_relative_error < 1.0e-6
    assert "basis_rank_below_expected" not in result.warnings
    assert "g_ill_conditioned" not in result.warnings
    assert "f_ill_conditioned" not in result.warnings
    assert int(basis.iloc[0]["positive_orca_mode_count"]) == 50
    assert int(basis.iloc[0]["positive_gf_eigenvalue_count"]) == 50
    assert "positive_orca_mode_count_below_expected_vibrational_rank" in result.warnings
    assert "positive_gf_eigenvalue_count_below_expected_vibrational_rank" in result.warnings
    assert "nonpositive_orca_modes_within_expected_vibrational_space" in result.warnings
    assert "nonpositive_gf_eigenvalues_within_expected_vibrational_space" in result.warnings
    assert result.orca_nonpositive_mode_count == 7
    assert result.orca_min_nonpositive_frequency_cm1 == pytest.approx(-49.58924532088554)
    assert result.gf_nonpositive_eigenvalue_count == 1
    assert int(basis.iloc[0]["orca_nonpositive_mode_count"]) == 7


def test_n_methylaniline_wilson_gf_large_system_basis_improves_f_condition():
    hess, internals, B, selected_idx, rank, _ = _pipeline_basis("N-methylaniline.hess")

    result = wilson_gf_diagonalization(hess, internals, B, selected_idx)
    basis = build_wilson_gf_basis_diagnostics_dataframe(result)

    assert rank == 45
    assert result.expected_vibrational_rank == 45
    assert result.internal_basis_size == 45
    assert tuple(selected_idx) != result.basis_indices
    assert result.validation_status == "PASS"
    assert result.g_rank == 45
    assert result.f_rank == 45
    assert result.g_condition < 1.0e6
    assert result.f_condition < 1.0e6
    assert result.max_relative_error < 1.0e-6
    assert "g_ill_conditioned" not in result.warnings
    assert "f_ill_conditioned" not in result.warnings
    assert int(basis.iloc[0]["positive_orca_mode_count"]) == 44
    assert int(basis.iloc[0]["positive_gf_eigenvalue_count"]) == 44
    assert "positive_orca_mode_count_below_expected_vibrational_rank" in result.warnings
    assert "positive_gf_eigenvalue_count_below_expected_vibrational_rank" in result.warnings
    assert "nonpositive_orca_modes_within_expected_vibrational_space" in result.warnings
    assert "nonpositive_gf_eigenvalues_within_expected_vibrational_space" in result.warnings
    assert result.orca_nonpositive_mode_count == 7
    assert result.orca_min_nonpositive_frequency_cm1 == pytest.approx(-198.14705984767406)
    assert result.gf_nonpositive_eigenvalue_count == 1
    assert int(basis.iloc[0]["orca_nonpositive_mode_count"]) == 7


def test_pipeline_wilson_gf_validation_is_opt_in_for_h2o():
    default_outdir = ROOT / "outputs" / "pytest_wilson_gf_default_h2o"
    opt_in_outdir = ROOT / "outputs" / "pytest_wilson_gf_opt_in_h2o"
    for outdir in (default_outdir, opt_in_outdir):
        if outdir.exists():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

    hess_path = ROOT / "data" / "hess" / "H2O_freq.hess"
    default_tables = analyze_orca_ped_like([hess_path], default_outdir)

    assert "wilson_gf_validation" not in default_tables
    assert not list(default_outdir.glob("*__wilson_gf_validation.csv"))
    assert not list(default_outdir.glob("*__wilson_gf_ped_audit.csv"))
    assert not list(default_outdir.glob("*__wilson_gf_basis_diagnostics.csv"))

    opt_in_tables = analyze_orca_ped_like([hess_path], opt_in_outdir, wilson_gf_validation=True)

    assert {"wilson_gf_validation", "wilson_gf_ped_audit", "wilson_gf_basis_diagnostics"}.issubset(opt_in_tables)
    validation = opt_in_tables["wilson_gf_validation"]
    assert not validation.empty
    assert set(validation["validation_status"]) == {"PASS"}
    assert float(validation["max_relative_error"].max()) < 1.0e-4
    assert next(opt_in_outdir.glob("*__wilson_gf_validation.csv")).is_file()
    assert next(opt_in_outdir.glob("*__wilson_gf_ped_audit.csv")).is_file()
    assert next(opt_in_outdir.glob("*__wilson_gf_basis_diagnostics.csv")).is_file()


def test_pipeline_veda_like_ped_is_opt_in_for_h2o():
    default_outdir = ROOT / "outputs" / "pytest_veda_like_default_h2o"
    opt_in_outdir = ROOT / "outputs" / "pytest_veda_like_opt_in_h2o"
    for outdir in (default_outdir, opt_in_outdir):
        if outdir.exists():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

    hess_path = ROOT / "data" / "hess" / "H2O_freq.hess"
    default_tables = analyze_orca_ped_like([hess_path], default_outdir)

    assert "veda_like_ped_audit" not in default_tables
    assert not list(default_outdir.glob("*__veda_like_ped_audit.csv"))
    assert not list(default_outdir.glob("*__veda_like_metadata.json"))

    opt_in_tables = analyze_orca_ped_like([hess_path], opt_in_outdir, veda_like_ped=True)

    expected = {
        "veda_like_ped_audit",
        "veda_like_ped_matrix",
        "veda_like_basis_diagnostics",
        "veda_like_mode_correspondence",
    }
    assert expected.issubset(opt_in_tables)
    assert "wilson_gf_validation" not in opt_in_tables
    assert not opt_in_tables["veda_like_ped_audit"].empty
    matrix = opt_in_tables["veda_like_ped_matrix"]
    assert np.allclose(
        matrix.groupby("mode")["contribution_percent"].sum().to_numpy(dtype=float),
        np.full(matrix["mode"].nunique(), 100.0),
        atol=1.0e-6,
    )
    assert next(opt_in_outdir.glob("*__veda_like_ped_audit.csv")).is_file()
    assert next(opt_in_outdir.glob("*__veda_like_ped_matrix.csv")).is_file()
    assert next(opt_in_outdir.glob("*__veda_like_basis_diagnostics.csv")).is_file()
    assert next(opt_in_outdir.glob("*__veda_like_mode_correspondence.csv")).is_file()
    metadata_path = next(opt_in_outdir.glob("*__veda_like_metadata.json"))
    metadata = metadata_path.read_text(encoding="utf-8")
    assert "does not reproduce original VEDA" in metadata
    assert "normal_modes[:, mode]" in metadata
