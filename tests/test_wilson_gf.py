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
    WILSON_GF_VALIDATION_METHOD,
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
