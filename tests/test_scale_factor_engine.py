from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scale_factor_engine import (  # noqa: E402
    apply_piecewise_scale,
    apply_scale,
    bootstrap_constant_scale,
    compare_scaling_models,
    fit_constant_scale,
    fit_constant_scale_robust,
    fit_inverse_frequency_scale,
    fit_piecewise_constant_scale,
    fit_power_law_frequency_dependent,
    fit_linear_frequency_dependent,
    metrics,
    predict_linear_frequency_dependent,
    predict_power_law_frequency_dependent,
)


def test_constant_scale_and_metrics_basic_case():
    omega = np.array([1000.0, 1500.0, 3000.0])
    nu = 0.96 * omega
    k = fit_constant_scale(omega, nu)
    pred = apply_scale(omega, k)
    m = metrics(pred, nu)
    assert abs(k - 0.96) < 1e-10
    assert m.mae < 1e-10
    assert m.rmse < 1e-10


def test_weighted_and_robust_scale_downweight_outlier():
    omega = np.array([1000.0, 1500.0, 1700.0, 3000.0])
    nu = np.array([960.0, 1440.0, 1632.0, 2400.0])  # last point is bad outlier
    plain = fit_constant_scale(omega, nu)
    weighted = fit_constant_scale(omega, nu, weights=np.array([1.0, 1.0, 1.0, 0.1]))
    huber = fit_constant_scale_robust(omega, nu, delta=50.0)
    assert abs(weighted - 0.96) < abs(plain - 0.96)
    assert abs(huber - 0.96) < abs(plain - 0.96)


def test_inverse_frequency_scale_exact_case():
    omega = np.array([200.0, 450.0, 800.0])
    nu = 0.92 * omega
    k = fit_inverse_frequency_scale(omega, nu)
    assert abs(k - 0.92) < 1e-10


def test_piecewise_scale_recovers_region_factors():
    omega = np.array([600.0, 900.0, 1300.0, 1800.0, 3100.0, 3400.0])
    nu = np.array([0.94 * 600.0, 0.94 * 900.0, 0.97 * 1300.0, 0.97 * 1800.0, 0.99 * 3100.0, 0.99 * 3400.0])
    regions = [(0.0, 1200.0), (1200.0, 3000.0), (3000.0, float("inf"))]
    fitted = fit_piecewise_constant_scale(omega, nu, regions)
    scaled = apply_piecewise_scale(omega, fitted)
    assert np.allclose(scaled, nu)
    assert abs(fitted[0]["k"] - 0.94) < 1e-10
    assert abs(fitted[1]["k"] - 0.97) < 1e-10
    assert abs(fitted[2]["k"] - 0.99) < 1e-10


def test_linear_and_power_law_models_reconstruct_training_data():
    omega = np.array([500.0, 900.0, 1500.0, 3000.0])

    pred_linear_true = 0.98 * omega + (-1.0e-5) * omega**2
    a_lin, b_lin = fit_linear_frequency_dependent(omega, pred_linear_true)
    pred_linear = predict_linear_frequency_dependent(omega, a_lin, b_lin)
    assert np.allclose(pred_linear, pred_linear_true)

    pred_power_true = omega * (1.21 * omega ** (-0.028))
    a_pow, b_pow = fit_power_law_frequency_dependent(omega, pred_power_true)
    pred_power = predict_power_law_frequency_dependent(omega, a_pow, b_pow)
    assert np.allclose(pred_power, pred_power_true)


def test_bootstrap_and_compare_models_return_expected_keys():
    omega = np.array([958.8, 1173.4, 1188.9, 1298.9, 1513.2, 1756.0, 3140.4, 3226.9])
    nu = np.array([953.0, 1078.0, 1035.0, 1358.0, 1460.0, 1597.0, 2933.0, 2858.0])
    weights = np.array([1.0, 0.5, 0.5, 1.0, 1.0, 1.0, 0.7, 0.7])

    boot = bootstrap_constant_scale(omega, nu, weights=weights, n_boot=200, seed=7)
    assert {"k_mean", "k_std", "k_p025", "k_p500", "k_p975"} <= set(boot)
    assert boot["k_std"] >= 0.0

    models = compare_scaling_models(omega, nu, weights=weights)
    assert {"global_ls", "global_weighted_ls", "global_huber", "piecewise_region", "power_law"} <= set(models)
    assert models["global_ls"]["metrics"]["mae"] >= 0.0
