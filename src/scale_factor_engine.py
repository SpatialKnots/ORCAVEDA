from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


ArrayLike1D = Sequence[float] | np.ndarray


@dataclass(frozen=True)
class ScaleFitMetrics:
    mae: float
    rmse: float
    median_abs: float
    bias: float
    max_abs: float


def _as_float_array(values: ArrayLike1D, *, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} is empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _validate_pair_arrays(omega: ArrayLike1D, nu: ArrayLike1D) -> tuple[np.ndarray, np.ndarray]:
    omega_arr = _as_float_array(omega, name="omega")
    nu_arr = _as_float_array(nu, name="nu")
    if omega_arr.shape != nu_arr.shape:
        raise ValueError(f"omega and nu shape mismatch: {omega_arr.shape} vs {nu_arr.shape}")
    return omega_arr, nu_arr


def _validate_weights(weights: ArrayLike1D | None, n: int) -> np.ndarray:
    if weights is None:
        return np.ones(n, dtype=float)
    w = _as_float_array(weights, name="weights")
    if w.shape != (n,):
        raise ValueError(f"weights shape mismatch: expected {(n,)}, got {w.shape}")
    if np.any(w < 0):
        raise ValueError("weights must be non-negative")
    if not np.any(w > 0):
        raise ValueError("weights must contain at least one positive value")
    return w


def metrics(pred: ArrayLike1D, exp: ArrayLike1D) -> ScaleFitMetrics:
    pred_arr, exp_arr = _validate_pair_arrays(pred, exp)
    err = pred_arr - exp_arr
    return ScaleFitMetrics(
        mae=float(np.mean(np.abs(err))),
        rmse=float(np.sqrt(np.mean(err**2))),
        median_abs=float(np.median(np.abs(err))),
        bias=float(np.mean(err)),
        max_abs=float(np.max(np.abs(err))),
    )


def apply_scale(omega: ArrayLike1D, k: float) -> np.ndarray:
    omega_arr = _as_float_array(omega, name="omega")
    return float(k) * omega_arr


def fit_constant_scale(omega: ArrayLike1D, nu: ArrayLike1D, weights: ArrayLike1D | None = None) -> float:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    w = _validate_weights(weights, len(omega_arr))
    denom = np.sum(w * omega_arr**2)
    if denom <= 0:
        raise ValueError("cannot fit constant scale: zero denominator")
    return float(np.sum(w * omega_arr * nu_arr) / denom)


def _huber_loss(residual: np.ndarray, delta: float) -> np.ndarray:
    abs_r = np.abs(residual)
    return np.where(abs_r <= delta, 0.5 * residual**2, delta * (abs_r - 0.5 * delta))


def fit_constant_scale_robust(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    *,
    weights: ArrayLike1D | None = None,
    loss: str = "huber",
    delta: float = 20.0,
    initial_scale: float | None = None,
    search_half_width: float = 0.08,
    refinement_steps: int = 7,
    grid_size: int = 401,
) -> float:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    w = _validate_weights(weights, len(omega_arr))
    if loss not in {"huber", "mae", "mse"}:
        raise ValueError(f"unknown robust loss: {loss}")
    if delta <= 0:
        raise ValueError("delta must be positive")

    center = float(initial_scale) if initial_scale is not None else fit_constant_scale(omega_arr, nu_arr, weights=w)

    def objective(k: float) -> float:
        residual = k * omega_arr - nu_arr
        if loss == "mae":
            return float(np.sum(w * np.abs(residual)))
        if loss == "mse":
            return float(np.sum(w * residual**2))
        return float(np.sum(w * _huber_loss(residual, delta)))

    left = center - float(search_half_width)
    right = center + float(search_half_width)
    best_k = center
    best_obj = objective(best_k)
    for _ in range(refinement_steps):
        grid = np.linspace(left, right, int(grid_size))
        values = np.array([objective(float(k)) for k in grid], dtype=float)
        idx = int(np.argmin(values))
        best_k = float(grid[idx])
        best_obj = float(values[idx])
        span = (right - left) / 8.0
        left = best_k - span
        right = best_k + span
    return best_k


def fit_inverse_frequency_scale(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    weights: ArrayLike1D | None = None,
) -> float:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    if np.any(omega_arr <= 0) or np.any(nu_arr <= 0):
        raise ValueError("inverse-frequency fitting requires strictly positive frequencies")
    w = _validate_weights(weights, len(omega_arr))
    a = 1.0 / omega_arr
    b = 1.0 / nu_arr
    denom = np.sum(w * a**2)
    if denom <= 0:
        raise ValueError("cannot fit inverse-frequency scale: zero denominator")
    p = np.sum(w * a * b) / denom
    return float(1.0 / p)


def fit_piecewise_constant_scale(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    regions: Iterable[tuple[float, float]],
    *,
    weights: ArrayLike1D | None = None,
) -> list[dict[str, object]]:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    w = _validate_weights(weights, len(omega_arr))
    result: list[dict[str, object]] = []
    for lo, hi in regions:
        mask = (omega_arr >= float(lo)) & (omega_arr < float(hi))
        count = int(np.sum(mask))
        if count < 2:
            result.append({"range": (float(lo), float(hi)), "k": float("nan"), "n": count})
            continue
        k = fit_constant_scale(omega_arr[mask], nu_arr[mask], weights=w[mask])
        result.append({"range": (float(lo), float(hi)), "k": float(k), "n": count})
    return result


def apply_piecewise_scale(omega: ArrayLike1D, piecewise: Sequence[dict[str, object]]) -> np.ndarray:
    omega_arr = _as_float_array(omega, name="omega")
    scaled = np.full_like(omega_arr, np.nan, dtype=float)
    for item in piecewise:
        lo, hi = item["range"]
        k = float(item["k"])
        if not np.isfinite(k):
            continue
        mask = (omega_arr >= float(lo)) & (omega_arr < float(hi))
        scaled[mask] = k * omega_arr[mask]
    return scaled


def fit_linear_frequency_dependent(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    weights: ArrayLike1D | None = None,
) -> tuple[float, float]:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    w = _validate_weights(weights, len(omega_arr))
    X = np.column_stack([omega_arr, omega_arr**2])
    sqrt_w = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(X * sqrt_w, nu_arr * np.sqrt(w), rcond=None)
    return float(coef[0]), float(coef[1])


def predict_linear_frequency_dependent(omega: ArrayLike1D, a: float, b: float) -> np.ndarray:
    omega_arr = _as_float_array(omega, name="omega")
    return float(a) * omega_arr + float(b) * omega_arr**2


def fit_power_law_frequency_dependent(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    weights: ArrayLike1D | None = None,
) -> tuple[float, float]:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    if np.any(omega_arr <= 0) or np.any(nu_arr <= 0):
        raise ValueError("power-law fitting requires strictly positive frequencies")
    w = _validate_weights(weights, len(omega_arr))
    y = np.log(nu_arr / omega_arr)
    X = np.column_stack([np.ones_like(omega_arr), np.log(omega_arr)])
    sqrt_w = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(X * sqrt_w, y * np.sqrt(w), rcond=None)
    log_a, b = coef
    return float(np.exp(log_a)), float(b)


def predict_power_law_frequency_dependent(omega: ArrayLike1D, a: float, b: float) -> np.ndarray:
    omega_arr = _as_float_array(omega, name="omega")
    if np.any(omega_arr <= 0):
        raise ValueError("power-law prediction requires strictly positive frequencies")
    return omega_arr * (float(a) * omega_arr ** float(b))


def bootstrap_constant_scale(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    *,
    weights: ArrayLike1D | None = None,
    n_boot: int = 500,
    seed: int = 1,
) -> dict[str, float]:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    w = _validate_weights(weights, len(omega_arr))
    rng = np.random.default_rng(seed)
    n = len(omega_arr)
    ks = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        ks.append(fit_constant_scale(omega_arr[idx], nu_arr[idx], weights=w[idx]))
    ks_arr = np.asarray(ks, dtype=float)
    return {
        "k_mean": float(np.mean(ks_arr)),
        "k_std": float(np.std(ks_arr, ddof=1)),
        "k_p025": float(np.percentile(ks_arr, 2.5)),
        "k_p500": float(np.percentile(ks_arr, 50.0)),
        "k_p975": float(np.percentile(ks_arr, 97.5)),
    }


def compare_scaling_models(
    omega: ArrayLike1D,
    nu: ArrayLike1D,
    *,
    weights: ArrayLike1D | None = None,
    piecewise_regions: Iterable[tuple[float, float]] | None = None,
) -> dict[str, dict[str, object]]:
    omega_arr, nu_arr = _validate_pair_arrays(omega, nu)
    w = _validate_weights(weights, len(omega_arr))
    models: dict[str, dict[str, object]] = {}

    k = fit_constant_scale(omega_arr, nu_arr, weights=w)
    pred = apply_scale(omega_arr, k)
    models["global_ls"] = {"params": {"k": k}, "metrics": metrics(pred, nu_arr).__dict__}

    k_weighted = fit_constant_scale(omega_arr, nu_arr, weights=w)
    pred_weighted = apply_scale(omega_arr, k_weighted)
    models["global_weighted_ls"] = {"params": {"k": k_weighted}, "metrics": metrics(pred_weighted, nu_arr).__dict__}

    k_huber = fit_constant_scale_robust(omega_arr, nu_arr, weights=w, loss="huber")
    pred_huber = apply_scale(omega_arr, k_huber)
    models["global_huber"] = {"params": {"k": k_huber}, "metrics": metrics(pred_huber, nu_arr).__dict__}

    regions = list(piecewise_regions) if piecewise_regions is not None else [(0.0, 1200.0), (1200.0, 3000.0), (3000.0, float("inf"))]
    piecewise = fit_piecewise_constant_scale(omega_arr, nu_arr, regions, weights=w)
    piecewise_pred = apply_piecewise_scale(omega_arr, piecewise)
    piecewise_mask = np.isfinite(piecewise_pred)
    models["piecewise_region"] = {
        "params": {"regions": piecewise},
        "metrics": metrics(piecewise_pred[piecewise_mask], nu_arr[piecewise_mask]).__dict__ if np.any(piecewise_mask) else None,
    }

    a_pow, b_pow = fit_power_law_frequency_dependent(omega_arr, nu_arr, weights=w)
    pred_pow = predict_power_law_frequency_dependent(omega_arr, a_pow, b_pow)
    models["power_law"] = {"params": {"a": a_pow, "b": b_pow}, "metrics": metrics(pred_pow, nu_arr).__dict__}

    return models
