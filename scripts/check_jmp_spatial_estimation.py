#!/usr/bin/env python3
"""Independent full-covariance SciPy oracle against the production public C ABI.

The production estimator profiles total variance and optimizes range/fraction
with bounded Nelder-Mead. This oracle instead optimizes all physical covariance
parameters with analytic likelihood gradients and L-BFGS-B.
"""
import argparse
import ctypes
import datetime
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg, optimize
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
REFERENCES = [
    "https://www.jmp.com/support/help/en/19.2/jmp/statistical-details-for-spatial-and-temporal-variability.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/launch-the-mixed-model-personality.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/repeated-effects-covariance-parameter-estimates.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/random-effects-covariance-parameter-estimates.shtml",
]


class Abi:
    def __init__(self, path):
        self.path = path
        self.library = ctypes.CDLL(str(path))
        double_p = ctypes.POINTER(ctypes.c_double)
        self.library.ds_analytics_run.argtypes = [
            ctypes.c_char_p, ctypes.c_size_t, double_p, ctypes.c_size_t,
            ctypes.c_size_t, double_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_size_t),
        ]
        self.library.ds_analytics_run.restype = ctypes.c_void_p
        self.library.ds_analytics_free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        self.library.ds_analytics_free.restype = None
        self.calls = 0

    def __call__(self, operation, data, args=()):
        data = np.ascontiguousarray(data, dtype=np.float64)
        args = np.ascontiguousarray(args, dtype=np.float64)
        assert data.ndim == 2
        operation = operation.encode()
        size = ctypes.c_size_t()
        result = self.library.ds_analytics_run(
            operation, len(operation), data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            *data.shape, args.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            len(args), ctypes.byref(size),
        )
        self.calls += 1
        if not result:
            assert size.value == 0
            return None
        try:
            raw = ctypes.string_at(result, size.value)
            return json.loads(raw, parse_constant=lambda v: (_ for _ in ()).throw(AssertionError(v)))
        finally:
            self.library.ds_analytics_free(result, size.value)


def correlation(distance, kind, length):
    z = distance / length
    if kind == 0:
        return np.where(z < 1, 1 - 1.5 * z + 0.5 * z**3, 0.0)
    if kind == 2:
        return np.exp(-z**2)
    return np.exp(-z)


def dense(data, kind, reml, theta, nugget=True):
    """Raw-unit covariance likelihood, QR fixed effects and analytic score."""
    n = len(data)
    location, feature, subject, y = data[:, :1], data[:, 1:2], data[:, 2], data[:, 3]
    x = np.column_stack([np.ones(n), feature])
    same = subject[:, None] == subject[None, :]
    h = cdist(location, location)
    sill = np.exp(theta[0])
    noise = np.exp(theta[1]) if nugget else 0.0
    length = np.exp(theta[-1])
    r = correlation(h, kind, length) * same
    covariance = sill * r + noise * np.eye(n)
    chol = linalg.cholesky(covariance, lower=True)
    wx = linalg.solve_triangular(chol, x, lower=True)
    wy = linalg.solve_triangular(chol, y, lower=True)
    q, triangular = linalg.qr(wx, mode="economic")
    beta = linalg.solve_triangular(triangular, q.T @ wy)
    residual = y - x @ beta
    inverse = linalg.cho_solve((chol, True), np.eye(n))
    vi_x = inverse @ x
    coefficient_cov = linalg.inv(x.T @ vi_x)
    projection = inverse - vi_x @ coefficient_cov @ vi_x.T
    alpha = inverse @ residual
    df = n - x.shape[1] if reml else n
    value = 0.5 * (
        df * np.log(2 * np.pi) + 2 * np.log(np.diag(chol)).sum()
        + residual @ alpha
        + (2 * np.log(np.abs(np.diag(triangular))).sum() if reml else 0)
    )
    z = h / length
    if kind == 0:
        dr = np.where(z < 1, 1.5 * z * (1 - z**2), 0.0)
    elif kind == 2:
        dr = 2 * z**2 * r
    else:
        dr = z * r
    derivatives = [sill * r]
    if nugget:
        derivatives.append(noise * np.eye(n))
    derivatives.append(sill * dr * same)
    score_matrix = (projection if reml else inverse) - np.outer(alpha, alpha)
    gradient = np.array([0.5 * np.sum(score_matrix * v) for v in derivatives])
    return value, gradient, beta, coefficient_cov, covariance, alpha


def reference_fit(data, kind, reml, nugget=True, bounds=(0.01, 30.0)):
    variance = np.var(data[:, -1])
    var_bounds = (np.log(variance) - 20, np.log(variance) + 20)
    parameter_bounds = [var_bounds] * (2 if nugget else 1) + [np.log(bounds)]
    solutions = []
    def objective(theta):
        try:
            return dense(data, kind, reml, theta, nugget)[:2]
        except (ValueError, linalg.LinAlgError):
            # A trial point can imply a numerically singular Gaussian kernel;
            # it is outside this numerical likelihood domain, not a fitted model.
            return 1e100, np.zeros(len(theta))
    for length in [0.08, 0.3, 1.5, 6.0, 20.0]:
        for fraction in ([0.1, 0.65] if nugget else [0.0]):
            initial = [np.log(variance * (1 - fraction))]
            if nugget:
                initial.append(np.log(variance * fraction))
            initial.append(np.log(length))
            solutions.append(optimize.minimize(
                objective,
                initial, jac=True, method="L-BFGS-B", bounds=parameter_bounds,
                options={"ftol": 2e-14, "gtol": 2e-8, "maxiter": 2000},
            ))
    result = min(solutions, key=lambda r: r.fun)
    return result, dense(data, kind, reml, result.x, nugget)


def raw_information(data, kind, reml, fit, nugget=True):
    # Central differences of an analytic full-covariance likelihood gradient,
    # unlike production's second differences of a profiled scalar likelihood.
    theta = np.log([fit["partial_sill"], *([fit["nugget"]] if nugget else []), fit["correlation_length"]])
    columns = []
    step = 1e-5
    for j in range(len(theta)):
        plus, minus = theta.copy(), theta.copy()
        plus[j] += step
        minus[j] -= step
        columns.append((dense(data, kind, reml, plus, nugget)[1] - dense(data, kind, reml, minus, nugget)[1]) / (2 * step))
    info = np.column_stack(columns)
    info = (info + info.T) / 2
    jac = np.diag(np.exp(theta))
    return jac @ linalg.inv(info) @ jac


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    opts = parser.parse_args()
    run = Abi(opts.library)
    groups, cases = {}, []
    maximum_error = 0.0
    field_count = 0

    def equal(actual, expected, label, tolerance=3e-5):
        nonlocal maximum_error, field_count
        actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        assert actual.shape == expected.shape, (label, actual.shape, expected.shape)
        assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), label
        error = float(np.max(np.abs(actual - expected) / (1 + np.abs(expected)), initial=0))
        assert error <= tolerance, (label, error, actual, expected)
        maximum_error = max(maximum_error, error)
        field_count += actual.size
        groups[label] = groups.get(label, 0) + actual.size

    rng = np.random.default_rng(2910)
    location = np.arange(24, dtype=float) / 4
    subjects = (np.arange(24) >= 12).astype(float)
    features = np.sin(location)
    x = np.column_stack([np.ones(24), features])
    h = np.abs(location[:, None] - location[None, :])
    generating = 1.3 * np.exp(-h / 0.8) * (subjects[:, None] == subjects[None, :]) + 0.25 * np.eye(24)
    y = x @ np.array([2.0, -0.7]) + linalg.cholesky(generating, lower=True) @ rng.normal(size=24)
    data = np.column_stack([location, features, subjects, y])
    for kind in range(4):
        for reml in [False, True]:
            args = [1, 1, kind, int(reml), 1, 1500, 1e-9, 0.01, 30]
            fit = run("spatial_fit", data, args)
            assert fit is not None, (kind, reml)
            reference, ref = reference_fit(data, kind, reml)
            value, _, beta, beta_cov, covariance, alpha = ref
            equal(fit["log_likelihood"], -value, "full_covariance_likelihood", 3e-6)
            equal([fit["intercept"], *fit["coefficients"]], beta, "independently_optimized_trend", 5e-4)
            # Repeated optima/variance boundaries can make covariance parameters
            # nonidentifiable. Their returned diagnostics, not parameter equality,
            # govern inference; mean/covariance/prediction remain checked below.
            theta = np.log([fit["partial_sill"], max(fit["nugget"], 1e-200), fit["correlation_length"]])
            at_fit = dense(data, kind, reml, theta)
            equal([fit["intercept"], *fit["coefficients"]], at_fit[2], "profiled_trend")
            equal(fit["coefficient_covariance"], at_fit[3], "trend_coefficient_covariance")
            equal(fit["marginal_fitted"], x @ at_fit[2], "marginal_fit")
            process_covariance = at_fit[4] - fit["nugget"] * np.eye(len(data))
            equal(fit["conditional_fitted"], x @ at_fit[2] + process_covariance @ at_fit[5], "conditional_fit")
            total = fit["total_variance"]
            equal(total, fit["partial_sill"] + fit["nugget"], "total_partial_sill_nugget")
            if fit["covariance_parameters_identified"]:
                equal([fit["partial_sill"], fit["nugget"], fit["correlation_length"]], np.exp(reference.x), "estimated_covariance_parameters", 5e-4)
                vcov = raw_information(data, kind, reml, fit)
                equal(fit["covariance_parameter_covariance"], vcov, "observed_information_covariance", 3e-4)
                equal(fit["covariance_parameter_standard_errors"], np.sqrt(np.diag(vcov)), "observed_information_standard_errors", 2e-4)
                estimates = np.array([fit["partial_sill"], fit["nugget"], fit["correlation_length"]])
                se = np.sqrt(np.diag(vcov))
                equal(fit["covariance_parameter_intervals"], np.column_stack([estimates - 1.959963984540054 * se, estimates + 1.959963984540054 * se]), "explicit_Wald_intervals", 3e-4)
            else:
                assert fit["covariance_parameter_standard_errors"] is None
                assert fit["covariance_parameter_intervals"] is None
            query = np.array([[6.0, np.sin(6), 1], [0, 0, 0], [0, 0, 3]], dtype=float)
            qx = np.column_stack([np.ones(3), query[:, 1]])
            cross = fit["partial_sill"] * correlation(cdist(query[:, :1], data[:, :1]), kind, fit["correlation_length"])
            cross *= query[:, 2:3] == data[None, :, 2]
            vi_cross = linalg.solve(at_fit[4], cross.T, assume_a="pos")
            delta = qx - vi_cross.T @ x
            var = fit["partial_sill"] - np.sum(cross * vi_cross.T, axis=1) + np.einsum("ij,jk,ik->i", delta, at_fit[3], delta)
            prediction = run("spatial_fit_predict", query, fit["predict_parameters"])
            assert prediction is not None
            equal(prediction["predictions"], qx @ at_fit[2] + cross @ at_fit[5], "saved_universal_kriging_mean")
            equal(prediction["latent_variances"], var, "saved_latent_variance")
            equal(prediction["observation_variances"], var + fit["nugget"], "saved_observation_variance")
            if kind in [1, 2]:
                altered = data.copy()
                altered[:, 0] *= 13
                altered[:, -1] = 7 + 3 * altered[:, -1]
                changed = run("spatial_fit", altered, [*args[:7], 0.13, 390])
                assert changed is not None
                equal(changed["correlation_length"] / 13, fit["correlation_length"], "coordinate_unit_equivariance", 4e-4)
                equal(changed["total_variance"] / 9, total, "response_variance_equivariance", 4e-4)
                equal((changed["intercept"] - 7) / 3, fit["intercept"], "response_location_equivariance", 4e-4)
            cases.append({"kind": kind, "reml": reml, "input": data.tolist(), "args": args, "expected": {"intercept": float(beta[0]), "coefficients": beta[1:].tolist(), "partial_sill": float(np.exp(reference.x[0])), "nugget": float(np.exp(reference.x[1])), "correlation_length": float(np.exp(reference.x[2])), "log_likelihood": float(-value)}, "actual_summary": {k: fit[k] for k in ["intercept", "coefficients", "partial_sill", "nugget", "correlation_length", "log_likelihood"]}, "identified": fit["covariance_parameters_identified"], "reference_gradient_inf": float(np.max(np.abs(reference.jac))), "reference_success": bool(reference.success)})
    for kind in range(4):
        for reml in [False, True]:
            args = [1, 1, kind, int(reml), 0, 1500, 1e-9, 0.01, 30]
            fit = run("spatial_fit", data, args)
            assert fit is not None, (kind, reml, "no nugget")
            reference, ref = reference_fit(data, kind, reml, False)
            equal(fit["log_likelihood"], -ref[0], "no_nugget_full_likelihood", 3e-6)
            equal([fit["intercept"], *fit["coefficients"]], ref[2], "no_nugget_optimized_trend", 5e-4)
            assert fit["nugget"] == 0
            if fit["covariance_parameters_identified"]:
                equal([fit["partial_sill"], fit["correlation_length"]], np.exp(reference.x), "no_nugget_estimated_covariance", 5e-4)
                cov = raw_information(data, kind, reml, fit, False)
                equal(fit["covariance_parameter_covariance"], cov, "no_nugget_information_covariance", 3e-4)
                equal(fit["covariance_parameter_standard_errors"], np.sqrt(np.diag(cov)), "no_nugget_standard_errors", 3e-4)
            if kind == 3:
                equal(fit["kernel_parameter"], np.exp(-1 / fit["correlation_length"]), "Power_parameterization", 1e-12)
            cases.append({"kind": kind, "reml": reml, "input": data.tolist(), "args": args, "expected": {"intercept": float(ref[2][0]), "coefficients": ref[2][1:].tolist(), "partial_sill": float(np.exp(reference.x[0])), "nugget": 0, "correlation_length": float(np.exp(reference.x[1])), "log_likelihood": float(-ref[0])}, "identified": fit["covariance_parameters_identified"], "reference_gradient_inf": float(np.max(np.abs(reference.jac))), "reference_success": bool(reference.success)})
    limited = run("spatial_fit", data, [1, 1, 1, 0, 1, 1, 1e-9, 0.01, 30])
    assert limited is not None and not limited["optimization_converged"]
    assert limited["covariance_parameter_intervals"] is None
    # Independent closed-form variance degrees of freedom and unseen-subject SE.
    for reml, variance in [(0, 14 / 3), (1, 7)]:
        fit = run("spatial_fit", [[0, 1, 1], [0, 2, 2], [0, 3, 6]], [1, 0, 1, reml, 0, 100, 1e-8, 0.1, 2])
        assert fit is not None and not fit["covariance_parameters_identified"]
        equal(fit["total_variance"], variance, "closed_form_iid_variance")
        pred = run("spatial_fit_predict", [[0, 4]], fit["predict_parameters"])
        equal(pred["observation_variances"], [variance * 4 / 3], "new_subject_iid_prediction_variance")
    assert run("spatial_fit", [[0, 0, 2], [1, 0, 2], [2, 0, 2]], [1, 0, 1, 1, 1]) is None
    assert run("spatial_fit_predict", [[0, 0]], []) is None
    assert run("spatial_fit", [[0, 1, 0, 1], [1, 1, 0, 2], [2, 1, 0, 4]], [1, 1, 1, 1, 1]) is None
    assert run("spatial_fit", [[0, 0.5, 1], [1, 0, 2], [2, 0, 4]], [1, 0, 1, 1, 1]) is None
    assert run("spatial_fit", [[0, 0, np.nan], [1, 0, 2], [2, 0, 4]], [1, 0, 1, 1, 1]) is None
    saved = fit["predict_parameters"]
    for index, value in [(0, 2), (1, 0), (2, 9), (3, 97), (5, -1), (6, -1), (7, 0), (8, np.nan)]:
        bad = list(saved)
        bad[index] = value
        assert run("spatial_fit_predict", [[0, 4]], bad) is None
    assert run("spatial_fit_predict", [[0, 4]], saved[:-1]) is None
    assert run("spatial_fit_predict", [[0, 4.5]], saved) is None
    record = {
        "schema_version": 1, "status": "PASS", "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "layer": "production ds_analytics_run/ds_analytics_free public C ABI from actual production static archive; native SQL recorded separately",
        "seed": 2910, "calls": run.calls, "numeric_fields": field_count, "maximum_scaled_error": maximum_error,
        "groups": groups, "fixtures": cases, "official_references": REFERENCES,
        "independent_method": "SciPy L-BFGS-B all-physical-covariance-parameter optimization with analytic score, dense Cholesky+QR; information from differenced analytic score",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "library_sha256": hashlib.sha256(opts.library.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((ROOT / "crates/analytics/src/spatial_estimation.rs").read_bytes()).hexdigest(),
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "limitations": ["bounded isotropic ML/REML and plug-in prediction only", "interior untruncated Wald approximation, not JMP bounded Satterthwaite intervals", "no anisotropy, G-side effects or small-sample df correction", "no global optimum certificate and no installed-JMP head-to-head comparison"],
    }
    opts.output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: record[k] for k in ["status", "calls", "numeric_fields", "maximum_scaled_error", "groups"]}, indent=2))


if __name__ == "__main__":
    main()
