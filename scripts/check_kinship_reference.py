#!/usr/bin/env python3
"""Independent dense-Cholesky/QR genomic profile likelihood references."""
import argparse
import datetime
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

import numpy as np
import scipy
from scipy.linalg import solve_triangular
from scipy.optimize import differential_evolution, minimize_scalar


def decimal_prediction(data, params):
    """Evaluate saved-model algebra independently, without binary64 intermediates.

    Decimal.from_float preserves the actual supplied binary64 values.  The
    precision also retains the smallest residuals in the cancellation fixtures;
    no production exponent-product or compensated-sum algorithm is reused.
    """
    with localcontext() as context:
        context.prec = 1200
        value = lambda x: Decimal.from_float(float(x))
        p = int(params[2])
        saved = list(map(value, params))
        output, marginal = [], []
        for observation in data:
            row = list(map(value, observation))
            fixed = saved[6] + sum(
                ((row[j]-saved[8+j])/saved[8+2*p+j]-saved[8+p+j])*saved[8+3*p+j]
                for j in range(p)
            )
            random = sum(row[p+j]*weight/saved[7]
                         for j, weight in enumerate(saved[8+4*p:]))
            marginal.append(float(saved[3]+saved[5]*(saved[4]+fixed)))
            output.append(float(saved[3]+saved[5]*(saved[4]+fixed+random)))
        return output, marginal


def boundary_references(request):
    """Independent deterministic regressions, separate from the 180 random cases."""
    cases = fields = rejected = 0
    errors, reports = [], []

    def run(name, op, data, params=(), expected_error=False):
        nonlocal cases, rejected
        result = request(op, data, params, expected_error=expected_error)
        cases += 1
        rejected += int(expected_error)
        reports.append({"case": name, "operation": op, "expected_error": expected_error})
        return result

    def close(actual, expected, relative=2e-12, absolute=0.):
        nonlocal fields
        actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
        # Two representable ulps permit correct rounding at the subnormal edge.
        ulps = max((2*math.ulp(float(x)) for x in expected.ravel()), default=0.)
        np.testing.assert_allclose(actual, expected, rtol=relative, atol=max(absolute, ulps))
        fields += expected.size
        errors.append(float(np.max(np.abs(actual-expected)/(1+np.abs(expected)))))

    def score(name, query, params, expected=None):
        predictions, marginal = decimal_prediction(query, params)
        if expected is not None:
            close(predictions, expected)
        result = run(name, "kinship_predict", query, params)
        close(result["predictions"], predictions)
        close(result["marginal_predictions"], marginal)
        return result

    n = 6
    block = np.equal.outer(np.arange(n)//2, np.arange(n)//2).astype(float)
    z = np.array([1., -1., 0., 0., 0., 0.])/np.sqrt(2.)
    kernel = block-1e-12*np.outer(z, z)
    y = np.repeat(np.arange(3, dtype=float), 2)+1e-5*z
    result = run("roundoff_negative_eigenvalue_fit", "kinship_model",
                 np.column_stack((y, kernel)), [0, 1, .999999999999, 128, 1e-12])
    assert result["converged"] and not result["upper_bound_active"]
    normalized = kernel/(np.trace(kernel)/n)
    eigenvalues, eigenvectors = np.linalg.eigh(normalized)
    assert -1e-11 < eigenvalues[0] < 0.
    clipped = np.maximum(eigenvalues, 0.)
    h = result["normalized_genetic_fraction"]
    diagonal = 1-h+h*clipped
    rotated_y = eigenvectors.T@y
    rotated_one = eigenvectors.T@np.ones(n)
    beta = np.sum(rotated_one*rotated_y/diagonal)/np.sum(rotated_one**2/diagonal)
    expected = eigenvectors@(h*clipped*(rotated_y-rotated_one*beta)/diagonal)
    close(result["genetic_values"], expected, relative=2e-8, absolute=2e-9)
    # The negative direction is an exact null contrast after PSD clipping.
    close(np.dot(z, result["genetic_values"]), 0., absolute=2e-10)
    psd_kernel = (eigenvectors*clipped)@eigenvectors.T
    saved = result["predict_parameters"]
    close(saved[5]*(psd_kernel@np.array(saved[8:])), result["genetic_values"],
          relative=2e-8, absolute=2e-9)
    scored = score("roundoff_negative_eigenvalue_training_score", kernel, saved)
    close(scored["predictions"], result["conditional_fitted"], absolute=2e-9)

    # All validation fixtures have supported n=4, so matrix validation is reached.
    valid_kernel = np.equal.outer(np.arange(4)//2, np.arange(4)//2).astype(float)+.2*np.eye(4)
    run("supported_dimension_psd_control", "kinship_model",
        np.column_stack(([0., 2., 4., 6.], valid_kernel)), [0])
    asymmetric = valid_kernel.copy()
    asymmetric[0, 1] = np.nextafter(asymmetric[0, 1], np.inf)
    run("nonexact_asymmetry_rejected", "kinship_model",
        np.column_stack(([0., 2., 4., 6.], asymmetric)), [0], expected_error=True)
    indefinite = np.eye(4)
    indefinite[0, 1] = indefinite[1, 0] = 2.
    run("supported_dimension_indefinite_rejected", "kinship_model",
        np.column_stack(([0., 2., 4., 6.], indefinite)), [0], expected_error=True)
    negative_diagonal = np.eye(4)
    negative_diagonal[0, 0] = -1.
    run("supported_dimension_negative_diagonal_rejected", "kinship_model",
        np.column_stack(([0., 2., 4., 6.], negative_diagonal)), [0], expected_error=True)

    tiny = 1e-320
    x = np.tile([0., tiny], 3)
    for name, last in [("subnormal_feature_nonzero_slope", 1.+1e-14),
                       ("subnormal_feature_zero_slope", 1.)]:
        feature_y = [0., 2., 2., 0., 1., last]
        fit = run(name+"_fit", "kinship_model", np.column_stack((feature_y, x, block)), [1, 1])
        assert fit["converged"] and fit["normalized_genetic_fraction"] == 0.
        predicted = score(name+"_score", [[1., 0., 0., 0., 0., 0., 0.]], fit["predict_parameters"])
        if last == 1.:
            close(predicted["predictions"], [1.])
        else:
            assert predicted["predictions"][0] > 1e300

    fit = run("subnormal_kernel_zero_genetic_fit", "kinship_model",
              np.column_stack(([0., 2., 0., 2., 0., 2.], block*tiny)), [0, 1])
    assert fit["converged"] and fit["normalized_genetic_fraction"] == 0.
    score("subnormal_kernel_zero_genetic_score", [[1.]*6], fit["predict_parameters"], [1.])

    # These valid version-1 models isolate scoring arithmetic from estimation.
    # Four training weights suffice; joint query PSD cannot be checked by this API.
    scalar_cases = [
        ("ratio_overflows_before_finite_physical_product", 1e-100, 1e-320,
         [1e-200, 0., 0., 0.], [1., 0., 0., 0.]),
        ("early_product_underflows_before_finite_ratio", 1e-200, 1e-100,
         [1e-200, 0., 0., 0.], [1e150, 0., 0., 0.]),
        ("small_finite_residual_of_large_kernel_terms", 1e-100, 1e-320,
         [1., 1e-300, -1., 0.], [1., 1., 1., 0.]),
        ("finite_sum_of_unrepresentable_opposite_physical_terms", 1e100, 1e-320,
         [1., 1e-300, -1., 0.], [1., 1., 1., 0.]),
        ("subnormal_saved_weight_survives_unit_products", 1e200, 1e-200,
         [float.fromhex("0x0.0000000000001p-1022"), 0., 0., 0.], [1e100, 0., 0., 0.]),
    ]
    for name, response_scale, kernel_scale, weights, query in scalar_cases:
        params = [1., 4., 0., 0., 0., response_scale, 0., kernel_scale, *weights]
        score(name, [query], params)

    # Include cancellation among mean and feature terms, not only random terms.
    params = [1., 4., 1., 1e300, 0., 1e300, -1., 1.,
              0., 0., 1., 1e-300, 0., 0., 0., 0.]
    score("mean_feature_physical_cancellation", [[1., 0., 0., 0., 0.]], params)
    params = [1., 4., 1., 1., 0., 1., 0., 1.,
              0., 0., tiny, 0., 0., 0., 0., 0.]
    score("zero_feature_slope_with_overflowing_standardized_query",
          [[1., 0., 0., 0., 0.]], params, [1.])

    return {"cases": cases, "numeric_fields": fields, "rejected_inputs": rejected,
            "maximum_scaled_error": max(errors),
            "reference_methods": "NumPy symmetric eigendecomposition with independently clipped PSD kernel, direct covariance spectral GLS/BLUP; Decimal.from_float at 1200 digits for exact supplied-binary64 saved-model algebra and finite cancellation; valid-dimension matrix rejection checks",
            "case_catalog": reports}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wrapper", required=True, type=Path)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    process = subprocess.Popen([str(args.wrapper.resolve())], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    rng = np.random.default_rng(293015)
    cases = fields = 0
    errors = []
    flags = {"lower_boundary": 0, "upper_boundary": 0, "interior": 0}

    def request(op, data, params=(), expected_error=False):
        data = np.asarray(data)
        n, p = data.shape
        process.stdin.write(" ".join([op, str(n), str(p), str(len(params)), *map(repr, map(float, data.ravel())), *map(repr, map(float, params))])+"\n")
        process.stdin.flush()
        line = process.stdout.readline()
        assert line, (op, "wrapper closed its output", data, params)
        if expected_error:
            assert line.startswith("ERROR:"), (op, line, data, params)
            return line.removeprefix("ERROR:").strip()
        assert not line.startswith("ERROR:"), (op, line, data, params)
        def invalid(value):
            raise AssertionError(value)
        return json.loads(line, parse_constant=invalid)

    def run(op, data, params=()):
        nonlocal cases
        result = request(op, data, params)
        cases += 1
        return result

    def close(a, b, relative=3e-5, absolute=1e-8):
        nonlocal fields
        np.testing.assert_allclose(a, b, rtol=relative, atol=absolute)
        aa, bb = np.asarray(a), np.asarray(b)
        fields += bb.size
        errors.append(float(np.max(np.abs(aa-bb)/(1+np.abs(bb)))))

    try:
        for trial in range(45):
            n, p = int(rng.integers(7, 22)), int(rng.integers(0, 4))
            markers = int(rng.integers(3, 15))
            Z = rng.normal(size=(n+3, markers))
            kernel_scale = float(10**rng.uniform(-3, 3))
            K = (Z[:n]@Z[:n].T/markers+.2*np.eye(n))*kernel_scale
            cross = Z[n:]@Z[:n].T/markers*kernel_scale
            X = np.column_stack((np.ones(n+3), rng.normal(size=(n+3, p))))
            beta = rng.normal(size=p+1)
            y = X[:n]@beta+np.linalg.cholesky(K/kernel_scale+float(10**rng.uniform(-2, 2))*np.eye(n))@rng.normal(size=n)
            normalized = K/np.trace(K)*n
            upper = .99999999
            data = np.column_stack((y, X[:n, 1:], K))
            query = np.column_stack((X[n:, 1:], cross))
            for reml in [0, 1]:
                df = n-(p+1 if reml else 0)
                def profile(h, detail=False):
                    covariance = (1-h)*np.eye(n)+h*normalized
                    L = np.linalg.cholesky(covariance)
                    white_x = solve_triangular(L, X[:n], lower=True)
                    white_y = solve_triangular(L, y, lower=True)
                    q, r = np.linalg.qr(white_x, mode="reduced")
                    coef = np.linalg.solve(r, q.T@white_y)
                    residual = white_y-white_x@coef
                    variance = residual@residual/df
                    value = .5*(df*(np.log(2*np.pi)+1+np.log(variance))+
                                 2*np.log(np.diag(L)).sum()+
                                 (2*np.log(np.abs(np.diag(r))).sum() if reml else 0))
                    return (value, coef, variance) if detail else value
                # An independent fine mesh brackets all resolved local minima;
                # bounded Brent uses a different solve and stopping mechanism.
                mesh = np.linspace(0, upper, 301)
                values = np.array([profile(h) for h in mesh])
                candidates = [(values[0], 0.), (values[-1], upper)]
                for i in range(1, len(mesh)-1):
                    if values[i] <= min(values[i-1], values[i+1]):
                        oracle = minimize_scalar(profile, bounds=(mesh[i-1], mesh[i+1]), method="bounded", options={"xatol": 1e-12})
                        assert oracle.success
                        candidates.append((oracle.fun, oracle.x))
                # A separate population search also resolves narrow extrema
                # between the last mesh point and the parameter boundary.
                global_check=differential_evolution(lambda v:profile(v[0]),[(0.,upper)],rng=np.random.default_rng(293015+2*trial+reml),tol=1e-11,atol=1e-12,polish=True)
                assert global_check.success,global_check.message
                candidates.append((global_check.fun,float(global_check.x[0])))
                objective, h = min(candidates)
                _, coefficient, variance = profile(h, True)
                genetic, residual = variance*h/(np.trace(K)/n), variance*(1-h)
                result = run("kinship_model", data, [p, reml, upper, 128, 1e-11])
                if abs(result["log_likelihood"]+objective)>1e-8+2e-8*abs(objective):
                    diagnostic={"trial":trial,"reml":reml,"oracle_h":h,"returned_h":result["normalized_genetic_fraction"],"returned_ll":result["log_likelihood"],"oracle_ll":-objective,"dense_ll_at_returned_h":-profile(result["normalized_genetic_fraction"]),"data":data.tolist(),"args":[p,reml,upper,128,1e-11]}
                    Path(str(args.wrapper)+"-failure.json").write_text(json.dumps(diagnostic,indent=2)+"\n")
                    print(json.dumps({k:v for k,v in diagnostic.items() if k!="data"}),flush=True)
                close(result["log_likelihood"], -objective, relative=2e-8)
                close([result["intercept"], *result["coefficients"]], coefficient)
                close(result["genetic_variance"], genetic, absolute=2e-7)
                close(result["residual_variance"], residual, absolute=2e-7)
                # Verify BLUP and original-unit covariance against a direct solve
                # at the actual returned estimates, separate from optimization.
                g, e = result["genetic_variance"], result["residual_variance"]
                V = g*K+e*np.eye(n)
                direct_beta = np.linalg.solve(X[:n].T@np.linalg.solve(V, X[:n]), X[:n].T@np.linalg.solve(V, y))
                alpha = np.linalg.solve(V, y-X[:n]@direct_beta)
                close(result["genetic_values"], g*K@alpha, relative=2e-7)
                close(result["conditional_fitted"], X[:n]@direct_beta+g*K@alpha, relative=2e-7)
                predicted = run("kinship_predict", query, result["predict_parameters"])
                close(predicted["predictions"], X[n:]@direct_beta+g*cross@alpha, relative=2e-7)
                key = "upper_boundary" if result["upper_bound_active"] else "lower_boundary" if result["boundary"] else "interior"
                flags[key] += 1
                if result["upper_bound_active"]:
                    assert not result["converged"]
        boundary_results = boundary_references(request)
        repo = Path(__file__).resolve().parents[1]
        paths = ["crates/analytics/src/kinship.rs", "crates/analytics/tests/candidates/kinship.rs", "scripts/check_kinship_reference.py"]
        if (repo/"crates/analytics/tests/kinship.rs").exists():
            paths.append("crates/analytics/tests/kinship.rs")
        evidence = {"schema_version": 1, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    "status": "PASS", "seed": 293015, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                    "boundary": "Independent oracle uses an isolated wrapper importing exact production kinship.rs; separate witnessed public execute16tests GREEN. C ABI/native SQL rebuild remains a separate pending gate.",
                    "cases": cases, "numeric_fields": fields, "maximum_scaled_error": max(errors), "fit_boundary_diagnostics": flags,
                    "boundary_regressions": boundary_results,
                    "source_sha256": {s: hashlib.sha256((repo/s).read_bytes()).hexdigest() for s in paths},
                    "reference_methods": "Independent full dense covariance Cholesky whitening, NumPy QR, 301-point+bounded Brent plus independently seeded differential-evolution likelihood search, direct mixed-model covariance equations and cross-kernel BLUP",
                    "limitations": "45 randomized datasets, both ML/REML, plus separate deterministic arithmetic/PSD boundary regressions; no installed JMP oracle, no global optimum theorem, no GWAS or variance-component uncertainty. Arbitrary saved-parameter probes validate scoring algebra, not estimation or a jointly valid query/training covariance.",
                    "harness_python": Path(__file__).read_text(),
                    "wrapper_source": Path(str(args.wrapper)+".rs").read_text() if Path(str(args.wrapper)+".rs").exists() else None}
        if args.evidence:
            if args.evidence.exists():
                previous = json.loads(args.evidence.read_text())
                if "TDD" in previous:
                    evidence["TDD"] = previous["TDD"]
            evidence["api_sha256"] = hashlib.sha256((repo/"tasks/jmp-kinship-api.json").read_bytes()).hexdigest()
            evidence["registry_sha256"] = hashlib.sha256((repo/"tasks/jmp-kinship-registry.inc").read_bytes()).hexdigest()
            args.evidence.write_text(json.dumps(evidence, indent=2)+"\n")
        print(json.dumps({k: v for k, v in evidence.items() if k not in ["source_sha256", "harness_python", "wrapper_source"]}, indent=2))
    finally:
        process.stdin.close()
        process.wait(timeout=5)


if __name__ == "__main__":
    main()
