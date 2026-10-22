#!/usr/bin/env python3
"""Independent ordinal ML/observed-information oracle at actual C ABI and SQL.

This checks a bounded native ML convention, not default JMP Lasso or its
undocumented sandwich details. SciPy optimization and original-unit numerical
Hessians are independent of production's analytic constrained Newton method.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import scipy
from scipy import optimize, special, stats

from check_jmp_censored import Abi
from check_jmp_continuation import Runner, expression, strict_json

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = [
    "https://www.jmp.com/support/help/en/19.2/jmp/specify-a-distribution.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/ordinal-responses.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/parameter-estimates-for-original-predictors.shtml",
]


def fixtures():
    yield "intercept_exact", [[30, 1], [10, 1], [20, 2]]
    yield "intercept_40_35_25", [[10, 40], [20, 35], [30, 25]]
    grouped = [[x, label, count] for x, counts in [(-1, [2, 8, 10]), (0, [5, 10, 5]), (1, [10, 8, 2])]
               for label, count in zip([10, 20, 30], counts)]
    yield "known_common_slope", grouped
    yield "shifted_predictor", [[10+3*x, label, count] for x, label, count in grouped]
    groups = [(-1, .4, [4, 7, 5, 2]), (-.6, -.7, [3, 6, 8, 3]),
              (-.1, .2, [6, 7, 6, 5]), (.3, .9, [8, 7, 5, 4]),
              (.8, -.4, [10, 9, 6, 3]), (1.2, .6, [13, 7, 5, 2])]
    yield "four_levels_two_features", [[x, z, label, n] for x, z, counts in groups
                                      for label, n in zip([-5, 0, 20, 90], counts)]


def parts(data):
    a = np.asarray(data, dtype=float)
    a = a[a[:, -1] > 0]
    labels = np.unique(a[:, -2])
    y = np.searchsorted(labels, a[:, -2])
    return a[:, :-2], y, a[:, -1], labels


def probabilities(theta, x, k):
    cuts, beta = theta[:k-1], theta[k-1:]
    cdf = special.expit(cuts[None, :]+(x@beta)[:, None])
    return np.diff(np.column_stack([np.zeros(len(x)), cdf, np.ones(len(x))]), axis=1)


def likelihood(theta, x, y, w, k):
    if np.any(np.diff(theta[:k-1]) <= 0):
        return 1e100
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        ps = probabilities(theta, x, k)[np.arange(len(y)), y]
        value = -float(w@np.log(ps))
    return value if math.isfinite(value) else 1e100


def hessian(fun, theta):
    # Independent original-unit Richardson central differences, including every
    # threshold/slope covariance entry; no production Hessian/score is reused.
    d = len(theta)
    result = np.empty((d, d))
    f0 = fun(theta)
    for i in range(d):
        for j in range(i+1):
            values = []
            for divisor in [1, 2]:
                hi, hj = 4e-4*(1+abs(theta[i]))/divisor, 4e-4*(1+abs(theta[j]))/divisor
                ei, ej = np.eye(d)[i]*hi, np.eye(d)[j]*hj
                if i == j:
                    value = (fun(theta+ei)-2*f0+fun(theta-ei))/hi**2
                else:
                    value = (fun(theta+ei+ej)-fun(theta+ei-ej)-fun(theta-ei+ej)+fun(theta-ei-ej))/(4*hi*hj)
                values.append(value)
            result[i, j] = result[j, i] = (4*values[1]-values[0])/3
    return result


def reference(data):
    x, y, w, labels = parts(data)
    k = len(labels)
    counts = np.bincount(y, weights=w, minlength=k)
    cumulative = np.cumsum(counts)[:-1]/sum(w)
    cuts = special.logit(cumulative)
    fun = lambda theta: likelihood(theta, x, y, w, k)
    if x.shape[1] == 0:
        theta = cuts
        pdf = cumulative*(1-cumulative)
        cov = (np.minimum.outer(cumulative, cumulative)-np.outer(cumulative, cumulative))/sum(w)/np.outer(pdf, pdf)
        method = "Closed-form empirical multinomial cumulative logits; exact delta covariance at observed occurrence fractions."
        optimizer = None
    else:
        def decode(a):
            return np.r_[a[0], a[0]+np.cumsum(np.exp(a[1:k-1])), a[k-1:]]
        start = np.r_[cuts[0], np.log(np.diff(cuts)), np.zeros(x.shape[1])]
        candidates = []
        for offset in [-.3, 0., .3]:
            a = start.copy()
            a[0] += offset
            first = optimize.minimize(lambda a: fun(decode(a)), a, method="Powell",
                                      options={"xtol": 1e-11, "ftol": 1e-13, "maxiter": 1500})
            second = optimize.minimize(lambda a: fun(decode(a)), first.x, method="BFGS", jac="3-point",
                                       options={"gtol": 1e-7, "maxiter": 600})
            candidates.extend([first, second])
        best = min(candidates, key=lambda r: r.fun)
        theta = decode(best.x)
        info = hessian(fun, theta)
        assert np.min(np.linalg.eigvalsh(info)) > 0
        cov = np.linalg.inv(info)
        optimizer = {"success": bool(best.success), "message": str(best.message)}
        method = "SciPy Powell plus independent finite-difference BFGS on log-gap coordinates; raw-coordinate Richardson Hessian."
    transform = np.eye(len(theta))
    for j in range(1, k-1):
        transform[j, j-1] = -1
    return {"labels": labels.tolist(), "theta": theta.tolist(), "log_likelihood": -fun(theta),
            "covariance": cov.tolist(), "genreg_theta": (transform@theta).tolist(),
            "genreg_covariance": (transform@cov@transform.T).tolist(),
            "reference_method": method, "optimizer": optimizer}


class Checks:
    def __init__(self):
        self.count = 0
    def close(self, actual, expected, rtol=4e-5, atol=3e-6):
        actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        assert actual.shape == expected.shape, (actual.shape, expected.shape)
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)
        self.count += expected.size
    def require(self, predicate, detail):
        assert predicate, detail
        self.count += 1


def edge_contracts(calls, nullable, checks):
    tiny = [1, 0, 3, .95, 1, 10, 20, 30, -math.log(3), math.log(3), 1e-320, 0, 0, 1e-320]
    out = {"tiny_valid_covariance": {}, "asymmetric_covariance_rejected": {}, "nonconverged_inference": {}}
    for layer, call in calls.items():
        result = call("genreg_ordinal_predict", [[.43]], tiny)
        # sqrt before products retains a representable SD whose variance underflows.
        unit = math.sqrt(1e-320)*.1875
        checks.close(result["probability_standard_errors"], [[unit, math.sqrt(2)*unit, unit]], rtol=2e-12, atol=0)
        out["tiny_valid_covariance"][layer] = result
        rejected = []
        for covariance in [[1e-200, 1e-201, 0, 1e-200], [1.7e308, 1.1e308, 1e308, 1.7e308]]:
            malformed = tiny[:10]+covariance
            checks.require(nullable[layer]("genreg_ordinal_predict", [[.43]], malformed) is None,
                           (layer, "accepted asymmetric covariance", covariance))
            rejected.append({"covariance": covariance, "returned_null": True})
        out["asymmetric_covariance_rejected"][layer] = rejected
        data = next(d for name, d in fixtures() if name == "known_common_slope")
        fit = call("genreg_ordinal_fit", data, [0, 1e-7, .95])
        checks.require(not fit["converged"] and fit["covariance"] is None, (layer, fit))
        prediction = call("genreg_ordinal_predict", [[0, .43]], fit["predict_parameters"])
        checks.require(prediction["probability_standard_errors"] == [[None, None, None]], (layer, prediction))
        out["nonconverged_inference"][layer] = {"fit": fit, "prediction": prediction}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--duckdb")
    parser.add_argument("--extension", type=Path)
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not (args.library or args.duckdb or args.reference_only) or bool(args.duckdb) != bool(args.extension):
        parser.error("select actual library and/or matching duckdb+extension,or explicit reference-only")
    calls, nullable = {}, {}
    if args.library:
        abi = Abi(args.library)
        calls["public_c_abi"], nullable["public_c_abi"] = abi, abi.nullable
    if args.duckdb:
        runner = Runner(args.duckdb, args.extension)
        calls["native_duckdb_sql"] = lambda op, d, a: runner.one(expression(op, d, a))
        def sql_nullable(op, data, params):
            rows = runner.query("SELECT "+expression(op, data, params)+" AS result;")
            assert len(rows) == 1
            value = rows[0]["result"]
            return None if value is None else strict_json(value)
        nullable["native_duckdb_sql"] = sql_nullable
    checks = Checks()
    result = {"schema": "jmp19.2-ordinal-independent-reference-v1", "primary_sources": PRIMARY,
              "scope": "Native unpenalized ML/observed-information alternative; not JMP default Lasso or verified sandwich.",
              "reference_versions": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__},
              "platform": platform.platform(), "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
              "boundaries": list(calls), "reference_only": not bool(calls), "cases": []}
    for name, data in fixtures():
        ref = reference(data)
        case = {"name": name, "data": data, "reference": ref, "results": {}}
        for layer, call in calls.items():
            actual = call("genreg_ordinal_fit", data, [600, 1e-7, .95])
            theta = np.r_[actual["cutpoints"], actual["coefficients"]]
            checks.require(actual["converged"] and actual["inference_valid"], (name, layer, actual))
            checks.close(actual["ordered_labels"], ref["labels"], rtol=0, atol=0)
            checks.close(theta, ref["theta"])
            checks.close(actual["log_likelihood"], ref["log_likelihood"], rtol=1e-9, atol=1e-7)
            checks.close(actual["covariance"], ref["covariance"], rtol=7e-4, atol=5e-6)
            checks.close(actual["genreg_covariance"], ref["genreg_covariance"], rtol=7e-4, atol=5e-6)
            x, y, w, labels = parts(data)
            k, p = len(labels), x.shape[1]
            score = optimize._numdiff.approx_derivative(lambda t: likelihood(t, x, y, w, k)/sum(w), theta, method="3-point")
            checks.require(np.max(np.abs(score)) <= 2e-7, (name, layer, "independent score", score))
            cov, gcov = np.array(ref["covariance"]), np.array(ref["genreg_covariance"])
            z = stats.norm.ppf(.975)
            se = np.sqrt(np.diag(cov))
            checks.close(actual["standard_errors"], se, rtol=5e-4)
            checks.close(actual["confidence_intervals"], np.column_stack([np.array(ref["theta"])-z*se, np.array(ref["theta"])+z*se]), rtol=7e-4)
            checks.close(actual["genreg_intercepts"], ref["genreg_theta"][:k-1])
            gse = np.sqrt(np.diag(gcov))
            checks.close(actual["genreg_standard_errors"], gse, rtol=5e-4)
            checks.close(actual["genreg_confidence_intervals"], np.column_stack([np.array(ref["genreg_theta"])-z*gse, np.array(ref["genreg_theta"])+z*gse]), rtol=7e-4)
            query_x = np.array([[(-.7, .2)[j] for j in range(p)], [(0.2, -.3)[j] for j in range(p)], [(0.85, .6)[j] for j in range(p)]]).reshape(3, p)
            if name == "shifted_predictor":
                query_x = 10+3*query_x
            q = [.2, .43, .81]
            query = [list(row)+[prob] for row, prob in zip(query_x, q)]
            pred = call("genreg_ordinal_predict", query, actual["predict_parameters"])
            rt = np.array(ref["theta"])
            ps = probabilities(rt, query_x, k)
            cdf = np.column_stack([special.expit(rt[None, :k-1]+(query_x@rt[k-1:])[:, None]), np.ones(3)])
            checks.close(pred["probabilities"], ps)
            checks.close(pred["cumulative_probabilities"], cdf)
            checks.close(pred["linear_predictor"], query_x@rt[k-1:])
            checks.close(pred["most_likely"], labels[np.argmax(ps, axis=1)], rtol=0, atol=0)
            checks.close(pred["quantiles"], [labels[np.flatnonzero(c >= u)[0]] for c, u in zip(cdf, q)], rtol=0, atol=0)
            jac = optimize._numdiff.approx_derivative(lambda t: probabilities(t, query_x, k).flatten(), rt, method="3-point")
            wanted_se = np.sqrt(np.einsum("ij,jk,ik->i", jac, cov, jac)).reshape(3, k)
            checks.close(pred["probability_standard_errors"], wanted_se, rtol=7e-4)
            checks.close(pred["probability_confidence_intervals"], np.stack([ps-z*wanted_se, ps+z*wanted_se], axis=-1), rtol=1e-3)
            case["results"][layer] = {"fit": actual, "prediction": pred, "independent_average_nll_score": score.tolist()}
        result["cases"].append(case)
    result["edge_contracts"] = edge_contracts(calls, nullable, checks)
    result["assertions"] = checks.count
    result["assertion_count_method"] = "Numeric scalar comparisons plus explicit convergence,independent-score,NULL and inference predicates; all decoded JSON is recursively finite."
    for key, path in [("library", args.library), ("extension", args.extension)]:
        if path:
            result[key] = {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"cases": len(result["cases"]), "assertions": checks.count, "boundaries": list(calls), "output": str(args.output)}))


if __name__ == "__main__":
    main()
