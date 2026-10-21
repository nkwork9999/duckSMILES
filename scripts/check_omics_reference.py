#!/usr/bin/env python3
"""Independent normalization and constrained-likelihood checks through exported C ABI."""
import argparse
import ctypes
import datetime
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import xlogy
from scipy.stats import rankdata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    lib = ctypes.CDLL(str(args.library.resolve()))
    ptr = ctypes.POINTER(ctypes.c_double)
    lib.ds_analytics_run.argtypes = [ctypes.c_char_p, ctypes.c_size_t, ptr,
                                   ctypes.c_size_t, ctypes.c_size_t, ptr,
                                   ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
    lib.ds_analytics_run.restype = ctypes.c_void_p
    lib.ds_analytics_free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    rng = np.random.default_rng(938512)
    cases = fields = 0
    operations = set()
    errors = []

    def run(op, data, params=()):
        nonlocal cases
        operations.add(op)
        x = np.ascontiguousarray(data, dtype=float)
        p = np.ascontiguousarray(params, dtype=float)
        size = ctypes.c_size_t()
        name = op.encode()
        result = lib.ds_analytics_run(name, len(name), x.ctypes.data_as(ptr),
                                      x.shape[0], x.shape[1], p.ctypes.data_as(ptr),
                                      len(p), ctypes.byref(size))
        assert result, (op, x, p)
        try:
            value = ctypes.string_at(result, size.value).decode()
        finally:
            lib.ds_analytics_free(result, size.value)
        def invalid(token):
            raise AssertionError((op, token))
        cases += 1
        return json.loads(value, parse_constant=invalid)

    def close(a, b, rtol=2e-8, atol=2e-10):
        nonlocal fields
        np.testing.assert_allclose(a, b, rtol=rtol, atol=atol)
        aa, bb = np.asarray(a), np.asarray(b)
        errors.append(float(np.max(np.abs(aa-bb)/(1+np.abs(bb)))))
        fields += bb.size

    for _ in range(40):
        n, p = int(rng.integers(2, 12)), int(rng.integers(3, 15))
        counts = rng.integers(0, 100, size=(n, p)).astype(float)
        counts[:, 0] += 1
        lengths = np.exp(rng.uniform(3, 12, p))
        totals = counts.sum(axis=1, keepdims=True)
        for method, params, expected in [
            (1, [], 1e6*counts/totals),
            (2, lengths, 1e9*counts/totals/lengths),
            (3, lengths, 1e6*(counts/lengths)/(counts/lengths).sum(axis=1, keepdims=True)),
            (5, [], 100*counts/totals), (6, [], counts/totals),
        ]:
            close(run("assay_normalize", counts, [method, *params])["rows"], expected)
        for ddof in [0, 1]:
            expected = (counts-counts.mean(axis=1, keepdims=True))/counts.std(axis=1, keepdims=True, ddof=ddof)
            close(run("assay_normalize", counts, [0, ddof])["rows"], expected)
        logs = np.log(counts+0.25)
        close(run("assay_normalize", counts, [4, .25])["rows"], logs-logs.mean(axis=1, keepdims=True))

        # Empirical rank lookup uses an independent SciPy average-rank grouping.
        x = rng.integers(-3, 9, size=(n, p)).astype(float)
        reference = np.sort(x, axis=1).mean(axis=0)
        def mapped(row):
            result = np.empty(p)
            ranks = rankdata(row, method="average")
            for value in np.unique(row):
                mask = row == value
                count = int(mask.sum())
                start = int(ranks[mask][0]-(count-1)/2-1)
                result[mask] = reference[start:start+count].mean()
            return result
        fitted = run("quantile_normalize", x)
        close(fitted["reference_quantiles"], reference)
        close(fitted["rows"], np.array([mapped(row) for row in x]))
        queries = rng.integers(-8, 12, size=(3, p)).astype(float)
        close(run("quantile_normalize_predict", queries, reference)["rows"], np.array([mapped(row) for row in queries]))

        genotypes = rng.integers(0, 3, size=(n, p)).astype(float)
        genotypes[rng.random((n, p)) < .2] = -1
        genotypes[0] = 0
        genotypes[:, 0] = 1
        observed = genotypes >= 0
        frequencies = np.where(observed, genotypes, 0).sum(axis=0)/(2*observed.sum(axis=0))
        frequencies = np.clip(frequencies, 1e-8, 1-1e-8)
        model = run("marker_admixture", genotypes, [1, 1000, 1e-12, 42, 1])
        close(model["allele_frequencies"], frequencies[None, :])
        probs = np.broadcast_to(frequencies, genotypes.shape)
        ll = np.where(observed, xlogy(np.maximum(genotypes, 0), probs)+
                      xlogy(2-np.maximum(genotypes, 0), 1-probs)+
                      np.where(genotypes == 1, np.log(2), 0), 0).sum()
        close(model["log_likelihood"], ll)
        scored = run("admixture_predict", genotypes[:2], model["predict_parameters"])
        close(scored["ancestry"], [[1], [1]])

    # Fixed frequencies yield a concave simplex likelihood. SLSQP with analytic
    # derivatives is independent of production EM; compare attained likelihood,
    # not nonidentifiable population labels or stopping iterations.
    for k in [2, 3, 4]:
        for _ in range(25):
            p = 12
            frequency = rng.uniform(.03, .97, (k, p))
            x = rng.integers(0, 3, size=(3, p)).astype(float)
            x[rng.random(x.shape) < .1] = -1
            x[:, 0] = 1
            params = [k, p, 10000, 1e-12, 1e-8, *([1/k]*k), *frequency.ravel()]
            scored = run("admixture_predict", x, params)
            for row, q, actual_ll in zip(x, scored["ancestry"], scored["sample_log_likelihood"]):
                mask = row >= 0
                f, g = frequency[:, mask], row[mask]
                def objective(weights):
                    pr = weights@f
                    return -(xlogy(g, pr)+xlogy(2-g, 1-pr)+np.where(g == 1, np.log(2), 0)).sum()
                def gradient(weights):
                    pr = weights@f
                    return -f@(g/pr-(2-g)/(1-pr))
                oracle = minimize(objective, np.full(k, 1/k), jac=gradient,
                                  method="SLSQP", bounds=[(0, 1)]*k,
                                  constraints={"type": "eq", "fun": lambda w: w.sum()-1,
                                               "jac": lambda w: np.ones(k)},
                                  options={"ftol": 1e-12, "maxiter": 2000})
                assert oracle.success, oracle.message
                close(actual_ll, -objective(np.asarray(q)), atol=2e-10)
                close(actual_ll, -oracle.fun, rtol=1e-8, atol=2e-7)
                close(np.sum(q), 1.)
                assert np.min(q) >= 0 and np.max(q) <= 1

    # A known separated two-population global supremum, subject to explicit floor.
    x = [[0, 2], [0, 2], [1, 1], [2, 0], [2, 0]]
    model = run("marker_admixture", x, [2, 10000, 1e-12, 42, 3])
    close(model["log_likelihood"], -2*np.log(2)+16*np.log1p(-1e-8), atol=3e-9)
    history = np.asarray(model["log_likelihood_history"])
    assert np.min(np.diff(history)) >= -1e-10
    assert model["converged"]
    close(np.asarray(model["ancestry"]).sum(axis=1), np.ones(5))
    assert len(operations) == 5
    repo = Path(__file__).resolve().parents[1]
    paths = ["crates/analytics/src/omics.rs", "crates/analytics/tests/omics.rs", "scripts/check_omics_reference.py"]
    evidence = {
        "schema_version": 1, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "status": "PASS", "boundary": "Production exported ds_analytics_run/free C ABI; native DuckDB and installed JMP are separate gates",
        "seed": 938512, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "cases": cases, "numeric_fields": fields, "operations": sorted(operations),
        "maximum_scaled_error": max(errors),
        "source_sha256": {s: hashlib.sha256((repo/s).read_bytes()).hexdigest() for s in paths},
        "reference_methods": "NumPy counts/std/log centering/sorted quantiles, SciPy rankdata for tie groups, closed-form one-population binomial MLE, SLSQP convex fixed-frequency ancestry likelihood with analytic gradients for K2/3/4, analytic separated-population likelihood.",
        "limitations": "EM fit is nonconvex in joint Q/F; random general-case global optima, parameter SEs, and installed JMP matching are not established. Fixed-frequency query fit likelihood tolerance2e-7+1e-8relative, standard numeric tolerance2e-10+2e-8relative.",
        "replay": "python3 scripts/check_omics_reference.py --library /path/to/analytics.dylib --evidence /path/to/evidence.json",
        "harness_python": Path(__file__).read_text(),
    }
    if args.evidence:
        args.evidence.write_text(json.dumps(evidence, indent=2)+"\n")
    print(json.dumps({k: v for k, v in evidence.items() if k not in ["source_sha256", "harness_python"]}, indent=2))


if __name__ == "__main__":
    main()
