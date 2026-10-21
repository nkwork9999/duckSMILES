#!/usr/bin/env python3
"""Independent exhaustive itemset/rule and NumPy scale-reliability ABI checks."""
import argparse
import ctypes
import datetime
import hashlib
import itertools
import json
from pathlib import Path
import platform

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    lib = ctypes.CDLL(str(args.library.resolve()))
    ptr = ctypes.POINTER(ctypes.c_double)
    lib.ds_analytics_run.argtypes = [ctypes.c_char_p, ctypes.c_size_t, ptr, ctypes.c_size_t,
                                   ctypes.c_size_t, ptr, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
    lib.ds_analytics_run.restype = ctypes.c_void_p
    lib.ds_analytics_free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    cases = fields = 0
    errors = []
    operations = set()

    def run(op, data, params=()):
        nonlocal cases
        x, p = np.ascontiguousarray(data, dtype=float), np.ascontiguousarray(params, dtype=float)
        size = ctypes.c_size_t()
        name = op.encode()
        result = lib.ds_analytics_run(name, len(name), x.ctypes.data_as(ptr), x.shape[0],
                                      x.shape[1], p.ctypes.data_as(ptr), len(p), ctypes.byref(size))
        assert result, (op, x, p)
        try:
            raw = ctypes.string_at(result, size.value).decode()
        finally:
            lib.ds_analytics_free(result, size.value)
        def invalid(token):
            raise AssertionError((op, token))
        cases += 1
        operations.add(op)
        return json.loads(raw, parse_constant=invalid)

    def close(a, b):
        nonlocal fields
        np.testing.assert_allclose(a, b, rtol=2e-10, atol=2e-12)
        aa, bb = np.asarray(a), np.asarray(b)
        fields += bb.size
        errors.append(float(np.max(np.abs(aa-bb)/(1+np.abs(bb)))))

    rng = np.random.default_rng(804129)
    for _ in range(70):
        n, p = int(rng.integers(5, 35)), int(rng.integers(2, 8))
        x = rng.random((n, p)) < rng.uniform(.2, .8, p)
        support, conf, lift = float(rng.choice([.1, .2, .4])), float(rng.choice([0, .4, .7])), float(rng.choice([0, 1, 1.2]))
        size = min(4, p)
        result = run("association_rules", x, [support, conf, lift, 3, size, 100000])
        # Enumerate all subsets directly, independent of production Apriori.
        counts = {}
        for width in range(1, size+1):
            for subset in itertools.combinations(range(p), width):
                counts[subset] = int(x[:, subset].all(axis=1).sum())
        expected_sets = {s: c for s, c in counts.items() if c/n >= support}
        actual_sets = {tuple(s["items"]): s["count"] for s in result["frequent_itemsets"]}
        assert actual_sets == expected_sets
        for s in result["frequent_itemsets"]:
            close(s["support"], counts[tuple(s["items"])]/n)
        expected_rules = {}
        for allocation in itertools.product([0, 1, 2], repeat=p):
            lhs = tuple(i for i in range(p) if allocation[i] == 1)
            rhs = tuple(i for i in range(p) if allocation[i] == 2)
            if not lhs or not rhs or len(lhs) > 3 or len(lhs)+len(rhs) > size:
                continue
            joint = tuple(sorted(lhs+rhs))
            c, cx, cy = counts[joint], counts[lhs], counts[rhs]
            if not c or not cx or not cy:
                continue
            confidence, ratio = c/cx, c*n/(cx*cy)
            if c/n >= support and confidence >= conf and ratio >= lift:
                expected_rules[(lhs, rhs)] = (c, cx, cy, confidence, ratio)
        actual_rules = {(tuple(r["antecedent"]), tuple(r["consequent"])): r for r in result["rules"]}
        assert set(actual_rules) == set(expected_rules)
        assert result["rule_count"] == len(expected_rules) and not result["truncated"]
        for key, (c, cx, cy, confidence, ratio) in expected_rules.items():
            r = actual_rules[key]
            close([r["support"], r["confidence"], r["lift"], r["leverage"]],
                  [c/n, confidence, ratio, c/n-cx*cy/n**2])
            if c == cx:
                assert r["conviction"] is None
                assert r["conviction_infinite"] == (cy < n)
            else:
                close(r["conviction"], (1-cy/n)/(1-confidence))

        result = run("basket_affinity", x)
        incidence = x.astype(float)
        joint = incidence.T@incidence
        totals = incidence.sum(axis=0)
        close(result["item_counts"], totals)
        assert len(result["pairs"]) == p*(p-1)//2
        for r in result["pairs"]:
            i, j = r["item0"], r["item1"]
            c, a, b = joint[i, j], totals[i], totals[j]
            close([r["joint_count"], r["support"], r["leverage"]], [c, c/n, c/n-a*b/n**2])
            for key, numerator, denominator in [("confidence0_to1", c, a), ("confidence1_to0", c, b),
                                                 ("lift", c*n, a*b), ("jaccard", c, a+b-c)]:
                if denominator == 0:
                    assert r[key] is None
                else:
                    close(r[key], numerator/denominator)

    for _ in range(100):
        n, p = int(rng.integers(5, 60)), int(rng.integers(2, 9))
        latent = rng.normal(size=(n, 1))
        x = (latent*rng.normal(size=p)+rng.normal(size=(n, p)))*np.exp(rng.uniform(-4, 4, p))
        for standardized in [0, 1]:
            z = (x-x.mean(axis=0))/x.std(axis=0, ddof=1) if standardized else x
            covariance = np.cov(z, rowvar=False)
            alpha = p/(p-1)*(1-np.trace(covariance)/covariance.sum())
            result = run("item_reliability", x, [standardized])
            close(result["alpha"], alpha)
            for j in range(p):
                rest = np.delete(z, j, axis=1)
                if p == 2:
                    assert result["alpha_if_deleted"][j] is None
                else:
                    rest_cov = np.cov(rest, rowvar=False)
                    close(result["alpha_if_deleted"][j], (p-1)/(p-2)*(1-np.trace(rest_cov)/rest_cov.sum()))
                close(result["corrected_item_total_correlations"][j], np.corrcoef(z[:, j], rest.sum(axis=1))[0, 1])

    repo = Path(__file__).resolve().parents[1]
    paths = ["crates/analytics/src/market.rs", "crates/analytics/tests/market.rs", "scripts/check_market_reference.py"]
    evidence = {"schema_version": 1, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "status": "PASS", "seed": 804129, "python": platform.python_version(), "numpy": np.__version__,
                "boundary": "Production exported ds_analytics_run/free C ABI; actual native DuckDB and installed JMP separate gates",
                "cases": cases, "numeric_fields": fields, "operations": sorted(operations), "maximum_scaled_error": max(errors),
                "source_sha256": {s: hashlib.sha256((repo/s).read_bytes()).hexdigest() for s in paths},
                "reference_methods": "Exhaustive powerset+ternary assignment rule enumeration, independent NumPy incidence matrix products, covariance/correlation matrices for raw/standardized alpha and deleted items.",
                "tolerance": {"relative": 2e-10, "absolute": 2e-12},
                "limitations": "Binary complete baskets and complete finite scale items only; no IRT/latent response parameters, frequency weighting or market causal claim.",
                "replay": "python3 scripts/check_market_reference.py --library /path/to/analytics.dylib --evidence /path/to/evidence.json",
                "harness_python": Path(__file__).read_text()}
    if args.evidence:
        args.evidence.write_text(json.dumps(evidence, indent=2)+"\n")
    print(json.dumps({k: v for k, v in evidence.items() if k not in ["source_sha256", "harness_python"]}, indent=2))


if __name__ == "__main__":
    main()
