#!/usr/bin/env python3
"""Independent NumPy checks through the production analytics C ABI.

Build the dependency-free crate as a cdylib and supply --library. This is the
direct exported-ABI layer; it does not claim native DuckDB or JMP concordance.
"""
import argparse
import ctypes
import datetime
import hashlib
import json
from pathlib import Path
import platform

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    lib = ctypes.CDLL(str(args.library.resolve()))
    pointer = ctypes.POINTER(ctypes.c_double)
    lib.ds_analytics_run.argtypes = [ctypes.c_char_p, ctypes.c_size_t, pointer,
                                   ctypes.c_size_t, ctypes.c_size_t, pointer,
                                   ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
    lib.ds_analytics_run.restype = ctypes.c_void_p
    lib.ds_analytics_free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    fields = 0
    cases = 0
    operations = set()

    def run(op, data, params=()):
        operations.add(op)
        x = np.ascontiguousarray(data, dtype=np.float64)
        p = np.ascontiguousarray(params, dtype=np.float64)
        n = ctypes.c_size_t()
        name = op.encode()
        result = lib.ds_analytics_run(name, len(name), x.ctypes.data_as(pointer),
                                      x.shape[0], x.shape[1], p.ctypes.data_as(pointer),
                                      len(p), ctypes.byref(n))
        assert result, (op, data, params)
        try:
            output = ctypes.string_at(result, n.value).decode()
        finally:
            lib.ds_analytics_free(result, n.value)
        def invalid(token):
            raise AssertionError((op, token))
        return json.loads(output, parse_constant=invalid)

    def close(actual, expected):
        nonlocal fields
        np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-12)
        fields += np.asarray(expected).size

    rng = np.random.default_rng(487112)
    for trial in range(100):
        n = int(rng.integers(4, 40))
        y = rng.integers(0, 2, n)
        y[:2] = [0, 1]
        w = rng.uniform(.05, 4, n)
        s = rng.integers(-3, 4, n) / 4
        W, P, N = w.sum(), w[y == 1].sum(), w[y == 0].sum()
        data = np.column_stack((y, s, w))
        # Pairwise weighted ranking is independent of the production sort/batch
        # trapezoid implementation and gives half credit for equal scores.
        ps, ns = s[y == 1, None], s[y == 0][None, :]
        pw = w[y == 1, None] * w[y == 0][None, :]
        expected_auc = ((ps > ns) * pw + .5 * (ps == ns) * pw).sum() / (P * N)
        thresholds = np.unique(s)[::-1]
        tp = np.array([w[(y == 1) & (s >= t)].sum() for t in thresholds])
        fp = np.array([w[(y == 0) & (s >= t)].sum() for t in thresholds])
        recall = np.r_[0., tp / P]
        fpr = np.r_[0., fp / N]
        precision = np.r_[1., tp / (tp + fp)]
        roc = run("roc_curve", data)
        close(roc["auc"], expected_auc)
        close(roc["thresholds"][1:], thresholds)
        close(roc["true_positive_rate"], recall)
        close(roc["false_positive_rate"], fpr)
        assert roc["thresholds"][0] is None
        pr = run("precision_recall_curve", data)
        close(pr["precision"], precision)
        close(pr["recall"], recall)
        close(pr["average_precision"], (np.diff(recall) * precision[1:]).sum())
        close(pr["trapezoid_area"], np.trapezoid(precision, recall))
        lift = run("lift_curve", data)
        close(lift["selection_fraction"][1:], (tp + fp) / W)
        close(lift["cumulative_gain"], recall)
        close(lift["lift"][1:], tp / (tp + fp) / (P / W))
        assert lift["lift"][0] is None
        cases += 3

        probability = rng.uniform(0, 1, n)
        binary = np.column_stack((y, probability, w))
        threshold = float(rng.uniform(.1, .9))
        decision = probability >= threshold
        b = run("binary_metrics", binary, [threshold])
        counts = [w[(y == a) & (decision == b_)].sum()
                  for a, b_ in [(0, 0), (0, 1), (1, 0), (1, 1)]]
        close(b["counts"], counts)
        close(b["accuracy"], np.average(decision == y, weights=w))
        close(b["brier_score"], np.average((probability - y)**2, weights=w))
        close(b["log_loss"], np.average(-np.log(np.where(y == 1, probability, 1-probability)), weights=w))
        actual = rng.normal(size=n)
        prediction = rng.normal(size=n)
        reg = run("regression_metrics", np.column_stack((actual, prediction, w)))
        residuals = actual - prediction
        close(reg["rmse"], np.sqrt(np.average(residuals**2, weights=w)))
        close(reg["mae"], np.average(np.abs(residuals), weights=w))
        close(reg["bias"], np.average(residuals, weights=w))
        close(reg["r_squared"], 1-np.average(residuals**2, weights=w)/np.average((actual-np.average(actual, weights=w))**2, weights=w))
        close(reg["residuals"], residuals)
        cases += 2

        true_labels = rng.choice([-3, 2, 9], n)
        predicted_labels = rng.choice([-3, 2, 9], n)
        cm = run("confusion_matrix", np.column_stack((true_labels, predicted_labels, w)))
        labels = np.union1d(true_labels, predicted_labels)
        count = np.zeros((len(labels), len(labels)))
        np.add.at(count, (np.searchsorted(labels, true_labels), np.searchsorted(labels, predicted_labels)), w)
        close(cm["labels"], labels)
        close(cm["matrix"], count)
        close(cm["accuracy"], np.trace(count)/W)
        close(cm["macro_f1"], np.mean(2*np.diag(count)/(count.sum(axis=0)+count.sum(axis=1))))
        cases += 1

        bins = int(rng.integers(2, 15))
        cal = run("calibration_curve", binary, [bins])
        group = np.minimum((probability*bins).astype(int), bins-1)
        masses = np.bincount(group, weights=w, minlength=bins)
        close(cal["bin_weight"], masses)
        expected_error = 0.
        for j in range(bins):
            chosen = group == j
            if not chosen.any():
                assert cal["mean_probability"][j] is None
                assert cal["observed_fraction"][j] is None
            else:
                pbar = np.average(probability[chosen], weights=w[chosen])
                ybar = np.average(y[chosen], weights=w[chosen])
                close(cal["mean_probability"][j], pbar)
                close(cal["observed_fraction"][j], ybar)
                expected_error += masses[j]/W*abs(pbar-ybar)
        close(cal["expected_calibration_error"], expected_error)
        profit_matrix = rng.normal(size=(2, 2))
        profit = run("profit_decision", binary, profit_matrix.ravel())
        expected_profit = np.column_stack((1-probability, probability)) @ profit_matrix
        best = expected_profit.argmax(axis=1)
        close(profit["decisions"], best)
        close(profit["actual_profit"], profit_matrix[y, best])
        close(profit["expected_profit"], expected_profit[np.arange(n), best])
        close(profit["mean_actual_profit"], np.average(profit_matrix[y, best], weights=w))
        close(profit["mean_expected_profit"], np.average(expected_profit[np.arange(n), best], weights=w))
        cases += 2

        members = rng.normal(size=(n, 4))
        member_weights = rng.uniform(.01, 10, 4)
        en = run("ensemble_predict", members, member_weights)
        close(en["predictions"], np.average(members, axis=1, weights=member_weights))
        comparison = run("model_compare", np.column_stack((actual, members)))
        scores = np.sqrt(np.mean((actual[:, None]-members)**2, axis=0))
        close(comparison["scores"], scores)
        close(comparison["ranked_model_indices"], np.argsort(scores, kind="stable"))
        compare_binary = run("model_compare", np.column_stack((y, probability, 1-probability)), [1, 1])
        close(compare_binary["scores"], [np.mean((probability-y)**2), np.mean((1-probability-y)**2)])
        cases += 3

    # Splits are checked by semantic invariants, not by replicating SplitMix.
    for seed in range(20):
        labels = np.repeat([0, 1, 2], [7, 11, 14])
        split = run("holdout_split", labels[:, None], [.7, seed, 1])
        train, test = split["train_indices"], split["test_indices"]
        assert set(train).isdisjoint(test)
        assert sorted(train + test) == list(range(len(labels)))
        close(np.bincount(labels[train], minlength=3), [4, 7, 9])
        close(np.bincount(labels[test], minlength=3), [3, 4, 5])
        rolling = run("rolling_origin", np.arange(32)[:, None], [8, 3, 2, 6])
        assert all(max(a) < min(b) for a, b in zip(rolling["train_indices"], rolling["test_indices"]))
        assert all(len(a) == 6 and len(b) == 3 for a, b in zip(rolling["train_indices"], rolling["test_indices"]))
        cases += 2
    assert len(operations) == 12
    repo = Path(__file__).resolve().parents[1]
    sources = ["crates/analytics/src/evaluation.rs", "crates/analytics/tests/evaluation.rs",
               "scripts/check_evaluation_reference.py"]
    evidence = {
        "schema_version": 1,
        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "boundary": "Actual production exported ds_analytics_run/free C ABI; no DuckDB/JMP claim",
        "seed": 487112,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "cases": cases,
        "numeric_fields": fields,
        "operations": sorted(operations),
        "status": "PASS",
        "source_sha256": {s: hashlib.sha256((repo/s).read_bytes()).hexdigest() for s in sources},
        "reference_method": "Independent NumPy weighted aggregates, pairwise weighted ROC ranking with half-ties, threshold-mask curves, matrix profit and confusion evaluation; literal/count/invariant split checks.",
        "tolerance": {"relative": 2e-10, "absolute": 2e-12},
        "replay": "python3 scripts/check_evaluation_reference.py --library /path/to/analytics.dylib --evidence /path/to/evidence.json",
        "harness_python": Path(__file__).read_text(),
    }
    if args.evidence:
        args.evidence.write_text(json.dumps(evidence, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({k:v for k,v in evidence.items() if k not in ["harness_python","source_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
