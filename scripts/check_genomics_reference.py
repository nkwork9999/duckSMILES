#!/usr/bin/env python3
"""Independent NumPy/SciPy references against the production exported Rust ABI."""
import argparse
import ctypes
import datetime
import hashlib
import json
from pathlib import Path
import platform
import warnings

import numpy as np
import scipy
from scipy.spatial.distance import cdist
from scipy.stats import chi2, false_discovery_control


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    lib = ctypes.CDLL(str(args.library.resolve()))
    ptr = ctypes.POINTER(ctypes.c_double)
    lib.ds_analytics_run.argtypes = [ctypes.c_char_p, ctypes.c_size_t, ptr,
                                   ctypes.c_size_t, ctypes.c_size_t, ptr,
                                   ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
    lib.ds_analytics_run.restype = ctypes.c_void_p
    lib.ds_analytics_free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    fields, cases = 0, 0
    operations = set()

    def run(op, data, params=()):
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
        return json.loads(value, parse_constant=invalid)

    def close(a, b):
        nonlocal fields
        np.testing.assert_allclose(a, b, rtol=2e-8, atol=2e-10)
        fields += np.asarray(b).size

    rng = np.random.default_rng(293487)
    for trial in range(50):
        n, p = int(rng.integers(5, 25)), int(rng.integers(2, 8))
        original = rng.integers(0, 3, (n, p)).astype(float)
        original[rng.random((n, p)) < .2] = -1
        original[:2] = np.array([0., 2.])[:, None]
        observed = np.where(original == -1, np.nan, original)
        center = np.nanmean(observed, axis=0)
        counts = np.array([np.sum(original == i, axis=0) for i in range(3)]).T
        called = counts.sum(axis=1)
        allele = center / 2
        expected = called[:, None] * np.column_stack(((1-allele)**2, 2*allele*(1-allele), allele**2))
        chisq = ((counts-expected)**2/expected).sum(axis=1)
        stats = run("marker_statistics", original)
        close(stats["genotype_counts"], counts)
        close(stats["called_count"], called)
        close(stats["allele_frequency"], allele)
        close(stats["minor_allele_frequency"], np.minimum(allele, 1-allele))
        close(stats["observed_heterozygosity"], counts[:, 1]/called)
        close(stats["expected_heterozygosity"], 2*allele*(1-allele))
        close(stats["polymorphism_information_content"], 2*allele*(1-allele)-2*(allele*(1-allele))**2)
        close(stats["hwe_chi_square"], chisq)
        close(stats["hwe_p_value"], chi2.sf(chisq, 1))
        close(stats["hwe_bh_p_value"], false_discovery_control(chi2.sf(chisq, 1)))
        cases += 1

        for method in [0, 1]:
            replacement = center if method == 0 else counts.argmax(axis=1)
            imputed = run("marker_impute", original, [method])
            close(imputed["rows"], np.where(original == -1, replacement, original))
            close(imputed["predict_parameters"], replacement)
            query = np.array([[-1.]*p, [2.]*p])
            scored = run("marker_impute_predict", query, imputed["predict_parameters"])
            close(scored["rows"], np.where(query == -1, replacement, query))
            cases += 2

        completed = np.where(original == -1, center, original)
        variance = 2*allele*(1-allele)
        for method in range(3):
            if method == 0:
                Z, denominator = completed-center, variance.sum()
            elif method == 1:
                Z, denominator = (completed-center)/np.sqrt(variance), p
            else:
                Z = np.where(original == 0, -2*allele**2,
                             np.where(original == 1, variance,
                                      np.where(original == 2, -2*(1-allele)**2, 0)))
                denominator = (variance**2).sum()
            result = run("marker_relatedness", original, [method])
            close(result["matrix"], Z@Z.T/denominator)
            close(result["denominator"], denominator)
            cases += 1

        dense = rng.integers(0, 5, (n, p)).astype(float)
        dense[:, -1] = 2
        ranges = np.ptp(dense, axis=0)
        gower = (np.abs(dense[:, None, :]-dense[None, :, :])/np.where(ranges == 0, 1, ranges)).mean(axis=2)
        for method, metric in enumerate(["euclidean", "jaccard", "braycurtis", "gower", "hamming"]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                expected = gower if metric == "gower" else cdist(dense.astype(bool) if metric == "jaccard" else dense,
                                                               dense.astype(bool) if metric == "jaccard" else dense,
                                                               metric=metric)
            result = run("distance_matrix", dense, [method])
            close(result["matrix"], np.nan_to_num(expected))
            cases += 1

        # Independent dual mixed-model equations validate the production primal
        # marker-ridge solve, including an explicitly fitted GLS intercept.
        train, queries = n-2, 2
        markers = original.copy()
        y = rng.normal(size=train)
        data = np.column_stack((markers, np.r_[y, [999.]*queries]))
        means = np.nanmean(np.where(markers[:train] == -1, np.nan, markers[:train]), axis=0)
        freq = means/2
        var = 2*freq*(1-freq)
        for method in [0, 1]:
            scale = np.ones(p) if method == 0 else np.sqrt(var)
            Z = (np.where(markers == -1, means, markers)-means)/scale
            c = var.sum() if method == 0 else p
            K = Z[:train]@Z[:train].T/c
            variance_ratio = float(10**rng.uniform(-3, 3))
            A = K+variance_ratio*np.eye(train)
            inv_y, inv_one = np.linalg.solve(A, y), np.linalg.solve(A, np.ones(train))
            intercept = inv_y.sum()/inv_one.sum()
            alpha = np.linalg.solve(A, y-intercept)
            prediction = intercept+Z[train:]@Z[:train].T@alpha/c
            fitted = intercept+K@alpha
            effects = Z[:train].T@alpha/c/scale
            model = run("genomic_blup", data, [train, variance_ratio, method])
            close(model["fitted"], fitted)
            close(model["predictions"], prediction)
            close(model["marker_effects"], effects)
            scored = run("genomic_blup_predict", markers[train:], model["predict_parameters"])
            close(scored["predictions"], prediction)
            cases += 2

    # LD-kNN uses a literal perfect-linkage fixture; this rejects a mean/mode
    # substitution while avoiding a second copy of the neighbor algorithm.
    ld = run("marker_impute", [[0, 0], [1, 1], [2, 2], [-1, 2]], [2, 1, 1])
    close(ld["rows"], [[0, 0], [1, 1], [2, 2], [2, 2]])
    simulated = run("marker_simulate", [[.3]], [30000, 42, .4])
    frequency = np.bincount(np.asarray(simulated["rows"], dtype=int).ravel(), minlength=3)/30000
    assert np.max(np.abs(frequency-np.array([.574, .252, .174]))) < .02
    parental = run("marker_cross", [[0, 2, 1], [2, 0, 1]], [0, 1, 30000, 31])
    offspring = np.asarray(parental["rows"])
    assert np.all(offspring[:, :2] == 1)
    assert np.max(np.abs(np.bincount(offspring[:, 2].astype(int), minlength=3)/30000-[.25, .5, .25])) < .02
    cases += 3
    assert len(operations) == 9
    repo = Path(__file__).resolve().parents[1]
    paths = ["crates/analytics/src/genomics.rs", "crates/analytics/tests/genomics.rs", "scripts/check_genomics_reference.py"]
    evidence = {
        "schema_version": 1,
        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "boundary": "Production exported ds_analytics_run/free C ABI; native DuckDB and actual JMP are separate gates",
        "seed": 293487,
        "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "cases": cases, "numeric_fields": fields,
        "operations": sorted(operations), "status": "PASS",
        "source_sha256": {s: hashlib.sha256((repo/s).read_bytes()).hexdigest() for s in paths},
        "tolerance": {"relative": 2e-8, "absolute": 2e-10},
        "reference_methods": "NumPy genotype aggregates, SciPy chi2.sf/BH/cdist, matrix-product VanRaden/dominance kernels, independent dual GLS equations for primal GBLUP; fixed-linkage literal and deterministic-seed simulation frequency checks.",
        "replay": "python3 scripts/check_genomics_reference.py --library /path/to/analytics.dylib --evidence /path/to/evidence.json",
        "harness_python": Path(__file__).read_text(),
    }
    if args.evidence:
        args.evidence.write_text(json.dumps(evidence, indent=2)+"\n")
    print(json.dumps({k:v for k,v in evidence.items() if k not in ["source_sha256", "harness_python"]}, indent=2))


if __name__ == "__main__":
    main()
