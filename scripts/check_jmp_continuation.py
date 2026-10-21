#!/usr/bin/env python3
"""Check the 62 continuation methods through native DuckDB SQL -> C++ -> Rust.

Standard library only. Requires matching native DuckDB CLI and extension builds.
Numerical fixtures are hand-derived or independently checked reference literals;
catalog smoke, model reuse, malformed inputs and vectorized execution are separate
checks. Does not establish numerical equivalence to an installed JMP executable.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
API_FILES = {
    "latent": ("jmp-latent-api.json", 8),
    "evaluation": ("jmp-evaluation-api.json", 12),
    "nonlinear": ("jmp-nonlinear-api.json", 7),
    "spectral": ("jmp-spectral-api.json", 10),
    "latent_inference": ("jmp-latent-inference-api.json", 4),
    "robust_learning": ("jmp-robust-learning-api.json", 5),
    "genomics": ("jmp-genomics-api.json", 9),
    "reliability": ("jmp-reliability-api.json", 7),
}

def sql_list(x):
    if isinstance(x, (list, tuple)):
        return "[" + ",".join(map(sql_list, x)) + "]"
    return repr(float(x))

def expression(op, data, args=()):
    return f"chem_{op}({sql_list(data)}, {sql_list(args)})"

def strict_json(value):
    def reject(token):
        raise AssertionError(f"non-JSON constant {token}")
    result = json.loads(value, parse_constant=reject)
    def finite(x):
        if isinstance(x, float):
            assert math.isfinite(x), x
        elif isinstance(x, (tuple, list)):
            for child in x:
                finite(child)
        elif isinstance(x, dict):
            for child in x.values():
                finite(child)
    finite(result)
    assert isinstance(result, dict) and result.get("version") == 1, result
    return result

def equal(actual, wanted, label, relative=1e-7, absolute=1e-9):
    if isinstance(wanted, (tuple, list)):
        assert len(actual) == len(wanted), (label, actual, wanted)
        for i, (a, b) in enumerate(zip(actual, wanted)):
            equal(a, b, f"{label}[{i}]", relative, absolute)
    elif wanted is None or isinstance(wanted, (bool, str)):
        assert actual == wanted, (label, actual, wanted)
    else:
        assert actual is not None and math.isclose(actual, wanted, rel_tol=relative, abs_tol=absolute), (label, actual, wanted)

def field(result, key):
    for part in key.split("."):
        result = result[int(part)] if part.isdigit() else result[part]
    return result

class Runner:
    def __init__(self, cli, extension):
        self.cli = cli
        self.extension = str(extension.resolve()).replace("'", "''")
        self.queries = 0
        self.decoded = 0
    def query(self, sql):
        self.queries += 1
        process = subprocess.run([self.cli, "-unsigned", "-batch", "-bail", "-json", ":memory:"], input=f"LOAD '{self.extension}';\n{sql}\n", text=True, capture_output=True, timeout=180)
        assert process.returncode == 0, f"DuckDB failed ({process.returncode}):\n{process.stderr}"
        return json.loads(process.stdout)
    def many(self, expressions):
        sql = " UNION ALL ".join(f"SELECT {i} AS id, {e} AS result" for i, e in enumerate(expressions))
        rows = self.query(sql + " ORDER BY id;")
        assert len(rows) == len(expressions)
        out = []
        for i, row in enumerate(rows):
            assert row["id"] == i and row["result"] is not None, (i, expressions[i], row)
            out.append(strict_json(row["result"]))
            self.decoded += 1
        return out
    def one(self, expr):
        return self.many([expr])[0]

def inventory():
    operations = []
    for group, (filename, count) in API_FILES.items():
        api = json.loads((ROOT / "tasks" / filename).read_text())
        assert len(api["operations"]) == count, (filename, len(api["operations"]), count)
        for item in api["operations"]:
            assert item["sql_name"] == "chem_" + item["operation"], item
            assert item["sql_example"].startswith(item["sql_name"] + "("), item
            operations.append((group, item))
    assert len(operations) == len({o[1]["operation"] for o in operations}) == 62
    return operations

def numerical_cases():
    cases = []
    def add(op, data, args, expected, tolerance=1e-7):
        cases.append((op, expression(op, data, args), expected, tolerance))
    add("mds", [[0], [1], [2]], [1], {"eigenvalues.0": 2, "normalized_raw_stress": 0})
    add("correspondence_analysis", [[20, 10], [10, 20]], [1], {"inertias.0": 1 / 9})
    add("factor_analysis", [[1, .64, .64], [.64, 1, .64], [.64, .64, 1]], [1, 0, 500, 1e-10, 1, 1], {"communalities": [.64] * 3, "uniquenesses": [.36] * 3})
    add("canonical_correlation", [[-1, -2], [0, 0], [1, 2], [2, 4]], [1], {"canonical_correlations": [1]})
    add("gaussian_mixture", [[-1], [0], [1]], [1], {"component_means": [[0]], "component_covariances_flat": [[2 / 3]], "probabilities": [[1]] * 3})
    add("self_organizing_map", [[-1], [0], [1]], [1, 1], {"codebook": [[0]], "labels": [0] * 3, "quantization_error": math.sqrt(2 / 3)})
    scored = [[0, .1], [1, .4], [0, .4], [1, .8]]
    add("regression_metrics", [[1, 2, 1], [2, 2, 2], [3, 1, 1]], [], {"rmse": math.sqrt(1.25), "mae": .75, "r_squared": -1.5, "bias": .25})
    add("confusion_matrix", [[10, 10], [10, 20], [20, 20], [30, 20]], [], {"accuracy": .5, "balanced_accuracy": .5, "macro_f1": 7 / 18})
    add("binary_metrics", scored, [.4], {"counts": [1, 1, 0, 2], "brier_score": .1425, "log_loss": .4389051056530454})
    add("roc_curve", scored, [], {"auc": .875, "thresholds.0": None, "true_positive_rate": [0, .5, 1, 1]})
    add("precision_recall_curve", scored, [], {"average_precision": 5 / 6, "trapezoid_area": 11 / 12})
    add("lift_curve", scored, [], {"selection_fraction": [0, .25, .75, 1], "lift": [None, 2, 4 / 3, 1]})
    add("calibration_curve", [[0, 0], [1, .5], [1, 1]], [4], {"bin_weight": [1, 0, 1, 1], "expected_calibration_error": 1 / 6})
    add("profit_decision", [[0, .1], [1, .6], [1, .2]], [2, -1, -5, 10], {"decisions": [0, 1, 1], "total_actual_profit": 22, "mean_expected_profit": 2.7})
    add("ensemble_predict", [[1, 3], [3, 7]], [1, 3], {"predictions": [2.5, 6]})
    add("curve_fit", [[0, 2.03], [.5, 2.64], [1, 3.7], [1.5, 4.88], [2, 6.68], [2.5, 8.94]], [0, 500, 1e-11, 1.8, .4], {"coefficients": [2.0045698299903227, .5987386163457715], "sse": .011118711912985363, "converged": True})
    add("bspline_basis", [[0], [.5], [1], [2]], [3, 0, 0, 0, 0, 1, 2, 2, 2, 2], {"basis.0": [1, 0, 0, 0, 0], "basis.1": [.125, .59375, .25, .03125, 0], "basis.3": [0, 0, 0, 0, 1]})
    add("spectral_snv", [[1, 2, 3], [10, 12, 14]], [], {"rows": [[-1, 0, 1], [-1, 0, 1]]})
    add("spectral_msc", [[5, 7, 11, 19]], [1, 2, 4, 8], {"rows": [[1, 2, 4, 8]], "slopes": [2], "intercepts": [3]})
    add("savitzky_golay", [[0, 1, 4, 9, 16, 25, 36]], [5, 2, 1, .5], {"rows": [[0, 4, 8, 12, 16, 20, 24]]})
    add("curve_derivative", [[0, 1, 4, 9, 16]], [2, .5], {"rows": [[8] * 5]})
    add("spectral_detrend", [[2, 3, 6, 11, 18]], [2], {"rows": [[0] * 5]})
    add("haar_transform", [[1, 2, 3, 4]], [], {"rows": [[5, -2, -math.sqrt(.5), -math.sqrt(.5)]]})
    add("fourier_transform", [[1, 2, 3, 4]], [], {"rows": [[10, -2, -2, -2, 0, 2, 0, -2]]})
    add("dynamic_time_warping", [[1, 0], [2, 0], [3, 0], [1, 1], [1, 1], [2, 1], [3, 1]], [], {"distance": 0})
    add("functional_pca", [[1, 2, 3], [2, 4, 6], [3, 6, 9]], [1], {"eigenvalues": [9], "quadrature_weights": [.5, 1, .5], "reconstructed_rows": [[1, 2, 3], [2, 4, 6], [3, 6, 9]]})
    add("peak_find", [[0, 3, 1, 4, 0, 2, 0]], [], {"peak_indices": [[1, 3, 5]], "prominences": [[2, 4, 2]]})
    manova = [[g, a + 2 * g, b + g] for g in [0, 1] for a in [-1, 1] for b in [-1, 1]]
    add("manova", manova, [1], {"hypothesis_sscp": [[8, 4], [4, 2]], "error_sscp": [[8, 0], [0, 8]], "wilks_lambda": 4 / 9, "p_values": [.1316872427983539] * 4})
    add("canonical_correlation_inference", [[.8], [.6]], [100, 2, 2], {"wilks_lambdas": [.2304, .64], "degrees_of_freedom": [4, 1], "bartlett_chi_square.0": 141.6560507904787})
    add("multiple_correspondence_analysis", [[0, 0], [0, 0], [1, 1], [1, 1]], [2], {"inertias": [1, 0], "total_inertia": 1, "benzecri_inertias": [1, 0]})
    add("factor_analysis_ml", [[1, .64, .64], [.64, 1, .64], [.64, .64, 1]], [1, 0, 500, 1e-9, 1e-5, 1, 0], {"communalities": [.64] * 3, "uniquenesses": [.36] * 3}, 1e-5)
    add("svm_fit", [[-2, 0], [-1, 0], [1, 1], [2, 1]], [4, 0, 10, 1, 1e-9, 10000], {"dual_objective": 1.25, "primal_objective": 1.25, "equality_residual": 0, "converged": True})
    add("quantile_regression", [[0], [1], [2], [3], [100]], [.75, 1, 50000, 1e-9], {"intercept": 3, "objective": 74.25, "converged": True}, 2e-6)
    add("huber_regression", [[0, 1], [1, 3], [2, 5], [3, 27], [4, 9], [5, 11], [6, 13]], [1.345, 1, 500, 1e-10], {"intercept": 1 + 1.345 / 6, "coefficients": [2], "robust_scale": 1, "converged": True})
    add("marker_statistics", [[0, 0], [1, 1], [1, 2], [2, -1]], [], {"called_count": [4, 3], "allele_frequency": [.5, .5], "polymorphism_information_content": [.375, .375], "hwe_chi_square": [0, 1 / 3]})
    add("marker_impute", [[0, 0], [1, 1], [2, 2], [-1, 2]], [2, 1, 1], {"rows": [[0, 0], [1, 1], [2, 2], [2, 2]]})
    add("marker_simulate", [[0], [1]], [2,608688947055533,0], {"rows": [[0,2],[0,2]]})
    add("marker_relatedness", [[0, 0, 2], [1, 1, 2], [2, 2, 2]], [0], {"matrix": [[2, 0, -2], [0, 0, 0], [-2, 0, 2]]})
    add("distance_matrix", [[0, 1], [2, 1]], [0], {"matrix": [[0, 2], [2, 0]]})
    add("genomic_blup", [[0, 1], [1, 2], [2, 5], [1, 999]], [3, 1], {"marker_effects": [1.6], "predictions": [8 / 3], "fitted": [8 / 3 - 1.6, 8 / 3, 8 / 3 + 1.6]})
    lifetime = [[x,t,e] for x,t,e in zip([-1,-1,-1,-.5,-.5,0,0,.5,.5,1,1,1],[1.2,2.1,1.6,3.5,2.7,3.2,5.5,4.8,7.1,6.2,8.4,9],[1,1,0,1,0,1,1,0,1,1,0,1])]
    add("aft_exponential", [[1,1],[2,1],[3,0],[4,0]], [], {"coefficients": [math.log(5)], "log_likelihood": -2*math.log(5)-2, "standard_errors": [math.sqrt(.5)]}, 2e-6)
    add("aft_weibull", lifetime, [], {"coefficients": [1.5422076020448903,.7390674142280002], "log_sigma": -1.7040161775568994, "log_likelihood": -12.599177632359707}, 2e-6)
    add("aft_lognormal", lifetime, [], {"coefficients": [1.440541204073499,.7697624298743153], "log_sigma": -1.3532392372547124, "log_likelihood": -13.428415667651329}, 2e-6)
    add("aft_predict", [[0],[5],[10]], [0,0,0,0,0,math.log(5)], {"survival_probabilities": [1,math.exp(-1),math.exp(-2)]})
    add("nhpp_powerlaw", [[1],[2],[4],[7]], [10], {"shape": .771457807511085, "log_eta": .5056049501753199, "expected_events_at_end": 4, "log_likelihood": -7.518067436566543})
    add("nhpp_predict", [[0],[5],[10]], [2,math.log(5)], {"expected_cumulative_events": [0,1,4], "intensities": [0,.4,.8]})
    add("reliability_blocks", [[.8,.7,.9]], [2,2], {"reliabilities": [.902]})
    return cases

def saved_models():
    tests = []
    def add(name, fit, predictor, queries, expected):
        sql = f"(WITH fit AS (SELECT {fit} AS model) SELECT chem_{predictor}({sql_list(queries)},json_extract(model,'$.predict_parameters')::DOUBLE[]) FROM fit)"
        tests.append((name, sql, expected))
    add("gaussian mixture", expression("gaussian_mixture", [[-1], [0], [1]], [1]), "gaussian_mixture_predict", [[-2], [0], [2]], {"probabilities": [[1]] * 3})
    add("SOM", expression("self_organizing_map", [[-1], [0], [1]], [1, 1]), "som_predict", [[-2], [0], [2]], {"distances": [math.sqrt(6), 0, math.sqrt(6)]})
    add("curve", expression("curve_fit", [[0, 2.03], [.5, 2.64], [1, 3.7], [1.5, 4.88], [2, 6.68], [2.5, 8.94]], [0, 500, 1e-11, 1.8, .4]), "curve_predict", [[0], [1]], {"predictions": [2.0045698299903227, 3.647959993479304]})
    add("spline", expression("spline_fit", [[0, 1], [.25, 1.5], [.5, 2], [.75, 2.5], [1, 3]], [1, 0, 0, 0, .5, 1, 1]), "spline_predict", [[.25], [.75]], {"predictions": [1.5, 2.5]})
    add("LOESS", expression("loess", [[0, 1], [1, 3], [2, 5], [3, 7], [4, 9]], [1, 1, 0]), "loess_predict", [[.5], [2.5]], {"predictions": [2, 6]})
    add("SVM", expression("svm_fit", [[-2, 0], [-1, 0], [1, 1], [2, 1]], [4, 0, 10, 1, 1e-9]), "svm_predict", [[-3], [0], [3]], {"decision_values": [-3, 0, 3], "predictions": [0, 1, 1]})
    add("quantile", expression("quantile_regression", [[0, 0], [1, 1], [2, 4], [3, 3], [4, 4]], [.5, 1, 50000, 1e-9]), "robust_predict", [[2.5]], {"predictions": [2.5]})
    add("Huber", expression("huber_regression", [[0, 1], [1, 3], [2, 5], [3, 27], [4, 9], [5, 11], [6, 13]], [1.345, 1, 500, 1e-10]), "robust_predict", [[10]], {"predictions": [21 + 1.345 / 6]})
    add("marker imputation", expression("marker_impute", [[0, -1], [1, 1], [-1, 2]], [0]), "marker_impute_predict", [[-1, -1], [2, 0]], {"rows": [[.5, 1.5], [2, 0]]})
    add("GBLUP", expression("genomic_blup", [[0, 1], [1, 2], [2, 5], [1, 999]], [3, 1]), "genomic_blup_predict", [[0], [-1], [2]], {"predictions": [8 / 3 - 1.6, 8 / 3, 8 / 3 + 1.6]})
    aft = expression("aft_exponential", [[1,1],[2,1],[3,0],[4,0]], [])
    tests.append(("AFT survival", f"(WITH fit AS (SELECT {aft} AS model) SELECT chem_aft_predict([[0.0],[5.0],[10.0]],list_concat([0.0],json_extract(model,'$.predict_parameters')::DOUBLE[])) FROM fit)", {"survival_probabilities": [1,math.exp(-1),math.exp(-2)]}))
    return tests

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", required=True, type=Path)
    args = parser.parse_args()
    runner = Runner(args.duckdb, args.extension)
    operations = inventory()
    names = [item["sql_name"] for _, item in operations]
    literals = ",".join("'" + name + "'" for name in names)
    catalog = runner.query(f"SELECT function_name,description,examples,parameter_types FROM duckdb_functions() WHERE function_name IN ({literals}) ORDER BY function_name;")
    assert len(catalog) == len(names) == 62, ("continuation catalog mismatch", len(catalog), set(names) - {r["function_name"] for r in catalog})
    assert all(r["description"] and r["examples"] and r["parameter_types"] == ["DOUBLE[][]", "DOUBLE[]"] for r in catalog), catalog
    examples = [item["sql_example"] for _, item in operations]
    runner.many(examples)
    cases = numerical_cases()
    numeric_fields = 0
    for (op, _, expected, tol), result in zip(cases, runner.many([c[1] for c in cases])):
        for key, wanted in expected.items():
            equal(field(result, key), wanted, f"{op}.{key}", tol)
            numeric_fields += 1
    reuse = saved_models()
    for (name, _, expected), result in zip(reuse, runner.many([t[1] for t in reuse])):
        for key, wanted in expected.items():
            equal(field(result, key), wanted, f"saved {name}.{key}", 3e-6)
    malformed = [
        "chem_mds([[0.0],[1.0]], [5.0])", "chem_spectral_snv([[1.0,1.0]], [])",
        "chem_bspline_basis([[2.0]], [1.0,0.0,0.0,1.0,1.0])",
        "chem_svm_fit([[0.0,0.0],[1.0,0.0]], [2.0])",
        "chem_quantile_regression([[1.0],[2.0]], [1.0])",
        "chem_marker_statistics([[3.0]], [])", "chem_marker_impute([[-1.0],[-1.0]], [])",
        "chem_analyze('unknown_continuation_operation', [[1.0]], [])",
        "chem_loess([[0.0,1.0],[1.0,2.0]], [0.5,2.0])",
        "chem_spectral_snv(NULL::DOUBLE[][], [])", "chem_spectral_snv([[1.0,NULL]], [])",
        "chem_spectral_snv([[1.0,2.0],[3.0]], [])", "chem_spectral_snv([[1.0,2.0]], [NULL])",
        "chem_aft_exponential([[0.0,1.0],[1.0,1.0]], [])", "chem_reliability_blocks([[1.1]], [])", "chem_nhpp_powerlaw([[1.0],[1.0]], [1.0])",
    ]
    bad = runner.query(" UNION ALL ".join(f"SELECT {i} AS id, ({expr}) IS NULL AS valid" for i, expr in enumerate(malformed)) + " ORDER BY id;")
    assert len(bad) == len(malformed) and all(row["valid"] for row in bad), bad
    selected = runner.query("SELECT i, chem_analyze(CASE WHEN i%2=0 THEN 'spectral_snv' ELSE 'curve_derivative' END, CASE WHEN i%2=0 THEN [[i::DOUBLE,i+1.0,i+2.0]] ELSE [[i::DOUBLE,i+1.0,i+4.0]] END, []::DOUBLE[]) AS result FROM range(5000) t(i) WHERE i%997=0 ORDER BY i;")
    for row in selected:
        result = strict_json(row["result"])
        equal(result["rows"], [[-1, 0, 1]] if row["i"] % 2 == 0 else [[0, 2, 4]], "selected dynamic operation")
        runner.decoded += 1
    # An actual NULL input interleaved with nonconstant vector data must preserve NULL rows.
    null_rows = runner.query("SELECT i,chem_spectral_snv(CASE WHEN i%3=0 THEN NULL ELSE [[i::DOUBLE,i+1.0,i+2.0]] END,[]) AS result FROM range(8) t(i) WHERE i<>4 ORDER BY i;")
    for row in null_rows:
        if row["i"] % 3 == 0:
            assert row["result"] is None, row
        else:
            equal(strict_json(row["result"])["rows"], [[-1, 0, 1]], "interleaved NULL vector")
            runner.decoded += 1
    # Prefix-only supervision must not learn suffix targets/labels in SVM or GBLUP.
    for op, prefix, parameters, query in [
        ("svm_fit", [[-2, 0], [-1, 0], [1, 1], [2, 1]], [4, 0, 10], [3]),
        ("genomic_blup", [[0, 1], [1, 2], [2, 5]], [3, 1], [1]),
    ]:
        results = runner.many([expression(op, prefix + [query + [label]], parameters) for label in [999, -1e100]])
        assert results[0] == results[1], f"{op} used suffix labels"
    print(json.dumps({"status": "PASS", "layer": "native DuckDB SQL -> C++ -> Rust", "continuation_catalog_examples": len(examples), "operations_with_numerical_fixtures": len(cases), "numeric_fields_checked": numeric_fields, "sql_saved_model_roundtrips": len(reuse), "invalid_inputs_returning_sql_null": len(malformed), "selected_dynamic_rows": len(selected), "interleaved_null_rows": len(null_rows), "prefix_models_ignoring_suffix_labels": 2, "native_sql_queries": runner.queries, "strict_json_results": runner.decoded}, indent=2))

if __name__ == "__main__":
    main()
