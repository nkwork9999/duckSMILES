#!/usr/bin/env python3
"""Check JMP-subset analytics through native DuckDB SQL -> C++ -> Rust.

Only the Python standard library is required. Numerical literals below are
hand-derived or independent NumPy 2.4.2/SciPy 1.17.1 fixtures, never results
obtained from this extension. Catalog examples exercise all 95 new functions.
"""

import argparse
import json
import math
from pathlib import Path
import subprocess


GROUPS = {
    "inference": "one_sample_t paired_t pooled_t variance_test variance_ratio_test proportion_test contingency_test fisher_exact mann_whitney wilcoxon_signed_rank kruskal_wallis levene_test equivalence_t",
    "glm": "linear_model logistic poisson glm_predict elastic_net polynomial_features",
    "quality": "ir_chart xbar_r_chart xbar_s_chart p_chart np_chart c_chart u_chart ewma_chart cusum_chart moving_average_chart western_electric gage_rr type1_gage agreement_kappa acceptance_oc levey_jennings laney_p_chart laney_u_chart",
    "resampling": "bootstrap_mean permutation_test sample_rows simulate_normal power_z sample_size_mean sample_size_proportion reliability_demonstration odds_ratio relative_risk cohen_d",
    "distributions": "fit_normal fit_lognormal fit_exponential fit_weibull fit_gamma fit_poisson normality_test kde_gaussian box_cox survival_exponential survival_weibull logrank_test",
    "timeseries": "acf pacf difference moving_average exponential_smoothing holt ar_model periodogram seasonal_decompose ljung_box",
    "multivariate": "covariance_matrix correlation_matrix hierarchical_cluster mahalanobis discriminant_analysis gaussian_naive_bayes knn_predict dbscan",
    "mixed": "mixed_model mixed_predict",
    "survival": "cox_ph cox_predict relative_risk_predict",
    "arima": "arima",
    "learning": "decision_tree random_forest gradient_boosting mlp",
    "advanced_design": "box_behnken latin_hypercube mixture_lattice fractional_factorial orthogonal_screening randomize_design design_evaluate",
}
OPERATIONS = tuple(op for group in GROUPS.values() for op in group.split())
assert len(OPERATIONS) == len(set(OPERATIONS)) == 95


def sql_list(value):
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(sql_list(v) for v in value) + "]"
    return repr(float(value))


def expression(operation, data, args=()):
    return f"chem_{operation}({sql_list(data)}, {sql_list(args)})"


def strict_json(value):
    def invalid(token):
        raise AssertionError(f"invalid JSON numeric token: {token}")

    result = json.loads(value, parse_constant=invalid)

    def finite_tree(item):
        if isinstance(item, float):
            assert math.isfinite(item), item
        elif isinstance(item, list):
            for child in item:
                finite_tree(child)
        elif isinstance(item, dict):
            for child in item.values():
                finite_tree(child)

    finite_tree(result)
    assert result["version"] == 1, result
    return result


def equal(actual, expected, label, tolerance=1e-7, absolute=1e-9):
    if isinstance(expected, (list, tuple)):
        assert len(actual) == len(expected), (label, actual, expected)
        for index, (left, right) in enumerate(zip(actual, expected)):
            equal(left, right, f"{label}[{index}]", tolerance, absolute)
    elif expected is None or isinstance(expected, (str, bool)):
        assert actual == expected, (label, actual, expected)
    else:
        assert actual is not None and math.isclose(
            actual, expected, rel_tol=tolerance, abs_tol=absolute
        ), (label, actual, expected)


def field(result, path):
    for name in path.split("."):
        result = result[int(name)] if name.isdigit() else result[name]
    return result


class Runner:
    def __init__(self, duckdb, extension):
        self.duckdb = str(duckdb)
        self.extension = str(extension.resolve()).replace("'", "''")
        self.sql_queries = 0
        self.json_results = 0

    def query(self, sql):
        self.sql_queries += 1
        process = subprocess.run(
            [self.duckdb, "-unsigned", "-batch", "-bail", "-json", ":memory:"],
            input=f"LOAD '{self.extension}';\n{sql}\n",
            text=True,
            capture_output=True,
            timeout=120,
        )
        if process.returncode:
            raise AssertionError(
                f"native DuckDB failed ({process.returncode}):\n{process.stderr}"
            )
        return json.loads(process.stdout)

    def many(self, expressions):
        sql = " UNION ALL ".join(
            f"SELECT {i} AS id, {expr} AS result"
            for i, expr in enumerate(expressions)
        )
        rows = self.query(sql + " ORDER BY id;")
        assert len(rows) == len(expressions)
        results = []
        for i, row in enumerate(rows):
            assert row["id"] == i
            assert row["result"] is not None, expressions[i]
            results.append(strict_json(row["result"]))
            self.json_results += 1
        return results

    def one(self, sql_expression):
        return self.many([sql_expression])[0]


def numerical_cases():
    cases = []

    def add(op, data, args, expected, tolerance=1e-7):
        cases.append((op, expression(op, data, args), expected, tolerance))

    a = [[1], [2], [4], [5], [8]]
    paired = [[1, 2], [2, 3], [4, 3], [5, 6], [8, 9]]
    grouped = [[v[0], 0] for v in a] + [[v, 1] for v in [2, 3, 3, 6, 9]]
    three = grouped + [[v, 2] for v in [1, 1, 2, 2, 3]]
    add("one_sample_t", a, [3], {"t_statistic": 0.8164965809277261, "p_value": 0.460050752843961, "confidence_interval": [0.5995630967087218, 7.400436903291278]})
    add("paired_t", paired, [], {"t_statistic": -1.5, "p_value": 0.208, "mean_difference": -0.6})
    add("pooled_t", grouped, [], {"t_statistic": -0.33752637027780696, "p_value": 0.7444056354601131, "degrees_of_freedom": 8})
    add("variance_test", a, [2], {"chi_square_statistic": 7.5, "p_value": 0.22341858563208655})
    add("variance_ratio_test", grouped, [], {"f_statistic": 0.9036144578313252, "p_value": 0.9241155364991409})
    add("proportion_test", [[v] for v in [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]], [.2], {"p_value": 0.42957465600000017, "confidence_interval": [0.10779126740630102, 0.6032218525388546]})
    add("contingency_test", [[12, 5, 3], [4, 9, 7]], [], {"chi_square_statistic": 6.742857142857142, "p_value": 0.034340544364887336, "expected_counts": [[8, 7, 5], [8, 7, 5]]})
    add("fisher_exact", [[1, 9], [11, 3]], [], {"p_value": 0.0027594561852200836, "odds_ratio": 1 / 33})
    add("mann_whitney", grouped, [], {"u_statistic": 10.5, "p_value": 0.7525671000692723})
    add("wilcoxon_signed_rank", paired, [], {"rank_statistic": 3, "p_value": 0.23303798227390538})
    add("kruskal_wallis", three, [], {"h_statistic": 4.220664206642074, "p_value": 0.12119770958717949})
    add("levene_test", three, [], {"f_statistic": 1.053763440860215, "p_value": 0.37877574177837225})
    add("equivalence_t", a, [2, 6], {"p_value": 0.08890390417811062, "equivalent": False, "confidence_interval": [1.3890315818537937, 6.6109684181462045]})

    linear = [[0, 1], [1, 2], [2, 1], [3, 4], [4, 4]]
    add("linear_model", linear, [], {"intercept": .8, "coefficients": [.8], "sse": 2.8, "residual_df": 3, "coefficient_standard_errors": [.30550504633038933], "overall_p_value": .07909569938205474})
    logistic = [[0, 0], [0, 0], [0, 0], [0, 1], [1, 0], [1, 1], [1, 1], [1, 1]]
    add("logistic", logistic, [], {"intercept": -math.log(3), "coefficients": [2 * math.log(3)], "fitted_values": [.25] * 4 + [.75] * 4, "deviance": 8.997362313900933, "converged": True})
    add("poisson", [[0, 0], [0, 1], [0, 2], [0, 1], [1, 2], [1, 4], [1, 6], [1, 4]], [], {"intercept": 0, "coefficients": [math.log(4)], "deviance": 4.865581297297973, "fitted_values": [1] * 4 + [4] * 4, "converged": True})
    add("glm_predict", [[-1000], [1000]], [1, 0, 1], {"predictions": [0, 1]})
    add("elastic_net", [[-1, -1, .8], [-1, 1, 1.2], [1, -1, 4.8], [1, 1, 5.2]], [.5, 1], {"intercept": 3, "coefficients": [1.5, 0], "objective": .895, "converged": True})
    add("polynomial_features", [[2, 3]], [3, 1], {"n_output_features": 10, "transformed_rows": [[1, 2, 3, 4, 6, 9, 8, 12, 18, 27]]})

    add("ir_chart", [[1], [2], [3], [4]], [], {"center": 2.5, "sigma": math.sqrt(math.pi) / 2})
    subgroups = [[1, 3], [2, 4], [3, 5]]
    add("xbar_r_chart", subgroups, [], {"center": 3, "sigma": 2 / 1.128})
    add("xbar_s_chart", subgroups, [], {"center": 3, "sigma": math.sqrt(math.pi)})
    add("p_chart", [[2, 100], [4, 100]], [], {"center": .03, "upper.0": .03 + 3 * math.sqrt(.03 * .97 / 100)})
    add("np_chart", [[2, 100], [4, 100]], [], {"center": 3})
    add("c_chart", [[4], [6]], [], {"center": 5})
    add("u_chart", [[2, 2], [6, 2]], [], {"center": 2, "upper": [5, 5]})
    add("ewma_chart", [[1], [2], [3]], [.5, 3, 0, 1], {"statistics": [.5, 1.25, 2.125], "upper.0": 1.5})
    add("cusum_chart", [[1], [2], [3]], [.5, 2, 0, 1], {"positive": [.5, 2, 4.5], "signals": [0, 0, 1]})
    add("moving_average_chart", [[1], [2], [3]], [2, 3, 0, 1], {"statistics": [1, 1.5, 2.5], "upper.1": 3 / math.sqrt(2)})
    add("western_electric", [[2.1], [2.2], [0]], [0, 1], {"rule_masks": [0, 0, 2]})
    gage = [[1, 0, 0], [3, 0, 0], [2, 0, 1], [4, 0, 1], [5, 1, 0], [7, 1, 0], [6, 1, 1], [8, 1, 1]]
    add("gage_rr", gage, [], {"repeatability_variance": 2, "operator_variance": .5, "part_variance": 8, "gage_variance": 2.5})
    add("type1_gage", [[9], [10], [11]], [10, 6], {"cg": .2, "cgk": .2})
    add("agreement_kappa", [[20, 5], [5, 20]], [], {"kappa": .6})
    add("acceptance_oc", [[0], [.1], [1]], [10, 0], {"acceptance_probability": [1, .9 ** 10, 0]})
    add("levey_jennings", [[0], [3.1]], [0, 1], {"lower": [-3, -3], "upper": [3, 3], "signals": [0, 1]})
    add("laney_p_chart", [[1, 10], [9, 10]], [], {"center": .5, "sigma_z": .8 / math.sqrt(.025) * math.sqrt(math.pi) / 2})
    add("laney_u_chart", [[1, 1], [9, 1]], [], {"center": 5, "upper.0": 5 + 12 * math.sqrt(math.pi)})

    add("bootstrap_mean", [[3], [3]], [100, .95, 42], {"estimate": 3, "bootstrap_standard_error": 0, "confidence_interval": [3, 3]})
    add("permutation_test", [[0, 0], [1, 0], [2, 1], [3, 1]], [0], {"statistic": -2, "p_value": 1 / 3, "assignments": 6, "exact": True})
    add("sample_rows", [[1], [2], [3]], [3, 0, 42], {})
    add("simulate_normal", [[2, 0]], [4, 42], {"samples": [2] * 4, "shape": [1, 4]})
    add("power_z", [[0, 1]], [25, .05, 0], {"power": [.05]})
    add("sample_size_mean", [[1]], [.2, .95], {"sample_sizes": [97]})
    add("sample_size_proportion", [[.5]], [.05, .95], {"sample_sizes": [385]})
    add("reliability_demonstration", [[.9]], [.95], {"sample_sizes": [29]})
    add("odds_ratio", [[20, 80], [10, 90]], [], {"estimate": 2.25})
    add("relative_risk", [[20, 80], [10, 90]], [], {"estimate": 2})
    add("cohen_d", [[1, 0], [2, 0], [3, 0], [3, 1], [4, 1], [5, 1]], [], {"cohen_d": -2, "pooled_sd": 1})

    sample = [[v] for v in [.5, 1, 1.5, 2.5, 4, 6]]
    for op, expected in [
        ("fit_normal", {"mean": 2.5833333333333335, "sigma": 1.9021187017522212, "log_likelihood": -12.371441425760134}),
        ("fit_lognormal", {"log_mean": .6344437482950532, "log_sigma": .836813848860365, "log_likelihood": -11.251371872859794}),
        ("fit_exponential", {"rate": .3870967741935484, "log_likelihood": -11.694483328182875}),
        ("fit_gamma", {"shape": 1.7371600335026973, "scale": 1.4871015240458112}),
        ("fit_weibull", {"shape": 1.3833562354097264, "scale": 2.839505353472658}),
        ("normality_test", {"statistic": .6664994693339898, "p_value": .7165912141396371}),
    ]:
        add(op, sample, [], expected)
    add("fit_poisson", [[0], [1], [2], [2], [3], [7]], [], {"rate": 2.5, "log_likelihood": -12.958854213301034})
    add("kde_gaussian", [[0], [1]], [.5, -1, 0, 1], {"densities": [.05412479673895295, .4529332469146207, .4529332469146207]})
    add("box_cox", sample, [0], {"transformed": [math.log(v[0]) for v in sample], "lambda": 0})
    censored = [[row[0], event] for row, event in zip(sample, [1, 0, 1, 1, 0, 1])]
    add("survival_exponential", censored, [], {"rate": .25806451612903225, "rate_standard_error": .12903225806451613})
    add("survival_weibull", censored, [], {"shape": 1.34901382272634, "scale": 3.8051173489356405})
    add("logrank_test", [[1, 1, 0], [2, 1, 0], [2, 0, 0], [2, 1, 1], [3, 1, 1], [4, 0, 1]], [], {"expected_group0": 1.3, "variance": .61, "statistic": .8032786885245902, "p_value": .3701148999172469})

    series = [[v] for v in [2, 4, 3, 5, 7, 6, 8, 9]]
    add("acf", series, [3], {"autocorrelations": [1, .48214285714285715, .2619047619047619, .05357142857142857]})
    add("pacf", series, [3], {"partial_autocorrelations": [1, .48214285714285715, .038360337903337476, -.1126752080124592]})
    add("difference", series, [2, 2], {"values": [3, 0, -3, 2], "original_indices": [4, 5, 6, 7]})
    add("moving_average", series, [3], {"values": [3, 4, 5, 6, 7, 23 / 3]})
    add("exponential_smoothing", series, [.4, 2], {"forecasts": [7.4719488] * 2})
    add("holt", series, [.4, .3, 3, .9], {"forecasts": [9.330366976752705, 9.917344847849883, 10.445624931837344]})
    add("ar_model", series, [2, 2], {"intercept": 1.07183908045977, "coefficients": [.4511494252873563, .617816091954023], "innovation_variance": 1.1135057471264371})
    add("periodogram", series, [2], {"frequencies": [0, .0625, .125, .1875, .25], "power_spectral_density": [0, 49.14213562373095, 10, 20.85786437626905, 4]})
    seasonal = [[10 + .25 * i + [0, 1, -2, 1][i % 4]] for i in range(12)]
    add("seasonal_decompose", seasonal, [4], {"seasonal_pattern": [0, 1, -2, 1], "trend.0": None, "trend.2": 10.5})
    add("ljung_box", series, [4, 2], {"q_statistic": 5.839434186372962, "degrees_of_freedom": 2, "p_value": .0539489476665005})

    add("covariance_matrix", [[1, 2], [2, 1], [3, 6]], [], {"covariance": [[1, 2], [2, 7]]})
    add("correlation_matrix", [[1, 2], [2, 1], [3, 6]], [], {"correlation": [[1, 2 / math.sqrt(7)], [2 / math.sqrt(7), 1]]})
    add("hierarchical_cluster", [[0], [1], [4], [5]], [2, 3], {"labels": [0, 0, 1, 1], "merges.2.2": math.sqrt(32)})
    add("mahalanobis", [[-1, -1], [-1, 1], [1, -1], [1, 1], [2, 0]], [4, 0], {"squared_distances": [3], "reference_mean": [0, 0]})
    classification = [[-2, 10], [0, 10], [2, 20], [4, 20], [-1, 99], [3, 99]]
    add("discriminant_analysis", classification, [4, 1, 0], {"predictions": [10, 20], "posterior_probabilities.0.0": 1 / (1 + math.exp(-4))})
    add("gaussian_naive_bayes", classification, [4, 1e-9], {"predictions": [10, 20], "posterior_probabilities.0.0": 1 / (1 + math.exp(-8))})
    add("knn_predict", [[-1, 9], [1, 3], [0, 99]], [2, 2, 1], {"predictions": [3], "neighbor_indices": [[0, 1]]})
    add("dbscan", [[0], [.1], [.2], [2], [2.1], [2.2], [9]], [.21, 3], {"labels": [0, 0, 0, 1, 1, 1, -1], "core_indices": [0, 1, 2, 3, 4, 5], "n_clusters": 2})

    balanced = [[y, g] for g, values in enumerate([[1, 2, 1, 2], [4, 5, 4, 5], [7, 8, 7, 8]]) for y in values]
    add("mixed_model", balanced, [], {"intercept": 4.5, "residual_variance": 1 / 3, "random_intercept_variance": 107 / 12, "log_likelihood": -15.490540829595014, "aic": None, "converged": True}, 3e-6)
    add("mixed_predict", [[0], [1], [99]], [0, 0, 4.5, 0, 2, 0, -1, 1, 1], {"marginal_predictions": [4.5] * 3, "conditional_predictions": [3.5, 5.5, 4.5], "known_group": [True, True, False]})
    cox = [[0, 1, 1, 1], [1, 0, 1, 1], [0, 0, 2, 0], [2, 1, 3, 1], [1, 2, 3, 0], [2, 0, 4, 1], [0, 2, 4, 0], [1, 1, 5, 1], [2, 2, 6, 1], [0, 1, 6, 0], [2, 0, 7, 1], [1, 2, 8, 0]]
    add("cox_ph", cox, [0], {"coefficients": [.22017852171500554, -.5010300380019524], "partial_log_likelihood": -11.911499078252731, "baseline_survival.2": .682398239031761})
    add("cox_predict", [[0], [1], [2]], [.8, math.log(2), 0, 0], {"relative_risks": [1, 2, 4], "survival_probabilities": [.8, .64, .4096]})
    add("relative_risk_predict", [[0], [1], [2]], [math.log(2), 0, 0], {"relative_risks": [1, 2, 4]})
    # Differenced data are [2,-1,2,2,-1,2,1], mean=1, ML variance=12/7.
    add("arima", series, [0, 1, 0, 2], {"sigma_squared": 12 / 7, "forecasts": [10, 11], "converged": True})
    tree_rows = [[-2, -4], [-1, -2], [1, 2], [2, 4], [-1.5, 99], [1.5, 99]]
    add("decision_tree", tree_rows, [4, 1, 1, 0], {"predictions": [-3, 3]})
    forest = [[0, 0], [1, 0], [2, 1], [3, 1], [-10, 99], [10, 99]]
    add("random_forest", forest, [4, 50, 3, 1, 1, 42], {"predictions": [0, 1]})
    add("gradient_boosting", [[0, 3], [1, 3], [2, 3], [3, 3], [1.5, 99]], [4, 5, .1, 2, 1, 0, 42], {"predictions": [3], "training_loss": [0] * 5})
    neural = [[i, 2 * i] for i in range(4)] + [[i, 99] for i in range(4)]
    add("mlp", neural, [4, 8, 4000, .02, 0, 42, 1e-10], {"predictions": [0, 2, 4, 6]}, 1e-6)

    add("box_behnken", [[-1, 1]] * 3, [], {"run_count": 13})
    add("latin_hypercube", [[-1, 1], [0, 10]], [4, 42], {"run_count": 4})
    add("mixture_lattice", [[0, 1]] * 3, [2], {"run_count": 6, "rows": [[0, 0, 1], [0, .5, .5], [0, 1, 0], [.5, 0, .5], [.5, .5, 0], [1, 0, 0]]})
    add("fractional_factorial", [[-1, 1]] * 4, [3, 1, 7], {"run_count": 9})
    add("orthogonal_screening", [[-1, 1]] * 5, [], {"run_count": 8})
    add("randomize_design", [[i] for i in range(6)], [2, 42, 2], {"run_count": 12, "replicate_labels": [0] * 6 + [1] * 6})
    add("design_evaluate", [[-1, 0], [0, 1], [1, 0], [2, 1], [-2, -1]], [], {"rank": 3, "log_determinant": math.log(60), "information_matrix": [[5, 0, 1], [0, 10, 4], [1, 4, 3]]})
    assert {case[0] for case in cases} == set(OPERATIONS)
    return cases, logistic, balanced, cox


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", type=Path, required=True)
    parser.add_argument("--extension", type=Path, required=True)
    args = parser.parse_args()
    runner = Runner(args.duckdb, args.extension)

    names = ",".join("'chem_" + op + "'" for op in OPERATIONS)
    catalog = runner.query(
        "SELECT function_name,examples FROM duckdb_functions() "
        f"WHERE function_name IN ({names}) ORDER BY function_name;"
    )
    assert {row["function_name"] for row in catalog} == {"chem_" + op for op in OPERATIONS}, "missing or unexpected new function registration"
    assert len(catalog) == 95, "unexpected overload count"
    examples = []
    for row in catalog:
        assert row["examples"] and row["examples"][0].startswith(row["function_name"] + "("), row
        examples.append(row["examples"][0])
    runner.many(examples)

    cases, logistic, balanced, cox = numerical_cases()
    results = runner.many([case[1] for case in cases])
    by_op = {}
    numerical_fields = 0
    for (op, _, expected, tolerance), result in zip(cases, results):
        by_op[op] = result
        for key, value in expected.items():
            equal(field(result, key), value, f"{op}.{key}", tolerance)
            numerical_fields += 1

    # Independent structural mathematics for stochastic/design outputs.
    assert sorted(by_op["sample_rows"]["indices"]) == [0, 1, 2]
    for j in range(2):
        assert sorted(row[j] for row in by_op["latin_hypercube"]["stratum_indices"]) == [0, 1, 2, 3]
    for row in by_op["fractional_factorial"]["coded_rows"][:-1]:
        equal(row[3], row[0] * row[1] * row[2], "fractional generator D=ABC")
    screen = by_op["orthogonal_screening"]["coded_rows"]
    for j in range(5):
        equal(sum(row[j] for row in screen), 0, "orthogonal column mean")
        for k in range(j):
            equal(sum(row[j] * row[k] for row in screen), 0, "orthogonal columns")
    randomized = by_op["randomize_design"]
    assert sorted(randomized["source_indices"][:6]) == list(range(6))
    assert sorted(randomized["source_indices"][6:]) == list(range(6))
    for index, row in zip(randomized["source_indices"], randomized["rows"]):
        assert row == [index]
    equal(sum(by_op["design_evaluate"]["leverages"]), 3, "design leverage trace")

    # Models remain in SQL as JSON: fit -> json_extract -> native typed lists.
    reuse_queries = [
        (
            "WITH fit AS (SELECT " + expression("logistic", logistic) + "::JSON j) "
            "SELECT chem_glm_predict([[0.0],[1.0]],json_extract(j,'$.stable_predict_parameters')::DOUBLE[]) AS result FROM fit;",
            "predictions", [.25, .75],
        ),
        (
            "WITH fit AS (SELECT " + expression("mixed_model", balanced) + "::JSON j) "
            "SELECT chem_mixed_predict([[0.0],[2.0],[99.0]],json_extract(j,'$.predict_parameters')::DOUBLE[]) AS result FROM fit;",
            "conditional_predictions", [4.5 - 107 / 36, 4.5 + 107 / 36, 4.5],
        ),
        (
            "WITH fit AS (SELECT " + expression("cox_ph", cox, [0]) + "::JSON j) "
            "SELECT chem_cox_predict([[1.0,1.0]],"
            "[json_extract(j,'$.baseline_survival[2]')::DOUBLE] || "
            "json_extract(j,'$.coefficients')::DOUBLE[] || "
            "json_extract(j,'$.feature_origins')::DOUBLE[] || "
            "json_extract(j,'$.mean_offsets')::DOUBLE[]) AS result FROM fit;",
            "survival_probabilities", [.682398239031761],
        ),
    ]
    for sql, key, expected in reuse_queries:
        row = runner.query(sql)[0]
        assert row["result"] is not None, sql
        result = strict_json(row["result"])
        equal(result[key], expected, f"SQL saved-model reuse {key}", 3e-6)

    malformed = [
        "chem_linear_model(NULL::DOUBLE[][],[])" ,
        "chem_linear_model([[0.0,1.0],[1.0,NULL],[2.0,4.0]],[])" ,
        "chem_logistic([[0.0,0.0],[1.0,0.5]],[])" ,
        "chem_poisson([[0.0,0.5],[1.0,1.0]],[])" ,
        "chem_glm_predict([[1000.0]],[2.0,0.0,1.0])" ,
        "chem_covariance_matrix([[1.0],[2.0]],NULL::DOUBLE[])" ,
        "chem_covariance_matrix([[1.0],[2.0,3.0]],[])" ,
        "chem_mahalanobis([[0.0,0.0],[1.0,1.0],[2.0,2.0],[1.0,2.0]],[3.0,0.0])" ,
        "chem_mixed_model([[1.0,0.0],[2.0,0.0],[3.0,0.0]],[])" ,
        "chem_cox_predict([[1.0]],[1.5,1.0,0.0,0.0])" ,
        "chem_arima([[1.0],[1.0],[1.0],[1.0]],[0.0,0.0,0.0,1.0])" ,
        "chem_box_behnken([[-1.0,1.0],[-1.0,1.0]],[])" ,
        "chem_mixture_lattice([[0.8,1.0],[0.8,1.0]],[2.0])" ,
        "chem_analyze('not_a_real_operation',[[1.0]],[])" ,
    ]
    null_rows = runner.query(" UNION ALL ".join(f"SELECT {i} id,{expr} IS NULL AS invalid" for i, expr in enumerate(malformed)) + ";")
    assert len(null_rows) == len(malformed) and all(row["invalid"] for row in null_rows), null_rows

    # Filtered selection vectors, nested-list validity, and dynamic operation names.
    selected = runner.query(
        "SELECT i,chem_glm_predict(CASE WHEN i=1994 THEN [[NULL::DOUBLE]] ELSE [[i::DOUBLE]] END,[0.0,3.0,2.0]) AS result "
        "FROM range(5000) t(i) WHERE i%997=0 ORDER BY i;"
    )
    for row in selected:
        if row["i"] == 1994:
            assert row["result"] is None
        else:
            equal(strict_json(row["result"])["predictions"], [3 + 2 * row["i"]], "selected nested-vector prediction")
    dynamic = runner.query(
        "SELECT i,chem_analyze(CASE WHEN i%2=0 THEN 'one_sample_t' ELSE 'variance_test' END,"
        "[[i::DOUBLE],[i::DOUBLE+1.0],[i::DOUBLE+2.0]],[1.0]) AS result "
        "FROM range(5000) t(i) WHERE i%997=0 ORDER BY i;"
    )
    for row in dynamic:
        result = strict_json(row["result"])
        key = "t_statistic" if row["i"] % 2 == 0 else "chi_square_statistic"
        expected = row["i"] * math.sqrt(3) if row["i"] % 2 == 0 else 2
        equal(result[key], expected, "selected dynamic-operation vector")

    # Every supervised prefix model must ignore the suffix dummy target/label.
    prefix = [[-2, 0], [-1, 0], [1, 1], [2, 1]]
    suffix_a = [[-1.5, 999], [1.5, 999]]
    suffix_b = [[-1.5, -999.123], [1.5, 321.5]]
    prefix_models = [
        ("discriminant_analysis", [4, 0, .1]),
        ("gaussian_naive_bayes", [4, 1e-9]),
        ("knn_predict", [4, 3, 1]),
        ("decision_tree", [4, 3, 1, 1]),
        ("random_forest", [4, 20, 3, 1, 1, 42]),
        ("gradient_boosting", [4, 20, .2, 1, 1, 1, 42]),
        ("mlp", [4, 4, 100, .03, 1, 42]),
    ]
    suffix_checks = []
    for op, parameters in prefix_models:
        suffix_checks.extend([expression(op, prefix + suffix_a, parameters), expression(op, prefix + suffix_b, parameters)])
    suffix_results = runner.many(suffix_checks)
    for i, (op, _) in enumerate(prefix_models):
        assert suffix_results[2 * i] == suffix_results[2 * i + 1], f"{op} learned from suffix labels"

    # Regression tests for the seven independently reviewed numerical defects.
    tiny_ewma = runner.one(expression("ewma_chart", [[2], [4]], [1e-20, 3, 0, 1]))
    equal(tiny_ewma["upper"][0] / 3e-20, 1, "tiny EWMA limit")
    assert tiny_ewma["signals"] == [0, 1]
    tiny_tree = runner.one(expression("decision_tree", [[0, 0], [1, 1e-200], [0, 0], [1, 0]], [2, 2, 1, 0]))
    equal(tiny_tree["predictions"][1] / 1e-200, 1, "tiny target tree")
    boots = runner.many([expression("bootstrap_mean", [[origin], [origin + 2]], [1000, .95, 42]) for origin in [0, 1e16]])
    equal(boots[0]["bootstrap_standard_error"], boots[1]["bootstrap_standard_error"], "bootstrap offset-safe SE")
    tiny_perm = runner.one(expression("permutation_test", [[0, 0], [1e-15, 0], [2e-15, 1], [3e-15, 1]], [0]))
    equal(tiny_perm["p_value"], 1 / 3, "permutation units")
    tiny_plan = runner.one(expression("sample_size_mean", [[1e-200]], [1e-200, .95]))
    assert tiny_plan["sample_sizes"] == [4]
    shifted_neural = [[i, 1e16 + 2 * i] for i in range(4)] + [[i, 0] for i in range(4)]
    shifted = runner.one(expression("mlp", shifted_neural, [4, 8, 4000, .02, 0, 42, 1e-10]))
    assert shifted["predictions"] == [1e16, 1e16 + 2, 1e16 + 4, 1e16 + 6]
    final_loss = runner.one(expression("mlp", [[i, 2 * i] for i in range(4)] + [[i, 0] for i in range(4)], [4, 3, 1, .1, 0, 42, 1e-10]))
    expected_loss = sum(.5 * ((prediction - 2 * i) / final_loss["target_scale"]) ** 2 / 4 for i, prediction in enumerate(final_loss["predictions"]))
    equal(final_loss["training_loss"], expected_loss, "returned MLP final model loss", 1e-10)

    print(json.dumps({
        "status": "PASS",
        "layer": "native DuckDB SQL -> C++ -> Rust",
        "catalog_examples": len(examples),
        "operations_with_numeric_or_structural_oracles": len(cases),
        "numeric_fields_checked": numerical_fields,
        "sql_model_reuse": len(reuse_queries),
        "invalid_inputs_returning_sql_null": len(malformed),
        "prefix_models_ignoring_suffix_labels": len(prefix_models),
        "reviewed_numeric_regressions": 7,
        "selected_vector_rows": len(selected) + len(dynamic),
        "native_sql_queries": runner.sql_queries,
        "json_results": runner.json_results,
    }, indent=2))


if __name__ == "__main__":
    main()
