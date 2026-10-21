#!/usr/bin/env python3
"""Native SQL validation for mixed models, causal inference and omics operations."""
import argparse
import json
import math
from pathlib import Path

from check_jmp_continuation import Runner, equal, expression, field, sql_list, strict_json

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", required=True, type=Path)
    args = parser.parse_args()
    runner = Runner(args.duckdb, args.extension)
    apis = {name: json.loads((ROOT/"tasks"/f"jmp-{name}-api.json").read_text())
            for name in ["mixed-advanced", "glmm", "causal", "omics", "market"]}
    operations = [item for api in apis.values() for item in api["operations"]]
    assert len(operations) == len({item["operation"] for item in operations}) == 22
    examples = [item["sql_example"] for item in operations]
    results = runner.many(examples)
    by_op = {item["operation"]: result for item, result in zip(operations, results)}
    names=",".join("'"+item["sql_name"]+"'" for item in operations)
    catalog=runner.query(f"SELECT function_name,description,examples,parameter_types FROM duckdb_functions() WHERE function_name IN ({names});")
    assert len(catalog)==22 and all(r["description"] and r["examples"] and r["parameter_types"]==["DOUBLE[][]","DOUBLE[]"] for r in catalog)
    numeric_fields = 0
    checked = set()

    def check(op, result, expected, tol=1e-7):
        nonlocal numeric_fields
        for key, wanted in expected.items():
            equal(field(result, key), wanted, f"{op}.{key}", tol)
            numeric_fields += 1
        checked.add(op)

    for item in apis["mixed-advanced"]["operations"]:
        check(item["operation"], by_op[item["operation"]], item["known_value"], 1e-5)
    # Each expected value is a literal or an analytic hand calculation,
    # separate from the broad NumPy/SciPy C-ABI reference scripts.
    for op, expected in [
        ("assay_normalize", {"rows": [[400000, 600000], [400000, 600000]]}),
        ("quantile_normalize", {"reference_quantiles": [1.5, 3.5, 5.5], "rows": [[5.5, 1.5, 3.5], [3.5, 1.5, 5.5]]}),
        ("quantile_normalize_predict", {"rows": [[4, 4, 1.5]]}),
        ("marker_admixture", {"allele_frequencies": [[.5, .5]], "log_likelihood": -9*math.log(2)}),
        ("admixture_predict", {"ancestry": [[1, 0], [.5, .5]]}),
        ("propensity_fit", {"intercept": -math.log(3), "coefficients": [2*math.log(3)]}),
        ("propensity_predict", {"propensities": [.25, .75]}),
        ("causal_balance", {"weighted_smd": [-1.4/math.sqrt(2)], "effective_sample_sizes": [25/13]*2}),
        ("propensity_match", {"estimate": 13/3, "standard_error": None, "control_indices": [3, 5, 5]}),
        ("ipw_effect", {"estimate": 2.6, "standard_error": math.sqrt(1.2288)}),
        ("aipw_effect", {"estimate": 2.25, "standard_error": .25}),
        ("uplift_curve", {"auuc": 1.4375, "qini_area": .9375}),
    ]:
        check(op, by_op[op], expected, 2e-6)
    # GLMM metadata carries independently integrated fixture literals.
    for item in apis["glmm"]["operations"]+apis["market"]["operations"]:
        if item.get("known_value"):
            check(item["operation"], by_op[item["operation"]], item["known_value"], 2e-5)
    assert len(checked) == 22, checked

    model_tests = []
    def saved(fit_op, predict_op, queries, expected, key="predict_parameters"):
        fit = next(item["sql_example"] for item in operations if item["operation"] == fit_op)
        sql = f"(WITH fit AS (SELECT {fit} AS model) SELECT chem_{predict_op}({sql_list(queries)},json_extract(model,'$.{key}')::DOUBLE[]) FROM fit)"
        model_tests.append((predict_op, sql, expected))
    saved("quantile_normalize", "quantile_normalize_predict", [[4, 4, 0]], {"rows": [[4.5, 4.5, 1.5]]}, "reference_quantiles")
    saved("marker_admixture", "admixture_predict", [[0, 2]], {"ancestry": [[1]]})
    saved("propensity_fit", "propensity_predict", [[0], [1], [2]], {"propensities": [.25, .75, 27/28]})
    saved("random_slope_model", "random_slope_predict", [[0, 0], [.5, 999]], {"known_subjects": [True, False]})
    saved("glmm_binomial", "glmm_predict", [[0], [999]], {"conditional_predictions": [.5,.5], "posterior_predictions": [.5,.5], "population_predictions": [.5,.5]})
    saved("glmm_poisson", "glmm_predict", [[0], [999]], {"conditional_predictions": [2,2], "posterior_predictions": [2,2], "population_predictions": [2,2]})
    for op, sql, expected in model_tests:
        check(op, runner.one(sql), expected, 2e-6)

    malformed = [
        expression("repeated_covariance", [[0], [1]], [4, 1, 1]),
        expression("repeated_model", [[1, 0, 0], [2, 0, 0]], [4, 2, 300, 1e-8, 1, .5]),
        expression("random_slope_predict", [[0, 1]], []),
        expression("assay_normalize", [[0, 0]], [1]),
        expression("assay_normalize", [[0, 1]], [4]),
        expression("quantile_normalize_predict", [[1, 2]], [2, 1]),
        expression("marker_admixture", [[-1, 0], [-1, 2]], [2]),
        expression("admixture_predict", [[1, 2]], [2, 2, 10, 1e-8, 1e-8, 1, 0, .1, .9, .9, .1]),
        expression("ipw_effect", [[4, 1, 0], [1, 0, .5]]),
        expression("aipw_effect", [[1, 1, .5, 0, 1], [2, 1, .5, 0, 1]]),
        expression("uplift_curve", [[1, 2, .5, .5]]),
        "chem_assay_normalize([[1.0,NULL]], [1.0])",
        "chem_assay_normalize([[1.0,2.0],[3.0]], [1.0])",
        "chem_propensity_fit(NULL::DOUBLE[][], [2.0])",
        "chem_glmm_predict([[1.0]], NULL::DOUBLE[])",
    ]
    bad = runner.query(" UNION ALL ".join(f"SELECT {i} AS id,({expr}) IS NULL AS valid" for i, expr in enumerate(malformed))+" ORDER BY id;")
    assert len(bad) == len(malformed) and all(row["valid"] for row in bad), bad
    selected = runner.query("SELECT i,chem_assay_normalize([[i::DOUBLE+1.0,i::DOUBLE+1.0]], [1.0]) AS result FROM range(5000) t(i) WHERE i%997=0 ORDER BY i;")
    for row in selected:
        equal(strict_json(row["result"])["rows"], [[500000, 500000]], "selected normalized rows")
        runner.decoded += 1
    interleaved = runner.query("SELECT i,chem_assay_normalize(CASE WHEN i%3=0 THEN NULL ELSE [[i::DOUBLE+1.0,i::DOUBLE+1.0]] END,[1.0]) AS result FROM range(8) t(i) WHERE i<>4 ORDER BY i;")
    for row in interleaved:
        if row["i"] % 3 == 0:
            assert row["result"] is None
        else:
            equal(strict_json(row["result"])["rows"], [[500000, 500000]], "interleaved NULL")
            runner.decoded += 1
    prefix = [[0, 0]]*3+[[0, 1], [1, 0]]+[[1, 1]]*3
    poisoned = runner.many([expression("propensity_fit", prefix+[[2, label]], [8]) for label in [-1e100, 1e100]])
    assert poisoned[0] == poisoned[1], "propensity suffix label leaked"
    print(json.dumps({"status": "PASS", "layer": "native DuckDB SQL -> C++ -> Rust",
                      "model_catalog_examples": len(examples), "operations_with_numerical_fixtures": len(checked),
                      "numeric_fields_checked": numeric_fields, "sql_saved_model_roundtrips": len(model_tests),
                      "invalid_inputs_returning_sql_null": len(malformed), "selected_rows": len(selected),
                      "interleaved_null_rows": len(interleaved), "prefix_models_ignoring_suffix_labels": 1,
                      "native_sql_queries": runner.queries, "strict_json_results": runner.decoded}, indent=2))


if __name__ == "__main__":
    main()
