#!/usr/bin/env python3
"""Verify estimated spatial models through native DuckDB SQL -> C++ -> Rust.

Uses the independent full-covariance SciPy reference fixtures already recorded
by check_jmp_spatial_estimation.py. This runner needs only the Python stdlib.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from check_jmp_continuation import Runner, equal, expression, sql_list, strict_json

ROOT = Path(__file__).resolve().parents[1]
API = ROOT / "tasks/jmp-spatial-estimation-api.json"
REFERENCE = ROOT / "tasks/evidence/jmp-spatial-estimation-reference.json"
IID_DATA = [[0, 0, 1], [1, 0, 2], [2, 0, 6]]
IID_ARGS = [1, 0, 0, 0, 0, 100, 1e-8, 0.1, 0.5]
IID_CARRIER = [1, 1, 0, 3, 0, 14 / 3, 0, 0.1, 0, 0, 1, 1, 0, 2, 2, 0, 6]

# Independently derived by differentiating the analytic full-covariance score
# at the SciPy optimum (kind0 spherical ML, first reference fixture). Negative
# lower Wald limits are intentional: these are untruncated interior intervals,
# not JMP's bounded Satterthwaite intervals or nonnegative parameter estimates.
INTERIOR_ML_ORACLE = {
    "covariance_parameter_standard_errors": [
        0.2916658455073921, 0.24058105950879358, 0.39557936891583184,
    ],
    "covariance_parameter_covariance": [
        [0.08506896543554193, -0.047725267591548166, -0.0178989331635952],
        [-0.047725267591548166, 0.05787924619437368, 0.024478267994216867],
        [-0.0178989331635952, 0.024478267994216867, 0.1564830371118478],
    ],
    "covariance_parameter_intervals": [
        [-0.25748076212939247, 0.8858283433004317],
        [-0.07451191809121127, 0.8685485059082345],
        [0.09729609389398408, 1.6479387260982117],
    ],
}


def numeric_leaves(value):
    if isinstance(value, (list, tuple)):
        return sum(numeric_leaves(v) for v in value)
    return int(isinstance(value, (int, float)) and not isinstance(value, bool))


def replay(fit_expression, query):
    return (
        "(WITH fit AS (SELECT " + fit_expression + " AS model) "
        "SELECT chem_spatial_fit_predict(" + sql_list(query)
        + ",json_extract(model,'$.predict_parameters')::DOUBLE[]) FROM fit)"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    runner = Runner(args.duckdb, args.extension)
    api = json.loads(API.read_text())
    reference = json.loads(REFERENCE.read_text())
    assert reference["status"] == "PASS"
    operations = api["operations"]
    assert [o["operation"] for o in operations] == ["spatial_fit", "spatial_fit_predict"]
    names = [o["sql_name"] for o in operations]
    catalog = runner.query(
        "SELECT function_name,description,examples,parameters,parameter_types "
        "FROM duckdb_functions() WHERE function_name IN "
        "('chem_spatial_fit','chem_spatial_fit_predict') ORDER BY function_name;"
    )
    assert [r["function_name"] for r in catalog] == names, catalog
    for row, operation in zip(catalog, operations):
        assert row["description"] and row["examples"], row
        assert row["parameter_types"] == ["DOUBLE[][]", "DOUBLE[]"], row
        assert row["parameters"] == ["data", "parameters"], row
        assert operation["sql_example"] in row["examples"], row
    runner.many([o["sql_example"] for o in operations])

    fixtures = reference["fixtures"]
    assert len(fixtures) == 16
    numerical_fields = 0
    expressions = [expression("spatial_fit", f["input"], f["args"]) for f in fixtures]
    results = runner.many(expressions)
    for fixture, result in zip(fixtures, results):
        for key, wanted in fixture["expected"].items():
            equal(result[key], wanted, f"independent spatial fit {fixture['kind']}/{fixture['reml']}/{fixture['args'][4]}.{key}", 5e-4)
            numerical_fields += numeric_leaves(wanted)
        equal(result["estimation"], "REML" if fixture["reml"] else "ML", "estimator label")
        assert result["converged"] and result["covariance_parameters_identified"], result
        assert result["projected_score_passed"]
        assert len(result["predict_parameters"]) <= 16384
    for key, wanted in INTERIOR_ML_ORACLE.items():
        equal(results[0][key], wanted, "independent interior inference." + key, 3e-4)
        numerical_fields += numeric_leaves(wanted)
    assert results[0]["covariance_parameter_intervals"][0][0] < 0
    assert results[0]["covariance_parameter_intervals"][1][0] < 0
    assert results[0]["partial_sill"] > 0 and results[0]["nugget"] > 0
    assert results[0]["ci_method"] == "interior_delta_method_untruncated_normal_Wald_95pct_not_JMP_bounded_Satterthwaite"

    closed = runner.many([
        expression("spatial_fit", IID_DATA, IID_ARGS),
        expression("spatial_fit", IID_DATA, [1, 0, 0, 1, 0, 100, 1e-8, 0.1, 0.5]),
    ])
    for result, variance in zip(closed, [14 / 3, 7]):
        equal(result["intercept"], 3, "iid trend")
        equal(result["total_variance"], variance, "iid ML/REML residual df")
        assert not result["covariance_parameters_identified"]
        assert result["covariance_parameter_intervals"] is None
        numerical_fields += 2

    saved = runner.many([
        replay(expression("spatial_fit", IID_DATA, IID_ARGS), [[10, 0]]),
        replay(expression("spatial_fit", [[0, 1, 1], [0, 2, 2], [0, 3, 6]], [1, 0, 1, 1, 0, 100, 1e-8, 0.1, 2]), [[0, 4]]),
    ])
    for result, variance in zip(saved, [56 / 9, 28 / 3]):
        for key, wanted in {"predictions": [3], "latent_variances": [variance], "observation_variances": [variance]}.items():
            equal(result[key], wanted, "saved model SQL roundtrip." + key)
            numerical_fields += numeric_leaves(wanted)

    # Saved interior ML model: SQL extracts its full covariance carrier, rather
    # than supplying hand-copied model coefficients to the prediction function.
    interior_replay = runner.one(replay(expressions[0], [[6, -0.27941549819892586, 1], [0, 0, 3]]))
    expected_interior = {
        "predictions": [3.1422515283332073, 3.271121013218532],
        "latent_variances": [0.29385644683392065, 0.3629487965554875],
        "observation_variances": [0.6908747407424323, 0.7599670904639991],
    }
    for key, wanted in expected_interior.items():
        equal(interior_replay[key], wanted, "independent saved interior prediction." + key, 3e-5)
        numerical_fields += numeric_leaves(wanted)

    limited_args = list(fixtures[2]["args"])
    limited_args[5] = 1
    limited = runner.one(expression("spatial_fit", fixtures[2]["input"], limited_args))
    assert not limited["optimization_converged"] and not limited["converged"]
    assert limited["covariance_parameter_intervals"] is None

    malformed = [
        "chem_spatial_fit(NULL::DOUBLE[][],[1,0,0,0,0])",
        "chem_spatial_fit([[0,0,1],[1,0,2],[2,0,6]],NULL::DOUBLE[])",
        "chem_spatial_fit([[0,0,NULL],[1,0,2],[2,0,6]],[1,0,0,0,0])",
        "chem_spatial_fit([[0,0,1],[1,0,2],[2,0,6]],[1,0,0,0,NULL])",
        "chem_spatial_fit([[0,0,1],[1,0],[2,0,6]],[1,0,0,0,0])",
        "chem_spatial_fit([[0,0,2],[1,0,2],[2,0,2]],[1,0,1,1,1])",
        "chem_spatial_fit([[0,1,0,1],[1,1,0,2],[2,1,0,4]],[1,1,1,1,1])",
        "chem_spatial_fit([[0,0.5,1],[1,0,2],[2,0,6]],[1,0,0,0,0])",
        "chem_spatial_fit([[0,0,1],[1,0,2],[2,0,6]],[1,0,0,0,0,100,1e-8,0.1])",
        "chem_spatial_fit_predict([[0,0]],[])",
        "chem_spatial_fit_predict(NULL::DOUBLE[][],[1])",
        "chem_spatial_fit_predict([[0,0]],[2,1,0,3,0,4.666666666666667,0,0.1,0,0,1,1,0,2,2,0,6])",
        "chem_spatial_fit_predict([[0,0.5]],[1,1,0,3,0,4.666666666666667,0,0.1,0,0,1,1,0,2,2,0,6])",
        "chem_spatial_fit_predict([[0,0]],[1,1,0,3,0,-1,0,0.1,0,0,1,1,0,2,2,0,6])",
        "chem_spatial_fit_predict([[0,0]],[1,1,0,3,0,4.666666666666667,0,0.1,0,0,1,1,0,2,2,0])",
    ]
    invalid = runner.query(" UNION ALL ".join(
        f"SELECT {i} AS id,({sql}) IS NULL AS valid" for i, sql in enumerate(malformed)
    ) + " ORDER BY id;")
    assert len(invalid) == len(malformed) and all(r["valid"] for r in invalid), invalid

    fit_args = sql_list(IID_ARGS)
    carrier = sql_list(IID_CARRIER)
    # Filtered nonconstant rows plus a per-row operation/shape/carrier switch.
    selected = runner.query(
        "SELECT i,chem_analyze(CASE WHEN i%2=0 THEN 'spatial_fit' ELSE 'spatial_fit_predict' END,"
        "CASE WHEN i%2=0 THEN [[i::DOUBLE,0.,1.],[i+1.,0.,2.],[i+2.,0.,6.]] ELSE [[i+10.,0.]] END,"
        f"CASE WHEN i%2=0 THEN {fit_args} ELSE {carrier} END) AS result "
        "FROM range(5000)t(i) WHERE i%997=0 ORDER BY i;"
    )
    for row in selected:
        result = strict_json(row["result"])
        if row["i"] % 2 == 0:
            equal(result["total_variance"], 14 / 3, "selected dynamic fit")
        else:
            equal(result["predictions"], [3], "selected dynamic saved prediction")
            equal(result["observation_variances"], [56 / 9], "selected dynamic variance")
        runner.decoded += 1
    # A materialized relation and repeated hash-join matches exercise selected
    # nested-list rows through the same callback used for dictionary vectors.
    dictionary = runner.query(
        "WITH source AS MATERIALIZED (SELECT i,[[i::DOUBLE,0.,1.],[i+1.,0.,2.],[i+2.,0.,6.]] AS data FROM range(4)t(i)) "
        f"SELECT j,chem_spatial_fit(data,{fit_args}) AS result FROM source JOIN range(12)d(j) "
        "ON i=j%4 WHERE j%3<>0 ORDER BY j;"
    )
    for row in dictionary:
        result = strict_json(row["result"])
        equal(result["total_variance"], 14 / 3, "joined dictionary fit")
        runner.decoded += 1
    null_rows = runner.query(
        "SELECT i,chem_spatial_fit(CASE WHEN i%3=0 THEN NULL::DOUBLE[][] ELSE "
        "[[i::DOUBLE,0.,1.],[i+1.,0.,2.],[i+2.,0.,6.]] END,"
        f"{fit_args}) AS result FROM range(8)t(i) WHERE i<>4 ORDER BY i;"
    )
    for row in null_rows:
        if row["i"] % 3 == 0:
            assert row["result"] is None, row
        else:
            equal(strict_json(row["result"])["total_variance"], 14 / 3, "interleaved SQL NULL")
            runner.decoded += 1
    null_parameters = runner.query(
        "SELECT i,chem_spatial_fit_predict([[i+10.,0.]],CASE WHEN i%3=0 THEN NULL::DOUBLE[] ELSE "
        f"{carrier} END) AS result FROM range(8)t(i) WHERE i<>4 ORDER BY i;"
    )
    for row in null_parameters:
        if row["i"] % 3 == 0:
            assert row["result"] is None, row
        else:
            equal(strict_json(row["result"])["observation_variances"], [56 / 9], "interleaved saved-carrier NULL")
            runner.decoded += 1

    record = {
        "schema_version": 1, "status": "PASS",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "layer": "actual native DuckDB SQL -> extension C++ vector adapter -> Rust",
        "duckdb_version": runner.query("SELECT version() AS version;")[0]["version"],
        "catalog_functions": 2, "catalog_examples_executed": 2,
        "independent_ml_reml_fits": len(fixtures),
        "numerical_reference_scalar_fields": numerical_fields,
        "saved_model_sql_roundtrips": len(saved) + 1,
        "invalid_inputs_returning_sql_null": len(malformed),
        "selected_dynamic_rows": len(selected), "join_dictionary_rows": len(dictionary),
        "interleaved_null_rows": len(null_rows),
        "interleaved_null_parameter_rows": len(null_parameters),
        "nonconverged_fit_withholding_intervals": 1,
        "negative_untruncated_Wald_lower_bounds_checked": 2,
        "native_sql_queries": runner.queries, "strict_json_results": runner.decoded,
        "extension_sha256": hashlib.sha256(args.extension.read_bytes()).hexdigest(),
        "reference_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256((ROOT / "crates/analytics/src/spatial_estimation.rs").read_bytes()).hexdigest(),
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "limitations": ["bounded isotropic ML/REML only; no whole-family JMP parity", "Wald intervals are explicitly untruncated and not JMP bounded Satterthwaite intervals", "plug-in prediction omits covariance estimation uncertainty"],
    }
    if args.output:
        args.output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
