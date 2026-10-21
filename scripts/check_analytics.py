#!/usr/bin/env python3
"""Exercise actual SQL -> C++ -> Rust analytics and JSON model reuse.

python3 scripts/check_analytics.py --duckdb /path/to/duckdb --extension /path/to/ducksmiles.duckdb_extension
Only Python standard-library modules are required.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", type=Path, required=True)
    args = parser.parse_args()
    extension = str(args.extension.resolve()).replace("'", "''")
    count = 0

    def query(sql):
        result = subprocess.run(
            [args.duckdb, "-unsigned", "-batch", "-bail", "-json", ":memory:"],
            input=f"LOAD '{extension}';\n{sql}", text=True, capture_output=True,
            check=True, timeout=120,
        )
        return json.loads(result.stdout)

    def analyze(expression):
        nonlocal count
        value = query(f"SELECT {expression} AS result;")[0]["result"]
        assert value is not None, expression
        result = json.loads(value, parse_constant=lambda value: (_ for _ in ()).throw(AssertionError(value)))
        assert result["version"] == 1
        count += 1
        return result

    def close(actual, expected, tolerance=1e-8):
        assert math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance), (actual, expected)

    pca = analyze("chem_pca([[1.,2.],[2.,4.],[3.,6.]], [1.,0.])")
    close(pca["explained_variance"][0], 5.)
    close(pca["explained_variance_ratio"][0], 1.)
    km = analyze("chem_kmeans([[0.],[1.],[9.],[10.]], [2.,100.,42.])")
    close(km["inertia"], 1.)
    assert km["labels"][0] == km["labels"][1] != km["labels"][2] == km["labels"][3]
    for operation in ("ridge", "pls"):
        model = analyze(f"chem_{operation}([[0.,1.],[1.,3.],[2.,5.]], [{'0.' if operation == 'ridge' else '1.'}])")
        close(model["intercept"], 1.)
        close(model["coefficients"][0], 2.)
        parameters = json.dumps([model["intercept"], *model["coefficients"]])
        predictions = analyze(f"chem_predict([[3.],[4.]], {parameters})")["predictions"]
        for actual, expected in zip(predictions, [7., 9.]):
            close(actual, expected)
    cv = analyze("chem_cross_validate([[0.,1.,0.],[1.,3.,0.],[2.,5.,1.],[3.,7.,1.]], [0.,0.])")
    close(cv["rmse"], 0.)
    folds = analyze("chem_scaffold_folds(['c1ccccc1','Cc1ccccc1','C1CCCCC1','CCO','CCN'], [2.,42.])")["folds"]
    assert folds[0] == folds[1] and folds[3] == folds[4]
    assert set(folds) == {0, 1}
    # Mirrored stereochemistry must not split the same ring skeleton into
    # different evaluation groups. The untouched fused-ring stereocenters
    # survive native Murcko extraction unless grouping explicitly removes them.
    mirrored = analyze("chem_scaffold_folds(['C1C[C@H]2CCCC[C@@H]2C1','C1C[C@@H]2CCCC[C@H]2C1','c1ccccc1','C1CCCCC1'], [2.,42.])")
    assert mirrored["n_groups"] == 3
    assert mirrored["folds"][0] == mirrored["folds"][1]
    group = analyze("chem_group_folds([[1.],[1.],[2.],[3.]], [2.,42.])")
    assert group["folds"][0] == group["folds"][1]
    assert len(analyze("chem_factorial_design([[10.,20.],[1.,3.]], [1.])")["rows"]) == 5
    assert len(analyze("chem_central_composite_design([[10.,20.],[1.,3.]], [1.,1.])")["rows"]) == 9
    optimal = analyze("chem_d_optimal_design([[-1.],[-0.5],[0.5],[1.]], [2.,20.,42.])")
    assert sorted(optimal["indices"]) == [0, 3]
    surface = analyze("chem_response_surface([[-1.,2.],[0.,1.],[1.,4.]], [0.])")
    parameters = json.dumps(surface["model_parameters"])
    response = analyze(f"chem_response_surface_predict([[2.]], {parameters})")["predictions"][0]
    close(response, 11.)
    gp = analyze("chem_gp_predict([[0.,2.],[0.,0.],[100.,0.]], [1.,1.,1.,1.,0.,0.])")
    close(gp["variances"][0], .5)
    assert gp["ranked_candidate_indices"] == [1, 0]
    ei = analyze("chem_acquisition([[2.,1.]], [2.,0.,0.])")
    close(ei["expected_improvement"][0], 1 / math.sqrt(2 * math.pi))
    selected = analyze("chem_diversity_select([[0.],[1.],[5.],[10.]], [2.,1.])")
    assert selected["selected_indices"] == [3, 2]
    for name, data, parameters in (
        ("summarize", "[[1.],[2.],[3.]]", "[]::DOUBLE[]"),
        ("welch_t_test", "[[1.,0.],[2.,0.],[3.,0.],[2.,1.],[3.,1.],[4.,1.]]", "[]::DOUBLE[]"),
        ("one_way_anova", "[[1.,0.],[2.,0.],[3.,0.],[2.,1.],[3.,1.],[4.,1.]]", "[]::DOUBLE[]"),
        ("adjust_pvalues", "[[0.01],[0.04],[0.03]]", "[1.]"),
        ("process_capability", "[[1.],[2.],[3.]]", "[-1.,5.]"),
        ("kaplan_meier", "[[1.,1.],[2.,0.],[3.,1.]]", "[]::DOUBLE[]"),
        ("random_intercept", "[[1.,0.],[3.,0.],[5.,1.],[7.,1.]]", "[]::DOUBLE[]"),
    ):
        analyze(f"chem_{name}({data}, {parameters})")

    # Every invalid row must become SQL NULL, without poisoning neighboring rows.
    malformed = [
        "chem_pca(NULL::DOUBLE[][], [1.,0.])",
        "chem_pca([[1.,NULL],[2.,3.]], [1.,0.])",
        "chem_pca([[1.],[2.,3.]], [1.,0.])",
        "chem_pca([[1.],[2.]], NULL::DOUBLE[])",
        "chem_ridge([[1.,1.],[1.,2.]], [0.])",
        "chem_analyze('missing', [[1.]], []::DOUBLE[])",
        "chem_scaffold_folds(['invalid','CCO'], [2.,42.])",
    ]
    for expression in malformed:
        assert query(f"SELECT {expression} IS NULL AS valid;")[0]["valid"], expression
    result = query("SELECT i, chem_analyze('acquisition', [[i::DOUBLE,1.]], [0.,0.,0.]) AS result FROM range(5000) t(i) WHERE i % 997 = 0 ORDER BY i;")
    for row in result:
        ei = json.loads(row["result"])["expected_improvement"][0]
        close(ei, row["i"] if row["i"] else 1 / math.sqrt(2 * math.pi))
    # Actual molecular descriptors feed the matrix, in stable ID order.
    rows = query("WITH molecules(id,smiles) AS (VALUES (1,'CCO'),(2,'CCN'),(3,'CCC'),(4,'c1ccccc1')) SELECT chem_pca(list([mol_weight(smiles),logp_crippen(smiles)] ORDER BY id), [2.,1.]) AS result FROM molecules;")
    assert len(json.loads(rows[0]["result"])["scores"]) == 4
    print(f"PASS: {count} numerical/model SQL checks; {len(malformed)} invalid inputs; selected-vector rows; molecule descriptor integration")


if __name__ == "__main__":
    main()
