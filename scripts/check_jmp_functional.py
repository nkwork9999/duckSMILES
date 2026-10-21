#!/usr/bin/env python3
"""Validate functional preprocessing through native DuckDB SQL -> C++ -> Rust.

Standard library only. The eight documented examples have independently derived
known values. Additional checks cover saved references, validation-only queries,
malformed inputs and nonconstant/vectorized SQL execution. This checks the bounded
CPU method contracts, not numerical equivalence to an installed JMP executable.
"""
import argparse
import json
from pathlib import Path

from check_jmp_continuation import Runner, equal, expression, field, strict_json


ROOT = Path(__file__).resolve().parents[1]
API = ROOT / "tasks" / "jmp-functional-preprocess-api.json"


def inventory():
    items = json.loads(API.read_text())["operations"]
    assert len(items) == 8
    assert len({item["operation"] for item in items}) == 8
    for item in items:
        assert item["sql_name"] == "chem_" + item["operation"]
        fixture = item["fixture"]
        assert fixture["data"] and item["known_value"]
        assert item["sql_example"] == expression(
            item["operation"], fixture["data"], fixture["args"]
        ), item["operation"]
    return items


def saved_references():
    """Keep all saved x/y arrays in SQL; no Python carrier reconstruction."""
    grid_fit = expression("functional_grid", [[0, 0, 7], [1, 2, 7], [2, 4, 7]])
    grid = f"""(WITH fit AS (SELECT {grid_fit} AS model),
        saved AS (SELECT json_extract(model,'$.grid')::DOUBLE[] AS grid FROM fit)
        SELECT chem_functional_grid([[0.0,1.0,8.0,1.0],[2.0,5.0,8.0,1.0]],
            [len(grid)::DOUBLE,grid[1],grid[len(grid)],2.0]) FROM saved)"""
    snip_fit = expression(
        "functional_snip", [[0, 1, 7], [1, 1, 7], [2, 9, 7], [3, 1, 7], [4, 1, 7]], [1, 1, 0]
    )
    baseline = f"""(WITH fit AS (SELECT {snip_fit} AS model),
        saved AS (SELECT json_extract(model,'$.x[0]')::DOUBLE[] AS x,
            json_extract(model,'$.baselines[0]')::DOUBLE[] AS baseline FROM fit),
        packed AS (SELECT list_transform(range(1,len(x)+1),
            lambda i: [x[i],baseline[i],900.0,0.0]) AS reference FROM saved)
        SELECT chem_functional_baseline(list_concat(reference,
            [[0.0,4.0,800.0,0.0],[2.0,6.0,800.0,0.0],[4.0,8.0,800.0,0.0],
             [0.0,4.0,801.0,1.0],[2.0,6.0,801.0,1.0],[4.0,8.0,801.0,1.0]]),
            [900.0,2.0]) FROM packed)"""
    target_fit = expression("functional_grid", [[0, 1, 7], [1, 2, 7], [2, 3, 7]])
    golden = f"""(WITH fit AS (SELECT {target_fit} AS model),
        saved AS (SELECT json_extract(model,'$.grid')::DOUBLE[] AS x,
            json_extract(model,'$.values[0]')::DOUBLE[] AS y FROM fit),
        packed AS (SELECT list_concat(
            list_transform(range(1,len(x)+1), lambda i: [x[i],y[i],900.0,1.0]),
            list_transform(range(1,len(x)+1), lambda i: [x[i],y[i]+1.0,800.0,0.0])) AS curves
            FROM saved)
        SELECT chem_functional_golden_match(curves,[900.0,0.0,-1.0,0.0]) FROM packed)"""
    return [
        ("saved training grid / validation-only interpolation", grid, {
            "ids": [8], "roles": [1], "grid": [0, 1, 2], "values": [[1, 3, 5]],
            "grid_estimation": "supplied_bounds_and_size",
        }),
        ("saved SNIP baseline / training subtraction / validation bypass", baseline, {
            "ids": [800, 801], "roles": [0, 1], "baselines": [[1, 1, 1], [0, 0, 0]],
            "corrected": [[3, 5, 7], [4, 6, 8]], "baseline_applied": [True, False],
            "baseline_excluded": True,
        }),
        ("saved resampled reference / golden RMSE", golden, {
            "target_id": 900, "candidate_ids": [800], "candidate_roles": [0],
            "scores": [1], "ranking_ids": [800], "target_excluded": True,
            "model_fitted": False,
        }),
    ]


def invalid_inputs(items):
    invalid = []

    def add(name, op, data, args=()):
        invalid.append((name, expression(op, data, args)))

    # Each algorithm gets a distinct domain rejection, beyond generic null input.
    add("grid duplicate x", "functional_grid", [[0, 1, 1], [0, 2, 1]])
    add("automatic grid has no training functions", "functional_grid", [[0, 1, 1, 1], [1, 2, 1, 1]])
    add("bin reversed bounds", "functional_bin", [[0, 1, 1]], [2, 1, 0])
    add("unrepresentable bin edges", "functional_bin", [[1e16, 1, 1]], [4, 1e16, 1e16 + 2])
    add("thinning offset beyond stride", "functional_thin", [[0, 1, 1]], [2, 2, 0])
    add("SNIP loglog requires nonnegative y", "functional_snip", [[0, 1, 1], [1, -1, 1], [2, 1, 1]], [1, 1, 0])
    add("SNIP requires uniform training x", "functional_snip", [[0, 1, 1], [1, 2, 1], [3, 1, 1]], [1, 0, 0])
    add("arWLS needs three training points", "functional_arwls", [[0, 1, 1], [1, 2, 1]], [1])
    add("arWLS lambda cap", "functional_arwls", [[0, 1, 1], [1, 2, 1], [2, 1, 1]], [1e11])
    add("missing supplied baseline", "functional_baseline", [[0, 1, 1], [1, 2, 1]], [99])
    ragged = [[0, 1, 1], [1, 2, 1], [0, 1, 2], [1, 2, 2], [2, 3, 2]]
    add("DTW impossible window", "functional_dtw_matrix", ragged, [0])
    add("normalized DTW constant function", "functional_dtw_matrix", [[0, 1, 1], [1, 1, 1], [0, 2, 2], [1, 3, 2]], [-1, 1])
    add("missing golden reference", "functional_golden_match", ragged, [99])
    add("golden RMSE mismatched grids", "functional_golden_match", ragged, [1, 0])
    # Check SQL and nested-vector NULL validity masks through every overload.
    for item in items:
        fn = item["sql_name"]
        fixture = item["fixture"]
        valid = expression(item["operation"], fixture["data"], fixture["args"])
        args = valid[valid.rfind(", ") + 2:-1]
        data = valid[len(fn) + 1:valid.rfind(", ")]
        invalid.extend([
            (fn + " outer data NULL", f"{fn}(NULL::DOUBLE[][],{args})"),
            (fn + " args NULL", f"{fn}({data},NULL::DOUBLE[])"),
            (fn + " inner value NULL", f"{fn}([[0.0,NULL,1.0],[1.0,2.0,1.0]],{args})"),
        ])
    invalid.extend([
        ("ragged observation rows", "chem_functional_grid([[0.0,1.0,1.0],[1.0,2.0]],[])"),
        ("null observation row", "chem_functional_thin([NULL::DOUBLE[],[1.0,2.0,1.0]],[])"),
        ("null parameter element", "chem_functional_thin([[0.0,1.0,1.0]],[NULL::DOUBLE])"),
        ("NaN observation", "chem_functional_thin([[0.0,'NaN'::DOUBLE,1.0]],[])"),
        ("infinite observation", "chem_functional_thin([[0.0,'Infinity'::DOUBLE,1.0]],[])"),
    ])
    return invalid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True, type=Path)
    parser.add_argument("--extension", required=True, type=Path)
    options = parser.parse_args()
    runner = Runner(str(options.duckdb), options.extension)
    items = inventory()
    names = [item["sql_name"] for item in items]
    literals = ",".join("'" + name + "'" for name in names)
    catalog = runner.query(
        "SELECT function_name,description,examples,parameter_types FROM duckdb_functions() "
        f"WHERE function_name IN ({literals}) ORDER BY function_name;"
    )
    found = {row["function_name"] for row in catalog}
    assert len(catalog) == 8 and found == set(names), (
        "functional catalog mismatch", "expected 8", "found", len(catalog),
        "missing", sorted(set(names) - found),
    )
    by_name = {item["sql_name"]: item for item in items}
    for row in catalog:
        assert row["parameter_types"] == ["DOUBLE[][]", "DOUBLE[]"], row
        assert row["description"] and row["examples"], row
        assert by_name[row["function_name"]]["sql_example"] in row["examples"], row

    known_fields = 0
    for item, result in zip(items, runner.many([item["sql_example"] for item in items])):
        for path, wanted in item["known_value"].items():
            equal(field(result, path), wanted, item["operation"] + "." + path)
            known_fields += 1
    # Exact single-iteration arWLS fixture must retain the honest false flag.
    arwls = runner.one(by_name["chem_functional_arwls"]["sql_example"])
    equal(arwls["converged"], [False], "single iteration arWLS convergence")
    equal(arwls["weights"], [[1, 1, 1]], "returned baseline uses returned weights")

    saved = saved_references()
    saved_fields = 0
    for (label, _, expected), result in zip(saved, runner.many([case[1] for case in saved])):
        for path, wanted in expected.items():
            equal(field(result, path), wanted, label + "." + path)
            saved_fields += 1

    invalid = invalid_inputs(items)
    null_query = " UNION ALL ".join(
        f"SELECT {i} AS id, {sql} AS result" for i, (_, sql) in enumerate(invalid)
    ) + " ORDER BY id;"
    rejected = runner.query(null_query)
    assert len(rejected) == len(invalid)
    for i, row in enumerate(rejected):
        assert row["id"] == i and row["result"] is None, (invalid[i][0], row)

    dynamic = runner.query("""SELECT i,chem_analyze(
        CASE WHEN i%2=0 THEN 'functional_grid' ELSE 'functional_thin' END,
        [[0.0,i::DOUBLE,7.0],[1.0,i+1.0,7.0],[2.0,i+4.0,7.0]],[]::DOUBLE[]) AS result
        FROM range(5000)t(i) WHERE i%997=0 ORDER BY i;""")
    assert [row["i"] for row in dynamic] == [0, 997, 1994, 2991, 3988, 4985]
    for row in dynamic:
        result = strict_json(row["result"])
        runner.decoded += 1
        i = row["i"]
        equal(result["ids"], [7], "selected dynamic function id")
        if i % 2 == 0:
            equal(result["grid"], [0, 1, 2], "selected dynamic grid")
            equal(result["values"], [[i, i + 1, i + 4]], "selected dynamic values")
        else:
            equal(result["x"], [[0, 2]], "selected dynamic thinning x")
            equal(result["values"], [[i, i + 4]], "selected dynamic thinning values")

    interleaved = runner.query("""SELECT i,chem_functional_thin(
        CASE WHEN i%3=0 THEN NULL::DOUBLE[][] ELSE
        [[0.0,i::DOUBLE,7.0],[1.0,i+1.0,7.0],[2.0,i+2.0,7.0]] END,
        [2.0,0.0,0.0]) AS result FROM range(8)t(i) WHERE i<>4 ORDER BY i;""")
    assert [row["i"] for row in interleaved] == [0, 1, 2, 3, 5, 6, 7]
    for row in interleaved:
        i = row["i"]
        if i % 3 == 0:
            assert row["result"] is None, row
        else:
            result = strict_json(row["result"])
            runner.decoded += 1
            equal(result["values"], [[i, i + 2]], "interleaved NULL thinning")
            equal(result["original_row_indices"], [[0, 2]], "interleaved NULL row mapping")

    print(json.dumps({
        "status": "PASS",
        "layer": "native DuckDB SQL -> C++ -> Rust",
        "functional_catalog_examples": len(items),
        "api_known_fields_checked": known_fields,
        "arwls_nonconvergence_contract_checks": 2,
        "sql_saved_reference_roundtrips": len(saved),
        "saved_reference_fields_checked": saved_fields,
        "invalid_inputs_returning_sql_null": len(invalid),
        "selected_dynamic_rows": len(dynamic),
        "interleaved_null_rows": len(interleaved),
        "native_sql_queries": runner.queries,
        "strict_json_results": runner.decoded,
    }, indent=2))


if __name__ == "__main__":
    main()
