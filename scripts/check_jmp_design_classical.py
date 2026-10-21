#!/usr/bin/env python3
"""Independent counting checks for six classical DOE SQL functions (stdlib only).

Requires matching native CLI/extension. This verifies plan structure, incidence,
restricted randomization and catalog OA balance; it does not analyze responses.
"""
import argparse
from collections import Counter, defaultdict
import itertools
import json
from pathlib import Path

from check_jmp_continuation import Runner, equal, expression, strict_json

ROOT = Path(__file__).resolve().parents[1]


def blocks(result, values, k, replication, concurrence, resolvable=False):
    """Count directly from returned rows/IDs, not returned summary fields."""
    n = len(result["design"])
    assert n == len(result["treatment_ids"]) == len(result["block_ids"])
    actual = defaultdict(list)
    for row, treatment, block in zip(result["design"], result["treatment_ids"], result["block_ids"]):
        assert row == values[treatment]
        actual[block].append(treatment)
    assert len(actual) == result["n_blocks"]
    visits = Counter()
    pair = Counter()
    for treatments in actual.values():
        assert len(treatments) == len(set(treatments)) == k
        visits.update(treatments)
        pair.update(itertools.combinations(sorted(treatments), 2))
    v = len(values)
    assert visits == Counter({i: replication for i in range(v)}), visits
    assert all(pair[a, b] == concurrence for a, b in itertools.combinations(range(v), 2)), pair
    equal(result["replication"], [replication]*v, "reported replications")
    expected = [[replication if i == j else concurrence for j in range(v)] for i in range(v)]
    equal(result["pair_concurrence"], expected, "reported concurrence")
    equal(result["treatment_information"],
          [[(replication if i == j else 0)-expected[i][j]/k for j in range(v)] for i in range(v)],
          "block adjusted information")
    assert result["run_order"] == list(range(1, n+1))
    if resolvable:
        classes = defaultdict(list)
        for cls, treatment in zip(result["parallel_class_ids"], result["treatment_ids"]):
            classes[cls].append(treatment)
        assert len(classes) == replication
        assert all(sorted(x) == list(range(v)) for x in classes.values())


def restricted(result, factors, replicas, strip):
    p = len(factors)
    roles = [int(f[2]) for f in factors]
    a = 2**roles.count(0)
    b = 2**roles.count(1)
    expected = set(itertools.product(*[(f[0], f[1]) for f in factors]))
    rows = result["design"]
    by_replica = defaultdict(list)
    for rep, row in zip(result["replicate_ids"], rows):
        by_replica[rep].append(tuple(row))
    assert len(rows) == replicas*a*b and len(by_replica) == replicas
    assert all(len(x) == len(expected) and set(x) == expected for x in by_replica.values())
    assert result["inference_ready"] == (replicas > 1)
    if strip:
        row_ids, col_ids = result["row_strip_ids"], result["column_strip_ids"]
        assert len(set(zip(row_ids, col_ids))) == len(rows), "strips must cross"
        for ids, role in [(row_ids, 0), (col_ids, 1)]:
            settings = defaultdict(set)
            for ident, row in zip(ids, rows):
                settings[ident].add(tuple(x for x, r in zip(row, roles) if r == role))
            assert all(len(x) == 1 for x in settings.values()), "strip settings changed"
        for key, ids, length in [("stage1_order", row_ids, b), ("stage2_order", col_ids, a)]:
            order = result[key]
            assert sorted(order) == list(range(len(rows)))
            for j in range(0, len(order), length):
                assert len({ids[i] for i in order[j:j+length]}) == 1, key
        equal(result["replication_error_df"],
              [(replicas-1)*(a-1), (replicas-1)*(b-1), (replicas-1)*(a-1)*(b-1)], "strip df")
    else:
        ids = result["whole_plot_ids"]
        for j in range(0, len(rows), b):
            assert len(set(ids[j:j+b])) == 1
            assert len({tuple(x for x, role in zip(row, roles) if role == 0) for row in rows[j:j+b]}) == 1
        equal(result["replication_error_df"], [(replicas-1)*(a-1), a*(replicas-1)*(b-1)], "split df")


def oa(result, levels):
    rows = result["design"]
    n = len(rows)
    for j, level_count in enumerate(levels):
        assert Counter(row[j] for row in rows) == Counter({i: n//level_count for i in range(level_count)})
        for k in range(j):
            counts = Counter((row[j], row[k]) for row in rows)
            assert counts == Counter({(a, b): n//(level_count*levels[k])
                                      for a in range(level_count) for b in range(levels[k])})
    assert result["strength2"] and result["main_effect_rank"] == 1+sum(l-1 for l in levels)
    gram = result["main_effect_gram"]
    equal(gram, [[int(i == j) for j in range(len(gram))] for i in range(len(gram))], "OA main contrast Gram")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", type=Path, required=True)
    args = parser.parse_args()
    runner = Runner(args.duckdb, args.extension)
    api = json.loads((ROOT/"tasks/jmp-design-classical-api.json").read_text())
    items = api["operations"]
    assert len(items) == len({x["operation"] for x in items}) == 6
    results = runner.many([x["sql_example"] for x in items])
    for item, result in zip(items, results):
        for key, value in item["known_value"].items():
            equal(result[key], value, item["operation"]+"."+key)
    by_op = {item["operation"]: result for item, result in zip(items, results)}
    blocks(by_op["rcbd_design"], items[0]["fixture"]["data"], 4, 3, 3)
    blocks(by_op["balanced_incomplete_block"], items[1]["fixture"]["data"], 3, 3, 1)
    blocks(by_op["resolvable_incomplete_block"], items[2]["fixture"]["data"], 3, 4, 1, True)
    restricted(by_op["split_plot_design"], items[3]["fixture"]["data"], 2, False)
    restricted(by_op["strip_plot_design"], items[4]["fixture"]["data"], 2, True)
    oa(by_op["taguchi_orthogonal_array"], [2]+[3]*7)
    structural = 6
    for v, k, lam in [(3, 2, 3), (6, 3, 2), (8, 4, 3), (9, 4, 3)]:
        data = [[i] for i in range(v)]
        blocks(runner.one(expression("balanced_incomplete_block", data, [k, lam, 42])),
               data, k, lam*(v-1)//(k-1), lam)
        structural += 1
    for strip in [False, True]:
        op = "strip_plot_design" if strip else "split_plot_design"
        factors = [[-1, 1, 0], [3, 5, 1]]
        restricted(runner.one(expression(op, factors, [1, 19])), factors, 1, strip)
        structural += 1
    array = runner.one(expression("taguchi_orthogonal_array", [[2]]*3, [4, 7]))
    oa(array, [2]*3)
    equal(array["max_main_interaction_alias"], 1, "L4 known complete main/interaction alias")
    structural += 1
    # Seed-dependent data comes through selected dictionary vectors; NULL rows
    # are interleaved with non-NULL selected rows.
    selected = runner.query("SELECT i,chem_rcbd_design([[i::DOUBLE],[i::DOUBLE+1]], [2.,i::DOUBLE]) AS result FROM range(5000) t(i) WHERE i%997=0 ORDER BY i;")
    for row in selected:
        result = strict_json(row["result"])
        blocks(result, [[row["i"]], [row["i"]+1]], 2, 2, 2)
        runner.decoded += 1
    interleaved = runner.query("SELECT i,chem_rcbd_design(CASE WHEN i%3=0 THEN NULL ELSE [[0.],[1.]] END,[2.,i::DOUBLE]) AS result FROM range(9) t(i) WHERE i<>4 ORDER BY i;")
    for row in interleaved:
        if row["i"] % 3 == 0:
            assert row["result"] is None
        else:
            blocks(strict_json(row["result"]), [[0], [1]], 2, 2, 2)
            runner.decoded += 1
    malformed = [
        expression("rcbd_design", [[0], [0]], [2]),
        expression("rcbd_design", [[0], [1]], [2, .5]),
        expression("balanced_incomplete_block", [[i] for i in range(6)], [3, 1]),
        expression("balanced_incomplete_block", [[i] for i in range(4)], [4, 1]),
        expression("resolvable_incomplete_block", [[i] for i in range(16)]),
        expression("split_plot_design", [[0, 1, 0]]),
        expression("strip_plot_design", [[0, 1, 0], [2, 1, 1]]),
        expression("taguchi_orthogonal_array", [[2], [2]], [18]),
        "chem_rcbd_design(NULL::DOUBLE[][],[2.])",
        "chem_rcbd_design([[0.],[1.]],NULL::DOUBLE[])",
        "chem_rcbd_design([[0.],[NULL]],[2.])",
        "chem_rcbd_design([[0.],[1.,2.]],[2.])",
        "chem_rcbd_design([NULL::DOUBLE[],[1.]],[2.])",
        "chem_rcbd_design([[0.],[1.]],[2.,'NaN'::DOUBLE])",
    ]
    invalid = runner.query(" UNION ALL ".join(f"SELECT {i} AS id,({e}) IS NULL AS valid" for i, e in enumerate(malformed))+" ORDER BY id;")
    assert len(invalid) == len(malformed) and all(r["valid"] for r in invalid), invalid
    print(json.dumps(dict(status="PASS", layer="native DuckDB SQL -> C++ -> Rust",
                          catalog_examples=6, independent_structural_cases=structural,
                          malformed_inputs_returning_null=len(malformed), selected_rows=len(selected),
                          interleaved_rows=len(interleaved), native_sql_queries=runner.queries,
                          strict_json_results=runner.decoded), indent=2))


if __name__ == "__main__":
    main()
