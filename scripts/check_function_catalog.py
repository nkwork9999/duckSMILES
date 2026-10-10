#!/usr/bin/env python3
"""Check catalog coverage and execute every scalar example using a DuckDB CLI.

Usage: python3 scripts/check_function_catalog.py --duckdb /path/to/duckdb \
    --extension /path/to/ducksmiles.duckdb_extension
The CLI and extension must be built for the same DuckDB version/platform.
"""

import argparse
import json
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", required=True, type=Path)
    args = parser.parse_args()
    extension = str(args.extension.resolve()).replace("'", "''")
    load = f"LOAD '{extension}';\n"

    def query(sql):
        result = subprocess.run(
            [args.duckdb, "-unsigned", "-batch", "-bail", "-json", ":memory:"],
            input=sql,
            text=True,
            capture_output=True,
            check=True,
            timeout=180,
        )
        return json.loads(result.stdout)

    rows = query(
        "CREATE TEMP TABLE before_load AS "
        "SELECT function_name, function_type, parameter_types FROM duckdb_functions();\n"
        + load
        + "SELECT function_name, function_type, parameters, parameter_types, "
        "description, examples, categories FROM duckdb_functions() f "
        "WHERE NOT EXISTS (SELECT 1 FROM before_load b "
        "WHERE b.function_name = f.function_name AND b.function_type = f.function_type "
        "AND b.parameter_types = f.parameter_types) "
        "ORDER BY function_name, len(parameters);"
    )
    if not rows:
        raise AssertionError("No extension functions discovered")
    examples = []
    for row in rows:
        signature = (row["function_name"], row["parameter_types"])
        if row["function_type"] != "scalar":
            raise AssertionError(f"Add explicit handling for non-scalar function: {signature}")
        if not (row["description"] or "").strip():
            raise AssertionError(f"Missing description: {signature}")
        names = row["parameters"]
        if len(names) != len(row["parameter_types"]) or any(
            not name or re.fullmatch(r"(?:col|arg)\d+", name) for name in names
        ):
            raise AssertionError(f"Missing meaningful parameter names: {signature}: {names}")
        if "ducksmiles" not in (row["categories"] or []):
            raise AssertionError(f"Missing extension category: {signature}")
        if not row["examples"]:
            raise AssertionError(f"Missing examples: {signature}")
        for example in row["examples"]:
            if not example.startswith(row["function_name"] + "("):
                raise AssertionError(f"Expected a scalar expression for {signature}: {example}")
            examples.append((signature, example))

    # Each documented example must bind, execute and produce a non-NULL result.
    checks = "\nUNION ALL\n".join(
        f"SELECT {i} AS example_id, ({example}) IS NOT NULL AS valid"
        for i, (_, example) in enumerate(examples)
    )
    results = query(load + checks + ";")
    if len(results) != len(examples):
        raise AssertionError("Example execution returned an unexpected row count")
    for result in results:
        if not result["valid"]:
            raise AssertionError(f"Example returned NULL: {examples[result['example_id']]}")
    print(f"PASS: {len(rows)} overloads documented; {len(examples)} examples executed")


if __name__ == "__main__":
    main()
