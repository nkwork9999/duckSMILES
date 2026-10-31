#!/usr/bin/env python3
"""Check chemistry SQL regressions and independent numerical/composition oracles.

Requires an existing CLI matching the supplied extension. Installs nothing.
The decimal constants are fixed reference values, not read from production code.
"""

import argparse
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import random
import subprocess


# token: (element, isotope, neutral isotope mass, average/labelled weight)
REFERENCE = {
    "C": ("C", None, "12", "12.011"),
    "[13C]": ("C", 13, "13.00335484", "13.00335484"),
    "H": ("H", None, "1.007825032", "1.008"),
    "[2H]": ("H", 2, "2.014101778", "2.014101778"),
    "N": ("N", None, "14.003074", "14.007"),
    "[15N]": ("N", 15, "15.0001089", "15.0001089"),
    "O": ("O", None, "15.99491462", "15.999"),
    "[18O]": ("O", 18, "17.999161", "17.999161"),
    "Cl": ("Cl", None, "34.96885268", "35.453"),
    "[37Cl]": ("Cl", 37, "36.96590259", "36.96590259"),
    "He": ("He", None, "4.002603254", "4.003"),
    "Li": ("Li", None, "7.01600455", "6.941"),
}
SQL_FILES = [
    "chemistry_tools.test", "mol_structure_integrity.test",
    "formula_monoisotopic_mass.test", "mol_input_validation.test",
    "mass_context.test", "function_catalog.test",
]


def display(value):
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return str(value).lower()
    return str(value) if value != "" else "(empty)"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    extension = args.extension.resolve()
    load = "LOAD '" + str(extension).replace("'", "''") + "';\n"

    def query(sql):
        result = subprocess.run(
            [args.duckdb, "-unsigned", "-batch", "-bail", "-json", ":memory:"],
            input=load + sql, text=True, capture_output=True, timeout=180,
        )
        if result.returncode:
            raise AssertionError(result.stderr + "\n" + sql[:1000])
        return json.loads(result.stdout)

    root = Path(__file__).resolve().parents[1]
    sql_checks = 0
    for name in SQL_FILES:
        lines = (root / "test/sql" / name).read_text().splitlines()
        i = 0
        while i < len(lines):
            if not lines[i].startswith("query "):
                if lines[i].startswith("statement "):
                    raise AssertionError("Add explicit statement handling: " + name)
                i += 1
                continue
            header = lines[i]
            i += 1
            sql = []
            # Empty lines may be inside a quoted MOL/SDF fixture.
            while lines[i] != "----":
                sql.append(lines[i])
                i += 1
            i += 1
            expected = []
            while i < len(lines) and lines[i].strip():
                expected.append(lines[i])
                i += 1
            actual = ["\t".join(map(display, row.values())) for row in query("\n".join(sql))]
            if "rowsort" in header:
                actual.sort()
                expected.sort()
            assert actual == expected, (name, sql, actual, expected)
            sql_checks += 1

    rng = random.Random(20261031)
    formulas, expected = [], []
    for _ in range(1000):
        terms = [(rng.choice(list(REFERENCE)), rng.randrange(1, 40))
                 for _ in range(rng.randrange(1, 12))]
        counts = Counter()
        mono, weight = Decimal(0), Decimal(0)
        for token, count in terms:
            symbol, isotope, mass, average = REFERENCE[token]
            counts[symbol, isotope] += count
            mono += Decimal(mass) * count
            weight += Decimal(average) * count
        composition = [dict(element=symbol, isotope=isotope, count=counts[symbol, isotope])
                       for symbol, isotope in sorted(counts, key=lambda key: (key[0], key[1] or 0))]
        formulas.append("".join(token + str(count) for token, count in terms))
        expected.append((float(mono), float(weight), composition))
    values = ",".join(f"({i},'{formula}')" for i, formula in enumerate(formulas))
    rows = query("SELECT id, formula_monoisotopic_mass(f) AS mono, formula_weight(f) AS weight, "
                 "formula_composition_json(f) AS composition FROM (VALUES " + values + ") t(id,f) ORDER BY id;")
    assert len(rows) == len(expected)
    mono_errors, weight_errors = [], []
    for row, (mono, weight, composition) in zip(rows, expected):
        mono_errors.append(abs(row["mono"] - mono))
        weight_errors.append(abs(row["weight"] - weight))
        assert json.loads(row["composition"]) == composition, row
    assert max(mono_errors) < 1e-9, max(mono_errors)
    assert max(weight_errors) < 1e-9, max(weight_errors)

    # Construct formulas with independently chosen DBE, then derive H counts.
    dbe_cases = []
    for i in range(1000):
        c, n, o = rng.randrange(1, 41), rng.randrange(9), rng.randrange(9)
        dbe = rng.randrange(c + 1)
        x = rng.randrange(min(3, 2 * c + 2 + n - 2 * dbe) + 1)
        h = 2 * c + 2 + n - x - 2 * dbe
        tokens = [("[13C]" if i % 2 else "C", c), ("[2H]" if i % 3 == 0 else "H", h),
                  ("N", n), ("O", o), (rng.choice(["F", "Cl", "Br", "I"]), x)]
        dbe_cases.append(("".join(token + str(count) for token, count in tokens if count), dbe))
    values = ",".join(f"({i},'{f}')" for i, (f, _) in enumerate(dbe_cases))
    rows = query("SELECT id, formula_dbe(f) AS dbe FROM (VALUES " + values + ") t(id,f) ORDER BY id;")
    assert len(rows) == len(dbe_cases)
    for row, (_, expected_dbe) in zip(rows, dbe_cases):
        assert row["dbe"] == expected_dbe, row

    ppm_cases = []
    with localcontext() as context:
        context.prec = 100
        for _ in range(1000):
            reference = rng.uniform(1, 9) * 10.0 ** rng.randrange(-250, 251)
            observed = reference * (1.0 + rng.uniform(-100, 100) / 1e6)
            exact_reference = Decimal.from_float(reference)
            expected_ppm = float((Decimal.from_float(observed) - exact_reference) / exact_reference * 1_000_000)
            ppm_cases.append((observed, reference, expected_ppm))
    values = ",".join(f"({i},{obs!r}::DOUBLE,{ref!r}::DOUBLE)" for i, (obs, ref, _) in enumerate(ppm_cases))
    rows = query("SELECT id, mass_error_ppm(observed, reference) AS ppm FROM (VALUES " + values + ") t(id,observed,reference) ORDER BY id;")
    assert len(rows) == len(ppm_cases)
    ppm_errors = [abs(row["ppm"] - expected_ppm) for row, (_, _, expected_ppm) in zip(rows, ppm_cases)]
    assert max(ppm_errors) < 1e-9, max(ppm_errors)

    receipt = {
        "status": "passed", "duckdb_version": query("SELECT version() AS version;")[0]["version"],
        "sql_files": SQL_FILES, "sql_queries_passed": sql_checks,
        "independent_cases": {"monoisotopic_mass": 1000, "average_weight": 1000,
                              "composition": 1000, "dbe": 1000, "ppm": 1000},
        "maximum_absolute_errors": {"monoisotopic_mass_Da": max(mono_errors),
                                    "average_weight_Da": max(weight_errors), "ppm": max(ppm_errors)},
        "oracle_seed": 20261031,
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
    }
    if args.receipt:
        args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
