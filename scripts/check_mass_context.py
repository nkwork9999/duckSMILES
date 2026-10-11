#!/usr/bin/env python3
"""Exercise the built extension using an existing, version-matched DuckDB Python.

No package is installed or updated. Run with --extension /path/to/extension.
The two SQL regression files, catalog examples and independent decimal mass
arithmetic are all exercised through the real SQL/FFI/Rust boundary.
"""

import argparse
import ctypes
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import random
import re

import _duckdb
import duckdb


ATOMS = dict(zip(
    ["C", "H", "N", "O", "S", "P", "Cl", "Br", "F", "I", "Si", "Se", "Na", "K", "B", "As", "Fe", "D"],
    map(Decimal, ["12", "1.00782503207", "14.0030740048", "15.9949146196", "31.97207100", "30.97376163",
                  "34.96885268", "78.9183371", "18.99840322", "126.904473", "27.9769265325", "79.9165213",
                  "22.9897692809", "38.96370668", "11.0093054", "74.9215965", "55.9349375", "2.0141017778"]),
))


def decimal_formula(formula):
    return sum(ATOMS[element] * int(count or 1)
               for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula))


def decimal_adduct(neutral_mass, label):
    # Independent grammar-based oracle; production uses an explicit lookup.
    n, terms, z, sign = re.fullmatch(r"\[(\d*)M(.*?)\](\d*)([+-])", label).groups()
    delta = Decimal(0)
    for op, multiplicity, formula in re.findall(r"([+-])(\d*)([A-Za-z][A-Za-z0-9]*)", terms):
        delta += (1 if op == "+" else -1) * int(multiplicity or 1) * decimal_formula(formula)
    charge = int(z or 1) * (1 if sign == "+" else -1)
    delta -= Decimal(charge) * Decimal("0.00054857990")
    return (Decimal(str(neutral_mass)) * int(n or 1) + delta) / abs(charge)


def sql_value(v):
    if v is None:
        return "NULL"
    if isinstance(v, bool):
        return str(v).lower()
    return str(v) if v != "" else "(empty)"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    ctypes.CDLL(_duckdb.__file__, mode=ctypes.RTLD_GLOBAL)
    con = duckdb.connect(config={"allow_unsigned_extensions": True})
    path = args.extension.resolve()
    con.execute("LOAD '" + str(path).replace("'", "''") + "'")
    root = Path(__file__).resolve().parents[1]
    checks = 0
    for name in ["mass_context.test", "function_catalog.test"]:
        text = (root / "test/sql" / name).read_text()
        for block in re.split(r"\n\s*\n", text):
            lines = [line for line in block.splitlines() if not line.startswith("#")]
            if not lines or not lines[0].startswith("query "):
                continue
            split = lines.index("----")
            actual = ["\t".join(map(sql_value, row))
                      for row in con.execute("\n".join(lines[1:split])).fetchall()]
            expected = lines[split + 1:]
            if "rowsort" in lines[0]:
                actual.sort()
                expected.sort()
            assert actual == expected, (name, lines[1:split], actual, expected)
            checks += 1

    # Deterministic N-present and 18-element coverage, no molecular labels used.
    rng = random.Random(20261011)
    formulas = list(ATOMS) + ["C8H10N4O2", "C5H5N", "CH3COOH", "D2O"]
    for _ in range(1000):
        elements = rng.sample(list(ATOMS), rng.randrange(1, 8))
        if "N" not in elements:
            elements.append("N")
        formulas.append("".join(e + str(rng.randrange(1, 40)) for e in elements))
    maximum = 0.0
    for formula in formulas:
        actual = con.execute("SELECT formula_exact_mass(?)", [formula]).fetchone()[0]
        error = abs(actual - float(decimal_formula(formula)))
        maximum = max(maximum, error)
        assert error < 1e-9, (formula, actual)

    source = (root / "crates/smiles/src/mass_context.rs").read_text()
    labels = re.findall(r'"(\[[^"\n]+\][+-]|\[[^"\n]+\]\d[+-])"', source.split("pub const SUPPORTED_ADDUCTS:")[1].split("];", 1)[0])
    assert len(labels) == len(set(labels)) == 29
    for label in labels:
        for neutral in [50.0, 194.0803755791, 1000.0]:
            expected = float(decimal_adduct(neutral, label))
            actual = con.execute("SELECT adduct_precursor_mz(?, ?)", [neutral, label]).fetchone()[0]
            assert abs(actual - expected) < 1e-10, (label, actual, expected)
            restored = con.execute("SELECT adduct_neutral_mass(?, ?)", [expected, label]).fetchone()[0]
            assert abs(restored - neutral) < 1e-10, (label, restored, neutral)

    catalog = con.execute("SELECT function_name, parameters, description, examples FROM duckdb_functions() WHERE list_contains(categories, 'ducksmiles')").fetchall()
    examples = 0
    for name, parameters, description, samples in catalog:
        assert description and samples and all(parameters)
        for sample in samples:
            assert con.execute("SELECT (" + sample + ") IS NOT NULL").fetchone()[0], name
            examples += 1
    receipt = {
        "status": "passed", "duckdb_version": duckdb.__version__,
        "platform": con.execute("PRAGMA platform").fetchone()[0],
        "sql_regression_queries": checks, "formula_decimal_comparisons": len(formulas),
        "maximum_formula_absolute_error": maximum,
        "adduct_labels": len(labels), "adduct_decimal_cases_each_direction": 3 * len(labels),
        "catalog_overloads": len(catalog), "catalog_examples_executed": examples,
        "extension_size": path.stat().st_size,
        "extension_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    if args.receipt:
        args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
