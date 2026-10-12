#!/usr/bin/env python3
"""Native graph SQL/FFI smoke and regression checks; no training or real data."""
import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import re

import _duckdb
import duckdb


def sql_value(value):
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return str(value).lower()
    return str(value) if value != "" else "(empty)"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    path = args.extension.resolve()
    ctypes.CDLL(_duckdb.__file__, mode=ctypes.RTLD_GLOBAL)
    con = duckdb.connect(config={"allow_unsigned_extensions": True,
                                 "autoinstall_known_extensions": False,
                                 "autoload_known_extensions": False})
    con.execute("LOAD '" + str(path).replace("'", "''") + "'")
    root = Path(__file__).resolve().parents[1]
    checks = 0
    for name in ("native_graph.test", "function_catalog.test"):
        for block in re.split(r"\n\s*\n", (root / "test/sql" / name).read_text()):
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

    def graph(smiles):
        raw = con.execute("SELECT mol_graph_json(?)", [smiles]).fetchone()[0]
        return json.loads(raw) if raw is not None else None

    invalid = [None, "", "C=", "C()", "C..N", ".C", "C.", "C1CC", "C11", "C1C1",
               "[Xx]", "[*]", "[D]", "[99999C]", "[C:4294967296]", "C=1CCCCC1",
               "C\x00C", "炭素", "C" * 257, "C" * 4097]
    for text in invalid:
        assert graph(text) is None, text
    charged = graph("CC[NH3+]")
    neutral = graph("CCN")
    assert charged["atoms"][2]["formal_charge"] == 1
    assert charged["atoms"][2]["bracket_h_count"] == 3
    assert neutral["atoms"][2]["implicit_h_count"] == 2
    assert charged != neutral
    explicit = graph("[H]N([2H])[H]")
    assert [a["atomic_number"] for a in explicit["atoms"]] == [1, 7, 1, 1]
    assert explicit["atoms"][2]["isotope"] == 2
    assert explicit["atoms"][1]["explicit_h_neighbor_count"] == 3
    chiral = graph("[H][13C@@H:7](N)C(=O)O")["atoms"][1]
    assert (chiral["isotope"], chiral["atom_map"], chiral["chirality_token"]) == (13, 7, "@@")
    assert chiral["stereo_neighbor_order"] == [0, -1, 2, 3]
    stereo = graph("F/C=C\\F")
    assert [b["direction"] for b in stereo["bonds"]] == [1, 0, -1]
    assert stereo["input_smiles"] == "F/C=C\\F"
    assert graph("CCO")["atoms"][0]["symbol"] == "C"
    assert graph("OCC")["atoms"][0]["symbol"] == "O"
    assert all(a["aromatic"] for a in graph("C1=CC=CC=C1")["atoms"])
    assert graph("[NH4+].[Cl-]")["bonds"] == []
    assert graph("[Og]")["atoms"][0]["atomic_number"] == 118
    assert graph("C" * 256)["atom_count"] == 256
    expected_computed = {"native_aromaticity": True, **dict.fromkeys(
        ["cip", "ez", "gasteiger", "hybridization", "ring_basis", "radical_electrons"], False)}
    assert neutral["computed"] == expected_computed
    assert neutral["pretrained_tensor_compatible"] is False
    assert neutral["schema"] == "ducksmiles-native-graph-v1"
    assert neutral["backend"] == "ducksmiles-native"
    assert neutral["backend_version"] == "0.3.0"
    # More than one vector, nullable branches, and repeated constants.
    count = con.execute("""SELECT count(*) FILTER (WHERE mol_graph_json(
        CASE WHEN i % 3 = 0 THEN NULL WHEN i % 3 = 1 THEN 'CN' ELSE 'C=' END
        ) IS NOT NULL) FROM range(5000) t(i)""").fetchone()[0]
    assert count == 1667
    first = con.execute("SELECT mol_graph_json('N[C@@H](C)C(=O)O')").fetchone()[0]
    assert con.execute("SELECT bool_and(mol_graph_json('N[C@@H](C)C(=O)O') = ?) FROM range(2500)", [first]).fetchone()[0]
    catalog = con.execute("SELECT function_name, parameters, description, examples FROM duckdb_functions() WHERE list_contains(categories, 'ducksmiles')").fetchall()
    examples = 0
    for name, parameters, description, samples in catalog:
        assert description and samples and all(parameters), name
        for sample in samples:
            assert con.execute("SELECT (" + sample + ") IS NOT NULL").fetchone()[0], name
            examples += 1
    receipt = {"status": "passed", "scope": "synthetic_native_graph_sql_only",
               "duckdb_version": duckdb.__version__, "platform": con.execute("PRAGMA platform").fetchone()[0],
               "sql_regression_queries": checks, "invalid_cases": len(invalid),
               "catalog_overloads": len(catalog), "catalog_examples_executed": examples,
               "extension_size": path.stat().st_size,
               "extension_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
               "model_accuracy_measured": False, "pretrained_tensor_parity_verified": False}
    if args.receipt:
        args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
