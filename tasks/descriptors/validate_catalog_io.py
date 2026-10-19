"""Independent Python CSV check of actual Swift import/export and native values."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

source, output, catalog, report = map(Path, sys.argv[1:])
specs = json.loads(catalog.read_text())["descriptors"]
with source.open(newline="") as f:
    inputs = list(csv.DictReader(f))
with (output / "bulk-output.csv").open(newline="") as f:
    reader = csv.DictReader(f)
    headers = reader.fieldnames
    exported = list(reader)
results = json.loads((output / "bulk-native-results.json").read_text())
legacy = ["mw", "exact_mass", "logp", "tpsa", "hbd", "hba", "heavy_atoms",
          "bonds", "formal_charge", "fragments", "rotatable_bonds", "fraction_csp3"]
columns = legacy + [(s["id"] + "_stable") if s["id"] in ("Ipc", "AvgIpc") else s["id"] for s in specs[12:]]
assert headers == ["name", "smiles", "status", "formula"] + columns
assert len(inputs) == len(results) == len(exported) == 10_000
assert len({r["smiles"] for r in inputs}) == 10_000
numeric = invalid = 0
for index, (original, row, result) in enumerate(zip(inputs, exported, results)):
    assert result["index"] == index
    assert original["name"] == row["name"]
    assert original["smiles"] == row["smiles"]
    assert row["status"] == result["status"]
    assert row["formula"] == result.get("formula", "")
    if result["status"] != "ok":
        invalid += 1
        assert all(row[c] == "" for c in columns)
        continue
    values = result["descriptor_values"]
    assert len(values) == len(specs)
    for i, column in enumerate(columns):
        value = result[column] if i < 12 else values[i]
        if value is None:
            assert row[column] == "", (index, column)
        else:
            actual = float(row[column])
            assert math.isfinite(actual) and actual == value, (index, column)
            numeric += 1
        if i < 12:
            assert value == values[i], (index, column, "legacy/vector disagreement")
summary = dict(rows=len(inputs), distinct_raw_smiles=10_000,
               descriptor_columns=len(specs), csv_columns=len(headers),
               invalid_rows=invalid, numeric_cells_checked=numeric,
               source_text_and_order_preserved=True, invalid_values_blank=True,
               every_csv_value_matches_swift_result=True,
               stable_information_columns={"Ipc": "Ipc_stable", "AvgIpc": "AvgIpc_stable"},
               parser="Python standard-library csv, independent from Swift CompoundImport",
               input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
               csv_sha256=hashlib.sha256((output / "bulk-output.csv").read_bytes()).hexdigest())
report.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary))
