"""Quantify pinned numerical IPC vs independently exact graph definitions.

Graph-kernel probes, not native SMILES acceptance or a speed benchmark.
High-degree synthetic graphs need not be chemically valid carbon molecules.
Existing production descriptor tolerances are not modified by this report.
"""
import hashlib
import json
import math
import platform
import warnings
from pathlib import Path

import numpy
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

assert rdkit.__version__ == "2026.03.6"
assert numpy.__version__ == "2.5.4"
RDLogger.DisableLog("rdApp.*")
root = Path(__file__).resolve().parents[2]
source = root / "tasks/descriptors/exact-ipc-reference.json"
cases = json.loads(source.read_text())["cases"]
records = []
for case in cases:
    rw = Chem.RWMol()
    for _ in range(case["atoms"]):
        rw.AddAtom(Chem.Atom(6))
    for a, b in case["edges"]:
        rw.AddBond(a, b, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    validation = Chem.SanitizeMol(Chem.Mol(mol), catchErrors=True)
    record = dict(name=case["name"], atoms=case["atoms"], edges=len(case["edges"]),
                  carbon_graph_passes_sanitization=validation == Chem.SanitizeFlags.SANITIZE_NONE,
                  values={})
    for id, exact_key in [("Ipc", "exact_ipc"), ("AvgIpc", "exact_avg_ipc")]:
        expected = float(case[exact_key])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                actual = float(getattr(Descriptors, id)(mol))
                finite = math.isfinite(actual)
                error = abs(actual - expected) if finite else None
                record["values"][id] = dict(
                    exact_decimal=case[exact_key], pinned=actual if finite else None,
                    pinned_nonfinite=not finite, absolute_error=error,
                    relative_error=(error / abs(expected) if expected and error is not None else None),
                    within_report_threshold=finite and error <= max(1e-10, 1e-12 * abs(expected)),
                    warnings=[str(w.message) for w in caught])
            except Exception as error:
                record["values"][id] = dict(exact_decimal=case[exact_key], pinned_error=f"{type(error).__name__}: {error}",
                                             within_report_threshold=False)
    records.append(record)
report = dict(
    scope="58 declared graph-kernel probes; some high-degree synthetic graphs are not chemically valid carbon molecules; not native input parity or a speed benchmark",
    rdkit=rdkit.__version__, numpy=numpy.__version__, platform=platform.platform(),
    oracle_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    experimental_native_source_sha256=hashlib.sha256((root / "tasks/descriptors/experimental_exact_ipc.rs").read_bytes()).hexdigest(),
    existing_production_tolerances_changed=False,
    report_only_threshold="max(absolute1e-10, relative1e-12); not a shipped compatibility policy",
    cases=len(records), sanitized_carbon_graphs=sum(r["carbon_graph_passes_sanitization"] for r in records),
    above_report_threshold={id: sum(not r["values"][id]["within_report_threshold"] for r in records) for id in ["Ipc", "AvgIpc"]},
    records=records)
destination = root / "tasks/evidence/catalog289-sps/experimental-exact-ipc-pinned58.json"
destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
print(json.dumps({k: v for k, v in report.items() if k != "records"}, indent=2))
for record in records:
    if record["name"] in ["path-40", "path-80", "path-128", "path-256"]:
        print(json.dumps(record, allow_nan=False))
