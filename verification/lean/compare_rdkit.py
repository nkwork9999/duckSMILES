"""Run the registered RDKit APIs and compare their actual results in Lean."""
import argparse
import csv
import json
import math
import platform
import subprocess
import sys
import tempfile
from collections import Counter

import numpy as np
from rdkit.Chem import Graphs
from pathlib import Path

from reference import CATALOG, extra50, graph_observations, observations, oracle, parse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
from gen_rdkit_fixtures import MOLECULES

HYDROGEN_CASES = ["CCO", "[H]O[H]", "[2H]O[2H]", "[H]C([H])([H])[H]", "[H]OC", "[H]N(C)C"]
INFORMATION_CASES = ["C" * n for n in (20, 40, 64, 80)]
UNSUPPORTED = ["C1", ""]
DOMAIN_CASES = ["C |$label$|"]


def raw_value(mol, name, functions):
    try:
        value = float(functions[name](oracle.Chem.Mol(mol)))
    except ZeroDivisionError:
        if name != "SPS" or mol.GetNumHeavyAtoms():
            raise
        return None, "ZeroDivisionError"
    if not math.isfinite(value):
        return None, "nonfinite"
    return value, "defined"


def compare(specs, inputs, functions, folder, default_h=False, information=False):
    rows = ([{"smiles": s, "observations": graph_observations(parse(s))} for s in inputs]
            if information else observations(inputs))
    actual = []
    states = {s["id"]: Counter() for s in specs}
    for row in rows:
        source = row["smiles"]
        if row["observations"] is None:
            actual.append({"smiles": source, "values": None})
            continue
        mol = oracle.Chem.MolFromSmiles(source) if default_h else parse(source)
        profile = oracle.profile_reference(mol) if any(s["origin"] == "native-profile" for s in specs) else {}
        graph_indices = extra50(mol) if any(s["origin"] == "native-extra" for s in specs) else {}
        values = []
        for spec in specs:
            name = spec["id"]
            if name in profile:
                x, status = profile[name], "defined"
            elif name in graph_indices:
                x = graph_indices[name]
                status = "defined" if x is not None else "undefined formula"
            else:
                x, status = raw_value(mol, name, functions)
            states[name][status] += 1
            values.append(x)
        actual.append({"smiles": source, "values": values})
    paths = [folder / f"{k}.json" for k in ("catalog", "observations", "actual")]
    for path, obj in zip(paths, [{"descriptors": specs}, rows, actual], strict=True):
        path.write_text(json.dumps(obj, allow_nan=False))
    command = [HERE / ".lake/build/bin/descriptor_check", "audit", *paths]
    differences = json.loads(subprocess.check_output(command, text=True))
    return states, differences


def path_coefficients(smiles):
    # det(xI-A) for a path follows P_n = x P_(n-1) - P_(n-2).
    n = len(smiles)
    prev, current = [1], [1, 0]
    for _ in range(2, n + 1):
        prev, current = current, [a - b for a, b in zip(current + [0], [0, 0] + prev, strict=True)]
    mol = parse(smiles)
    adj = np.equal(oracle.Chem.GetDistanceMatrix(mol), 1)
    computed = Graphs.CharacteristicPolynomial(mol, adj)
    differences = [i for i, v in enumerate(current) if v != float(computed[i])]
    return {"atoms": n, "different_coefficients": len(differences), "exact_constant": current[-1],
            "rdkit_constant": float(computed[-1])}


def run(output, library=None, check=False):
    specs = json.loads(CATALOG.read_text())["descriptors"]
    direct = [s for s in specs if s["origin"] == "rdkit-standard"]
    functions = dict(oracle.Descriptors._descList)
    if set(functions) != {s["id"] for s in direct}:
        raise ValueError("registered scalar inventory changed")
    if oracle.rdkit.__version__ != "2026.03.6":
        raise ValueError("requires RDKit 2026.03.6")
    inputs = list(dict.fromkeys(MOLECULES + ["[H]", "[H][H]", "[2H]O[2H]", "[NH4+]", "[Xe]"] + UNSUPPORTED))
    domain_differences = []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        states, diff = compare(specs, inputs, functions, work)
        _, hdiff = compare(direct, HYDROGEN_CASES, functions, work, default_h=True)
        ipc = [s for s in direct if s["id"] in ("Ipc", "AvgIpc")]
        _, idiff = compare(ipc, INFORMATION_CASES, functions, work, information=True)
        if library:
            native = oracle.Native(library.resolve())
            indexes = [next(i for i, spec in enumerate(native.catalog) if spec["id"] == key["id"]) for key in ipc]
            rust = [{"smiles": s, "values": [row["descriptor_values"][i] for i in indexes]}
                    for s, row in zip(INFORMATION_CASES, native.analyze(INFORMATION_CASES), strict=True)]
            rust_path = work / "rust.json"
            rust_path.write_text(json.dumps(rust, allow_nan=False))
            subprocess.run([HERE / ".lake/build/bin/descriptor_check", "check", work / "catalog.json",
                            work / "observations.json", rust_path], check=True)
        # Re-evaluate the smallest path control with an intentionally changed
        # coordinate. The audit must name that exact molecule and column.
        actual_path = work / "actual.json"
        actual = json.loads(actual_path.read_text())
        actual[0]["values"][0] += 1
        actual_path.write_text(json.dumps(actual))
        mutated = json.loads(subprocess.check_output(
            [HERE / ".lake/build/bin/descriptor_check", "audit", work / "catalog.json",
             work / "observations.json", actual_path], text=True))
        assert any(x["id"] == ipc[0]["id"] and x["smiles"] == INFORMATION_CASES[0] for x in mutated)
        if library:
            catalog_path = work / "catalog.json"
            input_path = work / "observations.json"
            catalog_path.write_text(json.dumps({"descriptors": specs}))
            input_path.write_text(json.dumps(observations(DOMAIN_CASES), allow_nan=False))
            modeled = json.loads(subprocess.check_output(
                [HERE / ".lake/build/bin/descriptor_check", "reference", catalog_path, input_path], text=True))
            actual = native.analyze(DOMAIN_CASES)
            for s, result, model in zip(DOMAIN_CASES, actual, modeled, strict=True):
                if (model["values"] is not None) != (result["status"] == "ok"):
                    domain_differences.append({"smiles": s, "lean_columns": len(model["values"] or []),
                                               "rust_status": result["status"]})
    counts = Counter(d["id"] for d in diff)
    hydrogen = Counter(d["id"] for d in hdiff)
    chains = Counter(d["id"] for d in idiff)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "rdkit-columns.csv").open("w", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(["id", "RDKit relationship", "reference defined", "reference nonfinite", "API exceptions",
                         "formula undefined", "same-graph discrepancies", "default-H discrepancies (6 inputs)",
                         "IPC chain discrepancies (4 inputs)"])
        for spec in specs:
            k = spec["id"]
            st = states[k]
            rel = {"rdkit-standard": "registered scalar", "native-profile": "composed from RDKit atom/bond/ring APIs",
                   "native-extra": "independent formula on RDKit graph"}[spec["origin"]]
            writer.writerow([k, rel, st["defined"], st["nonfinite"], st["ZeroDivisionError"], st["undefined formula"],
                             counts[k], hydrogen[k] if k in functions else "", chains[k] if spec in ipc else ""])
    report = {"rdkit": oracle.rdkit.__version__, "numpy": np.__version__, "python": platform.python_version(),
              "platform": platform.platform(), "source_rows": len(inputs),
              "valid_rows": sum(parse(s) is not None for s in inputs), "registered_scalars": len(direct),
              "invalid_inputs": [{"smiles": s, "rdkit_default_atom_count": None if (mol := oracle.Chem.MolFromSmiles(s)) is None
                                  else mol.GetNumAtoms(), "ducksmiles_profile": None} for s in UNSUPPORTED],
              "parser_domain_differences": domain_differences,
              "same_graph_differences": diff, "default_hydrogen_inputs": HYDROGEN_CASES,
              "default_hydrogen_differences": hdiff, "ipc_chain_lengths": list(map(len, INFORMATION_CASES)),
              "ipc_chain_differences": idiff, "ipc_rust_checked": bool(library),
              "ipc_path_coefficients": [path_coefficients(s) for s in INFORMATION_CASES], "sps_zero_division": [s for s in inputs if parse(s) is not None and not parse(s).GetNumHeavyAtoms()]}
    (output / "rdkit-differences.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if check:
        if diff:
            raise ValueError(f"{len(diff)} unexpected discrepancies on the curated same-graph corpus")
        if library and domain_differences != [{"smiles": DOMAIN_CASES[0], "lean_columns": len(specs),
                                                "rust_status": "invalid_or_unsupported"}]:
            raise ValueError("documented parser domain changed")
        water = next(item for item in hdiff if item["smiles"] == "[H]O[H]" and item["id"] == "Chi0")
        if water["observed"] != 0 or not math.isclose(water["modeled"], 2 + 1/math.sqrt(2)):
            raise ValueError("explicit hydrogen control changed")
        if report["sps_zero_division"] != ["[H]", "[H][H]"] or states["SPS"]["ZeroDivisionError"] != 2:
            raise ValueError("SPS zero-heavy-atom contract changed")
    print(len(direct), "registered functions;", sum(counts.values()), "same-graph discrepancies;",
          len(hdiff), "default-H discrepancies;", len(idiff), "IPC chain discrepancies")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--library", type=Path, help="also check the Rust information values on the chain controls")
    parser.add_argument("--check", action="store_true", help="fail on unexpected curated-corpus or domain differences")
    args = parser.parse_args()
    run(args.output, args.library, args.check)
