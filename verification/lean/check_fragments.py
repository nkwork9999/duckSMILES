"""Exercise every compiled SMARTS on the pinned RDKit positive and negative fixtures."""
import json
import subprocess
import tempfile
from pathlib import Path

from reference import CATALOG, oracle, graph_observations

HERE = Path(__file__).resolve().parent


def main():
    specs = [s for s in json.loads(CATALOG.read_text())["descriptors"] if s["id"].startswith("fr_")]
    names = [s["id"] for s in specs]
    fixtures = HERE.parents[1] / "crates/smiles/tests/fixtures/descriptor_fragments.tsv"
    expected = {}
    for line in fixtures.read_text().splitlines():
        if line.startswith("#"):
            continue
        smiles, name, count = line.split("\t")
        row = expected.setdefault(smiles, {})
        if name in row:
            raise ValueError(f"duplicate {smiles}, {name}")
        row[name] = int(count)
    assert all(set(row) == set(names) for row in expected.values())
    assert all(any(row[name] > 0 for row in expected.values()) for name in names)
    inputs, outputs = [], []
    for smiles, row in expected.items():
        mol = oracle.Chem.MolFromSmiles(smiles, oracle.PARAMS)
        inputs.append({"smiles": smiles, "observations": graph_observations(mol)})
        outputs.append({"smiles": smiles, "values": [row[name] for name in names]})
    with tempfile.TemporaryDirectory() as tmp:
        paths = [Path(tmp) / (k + ".json") for k in ("catalog", "observations", "actual")]
        for path, data in zip(paths, [{"descriptors": specs}, inputs, outputs], strict=True):
            path.write_text(json.dumps(data, allow_nan=False))
        subprocess.run([HERE / ".lake/build/bin/descriptor_check", "check", *paths], check=True)


if __name__ == "__main__":
    main()
