"""Independent prepared-graph oracle, original notation plus atom permutations.

Never calls native acceptance. Includes all seven retained partial-AR residuals.
Use SanitizeMol (not public canonical Kekulize) to match the parser pipeline.
"""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from compare_catalog import PARAMS, profile_reference

assert rdkit.__version__ == "2026.03.6"
root = Path(__file__).resolve().parents[2]
old = root / "tasks/evidence/catalog289-sps/v32-error-cohort415/differences.jsonl.gz"
residuals = [r for r in map(json.loads, gzip.open(old, "rt"))
             if "status" not in r["differences"] and r["index"] != 6]
assert [r["index"] for r in residuals] == [47, 87, 183, 217, 286, 305, 406]
sources = [r["smiles"] for r in residuals]
sources += [
    "C1:C:C:C:1", "c1ccc1", "C1:C:C:C:C:C:C:C:1", "c1ccccccc1",
    "C1:C:C:C:C:C:1", "c1ccccc1", "C1=CC=CC=C1", "C1=CC#CC=C1",
    "C1:C:C:C:C:C#1", "C1:C:C:C:C:C=1", "C:C", "C:N", "C:O",
    "[nH]1cccccc1", "[nH]1cccc1", "N1C=CC=CC=C1", "[NH]1C=CC=CC=C1",
    "[n+]1ccccc1", "[as]1ccccc1", "[asH]1cccc1", "[pH]1cccccc1",
    "[c]1ccccc1", "[cH]1ccccc1", "c1cc[nH]c1", "o1cccc1", "s1cccc1",
    "c1ccc2ccccc2c1", "c1ccc2[nH]ccc2c1", "O=C1C=CC(=O)C=C1",
    "[bH-]1ccccc1", "[nH+]1ccccc1", "[o+]1ccccc1",
    "[13cH]1ccccc1", "c1ccccc1-c1ccccc1", "C1:C:C:C:C:C:1C2:C:C:C:C:C:2",
    "CCO", "[HH]", "N->[Cu+2]<-N", "N1CCN->[Cu]<-1",
]
# Retain prior H/aromatic edge controls, including original formula residuals.
for line in (root / "crates/smiles/tests/fixtures/valence_h_contexts.tsv").read_text().splitlines():
    if not line.startswith("#"):
        sources.append(line.split("\t")[0])
sources = list(dict.fromkeys(sources))
original_count = len(sources)
for source in list(sources):
    mol = Chem.MolFromSmiles(source, PARAMS)
    assert mol is not None, source
    for order in [list(reversed(range(mol.GetNumAtoms()))),
                  list(range(1, mol.GetNumAtoms())) + [0]]:
        permuted = Chem.RenumberAtoms(mol, order)
        sources.append(Chem.MolToSmiles(permuted, canonical=False))
        sources.append(Chem.MolToSmiles(permuted, canonical=False, kekuleSmiles=True))
sources = list(dict.fromkeys(sources))
lines = ["# RDKit2026.03.6 sanitized graph; removeHs=False; original sources retained; no native filter"]
for source in sources:
    mol = Chem.MolFromSmiles(source, PARAMS)
    assert mol is not None, source
    profile = profile_reference(mol)
    atoms = list(mol.GetAtoms())
    fields = [source, rdMolDescriptors.CalcMolFormula(mol)]
    for f in [lambda a: a.GetTotalNumHs(), lambda a: int(a.GetIsAromatic()),
              lambda a: a.GetNumExplicitHs(), lambda a: a.GetNumImplicitHs(),
              lambda a: a.GetTotalValence(), lambda a: a.GetNumRadicalElectrons()]:
        fields.append(",".join(str(f(a)) for a in atoms))
    bonds = sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    int(b.GetBondType()), int(b.GetIsAromatic())) for b in mol.GetBonds())
    fields.append(";".join(",".join(map(str, b)) for b in bonds))
    lines.append("\t".join(fields))
fixture = root / "crates/smiles/tests/fixtures/kekulize_contexts.tsv"
fixture.write_text("\n".join(lines) + "\n")
with (root / "tasks/descriptors/kekulize-controls.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["name", "smiles"])
    writer.writerows((f"kekulize-{i}", s) for i, s in enumerate(sources))
meta = dict(rdkit=rdkit.__version__, original_sources=original_count, sources=sources,
            source_count=len(sources), fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
            native_acceptance_filter=False, residual_original_cohort_indices=[r["index"] for r in residuals])
invalid_sources = ["[b-]1cccc1", "c1cccc1", "cC", "C1:C:N:C:C:1"]
for source in invalid_sources:
    assert Chem.MolFromSmiles(source, PARAMS) is None, source
meta["declared_invalid_sources"] = invalid_sources
(root / "crates/smiles/tests/fixtures/kekulize_invalid.txt").write_text("\n".join(invalid_sources) + "\n")
(root / "tasks/descriptors/kekulize-controls.json").write_text(json.dumps(meta, indent=2) + "\n")
print(original_count, "original sources;", len(sources), "total with notation/atom-order controls")
