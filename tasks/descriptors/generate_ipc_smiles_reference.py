"""Exact graph information from original SMILES, independent native acceptance.

Small molecular graphs use integer determinants/rational interpolation. Long
carbon paths use analytic signed coefficients. Entropy uses85-digit Decimal.
"""
from pathlib import Path
import json

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from generate_exact_ipc_reference import determinant_interpolation, path, entropy

assert rdkit.__version__ == "2026.03.6"
root = Path(__file__).resolve().parents[2]
sources = [
    "C", "CC", "CCC", "CCO", "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "[Na+].[Cl-]", "[H]O[H]",
    "[2H]N([2H])[2H]", "[NH3]->[Cu+2]<-[NH3]", "C1:C:C:C:1",
    "c1ccc1", "C1:C:C:C:C:C=1", "[nH]1cccccc1", "[n+]1ccccc1",
    "[as]1ccccc1", "[HH]", "[H][H]", "C.C", "C.CC", "C1CCCCCCCCC1",
    "C12C3C4C1C5C2C3C45", "C12C3C1C23", "[CH3:7][C@@H](O)C(=O)O",
    "[Fe]", "CN(=O)=O", "[NH4+]", "[H-][Na+]",
]
sources += ["C" * n for n in [20, 40, 80, 128, 200, 256]]
sources = list(dict.fromkeys(sources))
parameters = Chem.SmilesParserParams()
parameters.removeHs = False
lines = ["# Original SMILES; formula; exact IPC; exact AvgIPC; graph vertices; Decimal85 entropy; no native filter"]
for source in sources:
    molecule = Chem.MolFromSmiles(source, parameters)
    assert molecule is not None, source
    n = molecule.GetNumAtoms()
    edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds()]
    coefficients = path(n) if n > 24 and source == "C" * n else determinant_interpolation(n, edges)
    ipc, average = entropy(coefficients)
    lines.append("\t".join([source, rdMolDescriptors.CalcMolFormula(molecule), ipc, average, str(n)]))
destination = root / "crates/smiles/tests/fixtures/descriptor_ipc_smiles.tsv"
destination.write_text("\n".join(lines) + "\n")
(root / "tasks/descriptors/ipc-smiles-controls.json").write_text(json.dumps(dict(rdkit=rdkit.__version__,
    sources=sources, source_count=len(sources), native_acceptance_filter=False,
    reference="integer determinant interpolation for small graphs; analytic long paths; Decimal85 entropy",
    atoms_include_explicit_hydrogen_vertices=True, unweighted_dative_connectivity=True), indent=2) + "\n")
print(len(sources), "original SMILES exact IPC/AvgIPC controls")
