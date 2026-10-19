"""Pinned RDKit oracle; retain every declared SMILES without native filtering."""
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from compare_catalog import PARAMS, profile_reference

assert rdkit.__version__ == "2026.03.6"
sources = [
    "[H]", "[HH]", "[H][H]", "[2HH]", "[2H][H]", "[HH+].[Cl-]",
    "C[H]", "[NH4+]", "[nH]1cccc1",
    r"N[C@@H](CO)C1=N/C(=C\C2:C:C:C(O):C:C:2)C(=O)N1CC=O.[HH]",
]
lines = ["# RDKit 2026.03.6; original inputs; removeHs=False; no native filter"]
for source in sources:
    mol = Chem.MolFromSmiles(source, PARAMS)
    assert mol is not None, source
    profile = profile_reference(mol)
    values = [str(int(profile[key])) for key in
              ("ds_ExplicitHydrogens", "ds_ImplicitHydrogens", "ds_TotalHydrogens")]
    lines.append("\t".join([source, rdMolDescriptors.CalcMolFormula(mol), *values]))
destination = Path(__file__).resolve().parents[2] / "crates/smiles/tests/fixtures/descriptor_hydrogen_profile.tsv"
destination.write_text("\n".join(lines) + "\n")
print(f"{len(sources)} sources → {destination}")
