"""Build trusted observations and independent oracle outputs for molecular profiles."""
import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tasks/descriptors"))
import compare_catalog as oracle
from extra50_reference import reference as extra50

CATALOG = ROOT / "tasks/descriptors/phase1-catalog.json"
CORPUS = ["CCO", "C1CCCCC1", "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O",
          "[NH4+].[Cl-]", "[13CH3]O", "CC#N", "F/C=C/F", "C1", ""]


def generate(smiles, specs, patterns):
    mol = oracle.Chem.MolFromSmiles(smiles, oracle.PARAMS) if smiles.strip() else None
    if mol is None or not mol.GetNumAtoms():
        return None
    graph_profile = oracle.profile_reference(mol)
    extra = extra50(mol)
    values = {}
    matches = {}
    for spec in specs:
        name = spec["id"]
        if name in patterns:
            matches[name] = mol.GetSubstructMatches(patterns[name], uniquify=True, maxMatches=1000)
        if name == "SPS" and not mol.GetNumHeavyAtoms():
            try:
                oracle.FUNCTIONS[name](mol)
            except ZeroDivisionError:
                values[name] = None
                continue
            raise ValueError("unexpected change in pinned SPS domain")
        value = oracle.reference_value(spec, mol, graph_profile, extra)
        values[name] = value if math.isfinite(value) else None
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    rings = list(mol.GetRingInfo().BondRings())
    return {"atoms": [a.GetSymbol() for a in atoms], "matches": matches, "oracle": values,
            "charges": [a.GetFormalCharge() for a in atoms],
            "explicitH": [1 if a.GetSymbol() == "H" else a.GetNumExplicitHs() for a in atoms],
            "implicitH": [a.GetNumImplicitHs() for a in atoms],
            "aromaticAtoms": [a.GetIsAromatic() for a in atoms],
            "ringAtoms": [a.IsInRing() for a in atoms],
            "bondOrders": [str(b.GetBondType()) for b in bonds],
            "ringBonds": [b.IsInRing() for b in bonds],
            "componentSizes": [len(c) for c in oracle.Chem.GetMolFrags(mol)],
            "ringSizes": [len(r) for r in rings]}


def observations(smiles):
    if oracle.rdkit.__version__ != "2026.03.6":
        raise RuntimeError("requires rdkit==2026.3.6")
    specs = json.loads(CATALOG.read_text())["descriptors"]
    patterns = {p["id"]: oracle.Chem.MolFromSmarts(p["smarts"])
                for p in json.loads((CATALOG.parent / "fragment-patterns.json").read_text())["descriptors"]}
    if set(patterns) != {s["id"] for s in specs if s["id"].startswith("fr_")}:
        raise ValueError("fragment pattern inventory differs from catalog")
    return [{"smiles": s, "observations": generate(s, specs, patterns)} for s in smiles]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="+")
    print(json.dumps(observations(parser.parse_args().smiles), allow_nan=False))
