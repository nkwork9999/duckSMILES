"""Build trusted observations and independent oracle outputs for molecular profiles."""
import argparse
import json
import math
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tasks/descriptors"))
import compare_catalog as oracle
from extra50_reference import reference as extra50
from rdkit.Chem import QED, Crippen, rdMolDescriptors, rdPartialCharges
from rdkit.Chem.EState import EStateIndices

CATALOG = ROOT / "tasks/descriptors/phase1-catalog.json"
CORPUS = ["CCO", "C1CCCCC1", "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O",
          "[NH4+].[Cl-]", "[13CH3]O", "CC#N", "F/C=C/F", "C1", ""]


def parse(smiles):
    mol = oracle.Chem.MolFromSmiles(smiles, oracle.PARAMS) if smiles.strip() else None
    if mol is None or not mol.GetNumAtoms():
        return None
    return mol


def oracle_values(smiles, specs):
    mol = parse(smiles)
    if mol is None:
        return None
    graph_profile = oracle.profile_reference(mol)
    extra = extra50(mol)
    values = {}
    for spec in specs:
        name = spec["id"]
        if name == "SPS" and not mol.GetNumHeavyAtoms():
            try:
                oracle.FUNCTIONS[name](mol)
            except ZeroDivisionError:
                values[name] = None
                continue
            raise ValueError("unexpected change in pinned SPS domain")
        value = oracle.reference_value(spec, mol, graph_profile, extra)
        values[name] = value if math.isfinite(value) else None
    return [values[s["id"]] for s in specs]


def bits(value):
    return int.from_bytes(struct.pack(">d", value), "big")


def burden_observations(mol):
    heavy = oracle.Chem.RemoveAllHs(mol)
    try:
        rdPartialCharges.ComputeGasteigerCharges(heavy, throwOnParamFailure=True)
        parameterized = True
    except ValueError as error:
        if "parameters" not in str(error).lower():
            raise
        parameterized = False
    crippen = Crippen._GetAtomContribs(heavy, force=True)
    diagonal = [[bits(a.GetMass()), bits(a.GetDoubleProp("_GasteigerCharge") if parameterized else math.nan),
                 bits(lp), bits(mr)] for a, (lp, mr) in zip(heavy.GetAtoms(), crippen, strict=True)]
    return {"diagonal": diagonal, "parameterized": parameterized,
            "edges": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondType())] for b in heavy.GetBonds()]}


def atomic_observations(mol):
    area, h_area = rdMolDescriptors._CalcLabuteASAContribs(mol)
    crippen = Crippen._GetAtomContribs(mol, force=True)
    estate = EStateIndices(mol)
    rdPartialCharges.ComputeGasteigerCharges(mol)
    rows = [{"area": bits(a), "mr": bits(mr), "logp": bits(lp), "estate": bits(e),
             "partialCharge": bits(atom.GetDoubleProp("_GasteigerCharge"))}
            for atom, a, (lp, mr), e in zip(mol.GetAtoms(), area, crippen, estate, strict=True)]
    table = oracle.Chem.GetPeriodicTable()
    tpsa = rdMolDescriptors._CalcTPSAContribs(mol)
    alpha = [0.0] * mol.GetNumAtoms()
    rdMolDescriptors.CalcHallKierAlpha(mol, alpha)
    stereo = oracle.Chem.Mol(mol)
    oracle.Chem.FindPotentialStereoBonds(stereo)
    centers = {i for i, _ in oracle.Chem.FindMolChiralCenters(
        stereo, includeUnassigned=True, includeCIP=False, useLegacyImplementation=False)}
    for bond in stereo.GetBonds():
        if bond.GetBondType() == oracle.Chem.BondType.DOUBLE and bond.GetStereo() != oracle.Chem.BondStereo.STEREONONE:
            centers.update([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    props = []
    for i, atom in enumerate(mol.GetAtoms()):
        z, iso = atom.GetAtomicNum(), atom.GetIsotope()
        exact = table.GetMassForIsotope(z, iso) if iso else table.GetMostCommonIsotopeMass(z)
        props.append(dict(mass=bits(atom.GetMass()), exactMass=bits(exact),
                          attachedH=atom.GetTotalNumHs(), outer=table.GetNOuterElecs(z),
                          radical=atom.GetNumRadicalElectrons(), tpsa=bits(tpsa[i]), alpha=bits(alpha[i]),
                          hybrid=atom.GetHybridization().name, stereoPossible=bool(atom.HasProp("_ChiralityPossible")),
                          cipAssigned=bool(atom.HasProp("_CIPCode")), modernStereo=i in centers))
    expanded = Crippen._GetAtomContribs(oracle.Chem.AddHs(mol), force=True)
    return {"qedProperties": list(map(bits, QED.properties(mol))), "burden": burden_observations(mol), "morganInitial": list(rdMolDescriptors.GetConnectivityInvariants(mol)),
            "atomicContributions": rows, "hydrogenArea": bits(h_area), "atomProperties": props,
            "expandedLogP": [bits(lp) for lp, _ in expanded], "expandedMR": [bits(mr) for _, mr in expanded]}


def graph_observations(mol):
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    rings = list(mol.GetRingInfo().BondRings())
    return {"atoms": [a.GetSymbol() for a in atoms],
            "vertices": [
                {"z": a.GetAtomicNum(), "aromatic": a.GetIsAromatic(),
                 "hydrogens": a.GetTotalNumHs(includeNeighbors=True),
                 "charge": a.GetFormalCharge(), "degree": a.GetDegree(),
                 "totalDegree": a.GetTotalDegree(), "rings": mol.GetRingInfo().NumAtomRings(a.GetIdx()),
                 "valence": a.GetTotalValence()} for a in atoms],
            "edges": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondType())] for b in bonds],
            "charges": [a.GetFormalCharge() for a in atoms],
            "explicitH": [1 if a.GetSymbol() == "H" else a.GetNumExplicitHs() for a in atoms],
            "implicitH": [a.GetNumImplicitHs() for a in atoms],
            "aromaticAtoms": [a.GetIsAromatic() for a in atoms],
            "ringAtoms": [a.IsInRing() for a in atoms],
            "bondOrders": [str(b.GetBondType()) for b in bonds],
            "ringBonds": [b.IsInRing() for b in bonds],
            "componentSizes": [len(c) for c in oracle.Chem.GetMolFrags(mol)],
            "ringBondIndices": rings,
            "perceivedRings": [
                {"orders": [str(bonds[i].GetBondType()) for i in ring],
                 "elements": [a.GetSymbol() for i in ring for a in
                              (bonds[i].GetBeginAtom(), bonds[i].GetEndAtom())]}
                for ring in rings],
            "ringSizes": [len(r) for r in rings]}


def observations(smiles):
    if oracle.rdkit.__version__ != "2026.03.6":
        raise RuntimeError("requires rdkit==2026.3.6")
    result = []
    for source in smiles:
        mol = parse(source)
        data = graph_observations(mol) | atomic_observations(mol) if mol is not None else None
        result.append({"smiles": source, "observations": data})
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="+")
    print(json.dumps(observations(parser.parse_args().smiles), allow_nan=False))
