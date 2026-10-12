#!/usr/bin/env python3
"""Deterministically generate the checked-in RDKit parity fixtures."""
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, MACCSkeys, QED, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
DATA = Path(__file__).resolve().parents[1] / "crates" / "smiles" / "tests" / "data"

# Hand-curated coverage corpus.  Keep explicit strings (rather than fetching a
# database) so regeneration is offline, reviewable, and stable across machines.
MOLECULES = [
    # aliphatics, functional groups, rings, bridged/spiro/macrocycles
    "C", "CC", "CCC", "CCCC", "CC(C)C", "CC(C)(C)C", "CCCCCC", "CCCCCCCCCCCCCCCC(=O)O",
    "CO", "CCO", "CCCO", "CC(C)O", "OCCO", "OCCCO", "CCN", "CCCN", "CCN(CC)CC",
    "CC=O", "CC(=O)C", "CC(=O)O", "CCC(=O)O", "CC(=O)OC", "CCOC(=O)C", "CC(=O)N",
    "CC(=O)NC", "CC#N", "C=CC=C", "C#CC", "C1CC1", "C1CCC1", "C1CCCC1", "C1CCCCC1",
    "C1CCCCCC1", "C1CCCCCCC1", "C1CCCCCCCC1", "C1CCCCCCCCCC1", "C1CCCCCCCCCCC1",
    "C1CCNCC1", "C1COCCN1", "C1COCCO1", "C1CC2CCC1CC2", "C1C2CC3CC1CC(C2)C3",
    "C1CCC2(CC1)CCCCC2", "C1CC2CCC3CC(C1)C23", "C1CC2(CC1)CCC2", "C1CCC2(CC1)CCCC2",
    "O=C1CCCCC1", "C=C1CCCCC1", "O=C1NC(=O)NC(=O)C1",
    # aromatic and heteroaromatic systems
    "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "Clc1ccccc1", "Brc1ccccc1",
    "FC(F)(F)c1ccccc1", "c1ccncc1", "c1ncccc1", "c1ncncc1", "c1ccoc1", "c1ccsc1",
    "c1cc[nH]c1", "c1ncc[nH]1", "c1n[nH]cc1", "c1ocnc1", "c1scnc1", "c1nn[nH]n1",
    "c1ccc(cc1)-c1ccccc1", "c1ccc2ccccc2c1", "C1=CC2=CC=CC=C2C=C1",
    "c1ccc2[nH]ccc2c1", "c1ccc2ncccc2c1", "c1ccc2nccnc2c1", "c1ccc2[nH]cnc2c1",
    "c1ncnc2[nH]cnc12", "c1ccc2[nH]cnc2c1", "c1ccc2sccc2c1", "c1ccc2occc2c1",
    "c1ccc2c(c1)[nH]c1ccccc12", "c1ccc2cc3ccccc3cc2c1", "c1ccc2c(c1)oc1ccccc12",
    "c1ccc2c(c1)sc1ccccc12", "c1ccc2c(c1)ncc1ccccc12",
    # quinones and exocyclic multiple bonds
    "O=C1C=CC(=O)C=C1", "CC1=CC(=O)C=CC1=O", "O=C1C2=CC=CC=C2C(=O)c2ccccc12",
    "N=C1CCCCC1", "S=C1NCCN1", "O=C1C=CC2=CC=CC=C12", "O=S1(=O)CCCCC1",
    # charges, zwitterions, salts
    "[Na+].[Cl-]", "CC(=O)[O-].[Na+]", "C[N+](C)(C)C", "CC(=O)[O-]", "C[NH3+]",
    "[NH3+]CC(=O)[O-]", "C[NH2+]CC(=O)[O-]", "[O-]S(=O)(=O)[O-].[Na+].[Na+]",
    "Cl.CN1CCC[C@H]1c1cccnc1", "CC(=O)O.[K+]", "[Li+].[O-]C(=O)C",
    # S, P, halogens, nitro, sulfonamide
    "CS(=O)(=O)O", "CSC", "CS(C)=O", "CS(=O)C", "NS(=O)(=O)C", "O=S(=O)(N)c1ccccc1",
    "COP(=O)(O)OC", "OP(=O)(O)O", "CP(=O)(O)O", "P(Cl)(Cl)(Cl)(Cl)Cl", "CCF", "CCCl",
    "CCBr", "CCI", "ClC(Cl)(Cl)Cl", "FC(F)(F)F", "O=[N+]([O-])c1ccccc1", "C[N+](=O)[O-]",
    # carbohydrates and peptides
    "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "OC[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O",
    "O=C[C@H](O)[C@H](O)[C@H](O)CO", "OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O",
    "NCC(=O)O", "C[C@H](N)C(=O)O", "N[C@@H](Cc1ccccc1)C(=O)O",
    "NCC(=O)NCC(=O)O", "C[C@H](N)C(=O)NCC(=O)O", "N[C@@H](CC(=O)O)C(=O)NCC(=O)O",
    "CC(C)[C@H](N)C(=O)N[C@@H](C)C(=O)O",
    # stereochemistry
    "C[C@@H](N)C(=O)O", "N[C@H](Cc1ccccc1)C(=O)O", "F/C=C/F", "F/C=C\\F",
    "C/C=C/C", "C/C=C\\C", "N[C@@H](CO)C(=O)O", "N[C@H](CO)C(=O)O",
    "C[C@H](O)[C@@H](O)C", "C[C@@H](O)[C@@H](O)C",
    # named drugs (>=30)
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",  # caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
    "CC(=O)Nc1ccc(O)cc1",  # acetaminophen
    "CN1CCC[C@H]1c1cccnc1",  # nicotine
    "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",  # salbutamol
    "CCN(CC)CCNC(=O)c1ccc(N)cc1",  # procainamide
    "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",  # atenolol
    "CC(=O)N1CCN(CC1)c1ccccc1",  # phenacetyl-like
    "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",  # diazepam
    "COc1ccc2cc(ccc2c1)C(C)C(=O)O",  # naproxen core
    "OC(=O)c1ccccc1O", "NC(=O)c1ccccc1", "CCOC(=O)c1ccccc1",
    "CN(C)CCOC(c1ccccc1)c1ccccc1",  # diphenhydramine
    "CNCCC(c1ccccc1)c1ccccc1",  # methadone fragment
    "COc1ccc(CCN(C)C)cc1OC",  # verapamil fragment
    "CN1CCN(CC1)c1ncccn1",  # buspirone fragment
    "CCOC(=O)N1CCC(CC1)N2CCC(CC2)OC", "COc1ccc(CN(C)C)cc1",
    "CC(C)c1ccc(cc1)[C@@H](C)C(=O)O", "COc1ccc(cc1OC)C(=O)NCCCN(C)C",
    "Clc1ccc(cc1)C(c1ccccc1)C(Cl)Cl",  # dicofol-like
    "CC1=C(C(=O)N(N1)C)c1ccccc1",  # antipyrine
    "CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1",  # atropine-like
    "CC(C)(C)c1ccc(O)cc1", "CC(C)(C)c1ccc(cc1)C(O)=O",
    "COc1ccccc1O", "CCOc1ccc(NC(C)=O)cc1", "CCN1C(=O)NC(C)c2ccccc21",
    "NCCc1ccc(O)c(O)c1",  # dopamine
    "CNCCc1ccc(O)c(O)c1",  # epinephrine fragment
    "COc1ccc(CCN)cc1OC", "CNC1(C)C2CCC1(C)C(O)C2", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
    "CC12CCC3C(C1CCC2=O)CCC4=CC(=O)CCC34C", "COc1ccc2[nH]cc(CCN(C)C)c2c1",
    "CN1CCC23C4=C5OCOC5=CC=C4C1CC2C=C3", "CCN(CC)C(=O)c1ccc(N)cc1",
    "CCOC(=O)c1ccc(N)cc1", "COc1ccc(C(=O)O)cc1", "CC(C)NCC(O)c1ccc(O)cc1",
    # more structural diversity
    "c1ccc(cc1)S(=O)(=O)Nc1ccccc1", "CCOC(=O)C1=CC=CC=C1", "CC(C)(C)OC(=O)NCC",
    "O=C(O)C(O)(P(=O)(O)O)P(=O)(O)O", "N=C(N)N", "NC(N)=O", "O=C(N)N1CCCC1",
    "C1=NC2=C(N1)N=CN2", "c1cc2ccc3cccc4ccc(c1)c2c34", "C1CCC2CCCCC2C1",
    "c1ccc(cc1)CCc1ccccc1", "OCC1OC(O)C(O)C(O)C1O", "C1CCCCC1.C1=CC=CC=C1",
    "CCCCC(=O)O", "CCCCOC", "CC(C)CC", "CC(C)C(C)C", "CCOC", "CCS", "CCP",
    "NCCN", "NCCO", "OCCS", "C1CNCCN1", "C1CSCCN1", "C1CCOC1", "C1CCSC1",
    "c1nccs1", "c1ncco1", "c1cn[nH]c1", "c1nnc[nH]1", "c1ccc2[nH]ncc2c1",
    "O=C(O)c1ncccc1", "NC(=O)c1ncccc1", "COC(=O)c1ccncc1", "CC(=O)c1ccccc1",
    "O=C(Nc1ccccc1)c1ccccc1", "CCS(=O)(=O)N", "N#CC1CCCCC1", "C1=CCCCC1",
]

HEADER = ["smiles", "num_atoms", "num_bonds", "formula", "mol_weight", "exact_mass", "tpsa",
          "logp", "mol_mr", "qed", "hbd", "hba", "rotb", "heteroatoms", "fcsp3", "ring_count",
          "arom_rings", "aliph_rings", "satur_rings", "arom_carbo", "arom_hetero", "satur_carbo",
          "satur_hetero", "aliph_carbo", "aliph_hetero", "lipinski_viol", "structural_alerts",
          "canonical", "murcko"]

# The public duckSMILES alert API documents this curated SMARTS contract. RDKit
# is used as the independent substructure matcher for the reference counts.
ALERT_SMARTS = [
    "[CX3](=[OX1])[F,Cl,Br,I]", "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]", "[CX4][Cl,Br,I]",
    "[CX3H1](=O)[#6]", "[CX3]=[CX3][CX3]=[OX1]", "C#C[CX3]=[OX1]", "C1OC1", "C1CN1",
    "C1CC(=O)N1", "[NX2]=C=[OX1]", "[NX2]=C=[SX1]", "[SX2]C#N", "[NX2]=[NX2]=[NX1]",
    "[NX2]=[NX2]", "[NX2]=O", "[NX3](=O)=O", "[a][NX3](=O)=O", "[NX3][NX3]",
    "[NX3][OX2H]", "[OX2][OX2]", "[SX2][SX2]", "[SX2H]", "[CX3H1]=O",
    "[CX3;!R]=[NX2;!R]", "[CX3](=[OX1])O[CX3](=[OX1])", "[#6]S(=O)(=O)O[#6]",
    "P(=O)([OX2][#6])[OX2][#6]", "[NX4+]", "[a][NX2]=O", "[#6](=O)[#6](=O)",
    "[CX2]#[CX2]", "[Cl,Br,I][c]",
]

def bits(fp):
    return ",".join(map(str, fp.GetOnBits()))

def lipinski(m):
    return sum((Descriptors.MolWt(m) > 500, Crippen.MolLogP(m) > 5,
                rdMolDescriptors.CalcNumHBD(m) > 5, rdMolDescriptors.CalcNumHBA(m) > 10))

def row(s, m, alerts):
    vals = [s, m.GetNumAtoms(), m.GetNumBonds(), rdMolDescriptors.CalcMolFormula(m),
            f"{Descriptors.MolWt(m):.4f}", f"{Descriptors.ExactMolWt(m):.4f}",
            f"{rdMolDescriptors.CalcTPSA(m):.4f}", f"{Crippen.MolLogP(m):.4f}",
            f"{Crippen.MolMR(m):.4f}", f"{QED.qed(m):.4f}", rdMolDescriptors.CalcNumHBD(m),
            rdMolDescriptors.CalcNumHBA(m), rdMolDescriptors.CalcNumRotatableBonds(m),
            rdMolDescriptors.CalcNumHeteroatoms(m), f"{rdMolDescriptors.CalcFractionCSP3(m):.4f}",
            rdMolDescriptors.CalcNumRings(m), rdMolDescriptors.CalcNumAromaticRings(m),
            rdMolDescriptors.CalcNumAliphaticRings(m), rdMolDescriptors.CalcNumSaturatedRings(m),
            rdMolDescriptors.CalcNumAromaticCarbocycles(m), rdMolDescriptors.CalcNumAromaticHeterocycles(m),
            rdMolDescriptors.CalcNumSaturatedCarbocycles(m), rdMolDescriptors.CalcNumSaturatedHeterocycles(m),
            rdMolDescriptors.CalcNumAliphaticCarbocycles(m), rdMolDescriptors.CalcNumAliphaticHeterocycles(m),
            lipinski(m), sum(m.HasSubstructMatch(q) for q in alerts), Chem.MolToSmiles(m),
            Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))]
    return list(map(str, vals))

def variants(m, seed):
    out = [Chem.MolToSmiles(m)]
    k = Chem.Mol(m)
    try:
        Chem.Kekulize(k, clearAromaticFlags=True)
        out.append(Chem.MolToSmiles(k, kekuleSmiles=True))
    except Exception:
        pass
    out.extend(Chem.MolToRandomSmilesVect(m, 6, randomSeed=0xD5A000 + seed))
    return list(dict.fromkeys(out))

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    alerts = [Chem.MolFromSmarts(s) for s in ALERT_SMARTS]
    if any(q is None for q in alerts):
        raise AssertionError("an alert SMARTS did not compile")
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    rows, var_rows, mac_rows, mor_rows = [], [], [], []
    seen = set()
    for source in MOLECULES:
        m = Chem.MolFromSmiles(source)
        if m is None:
            raise ValueError(f"unparseable corpus entry: {source}")
        canonical = Chem.MolToSmiles(m)
        if canonical in seen:
            continue
        seen.add(canonical)
        idx = len(rows)
        rows.append(row(source, m, alerts))
        var_rows.extend((str(idx), v) for v in variants(m, idx))
        mac_rows.append((source, bits(MACCSkeys.GenMACCSKeys(m))))
        mor_rows.append((source, bits(morgan.GetFingerprint(m))))
    if len(rows) < 200:
        raise AssertionError(f"corpus has only {len(rows)} unique molecules")
    outputs = {
        "rdkit_descriptors.tsv": (HEADER, rows),
        "smiles_variants.tsv": (["id", "smiles"], var_rows),
        "rdkit_maccs.tsv": (["smiles", "on_bits"], mac_rows),
        "rdkit_morgan.tsv": (["smiles", "on_bits"], mor_rows),
    }
    for name, (header, records) in outputs.items():
        with (DATA / name).open("w", encoding="utf-8", newline="\n") as f:
            f.write("\t".join(header) + "\n")
            for record in records:
                f.write("\t".join(record) + "\n")
    print(f"molecules: {len(rows)}, variants: {len(var_rows)}, rdkit {Chem.rdBase.rdkitVersion}")

if __name__ == "__main__":
    main()
