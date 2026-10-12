#!/usr/bin/env python3
"""RDKit 2025.09.6 oracle for the 21 structural descriptors. No runtime dependency."""
import csv
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors
from gen_rdkit_fixtures import MOLECULES, variants

ROOT = Path(__file__).resolve().parents[1]
APIS = ["mol_num_fragments","mol_formal_charge","mol_num_explicit_h","mol_num_implicit_h","mol_num_total_h","mol_num_single_bonds","mol_num_double_bonds","mol_num_triple_bonds","mol_num_aromatic_bonds","mol_num_ring_atoms","mol_num_ring_bonds","mol_largest_ring_size","mol_num_aromatic_atoms","mol_num_carbons","mol_num_nitrogens","mol_num_oxygens","mol_num_halogens","mol_heteroatom_fraction","mol_aromatic_fraction","mol_heavy_atom_mass","mol_mean_degree"]
EXTRA = ["[H]", "[H][H]", "[H]O[H]", "[2H]O[2H]", "[NH4+]", "[Cl-]",
         "[13CH3][18OH]", "[15NH2]C", "[235U]", "[Xe]", "[SiH4]", "[SeH2]",
         "C12C3C4C1C5C2C3C45", "C1=CC=CC=CC=CC=C1",
         "C1=CC2=C(C=C1)C1=CC=CC=C21", "O=C1C=CC(=O)C2=C1OC=CO2",
         "[CH3]", "[CH2]", "[C]", "[O]", "[N]", "C1=C[N]C=C1"]

def values(s):
    params = Chem.SmilesParserParams()
    params.removeHs = False  # Preserve input H atoms and bond counts.
    m = Chem.MolFromSmiles(s, params)
    if m is None:
        raise ValueError(s)
    atoms, bonds = list(m.GetAtoms()), list(m.GetBonds())
    heavy = [a for a in atoms if a.GetAtomicNum() != 1]
    explicit = sum(1 if a.GetAtomicNum() == 1 else a.GetNumExplicitHs() for a in atoms)
    implicit = sum(a.GetNumImplicitHs() for a in heavy)
    n = len(heavy)
    aromatic = sum(a.GetIsAromatic() for a in atoms)
    counts = [sum(b.GetBondType() == t for b in bonds) for t in
              (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)]
    return [len(Chem.GetMolFrags(m)), Chem.GetFormalCharge(m), explicit, implicit, explicit+implicit,
            *counts, sum(a.IsInRing() for a in atoms), sum(b.IsInRing() for b in bonds),
            max(map(len,m.GetRingInfo().AtomRings()),default=0), aromatic,
            *[sum(a.GetAtomicNum()==z for a in atoms) for z in (6,7,8)],
            sum(a.GetAtomicNum() in (9,17,35,53) for a in atoms),
            sum(a.GetAtomicNum() not in (1,6) for a in atoms)/n if n else 0,
            aromatic/n if n else 0, Descriptors.HeavyAtomMolWt(m),
            2*sum(b.GetBeginAtom().GetAtomicNum()!=1 and b.GetEndAtom().GetAtomicNum()!=1 for b in bonds)/n if n else 0]

def main():
    assert rdBase.rdkitVersion == "2025.09.6", rdBase.rdkitVersion
    corpus = dict.fromkeys(MOLECULES + EXTRA)
    for i,s in enumerate(MOLECULES + EXTRA):
        for v in variants(Chem.MolFromSmiles(s),i):
            corpus[v] = None
    # Exercise every supported element and every isotope in the pinned oracle,
    # not only the few labels that originally exposed the mass bug.
    table = Chem.GetPeriodicTable()
    for z in range(1,119):
        symbol = table.GetElementSymbol(z)
        corpus[f'[{symbol}]'] = None
        for isotope in range(1,351):
            if table.GetMassForIsotope(z,isotope):
                corpus[f'[{isotope}{symbol}]'] = None
    out = ROOT/"crates/smiles/tests/data/rdkit_profile.tsv"
    with out.open("w",newline="") as f:
        writer = csv.writer(f,delimiter="\t",lineterminator="\n")
        writer.writerow(["smiles",*APIS])
        for s in corpus:
            writer.writerow([s,*values(s)])
    print(f"{len(corpus)} inputs x {len(APIS)} descriptors; RDKit {rdBase.rdkitVersion}")
if __name__ == "__main__":
    main()
