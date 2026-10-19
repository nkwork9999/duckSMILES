"""Pinned RDKit regression graphs for cyclic donation and porphyrin fusion.

Definitions: RDKit Release_2026_03_6 Code/GraphMol/Aromaticity.cpp.
Source cases: archived NCI comparison plus quinone/sulfoxide/saturated controls.
The fixture stores literal flags; the native implementation is never consulted.
"""
from pathlib import Path
import json,gzip
import rdkit
from rdkit import Chem
assert rdkit.__version__ == "2026.03.6"
root=Path(__file__).resolve().parents[2]
selected={380,383,1218,2444,2455,2614,3043,3044,3045,3078,3079,3577,4321}
smiles=[]
with gzip.open(root/'tasks/evidence/catalog139-nci/differences.jsonl.gz','rt') as f:
 for line in f:
  d=json.loads(line)
  if d['index'] in selected: smiles.append(d['smiles'])
smiles += ['O=C1C=CC(=O)C=C1','O=S1(=O)C2=CC=CC=C2OC3=C1C=CC=C3','C1=CC2=C(C=C1)CCC2','C1=CC2=CC=CC=C2C=C1','c1ccc2ccccc2c1']
# All nine new ChEBI regressions from the completed v10 100k comparison.
smiles += ['COCCO[C@@H]1C2=C(N[C@H](C3=CC=CC=C3)[C@H]1O)C1=NC(C)=C(C)N1C=C2', 'C=CC1=C(C)/C2=C/C3=N/C(=C4\\C5=C(C(C)=C(/C=C6\\N=C(/C=C/1N2)C(C)=C6CC)N5)C(O)=C4C(=O)OC)C(/C=C/C(=O)O)=C3C', 'N#CC1=C2CCCCC2=C(NCCN2CCOCC2)N2C3=CC=CC=C3N=C12', 'C1=CC=C(CNC2=C3CCCCC3=NC3=NC=NN32)N=C1', 'CSC1=NN2C=C3C(=O)CC(C4=CC=CC=C4)CC3=NC2=N1', 'COC1=NC2=NC(C3=NOC(C)=N3)=CN2C2=C1CCCC2', 'CCCCN1CC2=C(N=C3C=C(C4=CC=CO4)NN3C2=O)C1=O', 'CC1=NC2=NC=NN2C2=C1CCN2C(C)(C)C', 'CCN1CCN(C2=NC3=C(C4=NC=NN42)C2=C(CCC(C)C2)S3)CC1']
params=Chem.SmilesParserParams();params.removeHs=False
out=['# RDKit 2026.03.6 sanitized, removeHs=False; original atom and bond order; aromatic flags']
for s in smiles:
 m=Chem.MolFromSmiles(s,params);assert m
 atoms=''.join('1' if a.GetIsAromatic() else '0' for a in m.GetAtoms())
 bonds=','.join(f'{b.GetBeginAtomIdx()}-{b.GetEndAtomIdx()}' for b in m.GetBonds() if b.GetIsAromatic())
 # Normalize edge orientation and ordering; native ring-closure bond insertion order differs.
 edges=sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx())) for b in m.GetBonds() if b.GetIsAromatic())
 bonds=','.join(f'{a}-{b}' for a,b in edges)
 out.append(f'{s}\t{atoms}\t{bonds}')
(root/'crates/smiles/tests/fixtures/aromatic_fused.tsv').write_text('\n'.join(out)+'\n')
print(len(smiles),'graphs')
