"""Pinned independent Graph/Lipinski references; no native acceptance filter."""
import json
from pathlib import Path
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
assert rdkit.__version__ == '2026.03.6'
root=Path(__file__).resolve().parents[2]
ids=['AvgIpc','BalabanJ','BertzCT','Chi0','Chi0n','Chi0v','Chi1','Chi1n','Chi1v','Chi2n','Chi2v','Chi3n','Chi3v','Chi4n','Chi4v','HallKierAlpha','Ipc','Kappa1','Kappa2','Kappa3','Phi']
registered_ids=[i for i in ids if i not in ('Ipc','AvgIpc')]
rows=[dict(id=i,title=i,group='Topology',unit='',digits=4,origin='rdkit-standard',verification='pending-expanded-validation') for i in registered_ids]
smiles=list(dict.fromkeys(l.split('\t')[0] for l in (root/'crates/smiles/tests/fixtures/descriptor_surface.tsv').read_text().splitlines() if not l.startswith('#')))
smiles+=['C1CC1','C1CCC1','C1CCCC1','CC(C)C','C12CC1C2','N1CC1','C1=C=CC1','N#N','N#CC#N','NC=O','Nc1ccccc1','COc1ccccc1','OC=O','CC(=O)[O-]','[CH3]','[CH2]','[C]','[C-]#[O+]','[N]','[NH]','[O]','[S]','[P]','[Na]','[U]','[Cl]','[n+]1ccccc1','[NH2]c1ccccc1','O=[N+]([O-])O','CP1(C)=CC=CN=C1C','P=CC','S=CC','S=C=S','C[Se]C','[H]O[H]','[2H]O','[H]C([H])([H])[H]']
smiles+=['[BH4-]','[BH4-].[Na+]','[BH4-].[NH4+]','[LiH2-]','[NaH2-]','[KH2-]','[FrH2-]','[AlH4-]']
smiles+=['C'*n for n in [20,40,80,128,200,256]]
smiles=list(dict.fromkeys(smiles))
p=Chem.SmilesParserParams();p.removeHs=False
fn=dict(Descriptors._descList)
ipc=['# RDKit 2026.03.6; experimental IPC numerical stability controls']
contexts=['# RDKit 2026.03.6; SMILES, atoms (symbol,H,charge,aromatic,bracket), bonds (a,b,order)']
scalar=['# RDKit 2026.03.6; removeHs=False; SMILES, descriptor, expected']
atomic=['# RDKit 2026.03.6; SMILES, atom, hybridization, total valence, radical electrons, HallKier contribution']
paths=['# RDKit 2026.03.6; SMILES, number of bonds, unique paths (useHs=False)']
for s in smiles:
 m=Chem.MolFromSmiles(s,p);assert m is not None,s
 for i in ids:
  target=ipc if i in ('Ipc','AvgIpc') else scalar
  target.append(f'{s}\t{i}\t{fn[i](m):.17g}')
 atoms=';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())}' for a in m.GetAtoms())
 bonds=';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{b.GetBondTypeAsDouble()}' for b in m.GetBonds())
 contexts.append(f'{s}\t{atoms}\t{bonds}')
 contrib=[0.]*m.GetNumAtoms();rdMolDescriptors.CalcHallKierAlpha(m,contrib)
 for i,a in enumerate(m.GetAtoms()): atomic.append(f'{s}\t{i}\t{a.GetHybridization().name}\t{a.GetTotalValence()}\t{a.GetNumRadicalElectrons()}\t{contrib[i]:.17g}')
 for n in range(1,5): paths.append(f'{s}\t{n}\t{len(Chem.FindAllPathsOfLengthN(m,n))}')
(root/'crates/smiles/tests/fixtures/graph_contexts.tsv').write_text('\n'.join(contexts)+'\n')
(root/'tasks/descriptors/ipc-numerical-controls.tsv').write_text('\n'.join(ipc)+'\n')
(root/'crates/smiles/tests/fixtures/descriptor_graph.tsv').write_text('\n'.join(scalar)+'\n')
(root/'crates/smiles/tests/fixtures/graph_atomic.tsv').write_text('\n'.join(atomic)+'\n')
(root/'crates/smiles/tests/fixtures/graph_paths.tsv').write_text('\n'.join(paths)+'\n')
(root/'tasks/descriptors/graph-definitions.json').write_text(json.dumps(dict(rdkit=rdkit.__version__,descriptors=rows,molecules=len(smiles),scalar_values=len(scalar)-1,atomic_values=len(atomic)-1),indent=2)+'\n')
pt=Chem.GetPeriodicTable()
# Minimal valence and radical tables are chemistry facts, not per-molecule answers.
vals=[list(pt.GetValenceList(z)) for z in range(1,119)]
s='//! RDKit 2026.03.6 PeriodicTable valence facts. See RDKIT-DATA-LICENSE.txt.\n'
s+='pub(crate) const VALENCES: [&[i32]; 118] = [\n'+''.join('    &'+str(v)+',\n' for v in vals)+'];\n'
(root/'crates/smiles/src/graph_data.rs').write_text(s)
print(len(smiles),'graphs;',len(scalar)-1,'scalar references;',len(atomic)-1,'atomic references')
