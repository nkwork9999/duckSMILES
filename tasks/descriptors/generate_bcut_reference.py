"""Independent registered BCUT and symmetric eigensolver controls; pinned RDKit."""
from pathlib import Path
import json,math
import numpy as np
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')
root=Path(__file__).resolve().parents[2]
smiles=list(dict.fromkeys(l.split('\t')[0] for l in (root/'crates/smiles/tests/fixtures/charge_contexts.tsv').read_text().splitlines() if not l.startswith('#')))
smiles+=['[H][NH]','C[H-]C','[H-][Na+]','[H][13C]([H])([H])[H]','[18OH2]','[3H]O','[H]C([H])([H])[H]','[H][H].[H]','[H-]=C','[H-]=[CH2]']
smiles=list(dict.fromkeys(smiles));p=Chem.SmilesParserParams();p.removeHs=False
ids=[f'BCUT2D_{s}' for s in ('MWHI','MWLOW','CHGHI','CHGLO','LOGPHI','LOGPLOW','MRHI','MRLOW')]
fn=dict(Descriptors._descList);scalar=['# registered RDKit2026.03.6 wrappers; exceptions become nan; original SMILES/removeHs=False']
contexts=['# BCUT-specific RemoveAllHs sanitized graph; source SMILES, atoms(symbol,H,charge,aromatic,bracket,isotope), bonds(a,b,order)']
original=['# Original sanitized RDKit graph; same schema as bcut_contexts.tsv']
for s in smiles:
 m=Chem.MolFromSmiles(s,p);assert m is not None,s
 for id in ids:scalar.append(f'{s}\t{id}\t{fn[id](m):.17g}')
 atoms0=';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())},{a.GetIsotope()}' for a in m.GetAtoms())
 bonds0=';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{b.GetBondTypeAsDouble()}' for b in m.GetBonds())
 original.append(f'{s}\t{atoms0}\t{bonds0}')
 h=Chem.RemoveAllHs(m)
 atoms=';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())},{a.GetIsotope()}' for a in h.GetAtoms())
 bonds=';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{b.GetBondTypeAsDouble()}' for b in h.GetBonds())
 contexts.append(f'{s}\t{atoms}\t{bonds}')
(root/'crates/smiles/tests/fixtures/descriptor_bcut.tsv').write_text('\n'.join(scalar)+'\n')
(root/'crates/smiles/tests/fixtures/bcut_original_contexts.tsv').write_text('\n'.join(original)+'\n')
(root/'crates/smiles/tests/fixtures/bcut_contexts.tsv').write_text('\n'.join(contexts)+'\n')
rows=[dict(id=id,title=id,group='BCUT',unit='',digits=4,origin='rdkit-standard',verification='pending-expanded-validation') for id in ids]
(root/'tasks/descriptors/bcut-definitions.json').write_text(json.dumps(dict(rdkit=rdkit.__version__,descriptors=rows,molecules=len(smiles),scalar_values=len(scalar)-1,hydrogen_policy='RemoveAllHs including isotope/isolated H; recompute charges/Crippen on that graph'),indent=2)+'\n')
# Matrices are defined independently by formula, not the native solver. Native
# tests reconstruct only these input matrices and compare to NumPy eigvalsh.
def matrix(family,n,seed):
 a=np.zeros((n,n))
 for i in range(n):
  for j in range(n):
   if family=='diagonal':v=float((i*13+seed)%17-8) if i==j else 0.
   elif family=='rankone':v=5. if i==j else .001
   elif family=='path':v=12.011 if i==j else 1. if abs(i-j)==1 else 0.
   elif family=='star':v=0. if i==j else 1. if min(i,j)==0 else 0.
   elif family=='clustered':v=1.+(i%3)*1e-10 if i==j else 1e-12
   else:
    lo,hi=min(i,j),max(i,j)
    v=float((i*13+seed)%50+1) if i==j else (((lo*131+hi*17+seed*31)%101)-50)/250.
   a[i,j]=v
 return a
lines=['# family, n, seed, smallest, largest; NumPy2.5.4 eigvalsh independent oracle']
for family,ns in [('diagonal',[1,2,16,256]),('rankone',[2,17,128,256]),('path',[2,3,40,128,256]),('star',[2,17,128]),('clustered',[2,17,64]),('random',[2,3,8,32,64,128,256])]:
 for n in ns:
  seed=n+7;a=matrix(family,n,seed);e=np.linalg.eigvalsh(a)
  lines.append(f'{family}\t{n}\t{seed}\t{e[0]:.17g}\t{e[-1]:.17g}')
(root/'crates/smiles/tests/fixtures/symmetric_extrema.tsv').write_text('\n'.join(lines)+'\n')
print(len(smiles),'BCUT inputs;',len(scalar)-1,'scalar controls;',len(lines)-1,'independent eigensolver matrices')
