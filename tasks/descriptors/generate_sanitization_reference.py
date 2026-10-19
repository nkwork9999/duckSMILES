"""Pinned parser outcomes for explicit valence, charge, aromatic and dative cases.

Original strings, no native acceptance filtering. Includes charged hypervalent
P/S/As/Se and two-coordinate hydride exceptions; not a universal sanitizer oracle.
"""
import hashlib,json
from pathlib import Path
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolDescriptors
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')
params=Chem.SmilesParserParams();params.removeHs=False
sources=['C(C)(C)(C)(C)C','N(C)(C)(C)C','[CH5]','[NH5]','O(C)(C)C','F(F)F','[S-](F)(F)(F)(F)(F)F',
 '[P-](F)(F)(F)(F)(F)F','[As-](F)(F)(F)(F)(F)F','[Se-](F)(F)(F)(F)F',
 '[H-](C)C','[H](C)C','[H+](C)C','CN(=O)=O','[NH3]->[Cu+2]<-[NH3]',
 '[nH]1cccccc1','c1cc[nH]c1','[n+]1ccccc1','C1:C:C:C:1','C1:C:C:C:C:C#1','C1=CC#CC=C1']
periodic=Chem.GetPeriodicTable()
charges=[-8,-4,-3,-2,-1,0,1,2,3,4,8]
def charge(q):return (('+' if q>0 else '-')+str(abs(q))) if q else ''
for z in range(1,119):
 symbol=periodic.GetElementSymbol(z)
 for q in charges:
  for h in range(8):
   sources.append(f'[{symbol}H{h}{charge(q)}]')
for symbol in ['H','B','C','N','O','F','Si','P','S','Cl','As','Se','Br','I','Na','Mg','Al','Ca','Fe','Cu','Zn','Sn','Te','Xe']:
 for q in [-3,-2,-1,0,1,2,3]:
  for degree in range(9):
   for neighbor in ['F','C','[H]']:
    sources.append(f'[{symbol}{charge(q)}]'+f'({neighbor})'*degree)
sources=list(dict.fromkeys(sources))
root=Path(__file__).resolve().parents[2]
lines=['# Original SMILES; pinned sanitized acceptance; formula; removeHs=False; no native filter']
accepted=0
for source in sources:
 mol=Chem.MolFromSmiles(source,params)
 ok=mol is not None and mol.GetNumAtoms()>0
 accepted+=ok
 lines.append('\t'.join([source,str(int(ok)),rdMolDescriptors.CalcMolFormula(mol) if ok else '-']))
dest=root/'crates/smiles/tests/fixtures/sanitization_valence.tsv';dest.write_text('\n'.join(lines)+'\n')
report=dict(rdkit=rdkit.__version__,rows=len(sources),accepted=accepted,rejected=len(sources)-accepted,native_filter=False,removeHs=False,fixture_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),scope='charged isolated atoms plus attached explicit H across118 elements; finite main-group/metal ligand stars; aromatic/coordination regressions; strict property-cache acceptance, not Atom.hasValenceViolation')
(root/'tasks/descriptors/sanitization-valence-controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
