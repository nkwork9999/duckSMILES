"""Pinned all-catalog controls for workbench aromatic candidate preparation."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,rdMolDescriptors
from compare_catalog import profile_reference
assert rdkit.__version__=='2026.03.6';RDLogger.DisableLog('rdApp.*')
root=Path(__file__).resolve().parents[2]
sources=['[n+]1ccccc1','[N+]1=CC=CC=C1','c1cc[n+]cc1','Cc1cc[n+]cc1',
 '[n+]1ccncc1','[n+]1ccccccc1','[nH+]1ccccc1','C[n+]1ccccc1',
 '[O-][n+]1ccccc1','O=N1=CC=CC=C1','[c]1ccccc1',
 'CP1(C)=CC=CN=C1C','CP1(C)=CC=CC=C1','CP1(C)=CC=CC=N1',
 'O=S1C=CC=C1','O=S1(=O)C=CC=C1','[NH3+]','C[NH2+]','[N+2](C)(C)C','[NH4+]',
 'O=C1OC[C@H]2OC(=O)C3=C(OC(=O)C4=CC(O)=C(O)C(O)=C4C4=C(O)C(O)=C(O)C=C4C(=O)O3)[C@@H]2OC(=O)C2=CC(O)=C(O)C(O)=C2C2=C(O)C(O)=C(O)C=C12','O1C=CC=CC=CC=C1','S1C=CC=CC=CC=C1','O1C=CC=C1','S1C=CC=C1','[O+]1=CC=CC=C1']
sources += ['c1cc[asH]c1','C1=C[AsH]C=C1','[as]1ccccc1',
 'C1=CC=[As]C=C1','c1cc[AsH]c1','c1cc[se]c1','c1cc[te]c1','c1cc[pH]c1','c1cc[sH+]c1']
p=Chem.SmilesParserParams();p.removeHs=False
specs=json.loads((root/'tasks/descriptors/phase1-catalog.json').read_text())['descriptors'];func=dict(Descriptors._descList)
lines=['# Original source, registered/profile ID, expected value; RDKit2026.03.6 removeHs=False']
for source in sources:
 m=Chem.MolFromSmiles(source,p);assert m is not None,source
 profile=profile_reference(m, include_extra=True);lines.append(f'{source}\t@formula\t{rdMolDescriptors.CalcMolFormula(m)}')
 for s in specs:
  id=s['id'];x=func[id](m,force=False) if 'EStateIndex' in id else func[id](m) if s['origin']=='rdkit-standard' else profile[id]
  lines.append(f'{source}\t{id}\t{float(x):.17g}')
(root/'crates/smiles/tests/fixtures/descriptor_aromatic_candidate.tsv').write_text('\n'.join(lines)+'\n')
(root/'tasks/descriptors/aromatic-candidate-controls.json').write_text(json.dumps(dict(rdkit=rdkit.__version__,source_smiles=sources,molecules=len(sources),catalog_count=len(specs),scalar_values=len(specs)*len(sources),formula_values=len(sources),native_acceptance_filter=False,scope='Isolated even pyridinyl-cation radical rings, default-valence eligibility, and preserved aromatic/neutral-carbon-radical controls. Includes first-ring macrocycle neutralO/S>=9 donor exclusion and small-ring/chargedO counterexamples. Not a general kekulizer.'),indent=2)+'\n')
print(len(sources),'sources;',len(specs)*len(sources),'scalar values + formulas')
