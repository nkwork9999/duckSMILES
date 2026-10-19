"""Independent all-catalog controls for pinned N/P/halogen cleanup rules."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from compare_catalog import profile_reference

assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
root=Path(__file__).resolve().parents[2]
specs=json.loads((root/'tasks/descriptors/phase1-catalog.json').read_text())['descriptors']
functions=dict(Descriptors._descList)
sources=[
 'CN(=O)=O','O=N(=O)O','O=N(=O)OO','CN(C)(C)=O','[NH](=O)=C',
 'C1=CC=CN(=O)=C1','n1(=O)ccccc1','CN=N#N','NN=N#N','O=N#N',
 '[18O]=N(=O)O','O=N(=[18O])O','[2H]ON(=O)=O','[OH]N(=O)=O',
 'O=N(=O)c1ccccc1','O=N(=O)C1=CC=CC=C1',
 'C=P(=O)O','C=P(=O)C','CN=P(=O)O','N=P(=O)O','CP(=O)=NC',
 'C=P(=O)[O-]','[13CH2]=P(=O)O','CP(=O)=C','CP(=O)=CC',
 'CP(C)(C)=O','CP(=O)(O)O','NP(=O)=O',
 'O=ClO','O=BrO','O=IO','O=Cl(=O)O','O=Br(=O)O','O=I(=O)O',
 'O=Cl(=O)(=O)[O-]','O=Br(=O)(=O)[O-]','O=I(=O)(=O)[O-]',
 'O=[Cl+](O)O','O=[Br+](O)O','CI(=O)=O',
 'CCO.CN(=O)=O','CN(C)(C)=O.O=Cl(=O)(=O)[O-]',
 'C[N+](=O)[O-]','[O-][N+](=O)O','C[N+](C)(C)[O-]',
 '[H]C=P(=O)O','[2H]C=P(=O)O','O=[Br+]O','O=[Cl+]O',
]
sources=list(dict.fromkeys(sources))
params=Chem.SmilesParserParams();params.removeHs=False
lines=['# Registered RDKit2026.03.6 + independent profile oracle; original source input, descriptor/formula, expected']
for smiles in sources:
    mol=Chem.MolFromSmiles(smiles,params)
    assert mol is not None,smiles
    profile=profile_reference(mol, include_extra=True)
    lines.append(f'{smiles}\t@formula\t{rdMolDescriptors.CalcMolFormula(mol)}')
    for spec in specs:
        id=spec['id']
        value=functions[id](mol) if spec['origin']=='rdkit-standard' else profile[id]
        lines.append(f'{smiles}\t{id}\t{float(value):.17g}')
(root/'crates/smiles/tests/fixtures/descriptor_cleanup.tsv').write_text('\n'.join(lines)+'\n')
metadata=dict(rdkit=rdkit.__version__,molecules=len(sources),catalog_count=len(specs),scalar_values=len(sources)*len(specs),formula_values=len(sources),source_smiles=sources,policy=dict(removeHs=False,sanitize=True,standardize=False),scope='All listed inputs retained; no native acceptance or agreement filtering.')
(root/'tasks/descriptors/cleanup-controls.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(len(sources),'source inputs;',metadata['scalar_values'],'scalar expectations plus formula')
