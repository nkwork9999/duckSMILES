"""Pinned RDKit scalar fixtures and periodic-table facts; never needed by the app."""
import json
from pathlib import Path
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
assert rdkit.__version__ == '2026.03.6', rdkit.__version__
specs=json.loads(Path('tasks/descriptors/phase1-catalog.json').read_text())['descriptors']
functions=dict(Descriptors._descList)
surface_ids={s['id'] for s in json.loads(Path('tasks/descriptors/surface-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/graph-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/charge-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/bcut-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/structural-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/morgan-density-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/qed-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/sps-definitions.json').read_text())['descriptors']}
surface_ids |= {s['id'] for s in json.loads(Path('tasks/descriptors/ipc-definitions.json').read_text())['descriptors']}
smiles=['CCO','CC(=O)Oc1ccccc1C(=O)O','Cn1c(=O)c2c(ncn2C)n(C)c1=O','c1ccccc1','C1CCCCC1','C1=CCCCC1','c1ccc2ccccc2c1','c1ccncc1','c1cc[nH]c1','C12C3C4C1C5C2C3C45','CC(=O)N','N','O','[H]O[H]','[2H]N([2H])[2H]','[NH4+].[Cl-]','[Na+].[Cl-]','[13CH4]','[Fe+2]','[Hg]','CCOP(=O)(O)O','CS(=O)(=O)O','C(F)(Cl)(Br)I','[H][H]','[He]','[H]OC','[H]N(C)C','[H]OC(=O)C','[2H]O[2H]','[H]c1ccccc1']
# Hypervalent ring sulfur must not lend a p orbital to the fused system.
smiles += ['CCN1C2=C(C=CC=C2)[S](=O)C3=C1C=CC(=C3)[N+]([O-])=O',
           'O=S1(=O)C2=CC=CC=C2OC3=C1C=CC=C3',
           'O=S1C2=C(NC3=C1C=CC=C3)C=CC=C2']
params=Chem.SmilesParserParams();params.removeHs=False
smiles += ['C=1CCCCC1','C1CCCCC=1','C1=CCCCC=1','C=1CCCCC=1',
           'C%10CCCCC=%10','C1:C:C:C:C:C:1','C1:C:N:S:N:1',
           'C1:C:C:C2:C(:C:1)CCC2']
out=['# RDKit 2026.03.6; SMILES, registered descriptor ID, expected value']
for s in smiles:
    m=Chem.MolFromSmiles(s,params);assert m is not None,s
    for spec in specs:
        # The separate fragment fixture covers these same inputs plus positive
        # examples for all 85 queries; do not duplicate those cells here.
        if spec['origin']=='rdkit-standard' and not spec['id'].startswith('fr_') and spec['id'] not in surface_ids:
            out.append(f"{s}\t{spec['id']}\t{functions[spec['id']](m):.16g}")
Path('crates/smiles/tests/fixtures').mkdir(exist_ok=True)
Path('crates/smiles/tests/fixtures/descriptor_phase1.tsv').write_text('\n'.join(out)+'\n')
outer=[Chem.GetPeriodicTable().GetNOuterElecs(z) for z in range(1,119)]
p=Path('crates/smiles/src/descriptor_catalog.rs')
s=p.read_text().split('// RDKit 2026.03.6 PeriodicTable facts;')[0].rstrip()
p.write_text(s+'\n\n// RDKit 2026.03.6 PeriodicTable facts; Z=1..118.\npub(crate) const OUTER_ELECTRONS: [i32; 118] = '+str(outer)+';\n')
print(len(specs),'catalog entries;',len(out)-1,'fixed reference values')
