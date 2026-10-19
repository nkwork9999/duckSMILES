"""Independent bond-iteration provenance, including label reuse and dot closures."""
from pathlib import Path
import hashlib
import json
from rdkit import Chem, rdBase

assert rdBase.rdkitVersion == '2026.03.6'
sources = [
    'C1CCC2CCCCC2C1', 'C2CCC1CCCCC1C2', 'C9CCC3CCCCC3C9',
    'C1CC1C1CC1', 'C12.C1C2', 'C2CC2.C1CC1', 'C1CC1.C1CC1',
    'C%11CCC%20CCCCC%20C%11', 'C%20CCC%11CCCCC%11C%20',
    'C1CC1C2CC2C1CC1', 'CC(CO)CC', 'C1(CO)CCC1', 'CC(C1CC1)CC',
    'CC(OC1CC1)C2CC2', 'c2ccc1ccccc1c2', 'F[C@H]1CCC(O)CC1',
    'C/1CCCCCC=C1', 'C1CCCCCC=C/1', 'C-1CCCCC1', 'C1CCCCC-1',
    'N1->BCCC1', 'N1CCCB<-1', 'N->B', 'B<-N',
    'C/1CCCCCC=C=1', 'C=1CCCCCC=C/1', 'C/1CCCCCC=C/1', 'C-1CCCCC=1',
]
params = Chem.SmilesParserParams()
params.removeHs = False
lines = ['# source; ordered oriented begin-end endpoint pairs; ordered bond types; RDKit2026.03.6']
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None, source
    pairs = [f'{b.GetBeginAtomIdx()}-{b.GetEndAtomIdx()}' for b in mol.GetBonds()]
    types = [str(b.GetBondType()) for b in mol.GetBonds()]
    lines.append(source+'\t'+','.join(pairs)+'\t'+','.join(types))
root = Path(__file__).resolve().parents[2]
fixture = root/'crates/smiles/tests/fixtures/smiles_bond_iteration.tsv'
fixture.write_text('\n'.join(lines)+'\n')
(root/'tasks/descriptors/bond-iteration-controls.json').write_text(json.dumps(dict(
    rdkit=rdBase.rdkitVersion, sources=len(sources), native_filter=False,
    source_smiles=sources, fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest()),indent=2)+'\n')
print(len(sources), 'independent oriented bond-iteration controls')
