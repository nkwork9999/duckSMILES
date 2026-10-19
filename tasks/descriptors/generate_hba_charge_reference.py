"""Independent HBA charge controls, including every v30/full100k HBA failure.

No native acceptance or output is used to calculate expected values. Keep the
original failure cohort visible alongside independent charge/H boundary cases.
"""
import csv
import gzip
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

assert rdBase.rdkitVersion == '2026.03.6'
root = Path(__file__).resolve().parents[2]
with gzip.open(root/'tasks/evidence/catalog289-sps/v30-full100k/differences.jsonl.gz', 'rt') as f:
    records = [r for r in map(json.loads, f) if 'NumHAcceptors' in r.get('differences', {})]
assert len(records) == 13
sources = [r['smiles'] for r in records] + [
    '[O-2]', '[S-2]', '[O-]', '[S-]', '[SH-]', '[S-3]', '[OH+]C',
    '[OH+]c1ccccc1', '[OH2+]C', '[O+]C', 'C[S+]C', '[SH+]C', '[OH+]=C',
    'CO', 'CS', '[H][O+]C', '[H][S+]C', '[O+2]', '[S+2]', '[OH]', '[SH]', '[O]', '[S]',
]
sources = list(dict.fromkeys(sources))
params = Chem.SmilesParserParams()
params.removeHs = False
lines = ['# original SMILES; NumHAcceptors; pinned RDKit2026.03.6; removeHs=False']
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None, source
    lines.append(f'{source}\t{rdMolDescriptors.CalcNumHBA(mol)}')
with (root/'tasks/descriptors/hba-charge-controls.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=['name', 'smiles'])
    writer.writeheader()
    writer.writerows(dict(name=f'hba-charge-{i}', smiles=s) for i, s in enumerate(sources))
(root/'crates/smiles/tests/fixtures/descriptor_hba_charge.tsv').write_text('\n'.join(lines)+'\n')
print(len(sources), 'independent charge/H acceptor sources')
