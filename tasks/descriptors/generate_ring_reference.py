"""Pinned prepared-ring membership and sequence; no native acceptance filter."""
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase

assert rdBase.rdkitVersion == '2026.03.6'
root = Path(__file__).resolve().parents[2]
sources = json.loads((root/'tasks/descriptors/stereo-controls.json').read_text())['source_smiles']
sources += [r['smiles'] for r in csv.DictReader((root/'tasks/evidence/catalog289-sps/ring-count104.csv').open())]
sources = list(dict.fromkeys(sources))
params = Chem.SmilesParserParams()
params.removeHs = False
lines = ['# original SMILES; ordered ring keys (each sorted endpoint-pair set); RDKit2026.03.6 sanitize=True removeHs=False']
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None, source
    keys = []
    for ring in mol.GetRingInfo().BondRings():
        edges = []
        for bi in ring:
            bond = mol.GetBondWithIdx(bi)
            a, b = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            edges.append(f'{a}-{b}')
        keys.append(','.join(sorted(edges)))
    lines.append(source+'\t'+';'.join(keys))
fixture = root/'crates/smiles/tests/fixtures/prepared_rings.tsv'
fixture.write_text('\n'.join(lines)+'\n')
(root/'tasks/descriptors/ring-controls.json').write_text(json.dumps(dict(
    rdkit=rdBase.rdkitVersion, native_filter=False, sources=len(sources),
    source_smiles=sources, fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
    policy='Compare original vertex endpoint pairs, not parser-specific bond indices. Retain ring-list sequence; ignore cycle rotation/direction.'),indent=2)+'\n')
print(len(sources), 'independent prepared-ring references')
