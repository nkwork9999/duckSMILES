#!/usr/bin/env python3
"""Independent large-domain boundaries and a >256-atom all-descriptor control."""
import hashlib
import json
import math
from pathlib import Path

import compare_catalog as c
from exact_ipc_reference import information
from flint_ipc_reference import coefficients


def clique(n, atom='C'):
    labels = [[] for _ in range(n)]
    index = 0
    for i in range(n):
        for j in range(i+2, n):
            index += 1
            token = str(index) if index < 10 else '%' + str(index)
            labels[i].append(token)
            labels[j].append(token)
    return ''.join(atom + ''.join(row) for row in labels)


def main():
    root = Path(__file__).resolve().parents[2]
    directory = root / 'tasks/evidence/catalog289-sps/large-next-controls'
    directory.mkdir(exist_ok=True)
    assert c.rdkit.__version__ == '2026.03.6'
    native = c.Native(Path('/private/tmp/chemob-descriptor-eval/libducksmiles-v39.dylib'))
    k5 = clique(5)
    boundaries = [
        ('path257', 'C'*257, True), ('path640', 'C'*640, True),
        ('atom641', 'C'*641, False),
        ('bond1280', '.'.join([k5]*128), True),
        ('bond1281', '.'.join([k5]*126 + [clique(6, '[Fe]'), clique(4)]), False),
        ('bytes8320', '[C:999999999]'*640, False),
    ]
    lines = ['# name\toriginal source\twithin prepared budgets\tatoms\tbonds\tformula']
    for name, source, supported in boundaries:
        mol = c.Chem.MolFromSmiles(source, c.PARAMS)
        assert mol is not None, name
        allowed = len(source.encode()) <= 4096 and mol.GetNumAtoms() <= 640 and mol.GetNumBonds() <= 1280
        assert allowed == supported, (name, len(source.encode()), mol.GetNumAtoms(), mol.GetNumBonds())
        lines.append('\t'.join([name, source, str(int(supported)), str(mol.GetNumAtoms()),
                                str(mol.GetNumBonds()), c.rdMolDescriptors.CalcMolFormula(mol)]))
    domain = directory / 'large_domain.tsv'
    domain.write_text('\n'.join(lines)+'\n')
    source = 'C'*257
    mol = c.Chem.MolFromSmiles(source, c.PARAMS)
    profile, extra = c.profile_reference(mol), c.extra_reference(mol)
    adjacency = [[j for j in (i-1, i+1) if 0 <= j < 257] for i in range(257)]
    stable = information(coefficients(adjacency))
    values = []
    for spec in native.catalog:
        value = stable[spec['id']] if spec['id'] in c.STABLE_INFORMATION_IDS else c.reference_value(spec, mol, profile, extra)
        values.append(repr(value) if math.isfinite(value) else 'null')
    path = directory / 'descriptor_path257.tsv'
    path.write_text('# ORIGINAL source\tformula\t291 independent values\n' + '\t'.join([
        source, c.rdMolDescriptors.CalcMolFormula(mol), *values]) + '\n')
    summary = {'rdkit': c.rdkit.__version__, 'domain_controls': len(boundaries),
               'descriptor_controls': 1, 'descriptor_count': len(values),
               'coefficient_reference': 'python-flint0.9.0/FLINT3.6.0; matched58 analytic/determinant and four large path controls',
               'sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [domain, path]},
               'scope': 'Prepared API candidate domain only; raw graph sidecar remains256atoms/512bonds'}
    (directory/'manifest.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
