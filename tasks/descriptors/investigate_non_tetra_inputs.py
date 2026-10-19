#!/usr/bin/env python3
"""Read-only non-tetrahedral projection experiment; NOT production normalization."""
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

import compare_catalog as c


def main():
    root = Path(__file__).resolve().parents[2]
    evidence = root / 'tasks/evidence/catalog289-sps'
    lib = Path('/private/tmp/chemob-descriptor-eval/libducksmiles-v38.dylib')
    native = c.Native(lib)
    with gzip.open(evidence / 'v37-full100k/differences.jsonl.gz', 'rt') as f:
        source = [json.loads(line) for line in f]
    records = []
    for row in source:
        original = row['smiles']
        if not re.search(r'@(SP|TB|OH)\d*', original):
            continue
        candidate = re.sub(r'@(SP|TB|OH)\d*', '', original)
        mol = c.Chem.MolFromSmiles(original, c.PARAMS)
        assert mol is not None
        actual = native.analyze([candidate])[0]
        d = {'source_index': row['index'], 'name': row['name'],
             'original_smiles': original, 'candidate_smiles': candidate,
             'candidate_status': actual['status'], 'original_atoms': mol.GetNumAtoms(),
             'differences': {}, 'reference_errors': {}}
        if actual['status'] == 'ok':
            profile, extra = c.profile_reference(mol), c.extra_reference(mol)
            information = c.exact_information_reference(mol)
            if actual['formula'] != c.rdMolDescriptors.CalcMolFormula(mol):
                d['differences']['formula'] = [actual['formula'], c.rdMolDescriptors.CalcMolFormula(mol)]
            for spec, x in zip(native.catalog, actual['descriptor_values']):
                id = spec['id']
                try:
                    y = (float(information[id]) if id in c.STABLE_INFORMATION_IDS
                         else c.reference_value(spec, mol, profile, extra))
                except Exception as exc:
                    d['reference_errors'][id] = repr(exc)
                    continue
                y = y if math.isfinite(y) else None
                if (x is None) != (y is None) or (x is not None and y is not None
                    and abs(x-y) > c.allowed_error(id, y)):
                    d['differences'][id] = [x, y]
        records.append(d)
        print(row['index'], actual['status'], d['differences'], flush=True)
    summary = {'scope': 'Read-only stripping experiment against ORIGINAL non-tetrahedral references; no support or stereo-equivalence claim',
               'library_revision': native.revision, 'library_sha256': hashlib.sha256(lib.read_bytes()).hexdigest(),
               'rows': len(records), 'records': records}
    (evidence / 'v38-non-tetra-base-experiment.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
