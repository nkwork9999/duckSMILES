#!/usr/bin/env python3
"""Read-only experiment: compare a base-graph candidate to ORIGINAL CX references.

This does not change production parsing or establish support for CXSMILES.
It records every source string, transformation, reference error and difference.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

import compare_catalog as comparison


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--library', type=Path, required=True)
    ap.add_argument('--differences', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    native = comparison.Native(args.library)
    assert comparison.rdkit.__version__ == '2026.03.6'
    with gzip.open(args.differences, 'rt') as f:
        inputs = [json.loads(line) for line in f]
    records = []
    for record in inputs:
        source = record['smiles']
        if '|' not in source:
            continue
        assert source.count('|') == 2 and not source.rsplit('|', 1)[1].strip()
        candidate, annotation, _ = source.split('|')
        candidate = candidate.rstrip()
        original = comparison.Chem.MolFromSmiles(source, comparison.PARAMS)
        assert original is not None
        actual = native.analyze([candidate])[0]
        result = {
            'source_index': record['index'], 'name': record['name'],
            'original_smiles': source, 'candidate_smiles': candidate,
            'cx_annotation': annotation, 'candidate_status': actual['status'],
            'original_atom_count': original.GetNumAtoms(), 'differences': {},
            'reference_errors': {},
        }
        if actual['status'] == 'ok':
            formula = comparison.rdMolDescriptors.CalcMolFormula(original)
            if actual['formula'] != formula:
                result['differences']['formula'] = [actual['formula'], formula]
            profile = comparison.profile_reference(original)
            extra = comparison.extra_reference(original)
            info = comparison.exact_information_reference(original)
            for spec, value in zip(native.catalog, actual['descriptor_values']):
                id = spec['id']
                try:
                    expected = (float(info[id]) if id in comparison.STABLE_INFORMATION_IDS
                                else comparison.reference_value(spec, original, profile, extra))
                except Exception as exc:
                    result['reference_errors'][id] = repr(exc)
                    continue
                expected = expected if math.isfinite(expected) else None
                if (value is None) != (expected is None) or (
                    value is not None and expected is not None
                    and abs(value - expected) > comparison.allowed_error(id, expected)
                ):
                    result['differences'][id] = [value, expected]
        records.append(result)
        print(f"CX {len(records)} original={record['index']} candidate={actual['status']} "
              f"differences={len(result['differences'])} errors={len(result['reference_errors'])}",
              flush=True)
    summary = {
        'scope': 'Read-only base-graph experiment against ORIGINAL annotated molecule; no production support claim',
        'library_revision': native.revision,
        'library_sha256': hashlib.sha256(args.library.read_bytes()).hexdigest(),
        'difference_source_sha256': hashlib.sha256(args.differences.read_bytes()).hexdigest(),
        'rdkit': comparison.rdkit.__version__,
        'cx_rows': len(records),
        'candidate_accepted': sum(r['candidate_status'] == 'ok' for r in records),
        'all_291_and_formula_matching_original': sum(
            r['candidate_status'] == 'ok' and not r['differences'] and not r['reference_errors']
            for r in records),
        'records': records,
    }
    args.output.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'records'}))


if __name__ == '__main__':
    main()
