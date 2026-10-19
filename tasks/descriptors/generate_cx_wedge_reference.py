#!/usr/bin/env python3
"""Pinned ORIGINAL CX references, independent of the candidate implementation."""
import ctypes
import hashlib
import json
import math
from pathlib import Path

import compare_catalog as c


ROOT = Path(__file__).resolve().parents[2]


def main():
    native = c.Native(Path('/private/tmp/chemob-descriptor-eval/libducksmiles-v37.dylib'))
    experiment_path = ROOT / 'tasks/evidence/catalog289-sps/v37-cx-base-experiment.json'
    experiment = json.loads(experiment_path.read_text())
    population = [r['original_smiles'] for r in experiment['records']
                  if r['original_atom_count'] <= 256]
    assert len(population) == 107
    bases = [
        'CC(O)F', 'C[C@H](O)F', 'N[C@@H](C)C(=O)O', 'F/C=C/F', 'F/C=C\\F',
        'CC=C(C)F', 'C1CCC(C)CC1', 'C1CC2CCC1C2',
        'c1ccc(-c2ccccc2)cc1', 'CC1=CC=CC(I)=C1N1C(C)=CC=C1Br',
        'C[C@H](O)c1cccc(-c2ccccc2F)c1', 'CS(C)=O', 'C[N+](C)(C)C',
        '[H][C@](F)(Cl)Br',
    ]
    fn = native.lib.ds_mol_graph_json
    fn.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t]
    fn.restype = ctypes.c_int32
    sources = population.copy()
    for base in bases:
        data = base.encode()
        size = fn(data, len(data), None, 0)
        assert size > 0
        out = ctypes.create_string_buffer(size)
        assert fn(data, len(data), out, size) == size
        graph = json.loads(out.raw)
        # Endpoint/index enumeration is a stimulus only; all expected chemistry
        # below is computed from the ORIGINAL annotated string with pinned RDKit.
        for bond in graph['bonds']:
            for atom in [bond['begin'], bond['end']]:
                for direction in ['wU', 'wD']:
                    sources.append(f"{base} |{direction}:{atom}.{bond['index']}|")
    sources.extend([
        'CC1=CC=CC(I)=C1N1C(C)=CC=C1Br |wD:7.7,wU:8.9|',
        'CC1=CC=CC(I)=C1N1C(C)=CC=C1Br |wD:7.7,wU:8.14|',
        'CC(O)F |wU:1.0,1.1,wD:1.2|',
    ])
    assert len(sources) == len(set(sources))
    lines = ['# ORIGINAL CXSMILES\tformula\t291 declared-reference values (null=undefined)']
    reference_errors = []
    for i, source in enumerate(sources):
        mol = c.Chem.MolFromSmiles(source, c.PARAMS)
        assert mol is not None, source
        profile, extra = c.profile_reference(mol), c.extra_reference(mol)
        information = c.exact_information_reference(mol)
        values = []
        for spec in native.catalog:
            try:
                value = (float(information[spec['id']])
                         if spec['id'] in c.STABLE_INFORMATION_IDS
                         else c.reference_value(spec, mol, profile, extra))
            except Exception as exc:
                reference_errors.append([i, source, spec['id'], repr(exc)])
                continue
            values.append(repr(value) if math.isfinite(value) else 'null')
        assert len(values) == 291, reference_errors[-5:]
        lines.append('\t'.join([source, c.rdMolDescriptors.CalcMolFormula(mol), *values]))
    path = ROOT / 'crates/smiles/tests/fixtures/descriptor_cx_wedge.tsv'
    path.write_text('\n'.join(lines) + '\n')
    summary = {
        'rdkit': c.rdkit.__version__, 'population_rows': len(population),
        'synthetic_rows': len(sources) - len(population), 'all_rows': len(sources),
        'all_values_from_original_annotated_reference': True,
        'columns': [s['id'] for s in native.catalog], 'reference_errors': reference_errors,
        'fixture_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'experiment_sha256': hashlib.sha256(experiment_path.read_bytes()).hexdigest(),
        'stimulus_graph_library_sha256': experiment['library_sha256'],
        'scope': 'Coordinate-free wU/wD only; original source retained. Not general CX or atropisomer depiction support.',
    }
    (ROOT / 'tasks/evidence/catalog289-sps/cx-wedge-reference.json').write_text(
        json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'columns'}))


if __name__ == '__main__':
    main()
