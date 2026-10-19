#!/usr/bin/env python3
"""Prepare independent references for the next non-tetrahedral input stage.

Writes evidence only. It does not enable new production grammar or remove tags.
The two still-ambiguous coordination cases stay explicitly excluded from this
focused cohort; the whole unselected population acceptance gate retains them.
"""
import hashlib
import json
import math
from pathlib import Path

import compare_catalog as c


def main():
    root = Path(__file__).resolve().parents[2]
    evidence = root / 'tasks/evidence/catalog289-sps'
    native = c.Native(Path('/private/tmp/chemob-descriptor-eval/libducksmiles-v38.dylib'))
    assert c.rdkit.__version__ == '2026.03.6'
    experiment = json.loads((evidence / 'v38-non-tetra-base-experiment.json').read_text())
    sources = [r['original_smiles'] for r in experiment['records'] if r['candidate_status'] == 'ok']
    excluded = [r['source_index'] for r in experiment['records'] if r['candidate_status'] != 'ok']
    assert len(sources) == 34 and excluded == [18084, 53180]
    for geometry, maximum, center, ligands in [
        ('SP', 3, 'Pt', ['F', 'Cl', 'Br', 'I']),
        ('TB', 20, 'Fe', ['F', 'Cl', 'Br', 'I', 'N']),
        ('OH', 30, 'Co', ['F', 'Cl', 'Br', 'I', 'N', 'O']),
    ]:
        def structure(tag, degree, h='', central=center):
            selected = ligands[:degree]
            return selected[0] + '[' + central + '@' + tag + h + ']' + ''.join('('+x+')' for x in selected[1:])
        for permutation in range(maximum + 1):
            sources.append(structure(geometry + str(permutation), len(ligands)))
        sources.append(structure(geometry, len(ligands)))
        for degree in range(2, len(ligands)+1):
            sources.extend([structure(geometry+'1', degree), structure(geometry+'1', degree-1, 'H')])
    for center in ['C', 'N+', 'P', 'S', 'As', 'Se']:
        sources.append('F[' + center.split('+')[0] + '@SP1' + ('+' if '+' in center else '') + '](Cl)(Br)I')
    sources.extend([
        'F[Pt@SP-](Cl)(Br)I', 'F[Pt@SP1-2](Cl)(Br)I',
        'F[57Fe@TB3+2](Cl)(Br)(I)N', 'F[Co@OH+](Cl)(Br)(I)(N)O',
        'F[Pt@SP4](Cl)(Br)I', 'F[Fe@TB21](Cl)(Br)(I)N',
        'F[Co@OH31](Cl)(Br)(I)(N)O',
        'F[Pt@SP1](Cl)(Br)(I)N', 'F[Fe@TB1](Cl)(Br)(I)(N)O',
        'F[Co@OH1](F)(Cl)(Br)(I)(N)O',
        'F[Pt@SP1H4]', 'F[Fe@TB1H5]', 'F[Co@OH1H6]',
        'F[C@SP1H]', 'F[C@SP1H](Cl)Br', 'F[S@SP1H]',
        'F[P@TB1H]', 'F[Se@OH1H]',
    ])
    sources = list(dict.fromkeys(sources))
    accepted, rejected, rows = [], [], ['# Original non-tetrahedral source\tformula\t291 independent reference values']
    states = []
    for source in sources:
        mol = c.Chem.MolFromSmiles(source, c.PARAMS)
        if mol is None:
            rejected.append(source)
            continue
        profile, extra = c.profile_reference(mol), c.extra_reference(mol)
        info = c.exact_information_reference(mol)
        values = []
        for spec in native.catalog:
            value = (float(info[spec['id']]) if spec['id'] in c.STABLE_INFORMATION_IDS
                     else c.reference_value(spec, mol, profile, extra))
            values.append(repr(value) if math.isfinite(value) else 'null')
        rows.append('\t'.join([source, c.rdMolDescriptors.CalcMolFormula(mol), *values]))
        accepted.append(source)
        states.append({'smiles': source, 'atoms': [
            {'index': a.GetIdx(), 'tag': str(a.GetChiralTag()), 'hybridization': str(a.GetHybridization()),
             'degree': a.GetDegree(), 'total_degree': a.GetTotalDegree(),
             'explicit_h': a.GetNumExplicitHs(), 'implicit_h': a.GetNumImplicitHs(),
             'legacy_chirality_possible': bool(a.HasProp('_ChiralityPossible') and a.GetIntProp('_ChiralityPossible')),
             'cip_code': a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None,
             'chiral_permutation': a.GetIntProp('_chiralPermutation') if a.HasProp('_chiralPermutation') else None}
            for a in mol.GetAtoms()]})
    directory = evidence / 'non-tetra-next-controls'
    directory.mkdir(exist_ok=True)
    path = directory / 'descriptor_non_tetra.tsv'
    path.write_text('\n'.join(rows)+'\n')
    (directory / 'non_tetra_rejected.txt').write_text('\n'.join(rejected)+'\n')
    (directory / 'prepared-atom-states.json').write_text(json.dumps(states, indent=2)+'\n')
    summary = {'scope': 'Next-stage independent ORIGINAL non-tetrahedral references; production support not enabled',
               'rdkit': c.rdkit.__version__, 'population_rows': 34,
               'excluded_ambiguous_coordination_source_indices': excluded,
               'accepted_rows': len(accepted), 'rejected_reference_sources': rejected,
               'ids': [s['id'] for s in native.catalog],
               'fixture_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    (directory / 'manifest.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'ids'}))


if __name__ == '__main__':
    main()
