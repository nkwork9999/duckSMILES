#!/usr/bin/env python3
"""Independent original/permuted metal-cluster references, never native filtered.

RDKit2026.03.6 supplies raw canonical ranks and sanitized descriptor definitions.
Random graph controls test ranking only, without claiming chemical validity.
"""
import hashlib
import json
import math
import random
from pathlib import Path

import compare_catalog as c


def main():
    root = Path(__file__).resolve().parents[2]
    evidence = root / 'tasks/evidence/catalog289-sps'
    output = evidence / 'cluster-next-controls'
    output.mkdir(exist_ok=True)
    assert c.rdkit.__version__ == '2026.03.6'
    native = c.Native(Path('/private/tmp/chemob-descriptor-eval/libducksmiles-v41-factor2.dylib'))
    original = json.loads((evidence / 'remaining-metal5-diagnostic.json').read_text())
    params = c.Chem.SmilesParserParams()
    params.sanitize = False
    params.removeHs = False
    rng = random.Random(20261018)
    sources = []
    for row in original:
        mol = c.Chem.MolFromSmiles(row['smiles'], params)
        mol.UpdatePropertyCache(False)
        sources.append(row['smiles'])
        for _ in range(40):
            order = list(range(mol.GetNumAtoms()))
            rng.shuffle(order)
            permuted = c.Chem.RenumberAtoms(mol, order)
            sources.append(c.Chem.MolToSmiles(permuted, canonical=False, allHsExplicit=True))
    sources.extend(['[O-]([Fe])[Ni]', '[Cl]([Au])[Ag]', '[O]([Fe])[Ni]',
                    '[Cl]1[Au][Cl][Ag]1', '[Cl]1[Au]([Cl])[Cl][Au]1([Cl])',
                    '[37Cl]([Au:4])[Au:5]', '[OH]([Fe])[Ni]'])
    sources = list(dict.fromkeys(sources))
    rank_sources = [s for s in sources if '@' not in s]
    # Independently generated simple graphs include maps, isotope, H and charge
    # distinctions, isolated vertices, rings, fused cages and tied partitions.
    for trial in range(350):
        rw = c.Chem.RWMol()
        n = rng.randrange(1, 25)
        for i in range(n):
            a = c.Chem.Atom(rng.choice([6, 8, 16, 17, 26, 79]))
            a.SetNoImplicit(True)
            if trial % 3 == 0:
                a.SetFormalCharge(rng.choice([-1, 0, 0, 1]))
            if trial % 5 == 0:
                a.SetAtomMapNum(rng.randrange(3))
            if trial % 7 == 0:
                a.SetIsotope(rng.choice([0, 0, 13]))
            if trial % 11 == 0:
                a.SetNumExplicitHs(rng.randrange(2))
            rw.AddAtom(a)
        for i in range(n):
            for j in range(i):
                if rng.random() < .14:
                    rw.AddBond(i, j, c.Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(False)
        rank_sources.append(c.Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True))
    rank_lines = ['# raw original SMILES\tpinned full canonical ranks']
    for source in dict.fromkeys(rank_sources):
        mol = c.Chem.MolFromSmiles(source, params)
        assert mol is not None, source
        mol.UpdatePropertyCache(False)
        ranks = list(c.Chem.CanonicalRankAtoms(mol))
        rank_lines.append(source + '\t' + ','.join(map(str, ranks)))
    ranks_path = root / 'crates/smiles/tests/fixtures/cluster_ranks.tsv'
    ranks_path.write_text('\n'.join(rank_lines) + '\n')
    rows = ['# original source\tformula\t291 independent descriptor values']
    accepted, rejected = [], []
    for source in sources:
        mol = c.Chem.MolFromSmiles(source, c.PARAMS)
        if mol is None:
            rejected.append(source)
            continue
        profile, extra = c.profile_reference(mol), c.extra_reference(mol)
        values = []
        for spec in native.catalog:
            value = c.reference_value(spec, mol, profile, extra, information_backend='flint')
            values.append(repr(value) if math.isfinite(value) else 'null')
        rows.append('\t'.join([source, c.rdMolDescriptors.CalcMolFormula(mol), *values]))
        accepted.append(source)
    path = root / 'crates/smiles/tests/fixtures/descriptor_clusters.tsv'
    path.write_text('\n'.join(rows) + '\n')
    assert all(r['smiles'] in accepted for r in original)
    import csv
    with (output/'clusters.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'smiles'])
        writer.writerows((f'cluster-{i}', s) for i, s in enumerate(accepted))
    manifest = dict(rdkit=c.rdkit.__version__, original_source_indices=[r['index'] for r in original],
                    native_filter=False, rank_controls=len(rank_lines)-1,
                    descriptor_controls=len(accepted), rejected_by_reference=rejected,
                    catalog=[s['id'] for s in native.catalog],
                    sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [ranks_path, path]},
                    ranking_scope='All bracket, single bond, non-stereo; <=64 atoms. Random graphs are ranking-only controls.')
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k != 'catalog'}))


if __name__ == '__main__':
    main()
