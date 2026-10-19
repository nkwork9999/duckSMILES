#!/usr/bin/env python3
"""Freeze a deterministic distinct-graph fixture from ChEBI's structure registry.

The source CANONICAL_SMILES is preserved verbatim. RDKit canonical isomeric
SMILES is used only as an identity key; no neutralization, salt/tautomer removal,
stereo merging or selection by native acceptance occurs.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import rdkit
from rdkit import Chem, RDLogger


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--count', type=int, default=100_000)
    args = ap.parse_args()
    assert rdkit.__version__ == '2026.03.6'
    assert args.count >= 10_000
    RDLogger.DisableLog('rdApp.*')
    params = Chem.SmilesParserParams()
    params.removeHs = False
    candidates = []
    source_rows = 0
    with gzip.open(args.source, 'rt', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            source_rows += 1
            layers = json.loads(row['layers'])
            smiles = layers.get('CANONICAL_SMILES') if isinstance(layers, dict) else None
            if smiles:
                # SHA ordering avoids selecting the first deposition block.
                candidates.append((hashlib.sha256(smiles.encode()).hexdigest(), row['structure_id'], smiles))
    candidates.sort()
    seen = set()
    selected = []
    scanned = duplicates = invalid = generic = 0
    for _, id, smiles in candidates:
        scanned += 1
        mol = Chem.MolFromSmiles(smiles, params)
        if mol is None or mol.GetNumAtoms() == 0:
            invalid += 1
            continue
        # R groups and polymer attachment points describe families, not a
        # definite input compound. Keep them in a separate stress corpus.
        if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
            generic += 1
            continue
        key = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if key in seen:
            duplicates += 1
            continue
        seen.add(key)
        selected.append(dict(name=f'ChEBI-structure-{id}', smiles=smiles))
        if len(selected) % 10_000 == 0:
            print(f'selected {len(selected)} distinct graphs; scanned {scanned}', flush=True)
        if len(selected) == args.count:
            break
    assert len(selected) == args.count, (len(selected), args.count)
    args.output.mkdir(parents=True, exist_ok=True)
    files = {}
    for count in (10_000, args.count):
        path = args.output / f'chebi-{count}.csv'
        with path.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['name', 'smiles'])
            writer.writeheader()
            writer.writerows(selected[:count])
        files[path.name] = dict(rows=count, bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    manifest = dict(
        schema=1, source_url='https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/structure_registry.tsv.gz',
        source_path=str(args.source), source_bytes=args.source.stat().st_size,
        source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(), source_rows=source_rows,
        candidates_with_smiles=len(candidates), scanned=scanned, duplicate_graphs_skipped=duplicates,
        invalid_rdkit_skipped=invalid, generic_graphs_skipped=generic,
        distinct_graphs=len(selected), rdkit=rdkit.__version__,
        identity='canonical isomeric SMILES, removeHs=False; no standardization',
        selection='ascending SHA256(source SMILES), tie by structure_id; exclude dummy-atom generic structures, stop after requested distinct count; no native acceptance filter',
        names='ChEBI structure record IDs, not ChEBI entity accession IDs', files=files,
    )
    (args.output / 'selection.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest), flush=True)


if __name__ == '__main__':
    main()
