"""Non-shipped diagnostic of pinned FindRings selection on native raw candidates.

This does NOT implement RDKit's candidate discovery or prepared coordination.
Retain all104 selected original failures, including any raw-topology limitation.
"""
import ast
import csv
import json
import sys
from pathlib import Path

from rdkit import Chem, rdBase

assert rdBase.rdkitVersion == '2026.03.6'


def select(candidates):
    # FindRings.cpp removeExtraRings + symmetrizeSSSR, Release_2026_03_6.
    rings = sorted(candidates, key=len)
    keep = set()
    available = set(range(len(rings)))
    union = set()
    for i, ring in enumerate(rings):
        if ring <= union:
            available.discard(i)
        if i not in available:
            continue
        union |= ring
        keep.add(i)
        consider = {j for j in available if j > i and len(rings[j]) == len(ring)}
        while consider:
            best = min(consider, key=lambda j: (-len(rings[j] & union), j))
            consider.remove(best)
            available.discard(best)
            if not rings[best] <= union:
                keep.add(best)
                union |= rings[best]
    basis = [r for i, r in enumerate(rings) if i in keep]
    counts = {b: sum(b in r for r in basis) for b in union}
    result = list(basis)
    for i, extra in enumerate(rings):
        if i in keep:
            continue
        if any(len(r) == len(extra) and bool(r & extra)
               and all(counts[b] > 1 or b in extra for b in r) for r in basis):
            result.append(extra)
    return result


root = Path(__file__).resolve().parents[2]
source = root/'tasks/evidence/catalog289-sps/ring-count104.csv'
raw = Path(sys.argv[1])
destination = Path(sys.argv[2])
originals = list(csv.DictReader(source.open()))
native_rows = raw.read_text().splitlines()
assert len(originals) == len(native_rows) == 104
params = Chem.SmilesParserParams()
params.removeHs = False
records = []
for source_row, line in zip(originals, native_rows):
    name, n, edges, cycles = line.split('\t')
    assert name == source_row['name']
    if n == 'unsupported':
        records.append(dict(name=name, raw_unsupported=True))
        continue
    n = int(n)
    edges = [tuple(sorted(x)) for x in ast.literal_eval(edges)]
    candidates = [set(edges[i] for i in ring) for ring in json.loads(cycles)]
    mol = Chem.MolFromSmiles(source_row['smiles'], params)
    assert mol is not None
    reference_edges = [tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()]
    same_graph = n == mol.GetNumAtoms() and set(edges) == set(reference_edges)
    gold = [set(reference_edges[i] for i in ring) for ring in mol.GetRingInfo().BondRings()]
    result = select(candidates)
    key = lambda rs: {frozenset(r) for r in rs}
    records.append(dict(name=name, raw_unsupported=False, same_raw_vertex_edge_graph=same_graph,
                        old_count=len(candidates), selected_count=len(result), reference_count=len(gold),
                        count_matches=len(result) == len(gold), exact_sets_match=same_graph and key(result) == key(gold),
                        removed=len(candidates)-len(result)))
report = dict(scope='diagnostic raw-candidate proxy; no shipped change, no full RDKit candidate discovery',
              inputs=104, native_filter=False, raw_unsupported=sum(r['raw_unsupported'] for r in records),
              count_matches=sum(r.get('count_matches', False) for r in records),
              exact_sets_match=sum(r.get('exact_sets_match', False) for r in records), records=records)
destination.write_text(json.dumps(report, indent=2)+'\n')
print({k: v for k, v in report.items() if k != 'records'})
