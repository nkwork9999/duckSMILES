"""New stereocounts through C ABI; compare all previous286 to actual Swiftv26.

This does not claim a current Swift/App run or full288 independent chemical parity.
Every original input and unsupported status is retained, with no native filter.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import sys
import rdkit
from compare_catalog import Native, PARAMS, FUNCTIONS, Chem

assert rdkit.__version__ == '2026.03.6', rdkit.__version__
assert Chem.GetUseLegacyStereoPerception() is True

source, library, previous_swift, output = map(Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
inputs = list(csv.DictReader(source.open(newline='')))
old = json.loads(previous_swift.read_text())
assert len(inputs) == len(old) == 10_000
assert len({r['smiles'] for r in inputs}) == 10_000
native = Native(library)
indices = {s['id']: i for i,s in enumerate(native.catalog)}
ids = ['NumAtomStereoCenters','NumUnspecifiedAtomStereoCenters']
assert len(native.catalog) == 288 and all(key in indices for key in ids)
stats = {key:dict(compared=0,mismatches=0) for key in ids}
all_rows = []
differences = []
changed = []
matching = unsupported = invalid = 0
for start in range(0, len(inputs), 64):
    batch = inputs[start:start+64]
    actual = native.analyze([r['smiles'] for r in batch])
    for offset,(record,row) in enumerate(zip(batch,actual)):
        i = start+offset; row['index'] = i
        assert old[i]['index'] == i
        all_rows.append(row)
        common = dict(row)
        if 'descriptor_values' in common: common['descriptor_values'] = common['descriptor_values'][:286]
        if old[i] != common: changed.append(i)
        mol = Chem.MolFromSmiles(record['smiles'],PARAMS)
        detail = {}
        if mol is None or not mol.GetNumAtoms():
            invalid += 1; detail['status'] = 'independent parser rejected'
        else:
            expected = {key:int(FUNCTIONS[key](mol)) for key in ids}
            if row['status'] != 'ok':
                unsupported += 1; detail = dict(status='native invalid or unsupported',reference_counts=expected)
            else:
                for key in ids:
                    value = row['descriptor_values'][indices[key]]; stats[key]['compared'] += 1
                    if value != expected[key]:
                        stats[key]['mismatches'] += 1
                        detail[key] = dict(native=value,rdkit=expected[key])
                if not detail: matching += 1
        if detail: differences.append(dict(index=i,**record,differences=detail))
    if start % 2048 == 0: print(min(start+64,len(inputs)),'/',len(inputs),flush=True)
(output/'native-results.json').write_text(json.dumps(all_rows,separators=(',',':'),allow_nan=False))
with gzip.open(output/'differences.jsonl.gz','wt') as f:
    for row in differences:f.write(json.dumps(row)+'\n')
value_rows = sum(any(k in row['differences'] for k in ids) for row in differences)
report = dict(gate='failed' if differences or changed else 'passed',value_gate='failed' if value_rows else 'passed',
              revision=native.revision,rows=len(inputs),native_catalog_count=288,standard_count=214,additional_count=50,
              matching_new_counts=matching,unsupported=unsupported,reference_invalid=invalid,value_differing_rows=value_rows,
              previous286_and_base_changed_rows=len(changed),previous286_changed_indices=changed,metrics=stats,
              source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
              baseline_swift_sha256=hashlib.sha256(previous_swift.read_bytes()).hexdigest(),
              scope='C ABI only; new2 independent pinned2026.03.6 legacy counts; previous286/base exact delta to actualSwiftv26. No currentSwift/UI or all288 independent parity claim.')
(output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
if differences or changed: raise SystemExit(1)
