"""Validate extra50 in actual Swift rows; retain every unsupported input.

No native recomputation, selection by native acceptance or RDKit built-in
descriptor claim. RDKit parses the independent dense-formula oracle's graph.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import sys
import rdkit
from rdkit import Chem, RDLogger
from extra50_reference import reference, IDS

assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
source, results, catalog, output = map(Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
inputs = list(csv.DictReader(source.open(newline='')))
actual = json.loads(results.read_text())
specs = json.loads(catalog.read_text())['descriptors']
indices = {s['id']: i for i,s in enumerate(specs)}
selected = [s for s in specs if s['origin'] == 'native-extra']
assert set(IDS) == {s['id'] for s in selected} and len(selected) == 50
assert len(inputs) == len(actual) and [r['index'] for r in actual] == list(range(len(inputs)))
params = Chem.SmilesParserParams(); params.removeHs = False
metrics = {s['id']:dict(compared=0,matching_undefined=0,mismatches=0,max_abs_error=0,tolerance=s['tolerance']) for s in selected}
matching = unsupported = invalid = 0
differences = []
for i, (record,row) in enumerate(zip(inputs,actual)):
    mol = Chem.MolFromSmiles(record['smiles'], params)
    detail = {}
    if mol is None or not mol.GetNumAtoms():
        invalid += 1; detail['status'] = 'independent parser rejected'
    elif row['status'] != 'ok':
        unsupported += 1; detail['status'] = 'native invalid or unsupported'
    else:
        assert len(row['descriptor_values']) == len(specs)
        expected = reference(mol)
        for spec in selected:
            key = spec['id']; m = metrics[key]; m['compared'] += 1
            a,b = row['descriptor_values'][indices[key]], expected[key]
            assert a is None or math.isfinite(a)
            if a is None and b is None:
                m['matching_undefined'] += 1; continue
            error = abs(a-b) if a is not None and b is not None else None
            if error is not None:
                m['max_abs_error'] = max(m['max_abs_error'],error)
                if error <= spec['tolerance']: continue
            m['mismatches'] += 1
            detail[key] = dict(native=a,independent_formula=b,abs_error=error)
        if not detail: matching += 1
    if detail: differences.append(dict(index=i,**record,differences=detail))
    if (i+1) % 2048 == 0: print(i+1,'/',len(inputs),flush=True)
with gzip.open(output/'differences.jsonl.gz','wt') as f:
    for row in differences:f.write(json.dumps(row)+'\n')
value_differences = sum(any(k!='status' for k in d['differences']) for d in differences)
report = dict(gate='failed' if differences else 'passed',value_gate='failed' if value_differences else 'passed',
              rows=len(inputs),distinct_raw_smiles=len({r['smiles'] for r in inputs}),native_catalog_count=len(specs),
              additional_count=50,all_extra_values_matching_rows=matching,unsupported=unsupported,
              independent_parser_invalid=invalid,value_differing_rows=value_differences,metrics=metrics,
              source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              actual_swift_results_sha256=hashlib.sha256(results.read_bytes()).hexdigest(),
              reference='extra50_reference.py independent dense formulas; RDKit2026.03.6 parses original source only',
              scope='Only extra50, not all RDKit standards. Every original input retained; no native recomputation/filter.')
(output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='metrics'}))
if differences: raise SystemExit(1)
