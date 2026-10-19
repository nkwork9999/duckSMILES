#!/usr/bin/env python3
"""Independent CSV-module check of actual Swift Workbench output and native JSON."""
import csv,gzip,json,math,hashlib
from pathlib import Path
from compare_rdkit import FIELDS
p=Path(__file__).resolve().parents[1]/'evidence'/'rdkit'
with (p/'bulk-input.csv').open(newline='') as f: source=list(csv.DictReader(f))
with (p/'bulk-output.csv').open(newline='') as f: exported=list(csv.DictReader(f))
results=json.loads((p/'bulk-native-results.json').read_text())
baseline=json.loads(gzip.decompress((p/'baseline-native-results.json.gz').read_bytes()))
assert len(source)==len(exported)==len(results)==len(baseline)==10000
assert results==baseline, 'speed optimization must preserve every output value'
for index,(a,b,r) in enumerate(zip(source,exported,results)):
    assert r['index']==index
    assert a['name']==b['name'], index
    assert a['smiles'].strip()==b['smiles'], index # documented import whitespace trim
    assert b['status']==r['status'], index
    assert b['formula']==r.get('formula',''), index
    for key in FIELDS:
        if r['status']=='ok':
            assert math.isfinite(float(b[key])) and float(b[key])==r[key], (index,key)
        else:assert b[key]=='', (index,key)
summary=dict(rows=10000,all_12_numeric_columns_roundtrip=True,order_and_names_and_smiles_preserved=True,
    invalid_rows_retained=True,optimized_outputs_identical_to_baseline=True,
    input_sha256=hashlib.sha256((p/'bulk-input.csv').read_bytes()).hexdigest(),
    output_sha256=hashlib.sha256((p/'bulk-output.csv').read_bytes()).hexdigest())
(p/'bulk-io-independent-check.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
