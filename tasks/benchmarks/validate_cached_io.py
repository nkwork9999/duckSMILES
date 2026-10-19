#!/usr/bin/env python3
"""Compare cached Swift output against uncached C ABI output for every source row."""
import argparse,csv,ctypes,hashlib,json,math
from pathlib import Path
from compare_rdkit import FIELDS
ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('library');args=ap.parse_args()
p=Path(args.directory)
with (p/'bulk-input.csv').open(newline='') as f: source=list(csv.DictReader(f))
with (p/'bulk-output.csv').open(newline='') as f: exported=list(csv.DictReader(f))
actual=json.loads((p/'bulk-native-results.json').read_text())
lib=ctypes.CDLL(args.library);fn=lib.ds_workbench_analyze
fn.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_void_p,ctypes.c_size_t];fn.restype=ctypes.c_int32
raw='\n'.join(row['smiles'].strip() or ' ' for row in source).encode()
buffer=ctypes.create_string_buffer(16*1024*1024);n=fn(raw,len(raw),buffer,len(buffer));assert 0<n<len(buffer)
fresh=json.loads(ctypes.string_at(buffer,n))['rows']
assert len(source)==len(exported)==len(actual)==len(fresh)==10000
assert actual==fresh, 'cached output differs from uncached native calculation'
for index,(a,b,r) in enumerate(zip(source,exported,actual)):
    assert r['index']==index and a['name']==b['name'] and a['smiles'].strip()==b['smiles']
    assert b['status']==r['status'] and b['formula']==r.get('formula','')
    for key in FIELDS:
        if r['status']=='ok':assert math.isfinite(float(b[key])) and float(b[key])==r[key],(index,key)
        else:assert b[key]==''
result=dict(rows=10000,cached_equals_uncached=True,all_numeric_columns_roundtrip=True,order_names_and_smiles_preserved=True,invalid_rows_retained=True,output_sha256=hashlib.sha256((p/'bulk-output.csv').read_bytes()).hexdigest())
(p/'bulk-io-independent-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
