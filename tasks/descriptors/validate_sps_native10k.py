"""Independent pinned SPS on original10k; all previous288 exact C ABI delta.

Does not claim current Swift/UI or full289 independent chemical parity.
"""
from pathlib import Path
import csv, gzip, hashlib, json, sys
import rdkit
from compare_catalog import Native, PARAMS, FUNCTIONS, Chem

assert rdkit.__version__ == '2026.03.6'
assert Chem.GetUseLegacyStereoPerception() is True
source, library, previous_native, output = map(Path, sys.argv[1:])
output.mkdir(parents=True,exist_ok=True)
inputs=list(csv.DictReader(source.open(newline='')))
old=json.loads(previous_native.read_text())
assert len(inputs)==len(old)==10_000 and len({r['smiles'] for r in inputs})==10_000
native=Native(library)
assert len(native.catalog)==289 and native.catalog[-1]['id']=='SPS'
all_rows=[]; differences=[]; changes=[]
matched=unsupported=invalid=undefined=mismatches=0
for start in range(0,len(inputs),64):
    batch=inputs[start:start+64]
    for offset,(record,row) in enumerate(zip(batch,native.analyze([r['smiles'] for r in batch]))):
        i=start+offset;row['index']=i;all_rows.append(row)
        assert old[i]['index']==i
        common=dict(row)
        if 'descriptor_values' in common:common['descriptor_values']=common['descriptor_values'][:288]
        if common!=old[i]:changes.append(i)
        mol=Chem.MolFromSmiles(record['smiles'],PARAMS)
        if mol is None or not mol.GetNumAtoms():invalid+=1;detail=dict(status='independent parser rejected')
        else:
            if mol.GetNumHeavyAtoms()==0:
                try:FUNCTIONS['SPS'](mol)
                except ZeroDivisionError:expected=None;undefined+=1
                else:raise AssertionError('unexpected SPS no-heavy behavior')
            else:expected=float(FUNCTIONS['SPS'](mol))
            if row['status']!='ok':unsupported+=1;detail=dict(status='native invalid or unsupported',expected=expected)
            else:
                observed=row['descriptor_values'][-1]
                same=(observed is None and expected is None) or (observed is not None and expected is not None and abs(observed-expected)<=1e-9)
                detail={} if same else dict(SPS=dict(native=observed,rdkit=expected))
                if same:matched+=1
                else:mismatches+=1
        if detail:differences.append(dict(index=i,**record,differences=detail))
    if start%2048==0:print(min(start+64,len(inputs)), '/',len(inputs),flush=True)
(output/'native-results.json').write_text(json.dumps(all_rows,separators=(',',':'),allow_nan=False))
with gzip.open(output/'differences.jsonl.gz','wt') as f:
    for row in differences:f.write(json.dumps(row)+'\n')
report=dict(gate='failed' if differences or changes else 'passed',value_gate='failed' if mismatches else 'passed',
    revision=native.revision,rows=len(inputs),catalog_count=289,standard_count=215,additional_count=50,
    matching_sps=matched,unsupported=unsupported,reference_invalid=invalid,named_domain_undefined=undefined,
    value_differing_rows=mismatches,previous288_and_base_changed_indices=changes,
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
    baseline_sha256=hashlib.sha256(previous_native.read_bytes()).hexdigest(),
    scope='C ABI only. Independent pinned SPS; previous288/base exact delta to v28 native. No current Swift/UI or all289 independent parity claim.')
(output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
if differences or changes:raise SystemExit(1)
