"""Independent RDKit check of selected descriptors in actual Swift result rows.

No native recomputation or native-acceptance filtering. Undefined values follow
the same pinned policy as the full comparator. Source order is checked first.
"""
import csv,gzip,hashlib,json,math,sys
from pathlib import Path
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors
from compare_catalog import tolerance
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')
source,results,catalog,definitions,output=map(Path,sys.argv[1:]);output.mkdir(parents=True,exist_ok=True)
inputs=list(csv.DictReader(source.open(newline='')));actual=json.loads(results.read_text())
specs=json.loads(catalog.read_text())['descriptors'];selected=json.loads(definitions.read_text())['descriptors']
indices={s['id']:i for i,s in enumerate(specs)};assert len(indices)==len(specs)
assert len(inputs)==len(actual) and [r['index'] for r in actual]==list(range(len(inputs)))
assert all(s['origin']=='rdkit-standard' and s['id'] in indices for s in selected)
functions=dict(Descriptors._descList);params=Chem.SmilesParserParams();params.removeHs=False
metrics={s['id']:dict(compared=0,mismatches=0,matching_undefined=0,reference_errors=0,max_abs_error=0,tolerance=tolerance(s['id'])) for s in selected}
unsupported=reference_invalid=matching=0;diffrows=[]
for i,(record,row) in enumerate(zip(inputs,actual)):
    mol=Chem.MolFromSmiles(record['smiles'],params)
    if mol is None or not mol.GetNumAtoms():
        reference_invalid+=1;diffrows.append(dict(index=i,**record,differences={'status':'reference invalid'}));continue
    if row['status']!='ok':
        unsupported+=1;diffrows.append(dict(index=i,**record,differences={'status':'native invalid or unsupported'}));continue
    values=row['descriptor_values'];assert len(values)==len(specs)
    assert all(v is None or math.isfinite(v) for v in values)
    differences={}
    for spec in selected:
        id=spec['id'];metric=metrics[id];value=values[indices[id]]
        try: expected=float(functions[id](mol))
        except Exception as error:
            metric['reference_errors']+=1;differences[id]={'reference_error':str(error)};continue
        metric['compared']+=1
        if not math.isfinite(expected):
            if value is None:metric['matching_undefined']+=1;continue
            detail=dict(native=value,rdkit=None,reference_nonfinite=True)
        elif value is None:detail=dict(native=None,rdkit=expected)
        else:
            error=abs(value-expected);metric['max_abs_error']=max(metric['max_abs_error'],error)
            if error<=metric['tolerance']:continue
            detail=dict(native=value,rdkit=expected,abs_error=error)
        metric['mismatches']+=1;differences[id]=detail
    if differences:diffrows.append(dict(index=i,**record,differences=differences))
    else:matching+=1
with gzip.open(output/'differences.jsonl.gz','wt',encoding='utf-8') as f:
    for d in diffrows:f.write(json.dumps(d,ensure_ascii=False)+'\n')
value_differences=sum(any(k!='status' for k in d['differences']) for d in diffrows)
report=dict(gate='failed' if diffrows else 'passed',value_gate='failed' if value_differences else 'passed',
    rdkit_version=rdkit.__version__,rows=len(inputs),distinct_raw_smiles=len({r['smiles'] for r in inputs}),
    native_catalog_count=len(specs),compared_descriptor_count=len(selected),descriptor_ids=[s['id'] for s in selected],
    common_valid=len(inputs)-unsupported-reference_invalid,all_selected_values_matching_rows=matching,
    native_unsupported=unsupported,reference_invalid=reference_invalid,value_differing_rows=value_differences,
    differing_rows=len(diffrows),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    actual_swift_results_sha256=hashlib.sha256(results.read_bytes()).hexdigest(),metrics=metrics,
    scope='Only the listed descriptor IDs in actual Swift results; not an all-catalog or performance claim.',
    reference_policy=dict(removeHs=False,sanitize=True,standardize=False,undefined='reference nonfinite + native null matches; exceptions fail'))
(output/'comparison.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='metrics'},ensure_ascii=False),flush=True)
if diffrows:raise SystemExit(1)
