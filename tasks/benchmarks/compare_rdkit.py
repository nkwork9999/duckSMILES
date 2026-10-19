#!/usr/bin/env python3
"""Same-input single-threaded RDKit/native C ABI comparison; no dedup cache.
Run with the pinned RDKit env, --library <host dylib> --output <evidence dir>.
"""
import argparse, csv, ctypes, hashlib, io, json, math, platform, statistics, time
from pathlib import Path
import rdkit
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import rdMolDescriptors as D, Descriptors, Crippen

FIELDS = ['mw','exact_mass','heavy_atoms','bonds','formal_charge','fragments','hbd','hba','tpsa','logp','rotatable_bonds','fraction_csp3']
TOLERANCE = {k: 0 for k in FIELDS}
TOLERANCE.update(mw=0.005, exact_mass=1e-5, tpsa=1e-6, logp=1e-6, fraction_csp3=1e-9)
PARAMS = Chem.SmilesParserParams(); PARAMS.removeHs = False
RDLogger.DisableLog('rdApp.*')

def reference(s, index):
    mol = Chem.MolFromSmiles(s, PARAMS) if s.strip() else None
    if mol is None or mol.GetNumAtoms() == 0:
        return {'index':index, 'status':'invalid_or_unsupported'}
    return dict(index=index,status='ok',formula=D.CalcMolFormula(mol),mw=Descriptors.MolWt(mol),
        exact_mass=D.CalcExactMolWt(mol),heavy_atoms=mol.GetNumHeavyAtoms(),bonds=mol.GetNumBonds(),
        formal_charge=Chem.GetFormalCharge(mol),fragments=len(Chem.GetMolFrags(mol)),
        hbd=D.CalcNumHBD(mol),hba=D.CalcNumHBA(mol),tpsa=D.CalcTPSA(mol,includeSandP=False),
        logp=Crippen.MolLogP(mol),rotatable_bonds=D.CalcNumRotatableBonds(mol,D.NumRotatableBondsOptions.Strict),
        fraction_csp3=D.CalcFractionCSP3(mol))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--library',required=True);ap.add_argument('--output',required=True);ap.add_argument('--trials',type=int,default=3);ap.add_argument('--strict',action='store_true')
    args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    source=Path(RDConfig.RDDataDir)/'NCI'/'first_5K.smi'
    records=[line.split()[0] for line in source.read_text().splitlines() if line.strip()]
    lib=ctypes.CDLL(args.library);fn=lib.ds_workbench_analyze
    fn.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_void_p,ctypes.c_size_t];fn.restype=ctypes.c_int32
    buffer=ctypes.create_string_buffer(16*1024*1024)
    def native_bytes(items):
        raw='\n'.join(items).encode();n=fn(raw,len(raw),buffer,len(buffer))
        assert 0<n<=len(buffer), (n,len(items))
        return ctypes.string_at(buffer,n)
    def native(items):
        result=json.loads(native_bytes(items));assert result['schema']==1
        rows=result['rows']; assert len(rows)==len(items)
        assert [r['index'] for r in rows]==list(range(len(items)))
        for r in rows:
            assert r['status'] in ['ok','invalid_or_unsupported']
            if r['status']=='ok':
                assert set(r)==set(FIELDS+['index','status','formula'])
                assert all(isinstance(r[k],(int,float)) and math.isfinite(r[k]) for k in FIELDS)
            else: assert set(r)=={'index','status'}
        return rows
    edge=['','invalid','C(','C=','C..O','[Na+].[Cl-]','[NH4+]','[13CH4]','[H]O[H]','[O-]','N[C@@H](C)C(=O)O','炭素','C'*257,'C1CC','C1CCCCC1','CC(=O)N','O=S(=O)(O)O','[Fe+2]','C\\C=C/C','C=1CCCCC1']
    # Blank sentinels mirror the Swift encoder; exact source strings remain in reports.
    fixtures=records+edge
    actual=[];t=time.perf_counter()
    for start in range(0,len(fixtures),64):
        actual.extend(native([s if s else ' ' for s in fixtures[start:start+64]]))
        if start%1024==0: print(f'native correctness {start}/{len(fixtures)}',flush=True)
    print(f'native correctness done {time.perf_counter()-t:.2f}s',flush=True)
    expected=[reference(s,i) for i,s in enumerate(fixtures)]
    mismatches=[]; common=[];status_counts={};field_mismatches={k:0 for k in FIELDS+['formula']};max_abs={k:0 for k in FIELDS}
    full_matches=0
    for i,(s,a,b) in enumerate(zip(fixtures,actual,expected)):
        pair=f"native_{a['status']}__rdkit_{b['status']}";status_counts[pair]=status_counts.get(pair,0)+1
        if a['status']!=b['status']:
            mismatches.append(dict(source_index=i,smiles=s,difference='status',native=a,rdkit=b));continue
        if a['status']!='ok':continue
        if i<len(records):common.append(s)
        differences={}
        for k in FIELDS:
            error=abs(a[k]-b[k]);max_abs[k]=max(max_abs[k],error)
            if error>TOLERANCE[k]:differences[k]={'native':a[k],'rdkit':b[k],'abs_error':error};field_mismatches[k]+=1
        if a['formula']!=b['formula']:
            differences['formula']={'native':a['formula'],'rdkit':b['formula']};field_mismatches['formula']+=1
        if differences:mismatches.append(dict(source_index=i,smiles=s,difference='values',fields=differences))
        else:full_matches+=1
    (out/'rdkit-differences.json').write_text(json.dumps(mismatches,ensure_ascii=False,indent=2)+'\n')
    # Large ABI contract: the app's 10,000-row limit, mixed valid/invalid rows, row order.
    bulk=(records*2)+['invalid',' '];assert len(bulk)==10000
    bulk_results=native(bulk)
    assert bulk_results[-2:]==[{'index':9998,'status':'invalid_or_unsupported'},{'index':9999,'status':'invalid_or_unsupported'}]
    # Bounds must be rejected before any partial output. Empty batch is zero rows.
    assert json.loads(native_bytes([]))['rows']==[]
    for raw,code in [(b'CCO\n'*10000+b'CCO',-2),(b'x'*1048577,-1),(b'\xff',-1)]:
        sentinel=ctypes.create_string_buffer(b'KEEP',4);assert fn(raw,len(raw),sentinel,4)==code;assert sentinel.raw==b'KEEP'
    small=ctypes.create_string_buffer(b'KEEP',4);assert fn(b'CCO',3,small,4)>4;assert small.raw==b'KEEP'
    assert fn(None,1,None,0)==-1
    print('ABI 10,000-row and bounds checks passed',flush=True)
    # Build CSV fixture used by the actual Swift Workbench import/export harness.
    csv_buffer=io.StringIO(newline='');w=csv.writer(csv_buffer);w.writerow(['name','smiles'])
    for i,s in enumerate(bulk):
        name=f'NCI-{i}' if i%101 else f'名前,{i} "quoted"\nsecond line'
        w.writerow([name,s])
    # Replace the two intentional failures with multiline/empty fields.
    fixture_text=csv_buffer.getvalue();(out/'bulk-input.csv').write_text(fixture_text)
    # Alternate order, fresh Mol each row/trial, no cached duplicate molecules.
    timed_inputs=common[:1000];assert len(timed_inputs)==1000
    def native_timed():
        return [json.loads(native_bytes(timed_inputs[i:i+64])) for i in range(0,len(timed_inputs),64)]
    def rdkit_timed():
        return [json.loads(json.dumps({'schema':1,'rows':[reference(s,j) for j,s in enumerate(timed_inputs[i:i+64])]},separators=(',',':'),ensure_ascii=False)) for i in range(0,len(timed_inputs),64)]
    native_timed();rdkit_timed();times={'native_ms':[],'rdkit_ms':[]}
    for trial in range(args.trials):
        order=[('native_ms',native_timed),('rdkit_ms',rdkit_timed)]
        if trial%2:order.reverse()
        for key,run in order:
            t=time.perf_counter();result=run();elapsed=(time.perf_counter()-t)*1000
            assert sum(len(b['rows']) for b in result)==1000
            times[key].append(elapsed);print(f'trial {trial} {key}={elapsed:.2f}',flush=True)
    report=dict(rdkit_compatibility_gate='failed' if mismatches else 'passed',rdkit_version=rdkit.__version__,platform=platform.platform(),arch=platform.machine(),dataset=dict(source='RDKit Data/NCI/first_5K.smi',sha256=hashlib.sha256(source.read_bytes()).hexdigest(),rows=len(records),unique_smiles=len(set(records))),
        fixture_rows=len(fixtures),edge_rows=len(edge),status_counts=status_counts,fully_matching_valid_rows=full_matches,field_mismatches=field_mismatches,max_abs_error=max_abs,tolerance=TOLERANCE,
        abi_checks=dict(bulk_rows=10000,preserved_order=True,invalid_rows_retained=True,rejects_10001_rows=True,rejects_over_1MiB=True,rejects_non_utf8=True,short_buffer_untouched=True),
        timing=dict(rows=1000,distinct_smiles=len(set(timed_inputs)),trials=args.trials,batch_size=64,raw=times,median_ms={k:(statistics.median(v) if v else None) for k,v in times.items()},includes='parse + 12 descriptors + formula + JSON encode/decode; in-process, single thread, no molecule cache',common_valid_selection='first 1000 accepted by both; correctness report includes every source row',rdkit_removeHs=False,tpsa_includeSandP=False,rotatable_bonds='Strict'))
    (out/'rdkit-comparison.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False),flush=True)
    if args.strict and mismatches: raise SystemExit(1)
if __name__=='__main__':main()
