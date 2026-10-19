#!/usr/bin/env python3
"""After builds finish, rotate before/after/RDKit ordering on the exact same inputs."""
import argparse,ctypes,json,statistics,time,os,platform,hashlib
from pathlib import Path
from rdkit import RDConfig
from compare_rdkit import reference
ap=argparse.ArgumentParser();ap.add_argument('baseline');ap.add_argument('optimized');ap.add_argument('output');args=ap.parse_args()
records=[s.split()[0] for s in (Path(RDConfig.RDDataDir)/'NCI'/'first_5K.smi').read_text().splitlines() if s.strip()]
def backend(path):
    lib=ctypes.CDLL(path);f=lib.ds_workbench_analyze;f.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_void_p,ctypes.c_size_t];f.restype=ctypes.c_int32
    buffer=ctypes.create_string_buffer(64*1024)
    def call(chunk):
        b='\n'.join(chunk).encode();n=f(b,len(b),buffer,len(buffer));assert 0<n<=len(buffer)
        return json.loads(ctypes.string_at(buffer,n))['rows']
    return call
before=backend(args.baseline);after=backend(args.optimized)
inputs=[]
for s in records:
    a=after([s])[0];b=reference(s,0)
    if a['status']=='ok' and b['status']=='ok':inputs.append(s)
    if len(inputs)==1000:break
assert len(inputs)==1000
chunks=[inputs[i:i+64] for i in range(0,1000,64)]
def native_run(call):return [call(chunk) for chunk in chunks]
def rdkit_run():
    return [json.loads(json.dumps([reference(s,j) for j,s in enumerate(chunk)],separators=(',',':'))) for chunk in chunks]
runners={'before_ms':lambda:native_run(before),'after_ms':lambda:native_run(after),'rdkit_ms':rdkit_run}
# Warm every path without retaining parsed molecules between rows or trials.
for run in runners.values():run()
times={key:[] for key in runners}
for trial in range(3):
    keys=list(runners);keys=keys[trial:]+keys[:trial]
    for key in keys:
        t=time.perf_counter();rows=runners[key]();elapsed=(time.perf_counter()-t)*1000
        assert sum(map(len,rows))==1000
        times[key].append(elapsed);print(trial,key,round(elapsed,2),flush=True)
medians={key:statistics.median(v) for key,v in times.items()}
result=dict(platform=platform.platform(),rows=1000,unique_smiles=len(set(inputs)),batch_size=64,trials=3,
    includes='SMILES parsing + 12 descriptors + formula + JSON serialization/deserialization',
    execution='single-threaded sequential paths; rotating order; fresh molecules; no duplicate cache',
    raw_ms=times,median_ms=medians,optimization_speedup=medians['before_ms']/medians['after_ms'],
    after_over_rdkit=medians['after_ms']/medians['rdkit_ms'],load_average=os.getloadavg(),
    input_sha256=hashlib.sha256('\n'.join(inputs).encode()).hexdigest())
Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
