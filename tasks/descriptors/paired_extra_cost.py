"""Measure eager-extra50 cost on fixed controls, including C ABI JSON I/O.

Old236 versus new286 is a feature-overhead experiment, not an equal-workload
RDKit comparison or physical-device measurement. Validate every old field first.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import platform
import statistics
import sys
import time
from compare_catalog import Native

source, old_library, new_library, output = map(Path, sys.argv[1:])
inputs = [r['smiles'] for r in csv.DictReader(source.open(newline=''))]
assert len(inputs) == len(set(inputs)) == 413
engines = dict(old236=Native(old_library), new286=Native(new_library))
assert len(engines['old236'].catalog) == 236 and len(engines['new286'].catalog) == 286
assert [s['id'] for s in engines['old236'].catalog] == [s['id'] for s in engines['new286'].catalog[:236]]

def run(name):
    rows = []
    for start in range(0, len(inputs), 64):
        rows.extend(engines[name].analyze(inputs[start:start+64]))
    return rows

baseline = {name: run(name) for name in engines}
for old,new in zip(baseline['old236'],baseline['new286']):
    assert old['status'] == new['status'] == 'ok'
    common = dict(new); common['descriptor_values'] = common['descriptor_values'][:236]
    assert old == common
times = {name:[] for name in engines}
load = os.getloadavg()
for trial in range(5):
    for name in (['old236','new286'] if trial % 2 == 0 else ['new286','old236']):
        start = time.perf_counter(); rows = run(name); elapsed = time.perf_counter()-start
        assert rows == baseline[name]
        times[name].append(elapsed)
medians = {name:statistics.median(values) for name,values in times.items()}
report = dict(platform=platform.platform(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              library_sha256={name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in [('old236',old_library),('new286',new_library)]},
              revisions={name:e.revision for name,e in engines.items()},rows=len(inputs),trials=5,warmup_each=1,
              seconds=times,median_seconds=medians,extra_overhead_ratio=medians['new286']/medians['old236']-1,
              previous236_and_base_exact_match=True,load_before=load,load_after=os.getloadavg(),
              scope='Fixed413regressioncontrols; fresh parse/calculate/JSON each run, no result cache. Different output counts measure eager50 feature overhead, not equal workloads. No RDKit/native/device/UI speed claim. Host background load uncontrolled.')
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
