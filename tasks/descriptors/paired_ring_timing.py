#!/usr/bin/env python3
"""Paired timing for ring reuse, with identical catalog-ordered outputs."""
import csv
import hashlib
import json
import platform
import statistics
import sys
import time
from pathlib import Path
from compare_catalog import Native

source, before_path, after_path, output = map(Path, sys.argv[1:])
with source.open(newline='') as f:
    inputs = [r['smiles'] for r in csv.DictReader(f)][:1000]
assert len(inputs) == len(set(inputs)) == 1000
engines = {'before_ms': Native(before_path), 'after_ms': Native(after_path)}
assert engines['before_ms'].catalog == engines['after_ms'].catalog


def run(engine):
    result = []
    for start in range(0, len(inputs), 64):
        result.extend(engine.analyze(inputs[start:start+64]))
    return result


baseline = run(engines['before_ms'])
assert run(engines['after_ms']) == baseline, 'refactor changed descriptor/status output'
times = {key: [] for key in engines}
for trial in range(3):
    order = list(engines)
    if trial % 2:
        order.reverse()
    for key in order:
        start = time.perf_counter()
        rows = run(engines[key])
        elapsed = (time.perf_counter() - start) * 1000
        assert rows == baseline
        times[key].append(elapsed)
        print(trial, key, elapsed, flush=True)
report = dict(platform=platform.platform(), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              before_library_sha256=hashlib.sha256(before_path.read_bytes()).hexdigest(),
              after_library_sha256=hashlib.sha256(after_path.read_bytes()).hexdigest(),
              rows=1000, distinct_source_smiles=1000, valid_rows=sum(r['status']=='ok' for r in baseline),
              descriptor_count=len(engines['before_ms'].catalog),
              batch_size=64, trials=3, output_equal=True, raw=times,
              median_ms={key:statistics.median(values) for key,values in times.items()},
              workload=f"fresh parse + same {len(engines['before_ms'].catalog)} descriptors + formula + C ABI JSON encode/decode, in process, one thread; no molecule-result cache",
              caveat='macOS host, not iPhone. Inspect load context in the report before interpreting small differences.')
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report), flush=True)
