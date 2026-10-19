"""Unfiltered original-input before/after native I/O timing; exact-output gate.

Sum per-batch C ABI + JSON decoding durations. Checksumming is outside the
timer. This does not compare RDKit, Swift caching, UI, or a physical device.
"""
import argparse
import csv
import ctypes
import hashlib
import json
import os
import platform
import statistics
import time
from collections import Counter
from pathlib import Path

from compare_catalog import Native


class MeasuredNative(Native):
    def __init__(self, path):
        super().__init__(path)
        raw = self.fn
        self.address = ctypes.cast(raw, ctypes.c_void_p).value
        self.stats = {}

        def counted(source, length, output, capacity):
            result = raw(source, length, output, capacity)
            self.stats['calls'] += 1
            self.stats['short_buffer_calls'] += int(result > capacity)
            self.stats['max_required_bytes'] = max(self.stats['max_required_bytes'], result)
            return result
        self.fn = counted

    def run(self, inputs):
        self.stats = dict(calls=0, short_buffer_calls=0, max_required_bytes=0)
        elapsed = 0.0
        hashes = []
        statuses = Counter()
        for start in range(0, len(inputs), 64):
            batch = inputs[start:start+64]
            clock = time.perf_counter()
            rows = self.analyze(batch)
            elapsed += time.perf_counter() - clock
            for row in rows:
                row['index'] += start
                statuses[row['status']] += 1
                data = json.dumps(row, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
                hashes.append(hashlib.sha256(data).hexdigest())
        return elapsed, hashes, dict(statuses), dict(self.stats)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trials', type=int, default=3)
    args = parser.parse_args()
    records = list(csv.DictReader(args.input.open(newline='')))
    inputs = [r['smiles'] for r in records]
    assert inputs and len(inputs) == len(set(inputs)), 'retain distinct original inputs'
    assert args.trials >= 3
    engines = {name: MeasuredNative(path) for name, path in [('before', args.before), ('after', args.after)]}
    assert engines['before'].address != engines['after'].address, 'distinct loaded native entry points required'
    assert engines['before'].catalog == engines['after'].catalog
    report = dict(
        gate='pending', scope='native C ABI+JSON decode per batch; checksumming excluded',
        input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(), rows=len(inputs),
        distinct_raw_smiles=len(set(inputs)), native_filter=False, batch_size=64,
        platform=platform.platform(), python=platform.python_version(),
        load_average_before=os.getloadavg(), trials=args.trials, warmup_per_library=1,
        libraries={name: dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             revision=engines[name].revision, entry_address=engines[name].address)
                   for name, path in [('before', args.before), ('after', args.after)]},
        caveat='No result cache. Includes native calculations for the full catalog, serialization, FFI and decoding. Unsupported original rows retained. Not a RDKit, UI, Swift-cache or physical-device benchmark; background load is uncontrolled.')
    baselines = {}
    for name, engine in engines.items():
        _, hashes, statuses, stats = engine.run(inputs)
        baselines[name] = hashes
        print('warmup', name, statuses, stats, flush=True)
        report.setdefault('status_counts', {})[name] = statuses
    changed = [i for i, (a, b) in enumerate(zip(baselines['before'], baselines['after'])) if a != b]
    if changed:
        report.update(gate='failed', changed_rows=changed)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        raise SystemExit('exact output differs: no speed result')
    results = {name: [] for name in engines}
    for trial in range(args.trials):
        order = ['before', 'after'] if trial % 2 == 0 else ['after', 'before']
        for name in order:
            elapsed, hashes, statuses, stats = engines[name].run(inputs)
            assert hashes == baselines[name], (name, trial, 'output changed')
            assert statuses == report['status_counts'][name]
            results[name].append(dict(seconds=elapsed, **stats))
            print(trial, name, round(elapsed, 6), stats, flush=True)
    report.update(gate='passed', changed_rows=[], measurements=results,
                  row_checksum_sha256=hashlib.sha256(''.join(baselines['before']).encode()).hexdigest(),
                  median_seconds={name: statistics.median(v['seconds'] for v in rows) for name, rows in results.items()},
                  load_average_after=os.getloadavg())
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ['gate', 'rows', 'median_seconds']}), flush=True)


if __name__ == '__main__':
    main()
