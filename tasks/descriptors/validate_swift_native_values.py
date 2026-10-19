"""Check the actual Swift decoder/row mapping against two frozen C ABI builds.

Original inputs, order and unsupported rows are retained. Compare binary64
values exactly: the descriptor contract is f64, including large JSON integers.
This is an I/O regression check, not an independent chemical reference or timing.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from compare_catalog import Native


def equal_value(a, b):
    if a is None or b is None:
        return a is None and b is None
    return float(a) == float(b)


def main():
    p = argparse.ArgumentParser()
    for name in ['input', 'actual', 'previous', 'current', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--expected-newly-accepted', type=Path,
                   help='Explicit JSON list of source indices whose status may change from unsupported to ok')
    args = p.parse_args()
    inputs = list(csv.DictReader(args.input.open(newline='')))
    actual = json.loads(args.actual.read_text())
    assert len(inputs) == len(actual) == 10_000
    assert len({r['smiles'] for r in inputs}) == len(inputs)
    assert [r['index'] for r in actual] == list(range(len(inputs)))
    previous, current = Native(args.previous), Native(args.current)
    previous_count = len(previous.catalog)
    assert 0 < previous_count <= len(current.catalog)
    assert [s['id'] for s in previous.catalog] == [s['id'] for s in current.catalog[:previous_count]]
    expected_new = json.loads(args.expected_newly_accepted.read_text()) if args.expected_newly_accepted else []
    assert expected_new == sorted(set(expected_new)) and all(0 <= i < len(inputs) for i in expected_new)
    new_indices = []
    report = dict(gate='pending', rows=len(inputs), distinct_raw_smiles=len(inputs),
                  source_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  actual_swift_rows_sha256=hashlib.sha256(args.actual.read_bytes()).hexdigest(),
                  scope='exact binary64 values, formula, status and source positions; no tolerance; no speed or independent-chemistry claim',
                  libraries={name: dict(revision=engine.revision, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                             for name, engine, path in [('previous', previous, args.previous), ('current', current, args.current)]})
    numeric_cells = 0
    unsupported = 0
    for start in range(0, len(inputs), 64):
        sources = [r['smiles'] for r in inputs[start:start+64]]
        before, now = previous.analyze(sources), current.analyze(sources)
        for offset, (old, native, swift) in enumerate(zip(before, now, actual[start:start+64])):
            i = start + offset
            assert old['index'] == native['index'] == offset and swift['index'] == i
            assert native['status'] == swift['status'], i
            if old['status'] != native['status']:
                assert i in expected_new and old['status'] == 'invalid_or_unsupported' and native['status'] == 'ok', i
                new_indices.append(i)
            if native['status'] != 'ok':
                unsupported += 1
                continue
            assert native['formula'] == swift['formula'], i
            if old['status'] == 'ok':
                assert old['formula'] == native['formula'], i
            assert len(swift['descriptor_values']) == len(current.catalog)
            for j, value in enumerate(native['descriptor_values']):
                assert equal_value(value, swift['descriptor_values'][j]), (i, current.catalog[j]['id'], 'Swift decode')
                numeric_cells += value is not None
                if old['status'] == 'ok' and j < previous_count:
                    assert equal_value(old['descriptor_values'][j], value), (i, current.catalog[j]['id'], 'prior value changed')
            for key in native.keys() - {'index', 'descriptor_values', 'formula', 'status'}:
                assert equal_value(native[key], swift[key]), (i, key)
                if old['status'] == 'ok':
                    assert equal_value(old[key], native[key]), (i, key)
        if start % 2048 == 0:
            print(min(start+64,len(inputs)), '/', len(inputs), flush=True)
    assert new_indices == expected_new, (new_indices, expected_new)
    report.update(gate='passed', unsupported_rows=unsupported,
                  every_swift_descriptor_matches_current_native=True,
                  previous_descriptor_count=previous_count,
                  previous_values_unchanged_on_common_accepted_rows=True,
                  newly_accepted_indices=new_indices,
                  finite_numeric_cells_checked=numeric_cells)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
