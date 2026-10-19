"""Whole-source comparison in disjoint, ordered ranges. No native filtering.

Workers use the same frozen library/oracles/policies. This is correctness work,
not a timing comparison. Original source SHA and global positions are retained.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def partitions(size, workers):
    assert size > 0 and workers > 0
    workers = min(size, workers)
    return [(size*i//workers, size*(i+1)//workers) for i in range(workers)]


def merge_metric(values):
    policies = ['reference_source', 'tolerance', 'relative_tolerance']
    result = {}
    for key in policies:
        assert all(v[key] == values[0][key] for v in values), ('policy drift', key)
        result[key] = values[0][key]
    keys = set().union(*(v.keys() for v in values)) - set(policies)
    for key in keys:
        if key == 'examples':
            result[key] = sorted((r for v in values for r in v.get(key, [])), key=lambda r:r['index'])[:5]
        elif key.startswith('max_'):
            result[key] = max(v.get(key, 0) for v in values)
        else:
            assert key in {'compared','mismatches','undefined_native','undefined_reference','matching_undefined','reference_errors','named_domain_undefined'}, key
            result[key] = sum(v.get(key, 0) for v in values)
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    for name in ['library','input','output']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--information-backend', choices=['python', 'flint'], default='python')
    args = p.parse_args()
    assert 1 <= args.workers <= 4
    records = list(csv.DictReader(args.input.open(newline='')))
    assert records and all('name' in r and 'smiles' in r for r in records)
    assert len({r['smiles'] for r in records}) == len(records), 'retain distinct original source'
    parts = partitions(len(records), args.workers)
    args.output.mkdir(parents=True, exist_ok=False)
    tools = [Path(__file__).resolve(), Path(__file__).with_name('compare_catalog.py'), Path(__file__).with_name('extra50_reference.py'), Path(__file__).with_name('exact_ipc_reference.py')]
    hashes = {str(f):sha(f) for f in tools}
    flint_oracle = Path(__file__).with_name('flint_ipc_reference.py')
    tools.append(flint_oracle)
    hashes[str(flint_oracle)] = sha(flint_oracle)
    source_hash, library_hash = sha(args.input), sha(args.library)
    jobs = []
    for i, (start, end) in enumerate(parts):
        out = args.output/f'part-{i}'
        log = (args.output/f'part-{i}.log').open('w')
        command = [sys.executable,str(tools[1]),'--library',str(args.library),'--input',str(args.input),'--output',str(out),'--offset',str(start),'--limit',str(end-start),'--strict']
        command += ['--information-backend', args.information_backend]
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        jobs.append((process,log,out,start,end))
        print('started',i,'source positions',start,end,'PID',process.pid,flush=True)
    reports = []
    for process,log,out,start,end in jobs:
        code = process.wait()
        log.close()
        assert code in (0,1), (out,code)
        report = json.loads((out/'comparison.json').read_text())
        assert code == (0 if report['gate']=='passed' else 1)
        d = report['dataset']
        assert (d['index_offset'],d['index_stop'],d['rows'],d['source_rows']) == (start,end,end-start,len(records))
        assert d['sha256']==source_hash and report['library_sha256']==library_hash
        reports.append(report)
        print('finished',out.name,report['gate'],report['common_valid'],'accepted',flush=True)
    assert {str(f):sha(f) for f in tools} == hashes, 'reference tools changed during run'
    assert sha(args.input)==source_hash and sha(args.library)==library_hash
    merged = dict(reports[0])
    invariant = ['revision','rdkit_version','platform','library_sha256','catalog_count','registered_standard_count','direct_pinned_standard_count','exact_definition_standard_count','profile_count','additional_count','reference_policy','caveat']
    for key in invariant:
        assert all(r[key]==reports[0][key] for r in reports), key
    merged['dataset'] = dict(path=str(args.input),sha256=source_hash,rows=len(records),distinct_raw_smiles=len(records),source_rows=len(records),index_offset=0,index_stop=len(records))
    sums = ['common_valid','all_values_matching_rows','standard_and_formula_matching_rows','differing_rows','direct_pinned_standard_and_formula_matching_rows','formula_mismatches']
    for key in sums:
        merged[key] = sum(r[key] for r in reports)
    statuses = Counter()
    for r in reports: statuses.update(r['status_counts'])
    assert sum(statuses.values())==len(records)
    merged['status_counts'] = dict(statuses)
    ids = list(reports[0]['metrics'])
    assert all(list(r['metrics'])==ids for r in reports)
    merged['metrics'] = {id:merge_metric([r['metrics'][id] for r in reports]) for id in ids}
    pinned = dict(reports[0]['pinned_information_compatibility'])
    pinned['differing_rows'] = sum(r['pinned_information_compatibility']['differing_rows'] for r in reports)
    pinned['gate'] = 'failed' if pinned['differing_rows'] else 'passed'
    pinned['metrics'] = {id:merge_metric([r['pinned_information_compatibility']['metrics'][id] for r in reports]) for id in pinned['metrics']}
    merged['pinned_information_compatibility'] = pinned
    for filename, expected in [('differences.jsonl.gz',merged['differing_rows']), ('pinned-information-differences.jsonl.gz',pinned['differing_rows'])]:
        count = 0
        previous = -1
        with gzip.open(args.output/filename,'wt') as output:
            for _,_,out,start,end in jobs:
                with gzip.open(out/filename,'rt') as source:
                    for line in source:
                        row=json.loads(line);index=row['index']
                        assert start <= index < end and index > previous
                        assert (row['name'],row['smiles'])==(records[index]['name'],records[index]['smiles'])
                        previous=index;count+=1;output.write(line)
        assert count==expected,(filename,count,expected)
    merged['gate'] = 'failed' if merged['differing_rows'] else 'passed'
    merged['accepted_value_gate'] = 'failed' if merged['formula_mismatches'] or any(m['mismatches'] or m['reference_errors'] for m in merged['metrics'].values()) else 'passed'
    merged['execution'] = dict(workers=len(parts),ranges=parts,all_source_positions_covered_once=True,native_filter=False,tool_sha256=hashes,timing_claim=False)
    (args.output/'comparison.json').write_text(json.dumps(merged,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:merged[k] for k in ['gate','accepted_value_gate','dataset','catalog_count','common_valid','all_values_matching_rows','differing_rows']}),flush=True)
    if merged['gate']!='passed': raise SystemExit(1)


if __name__=='__main__':
    main()
