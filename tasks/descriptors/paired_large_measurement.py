"""Alternating fresh C processes; one declared molecule, all291, no caches.

Run after owned compilers/audits finish. Other host activity is not controlled.
Every serialized result byte must agree before any speedup is reported.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    for key in ['before', 'after', 'driver', 'input', 'output']:
        p.add_argument('--'+key, type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    hashes = {k: sha(getattr(args, k)) for k in ['before','after','driver','input']}
    records = []
    reference = None
    for trial in range(3):
        for variant in (['before','after'] if trial % 2 == 0 else ['after','before']):
            output = args.output / f'{trial}-{variant}-actual.json'
            raw = subprocess.check_output([str(args.driver),str(getattr(args,variant)),str(args.input),str(output),'1'])
            result = json.loads(raw)
            actual = output.read_bytes()
            envelope = json.loads(actual)
            assert len(envelope['rows']) == 1
            row = envelope['rows'][0]
            assert row['status'] == 'ok' and len(row['descriptor_values']) == 291
            if reference is None:
                reference = actual
            assert actual == reference, 'Result bytes changed; reject timing claim'
            records.append(dict(trial=trial,variant=variant,**result))
            print(trial,variant,result['cold_seconds'],flush=True)
    assert hashes == {k: sha(getattr(args,k)) for k in hashes}, 'inputs changed'
    before = statistics.median(r['cold_seconds'] for r in records if r['variant']=='before')
    after = statistics.median(r['cold_seconds'] for r in records if r['variant']=='after')
    report = dict(scope='One original598-atom/632-bond molecule, macOS C ABI all291 plus JSON, fresh processes, loaded library, no cache/buffer retry; no physical iOS claim',
                  platform=platform.platform(),cpu=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),
                  owned_compilation_audit_overlap=False,other_host_activity_controlled=False,
                  sha256=hashes,records=records,all_six_result_files_byte_identical=True,
                  result_sha256=hashlib.sha256(reference).hexdigest(),
                  median_before_seconds=before,median_after_seconds=after,median_ratio=before/after,
                  before_max_process_rss_bytes=max(r['process_max_rss_bytes'] for r in records if r['variant']=='before'),
                  after_max_process_rss_bytes=max(r['process_max_rss_bytes'] for r in records if r['variant']=='after'))
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['records','sha256']}),flush=True)


if __name__ == '__main__':
    main()
