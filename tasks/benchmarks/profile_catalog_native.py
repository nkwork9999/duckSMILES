"""One uncached pass for external macOS sampling, explicitly not a benchmark."""
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "descriptors"))
from compare_catalog import Native

source, library, destination = map(Path, sys.argv[1:])
records = list(csv.DictReader(source.open(newline="")))
inputs = [r["smiles"] for r in records]
assert len(inputs) == len(set(inputs)) == 10_000
native = Native(library)
print(json.dumps(dict(pid=os.getpid(), revision=native.revision,
                     scope="sampling only; not elapsed-time evidence")), flush=True)
statuses = Counter()
checksum = hashlib.sha256()
for start in range(0, len(inputs), 64):
    rows = native.analyze(inputs[start:start + 64])
    for row in rows:
        row["index"] += start
        statuses[row["status"]] += 1
        checksum.update(json.dumps(row, sort_keys=True, separators=(",", ":")).encode())
report = dict(scope="sampling only; not a benchmark", rows=len(inputs),
              native_acceptance_filter=False, status_counts=dict(statuses),
              input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
              revision=native.revision, output_checksum_sha256=checksum.hexdigest())
destination.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report), flush=True)
