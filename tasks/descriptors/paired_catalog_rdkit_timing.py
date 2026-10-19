"""Fixed-input native/RDKit timing with the same catalog, formula and JSON I/O.

Reject any correctness/status difference before timing; never select by native
acceptance. This is a macOS host measurement, not a physical-device benchmark.
"""
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import sys
import time
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from compare_catalog import Native, PARAMS, FUNCTIONS, profile_reference, tolerance

source, library, destination = map(Path, sys.argv[1:])
assert rdkit.__version__ == '2026.03.6'
records = list(csv.DictReader(source.open(newline='')))
inputs = [r['smiles'] for r in records]
assert inputs and len(inputs) == len(set(inputs))
native = Native(library)
catalog = native.catalog
indices = {s['id']: i for i, s in enumerate(catalog)}
base_ids = dict(mw='MolWt', exact_mass='ExactMolWt', heavy_atoms='HeavyAtomCount',
                bonds='ds_BondCount', formal_charge='ds_FormalCharge',
                fragments='ds_FragmentCount', hbd='NumHDonors', hba='NumHAcceptors',
                tpsa='TPSA', logp='MolLogP', rotatable_bonds='NumRotatableBonds',
                fraction_csp3='FractionCSP3')


def run_native():
    rows = []
    for start in range(0, len(inputs), 64):
        chunk = native.analyze(inputs[start:start+64])
        for row in chunk:
            row['index'] += start
        rows.extend(chunk)
    return rows


def run_rdkit():
    rows = []
    for index, source in enumerate(inputs):
        mol = Chem.MolFromSmiles(source, PARAMS)
        assert mol is not None and mol.GetNumAtoms(), source
        profile = profile_reference(mol, include_extra=any(s["origin"]=="native-extra" for s in catalog))
        values = []
        for spec in catalog:
            key = spec['id']
            if spec['origin'] == 'rdkit-standard':
                if key == 'SPS' and mol.GetNumHeavyAtoms() == 0:
                    try:
                        FUNCTIONS[key](mol)
                    except ZeroDivisionError:
                        value = math.nan
                    else:
                        raise AssertionError('expected SPS no-heavy domain error')
                else:
                    value = FUNCTIONS[key](mol, force=False) if 'EStateIndex' in key else FUNCTIONS[key](mol)
            else:
                value = profile[key]
            value = float(value)
            values.append(value if math.isfinite(value) else None)
        row = dict(index=index, status='ok', formula=rdMolDescriptors.CalcMolFormula(mol),
                   descriptor_values=values)
        row.update({key: values[indices[id]] for key, id in base_ids.items()})
        rows.append(row)
    # Both paths include numeric/formula JSON encoding and decoding. Compact
    # Python JSON still represents integer-valued doubles differently from Rust.
    return json.loads(json.dumps(dict(schema=1, rows=rows), separators=(',', ':'), allow_nan=False))['rows']


baseline = run_native()
reference = run_rdkit()
assert len(baseline) == len(reference) == len(inputs)
for record, actual, expected in zip(records, baseline, reference):
    assert actual['status'] == expected['status'] == 'ok', record
    assert actual['formula'] == expected['formula'], record
    assert set(actual) == set(expected), record
    for key in base_ids:
        a, b = actual[key], expected[key]
        assert a == b or (a is not None and b is not None and abs(a-b) <= tolerance(base_ids[key])), (record,key,a,b)
    for spec, a, b in zip(catalog, actual['descriptor_values'], expected['descriptor_values']):
        assert a == b or (a is not None and b is not None and abs(a-b) <= tolerance(spec['id'])), (record,spec['id'],a,b)

functions = dict(native=run_native, rdkit=run_rdkit)
times = {name: [] for name in functions}
start_load = os.getloadavg()
for trial in range(5):
    order = ['native', 'rdkit'] if trial % 2 == 0 else ['rdkit', 'native']
    for name in order:
        start = time.perf_counter()
        rows = functions[name]()
        elapsed = time.perf_counter() - start
        assert rows == (baseline if name == 'native' else reference)
        times[name].append(elapsed)
        print(trial, name, round(elapsed, 6), flush=True)
report = dict(platform=platform.platform(), python=platform.python_version(), rdkit=rdkit.__version__,
              revision=native.revision, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
              rows=len(inputs), distinct_source_smiles=len(set(inputs)), catalog_count=len(catalog),
              standard_count=sum(s['origin']=='rdkit-standard' for s in catalog),
              profile_count=sum(s['origin']=='native-profile' for s in catalog),
              additional_count=sum(s['origin']=='native-extra' for s in catalog),
              warmup_per_engine=1, trials=5, batch_size=64, seconds=times,
              median_seconds={k:statistics.median(v) for k,v in times.items()},
              correctness_gate='all input rows/status/formula/catalog fields/base fields passed; no native filtering',
              load_average_before=start_load, load_average_after=os.getloadavg(),
              workload='Fresh parse per molecule/trial; same catalog+formula+base outputs; JSON encoding/decoding; no molecule-result cache. Reference is RDKit standard fields + Python24 profiles + independent NumPy extra50 formulas when present. Native C ABI chunks64; reference Python loop. In-process, caller sequential.',
              caveat='Fixed chemical regression controls, not representative10k throughput, UI rendering or physical iPhone speed. With extra50 this compares complete equivalent-output workflows, not RDKit built-in kernel speed. Host background load remains uncontrolled.')
destination.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)
