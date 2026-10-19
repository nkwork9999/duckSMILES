"""Equivalent 291-output workflow timing; not a pure RDKit kernel benchmark.

The reference workflow uses RDKit215, Python24 profiles, independent NumPy50
formulas, and FLINT2 stable information definitions. Original inputs are fixed
before either engine runs. Correctness is a prerequisite, never a row filter.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import time
from pathlib import Path

import compare_catalog as c
import flint
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    for name in ['input', 'library', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--trials', type=int, default=3)
    args = parser.parse_args()
    assert args.trials >= 3 and c.rdkit.__version__ == '2026.03.6'
    records = list(csv.DictReader(args.input.open(newline='')))
    inputs = [r['smiles'] for r in records]
    assert inputs and len(inputs) == len(set(inputs))
    engine = c.Native(args.library)
    assert len(engine.catalog) == 291
    catalog = engine.catalog
    index = {s['id']: i for i,s in enumerate(catalog)}
    base_ids = dict(mw='MolWt', exact_mass='ExactMolWt', logp='MolLogP', tpsa='TPSA',
                    hbd='NumHDonors', hba='NumHAcceptors', heavy_atoms='HeavyAtomCount',
                    bonds='ds_BondCount', formal_charge='ds_FormalCharge', fragments='ds_FragmentCount',
                    rotatable_bonds='NumRotatableBonds', fraction_csp3='FractionCSP3')

    def native():
        rows = []
        for start in range(0,len(inputs),64):
            batch = engine.analyze(inputs[start:start+64])
            for row in batch: row['index'] += start
            rows.extend(batch)
        return rows

    def reference():
        rows = []
        for start in range(0,len(inputs),64):
            batch = []
            for offset,source in enumerate(inputs[start:start+64]):
                mol = c.Chem.MolFromSmiles(source,c.PARAMS)
                assert mol is not None and mol.GetNumAtoms(),start+offset
                profile, extra = c.profile_reference(mol), c.extra_reference(mol)
                values = []
                for spec in catalog:
                    if spec['id'] == 'SPS' and not mol.GetNumHeavyAtoms():
                        try: c.FUNCTIONS['SPS'](mol)
                        except ZeroDivisionError: value = math.nan
                        else: raise AssertionError('SPS no-heavy domain changed')
                    else:
                        value = c.reference_value(spec,mol,profile,extra,information_backend='flint')
                    values.append(value if math.isfinite(value) else None)
                row = dict(index=start+offset,status='ok',formula=c.rdMolDescriptors.CalcMolFormula(mol),descriptor_values=values)
                row.update({key:values[index[id]] for key,id in base_ids.items()})
                batch.append(row)
            rows.extend(json.loads(json.dumps(batch,separators=(',',':'),allow_nan=False)))
        return rows

    functions = dict(native=native,reference_workflow=reference)
    baselines = {name:fn() for name,fn in functions.items()}
    assert len(baselines['native']) == len(baselines['reference_workflow']) == len(inputs)
    for i,(actual,expected) in enumerate(zip(baselines['native'],baselines['reference_workflow'])):
        assert actual['status'] == expected['status'] == 'ok',i
        assert actual['index'] == expected['index'] == i
        assert actual['formula'] == expected['formula'] and set(actual) == set(expected),i
        for spec,a,b in zip(catalog,actual['descriptor_values'],expected['descriptor_values']):
            assert (a is None and b is None) or (a is not None and b is not None and abs(a-b) <= c.allowed_error(spec['id'],b)),(i,spec['id'],a,b)
        for key,id in base_ids.items():
            a,b = actual[key],expected[key]
            assert a == b or (a is not None and b is not None and abs(a-b) <= c.allowed_error(id,b)),(i,key)
    measured = {name:[] for name in functions}
    load_before = os.getloadavg()
    for trial in range(args.trials):
        order = list(functions) if trial%2 == 0 else list(reversed(functions))
        for name in order:
            begin = time.perf_counter(); rows = functions[name](); elapsed = time.perf_counter()-begin
            assert rows == baselines[name],(trial,name,'nondeterministic result')
            measured[name].append(elapsed)
            print(trial,name,elapsed,flush=True)
    report = dict(gate='passed',platform=platform.platform(),python=platform.python_version(),rdkit=c.rdkit.__version__,flint=flint.__version__,numpy=np.__version__,
                  revision=engine.revision,rows=len(inputs),distinct_original_smiles=len(set(inputs)),native_filter=False,
                  source_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),catalog_count=len(catalog),
                  reference_components=dict(pinned_rdkit=215,stable_exact_information_flint=2,python_profiles=24,numpy_extra=50),
                  information_compatible_with_pinned_rdkit_numerics=False,
                  correctness='Every original input, status, formula, all291 and legacy base fields before timing; no exclusions',
                  batch_size=64,warmup_per_workflow=1,trials=args.trials,seconds=measured,
                  median_seconds={name:statistics.median(v) for name,v in measured.items()},
                  load_average_before=load_before,load_average_after=os.getloadavg(),
                  tool_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(c.__file__),Path(c.__file__).with_name('extra50_reference.py'),Path(c.__file__).with_name('exact_ipc_reference.py'),Path(c.__file__).with_name('flint_ipc_reference.py')]},
                  scope='Fresh parse and all291/formula/base output plus JSON encode/decode, no molecule-result cache. Native Rust C ABI vs RDKit/Python/NumPy/FLINT workflow, not pure RDKit kernel speed. Fixed source prefix, not random population or physical iOS/UI measurement. Other host activity uncontrolled.')
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    main()
