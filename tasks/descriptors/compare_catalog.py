#!/usr/bin/env python3
"""Compare every exported value against an explicitly declared reference.

This is a correctness harness, not a speed comparison: native/profile fields and
registered RDKit fields are reported separately. Input CSV is name,smiles; source
strings and order are preserved. It does not normalize molecules before analysis.
"""
import argparse
import csv
import ctypes
import gzip
import hashlib
import json
import math
import platform
from pathlib import Path

import rdkit
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from extra50_reference import reference as extra_reference
from exact_ipc_reference import reference as exact_information_reference

PARAMS = Chem.SmilesParserParams()
PARAMS.removeHs = False
RDLogger.DisableLog('rdApp.*')
FUNCTIONS = dict(Descriptors._descList)
STABLE_INFORMATION_IDS = frozenset(('Ipc', 'AvgIpc'))


class Native:
    def __init__(self, library):
        self.lib = ctypes.CDLL(str(library))
        self.fn = self.lib.ds_workbench_analyze
        self.fn.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t]
        self.fn.restype = ctypes.c_int32
        self.catalog = self.read_static('ds_workbench_catalog', json.loads)['descriptors']
        self.revision = self.read_static('ds_workbench_revision', lambda b: b.decode())
        assert len({s['id'] for s in self.catalog}) == len(self.catalog)

    def read_static(self, symbol, decode):
        fn = getattr(self.lib, symbol)
        fn.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        fn.restype = ctypes.c_int32
        n = fn(None, 0)
        assert 0 < n <= 1_048_576
        short = ctypes.create_string_buffer(b'KEEP', 4)
        assert fn(short, 4) == n and short.raw == b'KEEP'
        assert fn(None, 4) == -1
        out = ctypes.create_string_buffer(n)
        assert fn(out, n) == n
        return decode(out.raw)

    def analyze(self, items):
        data = '\n'.join(s if s else ' ' for s in items).encode()
        buf = ctypes.create_string_buffer(max(4096, len(items) * 4096))
        n = self.fn(data, len(data), buf, len(buf))
        assert 0 < n <= 64_000_000, (n, len(items))
        if n > len(buf):
            buf = ctypes.create_string_buffer(n)
            assert self.fn(data, len(data), buf, len(buf)) == n
        envelope = json.loads(buf.raw[:n])
        assert envelope['schema'] == 1
        rows = envelope['rows']
        assert [r['index'] for r in rows] == list(range(len(items)))
        for r in rows:
            assert r['status'] in ('ok', 'invalid_or_unsupported')
            if r['status'] == 'ok':
                v = r['descriptor_values']
                assert len(v) == len(self.catalog)
                assert all(x is None or (isinstance(x, (int, float)) and math.isfinite(x)) for x in v)
            else:
                assert set(r) == {'index', 'status'}
        return rows


def tolerance(id):
    if id in STABLE_INFORMATION_IDS:
        return 1e-10 if id == 'Ipc' else 1e-12
    if id in ('MolWt', 'HeavyAtomMolWt', 'ExactMolWt'):
        return 1e-8
    if id == 'FractionCSP3':
        return 1e-12
    return 1e-9 if id in ('SPS', 'qed', 'MolLogP', 'TPSA', 'MolMR','BalabanJ','BertzCT','HallKierAlpha','Kappa1','Kappa2','Kappa3','Phi','ds_Harary') or 'EStateIndex' in id or id == 'LabuteASA' or id.endswith('PartialCharge') or id.startswith(('ds_MDE','ds_GGI','ds_MWC','ds_SRW','BCUT2D_','FpDensityMorgan','PEOE_VSA','Chi','SMR_VSA', 'SlogP_VSA', 'EState_VSA', 'VSA_EState')) else 0


def relative_tolerance(id):
    return 1e-12 if id in STABLE_INFORMATION_IDS else 0


def allowed_error(id, expected):
    return max(tolerance(id), relative_tolerance(id) * abs(expected))


def reference_source(spec):
    if spec['id'] in STABLE_INFORMATION_IDS:
        return 'exact-integer-characteristic-polynomial'
    return {'rdkit-standard': 'RDKit2026.03.6',
            'native-profile': 'RDKit2026.03.6-graph-profile',
            'native-extra': 'independent-extra50-formulas'}[spec['origin']]


def reference_value(spec, mol, profile, extra, *, information_backend='python'):
    id = spec['id']
    if id in STABLE_INFORMATION_IDS:
        return float(exact_information_reference(mol, backend=information_backend)[id])
    if spec['origin'] == 'native-extra':
        return float(extra[id]) if extra[id] is not None else math.nan
    if 'EStateIndex' in id:
        return float(FUNCTIONS[id](mol, force=False))
    return float(FUNCTIONS[id](mol) if spec['origin'] == 'rdkit-standard' else profile[id])


def profile_reference(mol, include_extra=False):
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    counts = {e: sum(a.GetSymbol() == e for a in atoms) for e in ('C','N','O','S','P','F','Cl','Br','I')}
    values = {f'ds_Count{e}': n for e, n in counts.items()}
    values.update(
        ds_BondCount=len(bonds), ds_FormalCharge=Chem.GetFormalCharge(mol),
        ds_FragmentCount=len(Chem.GetMolFrags(mol)),
        ds_ExplicitHydrogens=sum(a.GetNumExplicitHs() for a in atoms) + sum(a.GetAtomicNum() == 1 for a in atoms),
        ds_ImplicitHydrogens=sum(a.GetNumImplicitHs() for a in atoms),
        ds_TotalHydrogens=sum(a.GetNumExplicitHs() + a.GetNumImplicitHs() for a in atoms) + sum(a.GetAtomicNum() == 1 for a in atoms),
        ds_HalogenCount=sum(counts[e] for e in ('F','Cl','Br','I')),
        ds_RingAtoms=sum(a.IsInRing() for a in atoms), ds_RingBonds=sum(b.IsInRing() for b in bonds),
        ds_LargestRingSize=max((len(r) for r in Chem.GetSymmSSSR(mol)), default=0),
        ds_AromaticAtoms=sum(a.GetIsAromatic() for a in atoms),
    )
    for id, kind in [('Single', Chem.BondType.SINGLE), ('Double', Chem.BondType.DOUBLE), ('Triple', Chem.BondType.TRIPLE), ('Aromatic', Chem.BondType.AROMATIC)]:
        values[f'ds_{id}Bonds'] = sum(b.GetBondType() == kind for b in bonds)
    if include_extra:
        values.update({key: (value if value is not None else math.nan) for key,value in extra_reference(mol).items()})
    return values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--library', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--input', type=Path, help='CSV with name,smiles columns; defaults to RDKit NCI5K')
    ap.add_argument('--limit', type=int)
    ap.add_argument('--offset', type=int, default=0, help='Original source position; retained in every difference record')
    ap.add_argument('--actual-rows', type=Path, help='Actual Swift JSON result array; validate it without recomputing native results')
    ap.add_argument('--strict', action='store_true')
    ap.add_argument('--information-backend', choices=['python', 'flint'], default='python',
                    help='Exact integer coefficient oracle; python uses FLINT beyond256 vertices/512 edges')
    args = ap.parse_args()
    assert rdkit.__version__ == '2026.03.6', rdkit.__version__
    args.output.mkdir(parents=True, exist_ok=True)
    if args.input:
        source = args.input
        with source.open(newline='') as f:
            inputs = list(csv.DictReader(f))
        assert all('smiles' in r and 'name' in r for r in inputs)
    else:
        source = Path(RDConfig.RDDataDir) / 'NCI' / 'first_5K.smi'
        inputs = [dict(name=f'NCI-{i}', smiles=line.split()[0]) for i, line in enumerate(source.read_text().splitlines()) if line.strip()]
    source_rows = len(inputs)
    assert 0 <= args.offset < source_rows
    assert args.limit is None or args.limit > 0
    inputs = inputs[args.offset:args.offset + args.limit if args.limit else None]
    native = Native(args.library)
    specs = native.catalog
    stored_rows = None
    if args.actual_rows:
        assert args.offset == 0 and args.limit is None, 'stored Swift rows require their complete original source'
        stored_rows = json.loads(args.actual_rows.read_text())
        assert len(stored_rows) == len(inputs)
        assert [r['index'] for r in stored_rows] == list(range(len(inputs)))
        for row in stored_rows:
            assert row['status'] in ('ok', 'invalid_or_unsupported')
            if row['status'] == 'ok':
                assert len(row['descriptor_values']) == len(specs)
                assert all(v is None or (isinstance(v, (int, float)) and math.isfinite(v)) for v in row['descriptor_values'])
    standard = [s for s in specs if s['origin'] == 'rdkit-standard']
    assert all(s['id'] in FUNCTIONS for s in standard)
    stats = {s['id']: dict(reference_source=reference_source(s), tolerance=tolerance(s['id']), relative_tolerance=relative_tolerance(s['id']), compared=0, mismatches=0, undefined_native=0, undefined_reference=0, matching_undefined=0, reference_errors=0, max_abs_error=0, examples=[]) for s in specs}
    assert all(s.get('tolerance', tolerance(s['id'])) == tolerance(s['id']) for s in specs), 'native tolerance differs from independent verifier policy'
    assert all(s.get('relative_tolerance', 0) == relative_tolerance(s['id']) for s in specs), 'native relative tolerance differs from independent verifier policy'
    stable = [s for s in standard if s['id'] in STABLE_INFORMATION_IDS]
    direct_standard = [s for s in standard if s['id'] not in STABLE_INFORMATION_IDS]
    pinned_compatibility = {s['id']: dict(reference_source='RDKit2026.03.6', tolerance=tolerance(s['id']), relative_tolerance=relative_tolerance(s['id']), compared=0, mismatches=0, reference_errors=0, max_abs_error=0, max_relative_error=0, examples=[]) for s in stable}
    pinned_differing_rows = 0
    direct_standard_matches = 0
    status = {}
    formula_mismatches = 0
    full_matches = 0
    standard_matches = 0
    common_valid = 0
    differences_count = 0
    with gzip.open(args.output / 'differences.jsonl.gz', 'wt') as differences, gzip.open(args.output / 'pinned-information-differences.jsonl.gz', 'wt') as pinned_differences:
        for start in range(0, len(inputs), 64):
            batch = inputs[start:start+64]
            rows = (stored_rows[start:start+64] if stored_rows is not None
                    else native.analyze([r['smiles'] for r in batch]))
            for offset, (record, actual) in enumerate(zip(batch, rows)):
                index = args.offset + start + offset
                mol = Chem.MolFromSmiles(record['smiles'], PARAMS) if record['smiles'].strip() else None
                valid = mol is not None and mol.GetNumAtoms() > 0
                status_key = f"native_{actual['status']}__rdkit_{'ok' if valid else 'invalid_or_unsupported'}"
                status[status_key] = status.get(status_key, 0) + 1
                diff = {}
                if actual['status'] != ('ok' if valid else 'invalid_or_unsupported'):
                    diff['status'] = status_key
                elif valid:
                    common_valid += 1
                    profile = profile_reference(mol)
                    extra = extra_reference(mol) if any(s['origin'] == 'native-extra' for s in specs) else {}
                    pinned_diff = {}
                    for spec, value in zip(specs, actual['descriptor_values']):
                        id = spec['id']
                        metric = stats[id]
                        try:
                            if id == 'SPS' and mol.GetNumHeavyAtoms() == 0:
                                try:
                                    FUNCTIONS[id](mol)
                                except ZeroDivisionError:
                                    metric['named_domain_undefined'] = metric.get('named_domain_undefined', 0) + 1
                                    expected = math.nan
                                else:
                                    raise AssertionError('SPS no-heavy domain no longer raises')
                            else:
                                expected = reference_value(spec, mol, profile, extra, information_backend=args.information_backend)
                        except Exception as error:
                            metric['reference_errors'] += 1
                            diff[id] = dict(reference_error=str(error), native=value)
                            continue
                        metric['compared'] += 1
                        if value is None:
                            metric['undefined_native'] += 1
                        if not math.isfinite(expected):
                            metric['undefined_reference'] += 1
                            if value is None:
                                metric['matching_undefined'] += 1
                                continue
                            error = None
                            mismatch = True
                            detail = dict(native=value, reference=None, reference_source=reference_source(spec), reference_nonfinite=True, abs_error=None)
                        else:
                            error = None if value is None else abs(value - expected)
                            if error is not None:
                                metric['max_abs_error'] = max(metric['max_abs_error'], error)
                            mismatch = error is None or error > allowed_error(id, expected)
                            detail = dict(native=value, reference=expected, reference_source=reference_source(spec), abs_error=error)
                        if mismatch:
                            metric['mismatches'] += 1
                            diff[id] = detail
                            if len(metric['examples']) < 5:
                                metric['examples'].append(dict(index=index, **record, **detail))
                    # Pinned numerical compatibility is independently measured.
                    # It never substitutes for the exact-definition reference.
                    actual_values = dict(zip((s['id'] for s in specs), actual['descriptor_values']))
                    for spec in stable:
                        id = spec['id']
                        metric = pinned_compatibility[id]
                        value = actual_values[id]
                        try:
                            expected = float(FUNCTIONS[id](mol))
                        except Exception as error:
                            metric['reference_errors'] += 1
                            pinned_diff[id] = dict(native=value, reference_error=str(error))
                            continue
                        metric['compared'] += 1
                        finite = math.isfinite(expected)
                        error = abs(value - expected) if value is not None and finite else None
                        if error is not None:
                            metric['max_abs_error'] = max(metric['max_abs_error'], error)
                            if expected:
                                metric['max_relative_error'] = max(metric['max_relative_error'], error / abs(expected))
                        mismatch = (value is not None) if not finite else (error is None or error > allowed_error(id, expected))
                        if mismatch:
                            metric['mismatches'] += 1
                            detail = dict(native=value, rdkit=expected if finite else None, reference_nonfinite=not finite, abs_error=error)
                            pinned_diff[id] = detail
                            if len(metric['examples']) < 5:
                                metric['examples'].append(dict(index=index, **record, **detail))
                    if pinned_diff:
                        pinned_differing_rows += 1
                        pinned_differences.write(json.dumps(dict(index=index, **record, differences=pinned_diff), ensure_ascii=False) + '\n')
                    formula = rdMolDescriptors.CalcMolFormula(mol)
                    if actual['formula'] != formula:
                        formula_mismatches += 1
                        diff['formula'] = dict(native=actual['formula'], rdkit=formula)
                    if not any(s['id'] in diff for s in standard) and 'formula' not in diff:
                        standard_matches += 1
                    if not any(s['id'] in diff for s in direct_standard) and 'formula' not in diff:
                        direct_standard_matches += 1
                    if not diff:
                        full_matches += 1
                if diff:
                    differences_count += 1
                    differences.write(json.dumps(dict(index=index, **record, differences=diff), ensure_ascii=False) + '\n')
            if start % 4096 == 0:
                print(f'compared {min(start+64,len(inputs))}/{len(inputs)}; differing rows {differences_count}', flush=True)
    report = dict(
        gate='failed' if differences_count else 'passed', revision=native.revision,
        rdkit_version=rdkit.__version__, platform=platform.platform(), library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),
        dataset=dict(path=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest(), rows=len(inputs), distinct_raw_smiles=len({r['smiles'] for r in inputs}), source_rows=source_rows, index_offset=args.offset, index_stop=args.offset+len(inputs)),
        actual_rows=(dict(path=str(args.actual_rows), sha256=hashlib.sha256(args.actual_rows.read_bytes()).hexdigest(), native_recomputed=False) if args.actual_rows else dict(source='native C ABI', native_recomputed=True)),
        catalog_count=len(specs), registered_standard_count=len(standard),
        direct_pinned_standard_count=len(direct_standard), exact_definition_standard_count=len(stable),
        profile_count=sum(s['origin']=='native-profile' for s in specs), additional_count=sum(s['origin']=='native-extra' for s in specs),
        reference_policy=dict(removeHs=False, sanitize=True, standardize=False, tpsa_includeSandP=False, rotatable_bonds='Strict', ring_basis='RDKit symmetric SSSR', standard='215 pinned numerical definitions; Ipc/AvgIpc use exact integer characteristic coefficients with independent Python full matrices and f64 entropy; each field declares its reference_source', additional='extra50_reference.py independent dense formulas; RDKit parses graph only', undefined_values='reference nonfinite + native null matches as undefined; SPS no-heavy ZeroDivisionError is explicitly asserted and counted as named_domain_undefined; all other computation exceptions remain errors'),
        status_counts=status, common_valid=common_valid, all_values_matching_rows=full_matches,
        standard_and_formula_matching_rows=standard_matches, differing_rows=differences_count,
        direct_pinned_standard_and_formula_matching_rows=direct_standard_matches,
        pinned_information_compatibility=dict(gate='failed' if pinned_differing_rows else 'passed', scope='common accepted rows only; separate from exact-definition correctness; same declared tolerances', differing_rows=pinned_differing_rows, metrics=pinned_compatibility, differences='pinned-information-differences.jsonl.gz'),
        formula_mismatches=formula_mismatches, metrics=stats,
        caveat='Agreement with explicitly declared mixed references, not blanket pinned RDKit parity or experimental property accuracy. Pinned Ipc/AvgIpc numerical compatibility is reported separately. No speed comparison in this run.',
    )
    import flint
    report['reference_policy']['information_coefficients'] = dict(
        backend=args.information_backend,
        python_small_domain_full_object_matrices=args.information_backend == 'python',
        flint_fallback_domain='>256 vertices or >512 edges',
        flint_version=flint.__version__, flint_library_version=flint.__FLINT_VERSION__,
        exact_integer=True, entropy='unchanged math.fsum binary64',
    )
    report['reference_policy']['standard'] = '215 pinned numerical definitions; Ipc/AvgIpc use exact integer characteristic coefficients and f64 entropy; coefficient backend declared separately; each field declares its reference_source'
    (args.output / 'comparison.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'metrics'}, ensure_ascii=False), flush=True)
    if args.strict and differences_count:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
