#!/usr/bin/env python3
"""Real SQL/FFI tests and independent scientific comparisons (dev dependencies only).

Run with pinned RDKit 2026.3.6, DuckDB matching the built extension, pymatgen,
scipy and chemicals. Example: python scripts/check_materials.py --extension ...
"""
import argparse
import ctypes
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import random
import re

import _duckdb
import duckdb
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors
from pymatgen.core import Composition
from chemicals.utils import mixing_logarithmic
from scipy.stats import linregress


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extension', type=Path, required=True)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.03.6', rdBase.rdkitVersion
    ctypes.CDLL(_duckdb.__file__, mode=ctypes.RTLD_GLOBAL)
    con = duckdb.connect(config={'allow_unsigned_extensions': True})
    extension = args.extension.resolve()
    con.execute("LOAD '" + str(extension).replace("'", "''") + "'")
    root = Path(__file__).resolve().parents[1]
    checks = 0
    def profile(name, data, params=None):
        nonlocal checks
        value = con.execute('SELECT ' + name + '(' + ('?,?' if params is not None else '?') + ')',
                            [data, params] if params is not None else [data]).fetchone()[0]
        assert value is not None, (name, data, params)
        checks += 1
        return json.loads(value)
    def close(actual, expected, tol=1e-10):
        nonlocal checks
        assert math.isclose(actual, expected, abs_tol=tol, rel_tol=tol), (actual, expected, tol)
        checks += 1
    for name in ['materials.test', 'function_catalog.test']:
        text = (root / 'test/sql' / name).read_text()
        for block in re.split(r'\n\s*\n', text):
            lines = [line for line in block.splitlines() if not line.startswith('#')]
            if lines and lines[0].startswith('require '):
                assert len(lines) == 1, 'sqllogictest requires an empty line after each require'
                continue
            if not lines or not lines[0].startswith('query '):
                continue
            split = lines.index('----')
            def cell(x):
                return 'NULL' if x is None else str(x).lower() if isinstance(x, bool) else str(x)
            actual = ['\t'.join(map(cell, row)) for row in con.execute('\n'.join(lines[1:split])).fetchall()]
            expected = lines[split+1:]
            if 'rowsort' in lines[0]:
                actual.sort(); expected.sort()
            assert actual == expected, (name, actual, expected)
            checks += 1
    solvents = ['O', 'CCO', 'CC(C)O', 'CS(C)=O', 'CN(C)C=O', 'CN1CCCC1=O',
                'O=C1OCCC1', 'CC#N', 'COc1ccccc1', 'Clc1ccccc1', 'C[NH3+]',
                'NC(N)=[NH2+]', 'CC(=O)O', 'CCOC(C)=O', 'ClCCl']
    largest_error = {}
    for smiles in solvents:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        actual = profile('solvent_descriptors_json', smiles)['descriptors']
        for key, value in actual.items():
            expected = getattr(Descriptors, key)(mol)
            assert value is not None, (smiles, key)
            close(value, expected, 1e-8 if key == 'MolWt' else 1e-9)
            largest_error[key] = max(largest_error.get(key, 0), abs(value-expected))
    formulas = ['CsPbI3', 'Cs0.5Rb0.5Pb(I0.8Br0.2)3', 'Ca(OH)2', 'BaTiO3',
                'SrTiO3', 'Fe2O3', 'Li0.5Na0.5Cl', 'Al2(SO4)3', 'CH6NPbI3']
    pt = Chem.GetPeriodicTable()
    for formula in formulas:
        expected = Composition(formula).get_el_amt_dict()
        total = sum(expected.values())
        j = profile('composition_descriptors_json', formula)
        actual = {v['element']: v for v in j['elements']}
        assert actual.keys() == expected.keys()
        for e, amount in expected.items():
            close(actual[e]['amount'], amount)
            close(actual[e]['atomic_fraction'], amount/total)
        close(j['formula_weight_g_mol'], sum(pt.GetAtomicWeight(e)*n for e,n in expected.items()))
        close(j['mean_atomic_number'], sum(pt.GetAtomicNumber(e)*n for e,n in expected.items())/total)
        for p, norm in zip(j['stoichiometry_p'], j['stoichiometry_norms']):
            close(norm, sum((n/total)**p for n in expected.values())**(1/p))
    rng = random.Random(4281)
    for _ in range(100):
        ra, rb, rx = [rng.uniform(0.5, 3.) for _ in range(3)]
        row = profile('chem_perovskite_descriptors', [[ra,rb,rx,1.,2.,-1.]], [])['rows'][0]
        close(row['goldschmidt_t'], (ra+rx)/math.sqrt(2)/(rb+rx))
        close(row['octahedral_mu'], rb/rx)
        assert row['charge_balanced']
        x = rng.random(); eta = [10**rng.uniform(-5,0), 10**rng.uniform(-5,0)]
        mixed = profile('chem_viscosity_mix', [[x,eta[0]],[1-x,eta[1]]], [])
        close(mixed['viscosity_Pa_s'], mixing_logarithmic([x,1-x],eta))
        g = rng.uniform(-2,2)
        mixed = profile('chem_viscosity_mix', [[x,eta[0]],[1-x,eta[1]]], [g])
        close(mixed['viscosity_Pa_s'], mixing_logarithmic([x,1-x],eta)*math.exp(x*(1-x)*g))
    temperatures = [280.+i*5 for i in range(20)]
    measured = [math.exp(-12.+1500/t+rng.gauss(0,0.02)) for t in temperatures]
    fit = profile('chem_viscosity_fit', list(map(list,zip(temperatures,measured))), [])
    oracle = linregress([1/t for t in temperatures], [math.log(v) for v in measured])
    close(fit['parameters'][0], oracle.intercept)
    close(fit['parameters'][1], oracle.slope)
    predicted = profile('chem_viscosity_arrhenius', [[t] for t in temperatures],fit['parameters'])
    for t, v, observed in zip(temperatures, predicted['viscosity_Pa_s'], measured):
        close(v, math.exp(oracle.intercept+oracle.slope/t))
    close(fit['rmse_log'], math.sqrt(sum((math.log(y/v))**2 for y,v in zip(measured,predicted['viscosity_Pa_s']))/len(measured)))
    # Data-dependent nested lists cross multiple DuckDB vector batches and dictionary selections.
    rows = con.execute("SELECT i, chem_viscosity_arrhenius([[280.0+i%20]], [-12.0,1500.0]) FROM range(5000) t(i) WHERE i%3<>0").fetchall()
    assert len(rows) == 3333
    for i, j in rows:
        close(json.loads(j)['viscosity_Pa_s'][0], math.exp(-12.+1500/(280+i%20)))
    # Execute every new API example, including the joined solvent/mixture feature pipeline.
    for block in re.findall(r'```sql\n(.*?)```', (root/'MATERIALS.md').read_text(), re.S):
        assert con.execute(block).fetchone()[0] is not None
        checks += 1
    versions = {p:importlib.metadata.version(p) for p in ['duckdb','rdkit','pymatgen','scipy','chemicals']}
    receipt = {'schema':1, 'checks':checks, 'versions':versions,
               'solvents':len(solvents), 'solvent_max_abs_errors':largest_error,
               'compositions':len(formulas), 'random_perovskite_and_mixture_cases':100,
               'extension_sha256':hashlib.sha256(extension.read_bytes()).hexdigest(),
               'scope':'native SQL, descriptors and correlation numerics; no experimental predictive-accuracy or full RDKit/structure parity claim'}
    if args.receipt:
        args.receipt.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
