"""Generate independent controls before registering/implementing extra50."""
from pathlib import Path
import hashlib
import json
import math
import random
import rdkit
from rdkit import Chem, RDLogger
from extra50_reference import IDS, MDE, reference

assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
params = Chem.SmilesParserParams(); params.removeHs = False
sources = ['C', '[Na+].[Cl-]', '[H][H]', '[2H]', 'CC', 'CCC', 'CCCC',
           'CC(C)C', 'CC(C)(C)C', 'CC(C)C(C)C', 'C1CC1', 'C1CCC1',
           'c1ccccc1', 'C1=CC=CC=C1', 'C1CC2CCC1C2', 'C12C3C4C1C5C2C3C45',
           'CC.O', 'CC.CC', '[13CH3]CO', '[2H]C([2H])([2H])C', '[H]C([H])([H])C',
           'NCCN', 'NCC(N)N', 'CN(C)CCN(C)C', 'CN(C)CCN(C)CCN(C)C',
           'OCCO', 'COCCOC', 'COCC(O)CO', 'OCC(O)COCCO',
           'C[N+](C)(C)CCO', 'C:C', 'C1:C:C:C:C:C#1', '[OH2][Mg+2]',
           '[NH3]->[Cu+2]<-[NH3]', 'C[H-]C', 'C.C', 'N.N', 'C'*64]
for path in [root/'tasks/evidence/catalog236-qed/controls413.csv']:
    import csv
    sources += [r['smiles'] for r in csv.DictReader(path.open(newline=''))]
# Seeded carbon trees and atom permutations test graph invariance. Selection
# uses the independent parser only; no native acceptance/results are consulted.
rng = random.Random(91376)
for size in range(5, 25):
    rw = Chem.RWMol()
    for _ in range(size): rw.AddAtom(Chem.Atom(6))
    degrees = [0]*size
    for atom in range(1, size):
        parent = rng.choice([i for i in range(atom) if degrees[i] < 4])
        rw.AddBond(parent, atom, Chem.BondType.SINGLE)
        degrees[parent] += 1; degrees[atom] += 1
    mol = rw.GetMol(); Chem.SanitizeMol(mol)
    sources.append(Chem.MolToSmiles(mol, canonical=False))
    order = list(range(size)); rng.shuffle(order)
    sources.append(Chem.MolToSmiles(Chem.RenumberAtoms(mol, order), canonical=False))
sources = list(dict.fromkeys(sources))
lines = ['# Original SMILES; extra50 ID; independent dense reference (None=null)']
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None and mol.GetNumAtoms(), source
    values = reference(mol)
    for key in IDS:
        value = values[key]
        assert value is None or math.isfinite(value), (source, key)
        lines.append(f'{source}\t{key}\t' + ('null' if value is None else f'{value:.17g}'))

# Hand-calculated graph controls independently check the matrix oracle.
ethane = reference(Chem.MolFromSmiles('CC'))
assert ethane['ds_Wiener'] == ethane['ds_Harary'] == 1
assert ethane['ds_DegreeDistance'] == 2 and ethane['ds_Gutman'] == 1
assert ethane['ds_EccentricConnectivity'] == 2 and ethane['ds_MDEC11'] == 1
assert ethane['ds_SRW3'] == 0 and abs(ethane['ds_SRW4']-math.log(3)) < 1e-14
chain = reference(Chem.MolFromSmiles('CCCC'))
assert chain['ds_Wiener'] == 10 and chain['ds_DegreeDistance'] == 28
assert chain['ds_Gutman'] == 19 and chain['ds_EccentricConnectivity'] == 14
assert abs(chain['ds_MDEC11'] - 1/3) < 1e-14
assert reference(Chem.MolFromSmiles('CC.CC'))['ds_Wiener'] == 2
assert all(v is None for v in reference(Chem.MolFromSmiles('[H][H]', params)).values())
fixture = root/'crates/smiles/tests/fixtures/descriptor_extra50.tsv'
fixture.write_text('\n'.join(lines)+'\n')
(root/'tasks/descriptors/extra50-controls.json').write_text(json.dumps(dict(
    rdkit=rdkit.__version__, source_count=len(sources), scalar_count=len(sources)*50,
    descriptor_ids=IDS, source_smiles=sources, native_acceptance_filter=False,
    fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),
    oracle='RDKit parsing only; independent NumPy dense powers/Floyd-Warshall, explicit formulas and hand graph checks'), indent=2)+'\n')
print(len(sources), 'sources;', len(sources)*50, 'independent extra50 cells')
