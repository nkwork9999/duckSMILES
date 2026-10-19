"""Fixed v42 residuals and atom-order controls; no tolerance or API changes."""
import gzip, json, math, random, hashlib, sys, os
from pathlib import Path
import numpy, rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from flint_ipc_reference import coefficients
from compare_catalog import Native, PARAMS

assert rdkit.__version__ == '2026.03.6' and numpy.__version__ == '2.5.4'
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
dest = root/'tasks/evidence/depiction-ipc-v43'
residual = root/'tasks/evidence/catalog289-sps/v42-full100k/pinned-information-differences.jsonl.gz'
rows = [json.loads(s) for s in gzip.open(residual,'rt')]
assert len(rows) == 5723
counts = {}; ranges = {}
for row in rows:
    for key, d in row['differences'].items():
        counts[key] = counts.get(key,0)+1
        if not d['reference_nonfinite']:
            relative = abs(d['native']-d['rdkit'])/max(abs(d['native']),1e-300)
            bucket = 'over100percent' if relative > 1 else 'over1percent' if relative > .01 else 'atMost1percent'
            ranges.setdefault(key,{}); ranges[key][bucket] = ranges[key].get(bucket,0)+1
native = Native(Path('/private/tmp/chemob-descriptor-eval/libducksmiles-v42-catalog.dylib'))
indices = {s['id']:i for i,s in enumerate(native.catalog)}
sources = [('path40','C'*40),('path80','C'*80),('original22',rows[0]['smiles'])]
results=[]
for name,source in sources:
    mol = Chem.MolFromSmiles(source,PARAMS)
    assert mol is not None
    order = list(range(mol.GetNumAtoms()))
    variants = [order, list(reversed(order))]
    random.Random(438).shuffle(order); variants.append(order)
    versions=[]
    for permutation in [None]+variants:
        m = Chem.Mol(mol) if permutation is None else Chem.RenumberAtoms(mol,permutation)
        # Serialize and reparse the same original graph; source identity verified.
        smi=source if permutation is None else Chem.MolToSmiles(m,canonical=False,isomericSmiles=True)
        m=Chem.MolFromSmiles(smi,PARAMS)
        assert Chem.MolToSmiles(m)==Chem.MolToSmiles(mol)
        adjacency = [sorted(n.GetIdx() for n in a.GetNeighbors()) for a in m.GetAtoms()]
        c = [abs(v) for v in coefficients(adjacency)]
        total = sum(c)
        avg = -sum((v/total)*math.log2(v/total) for v in c if v)
        actual = native.analyze([smi])[0]; assert actual['status']=='ok'
        expected = {'Ipc':total*avg,'AvgIpc':avg}
        stable = {k:actual['descriptor_values'][indices[k]] for k in expected}
        for k,v in expected.items():assert abs(stable[k]-v) <= max(1e-10,abs(v)*1e-12)
        versions.append(dict(smiles=smi,pinned=dict(Ipc=Descriptors.Ipc(m),AvgIpc=Descriptors.AvgIpc(m)),stable=stable,exact_reference=expected))
    results.append(dict(name=name,atoms=mol.GetNumAtoms(),variants=versions,
        pinned_Ipc_max_to_min=max(v['pinned']['Ipc'] for v in versions)/min(v['pinned']['Ipc'] for v in versions),
        stable_values_invariant=all(v['stable']==versions[0]['stable'] for v in versions)))
report = dict(rdkit=rdkit.__version__,numpy=numpy.__version__,residual_rows=len(rows),
    thread_environment={k:os.environ.get(k) for k in ["VECLIB_MAXIMUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS"]},
    residual_sha256=hashlib.sha256(residual.read_bytes()).hexdigest(),counts=counts,
    relative_to_stable_bins=ranges,atom_order_controls=results,
    compatibility='No portable pinned-numeric compatibility claim. Stable fields remain separate; do not silently ingest them as legacy RDKit features.',
    tolerance_changed=False,calculation_code_changed=False)
(dest/'information-order-diagnostic.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='atom_order_controls'},indent=2))
print([(r['name'],r['pinned_Ipc_max_to_min'],r['stable_values_invariant']) for r in results])
