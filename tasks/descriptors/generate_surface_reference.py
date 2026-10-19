"""Independent RDKit surface fixtures, atomic radii and bin definitions.

Run before generate_phase1_catalog.py and generate_reference_fixture.py.
The native implementation is never consulted for selection or expected values.
"""
import json
from pathlib import Path
import rdkit
from rdkit import Chem, RDConfig
from rdkit.Chem import Descriptors, Crippen, MolSurf, rdMolDescriptors
from rdkit.Chem.EState import EState_VSA
assert rdkit.__version__ == '2026.03.6'
root = Path(__file__).resolve().parents[2]
params = Chem.SmilesParserParams(); params.removeHs = False
families = [('SMR_VSA', MolSurf.mrBins, 'MR atom contribution', 'Å²', 'Surface area'),
            ('SlogP_VSA', MolSurf.logpBins, 'logP atom contribution', 'Å²', 'Surface area'),
            ('EState_VSA', EState_VSA.estateBins, 'EState index', 'Å²', 'Electronic surface'),
            ('VSA_EState', EState_VSA.vsaBins, 'Labute atomic area', '', 'Electronic surface')]
rows = [dict(id='LabuteASA', title='Labute ASA', group='Surface area', unit='Å²', digits=4,
             origin='rdkit-standard', verification='pending-expanded-validation',
             description='Labute approximate surface area, includeHs=True; atom areas plus the aggregate hydrogen correction.')]
for prefix, bins, axis, unit, group in families:
    for i in range(len(bins)+1):
        lower = bins[i-1] if i else None
        upper = bins[i] if i < len(bins) else None
        rows.append(dict(id=f'{prefix}{i+1}', title=f'{prefix.replace("_", " ")} {i+1}',
                         group=group, unit=unit, digits=4, origin='rdkit-standard',
                         verification='pending-expanded-validation', axis=axis,
                         lower_inclusive=lower, upper_exclusive=upper,
                         description=('Sum EState indices by atomic area' if prefix=='VSA_EState' else 'Sum atomic Labute areas') + '; lower-inclusive, upper-exclusive bins (lower <= value < upper).'))
assert len(rows) == 44
functions = dict(Descriptors._descList)
assert all(s['id'] in functions for s in rows)
smiles = []
for name in ['descriptor_phase1.tsv','descriptor_fragments.tsv','aromatic_fused.tsv']:
    for line in (root/'crates/smiles/tests/fixtures'/name).read_text().splitlines():
        if line.startswith('#'): continue
        s = line.split('\t')[0]
        if s not in smiles: smiles.append(s)
for s in ['[2H]O', '[H]O', '[2H]N', '[2H]C', '[2H]P', '[2H]S', '[H]N([H])', '[H]', '[Og]', '[Og][Og]', 'C[Si](C)(C)C', 'CI', 'C[Se]C', '[H][N+]([H])([H])[H]']:
    if s not in smiles: smiles.append(s)
seen_positive = set()
def record_positive(s):
    m = Chem.MolFromSmiles(s,params); assert m is not None,s
    for spec in rows:
        value = functions[spec['id']](m)
        if value > 0: seen_positive.add(spec['id'])
for s in smiles: record_positive(s)
# SMR_VSA8 is unreachable with the pinned standard Crippen table: no MR
# parameter is in [3.63,3.8). It remains a real registered zero-valued bin.
unreachable = {'SMR_VSA8', 'SlogP_VSA9'}
for line in (Path(RDConfig.RDDataDir)/'NCI/first_5K.smi').read_text().splitlines():
    if len(seen_positive | unreachable) == len(rows): break
    s = line.split()[0]
    if s in smiles: continue
    m = Chem.MolFromSmiles(s,params)
    if m is None or m.GetNumAtoms()>128 or any(a.GetAtomicNum()==0 for a in m.GetAtoms()): continue
    positive = {spec['id'] for spec in rows if functions[spec['id']](m)>0}
    if positive - seen_positive:
        smiles.append(s); seen_positive |= positive
missing = {r['id'] for r in rows} - seen_positive
assert missing == unreachable, missing
# Independently inspect the full reference table for both unreachable bins.
Crippen._Init()
assert not any(3.63 <= vals[3] < 3.8 for patterns in Crippen._smartsPatterns.values() for vals in patterns)
assert not any(0.3 <= vals[2] < 0.4 for patterns in Crippen._smartsPatterns.values() for vals in patterns)
scalar = ['# RDKit 2026.03.6; removeHs=False; SMILES, descriptor, expected; nan means undefined']
atomic = ['# RDKit 2026.03.6; SMILES, kind (area/hydrogen/logp/mr), atom index, expected']
for s in smiles:
    m = Chem.MolFromSmiles(s,params); assert m is not None
    for spec in rows: scalar.append(f"{s}\t{spec['id']}\t{functions[spec['id']](m):.16g}")
    area, h = rdMolDescriptors._CalcLabuteASAContribs(m)
    atomic.append(f'{s}\thydrogen\t0\t{h:.16g}')
    contributions = Crippen._GetAtomContribs(m,force=True)
    for i, (a, (lp,mr)) in enumerate(zip(area,contributions)):
        for key,value in [('area',a),('logp',lp),('mr',mr)]:
            atomic.append(f'{s}\t{key}\t{i}\t{value:.16g}')
(root/'crates/smiles/tests/fixtures/descriptor_surface.tsv').write_text('\n'.join(scalar)+'\n')
(root/'crates/smiles/tests/fixtures/surface_atomic.tsv').write_text('\n'.join(atomic)+'\n')
report = dict(rdkit=rdkit.__version__, descriptors=rows, molecules=len(smiles),
              scalar_reference_values=len(scalar)-1, atomic_reference_values=len(atomic)-1,
              positive_ids=sorted(seen_positive), unreachable_positive_bins=sorted(unreachable))
(root/'tasks/descriptors/surface-definitions.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
radii = [Chem.GetPeriodicTable().GetRb0(z) for z in range(1,119)]
source = '//! RDKit 2026.03.6 PeriodicTable Rb0 radii and default VSA bin boundaries.\n'
source += f'pub(crate) const RADII: [f64; 118] = {radii!r};\n'
for name,bins in [('MR_BINS',MolSurf.mrBins),('LOGP_BINS',MolSurf.logpBins),('ESTATE_BINS',EState_VSA.estateBins),('AREA_BINS',EState_VSA.vsaBins)]:
    source += f'pub(crate) const {name}: [f64; {len(bins)}] = {list(map(float,bins))!r};\n'
(root/'crates/smiles/src/surface_data.rs').write_text(source)
print(len(smiles),'molecules;',len(scalar)-1,'scalar refs;',len(atomic)-1,'atomic refs; missing positives',missing)
