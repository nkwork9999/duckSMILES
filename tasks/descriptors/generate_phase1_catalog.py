from pathlib import Path
import json
core=[('MolWt','Molecular weight','g/mol',2),('ExactMolWt','Exact mass','Da',5),('MolLogP','cLogP','',2),('TPSA','TPSA','Å²',2),('NumHDonors','H-bond donors','',0),('NumHAcceptors','H-bond acceptors','',0),('HeavyAtomCount','Heavy atoms','',0),('ds_BondCount','Bonds','',0),('ds_FormalCharge','Formal charge','',0),('ds_FragmentCount','Fragments','',0),('NumRotatableBonds','Rotatable bonds','',0),('FractionCSP3','Fraction sp³ carbons','',2)]
rings=[('RingCount','Rings'),('NumAromaticRings','Aromatic rings'),('NumAliphaticRings','Aliphatic rings'),('NumSaturatedRings','Saturated rings'),('NumAromaticHeterocycles','Aromatic heterocycles'),('NumAromaticCarbocycles','Aromatic carbocycles'),('NumAliphaticHeterocycles','Aliphatic heterocycles'),('NumAliphaticCarbocycles','Aliphatic carbocycles'),('NumSaturatedHeterocycles','Saturated heterocycles'),('NumSaturatedCarbocycles','Saturated carbocycles'),('NumHeterocycles','Heterocycles')]
composition=[('HeavyAtomMolWt','Heavy-atom mass','g/mol',2),('NumHeteroatoms','Heteroatoms','',0),('NOCount','N + O atoms','',0),('NHOHCount','NH + OH hydrogens','',0),('NumValenceElectrons','Valence electrons','',0)]
composition += [(f'ds_{key}',title,'',0) for key,title in [('ExplicitHydrogens','Explicit hydrogens'),('ImplicitHydrogens','Implicit hydrogens'),('TotalHydrogens','Total hydrogens')]]+[(f'ds_Count{e}',f'{e} atoms','',0) for e in ['C','N','O','S','P','F','Cl','Br','I']]+[('ds_HalogenCount','Halogens','',0)]
graph=[(f'ds_{key}',title,'',0) for key,title in [('SingleBonds','Single bonds'),('DoubleBonds','Double bonds'),('TripleBonds','Triple bonds'),('AromaticBonds','Aromatic bonds'),('RingAtoms','Ring atoms'),('RingBonds','Ring bonds'),('LargestRingSize','Largest ring size'),('AromaticAtoms','Aromatic atoms')]]+[('MolMR','Molar refractivity','',2)]
rows=[]
for group,items in [('Essentials',core),('Rings',[(k,v,'',0) for k,v in rings]),('Composition',composition),('Graph',graph)]:
 for id,title,unit,digits in items: rows.append(dict(id=id,title=title,group=group,unit=unit,digits=digits,origin='native-profile' if id.startswith('ds_') else 'rdkit-standard',verification='pending-expanded-validation'))
assert len(rows)==50
fragments = json.loads(Path('tasks/descriptors/fragment-patterns.json').read_text())['descriptors']
rows += fragments
assert len(rows) == 135
for id, title in [('MaxEStateIndex','Max EState index'),('MinEStateIndex','Min EState index'),
                  ('MaxAbsEStateIndex','Max |EState|'),('MinAbsEStateIndex','Min |EState|')]:
 rows.append(dict(id=id,title=title,group='Electronic',unit='',digits=4,
                  origin='rdkit-standard',verification='pending-expanded-validation'))
assert len(rows) == 139
rows += json.loads(Path('tasks/descriptors/surface-definitions.json').read_text())['descriptors']
assert len(rows) == 183
rows += json.loads(Path('tasks/descriptors/graph-definitions.json').read_text())['descriptors']
assert len(rows) == 202
rows += json.loads(Path('tasks/descriptors/charge-definitions.json').read_text())['descriptors']
assert len(rows) == 221
rows += json.loads(Path('tasks/descriptors/bcut-definitions.json').read_text())['descriptors']
assert len(rows) == 229
rows += json.loads(Path('tasks/descriptors/structural-definitions.json').read_text())['descriptors']
assert len(rows) == 232
rows += json.loads(Path('tasks/descriptors/morgan-density-definitions.json').read_text())['descriptors']
assert len(rows) == 235
rows += json.loads(Path('tasks/descriptors/qed-definitions.json').read_text())['descriptors']
assert len(rows) == 236
rows += json.loads(Path('tasks/descriptors/extra50-definitions.json').read_text())['descriptors']
assert len(rows) == 286
rows += json.loads(Path('tasks/descriptors/stereo-definitions.json').read_text())['descriptors']
assert len(rows) == 288
rows += json.loads(Path('tasks/descriptors/sps-definitions.json').read_text())['descriptors']
assert len(rows) == 289
rows += json.loads(Path('tasks/descriptors/ipc-definitions.json').read_text())['descriptors']
assert len(rows) == 291
for row in rows:
 if row['id'] not in ('Ipc', 'AvgIpc'):
  row['tolerance'] = 0.0 if row['digits']==0 else (1e-8 if row['id'] in ['MolWt','HeavyAtomMolWt','ExactMolWt'] else 1e-12 if row['id']=='FractionCSP3' else 1e-9)
 row['relative_tolerance'] = row.get('relative_tolerance', 0.0)
known=json.loads(Path('tasks/descriptors/rdkit-2026.03.6-inventory.json').read_text());ids={r['id'] for r in known['descriptors']}
assert all(r['id'] in ids for r in rows if r['origin']=='rdkit-standard')
Path('tasks/descriptors/phase1-catalog.json').write_text(json.dumps({'schema':1,'count':len(rows),'descriptors':rows},ensure_ascii=False,indent=2)+'\n')
# Serialized once per catalog call; byte literal keeps serialization dependency-free.
serialized=json.dumps({'schema':1,'descriptors':rows},ensure_ascii=False,separators=(',',':'))
s='//! Versioned descriptor catalog and parse-once profile calculation.\n'
s+='pub struct DescriptorSpec { pub id: &\'static str, pub group: &\'static str, pub tolerance: f64, pub relative_tolerance: f64 }\n'
s+='pub const SPECS: &[DescriptorSpec] = &[\n'
for r in rows:
 tol=r['tolerance']
 s+=f'    DescriptorSpec {{ id: "{r["id"]}", group: "{r["group"]}", tolerance: {float(tol)}, relative_tolerance: {float(r["relative_tolerance"])} }},\n'
s+='];\n\npub const CATALOG_JSON: &str = r###"'+serialized+'"###;\n'
Path('crates/smiles/src/descriptor_catalog.rs').write_text(s)
