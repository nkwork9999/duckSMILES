"""Independent directed-coordination graphs and all-catalog expectations."""
from pathlib import Path
import csv,json
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,rdMolDescriptors,rdqueries
from compare_catalog import profile_reference
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')
root=Path(__file__).resolve().parents[2];fixture=root/'crates/smiles/tests/fixtures'
rows=list(csv.DictReader((root/'tasks/evidence/catalog232-cleanup/changed-inputs/source.csv').open(newline='')))
ce=next(r['smiles'] for r in rows if r['name']=='ChEBI-structure-23350')
sources=[ce,'[NH3][Cu+2]','CN(C)(C)[Cu+2]','[Cu+2]N(C)(C)C',
 '[OH2][Mg+2]','[OH2][Mg+2]([OH2])([OH2])([OH2])([OH2])[OH2]',
 '[OH2][Zn+2]([OH2])([OH2])[OH2]','CN(C)[Cu+2]',
 '[Cu+2]N1=CC=CC=C1','C1=CC=N([Cu+2])C=C1',
 '[NH3][Co+3]([NH3])([NH3])([NH3])([NH3])[NH3]',
 '[NH3][Pt]([NH3])(Cl)Cl','[Cu+2]1[NH2]CC[NH2]1',
 '[CH3-][Na+]','[PH3][Cu+]','[BH3-][Na+]','[NH3][Ge]',
 '[NH3][Sb]','N.[Cu+2]','Cl[Cu+2]Cl',
 '[Fe]c1([H])ccccc1','[Fe]c1([2H])ccccc1','[Fe]C1([H])=CC=CC=C1']
chebi=Path('/private/tmp/chemob-descriptor-data/concrete/chebi-10000.csv')
chebi_rows=list(csv.DictReader(chebi.open(newline='')))
# Existing coordination error rows, selected from the prior independent audit,
# not from whether this implementation accepts or matches their values.
audit_indices=[39,82,3394,3573,4657,4775,5519,5671,6085,7684,8334,9156,9517]
sources.extend(chebi_rows[i]['smiles'] for i in audit_indices)
sources=list(dict.fromkeys(sources));p=Chem.SmilesParserParams();p.removeHs=False
rejected_sources=['[Cu+2]n1ccccc1','c1ccn([Cu+2])cc1']
for smiles in rejected_sources:
 assert Chem.MolFromSmiles(smiles,p) is None,smiles
specs=json.loads((root/'tasks/descriptors/phase1-catalog.json').read_text())['descriptors'];functions=dict(Descriptors._descList)
scalar=['# RDKit2026.03.6 registered/profile values and formula; original source retained']
contexts=['# source, atoms(symbol,H,charge,arom,noImplicit,isotope,totalValence,hybrid,radicals), bonds(a,b,RDKitType,conjugated), symmetric bond rings']
counts=[]
for smiles in sources:
 mol=Chem.MolFromSmiles(smiles,p);assert mol is not None,smiles
 profile=profile_reference(mol, include_extra=True)
 scalar.append(f'{smiles}\t@formula\t{rdMolDescriptors.CalcMolFormula(mol)}')
 for spec in specs:
  id=spec['id'];value=functions[id](mol) if spec['origin']=='rdkit-standard' else profile[id]
  scalar.append(f'{smiles}\t{id}\t{float(value):.17g}')
 atoms=';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())},{a.GetIsotope()},{a.GetTotalValence()},{int(a.GetHybridization())},{a.GetNumRadicalElectrons()}' for a in mol.GetAtoms())
 bonds=';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{int(b.GetBondType())},{int(b.GetIsConjugated())}' for b in mol.GetBonds())
 rings=';'.join(','.join(map(str,r)) for r in mol.GetRingInfo().BondRings())
 contexts.append(f'{smiles}\t{atoms}\t{bonds}\t{rings}')
 counts.append(sum(b.GetBondType()==Chem.BondType.DATIVE for b in mol.GetBonds()))
(fixture/'descriptor_coordination.tsv').write_text('\n'.join(scalar)+'\n')
(fixture/'coordination_contexts.tsv').write_text('\n'.join(contexts)+'\n')
query=rdqueries.MAtomQueryAtom()
meta=dict(rdkit=rdkit.__version__,catalog_count=len(specs),molecules=len(sources),scalar_values=len(specs)*len(sources),formula_values=len(sources),source_smiles=sources,dative_bond_counts=counts,metal_atomic_numbers=[z for z in range(1,119) if query.Match(Chem.Atom(z))],scope='Curated automatic-coordination and negative controls; no native acceptance/agreement filtering. Raw arrow notation and ambiguous multi-metal ranking require separate work.')
meta['reference_rejected_source_smiles']=rejected_sources
meta['prior_audit_source_indices']=audit_indices
meta['prior_audit_source_sha256']=__import__('hashlib').sha256(chebi.read_bytes()).hexdigest()
(root/'tasks/descriptors/coordination-controls.json').write_text(json.dumps(meta,indent=2)+'\n')
print(len(sources),'source inputs;',meta['scalar_values'],'scalar values;',counts,'dative bonds')
