"""Independent registered Morgan density and sparse-key controls (no native filtering)."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,rdMolDescriptors,rdFingerprintGenerator
assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
root=Path(__file__).resolve().parents[2];fixture=root/'crates/smiles/tests/fixtures'
sources=[]
for filename in ['structural_contexts.tsv','coordination_contexts.tsv','descriptor_cleanup.tsv']:
 sources.extend(line.split('\t')[0] for line in (fixture/filename).read_text().splitlines() if not line.startswith('#'))
sources += ['[H]','[2H]','[H][H]','[2H][H]','[H-].[Na+]','[H]C([H])([H])[H]',
 '[H]C([2H])([3H])[H]','[13CH4]','[14CH4]','[12CH4]','[0CH4]',
 '[15NH3]','[18OH2]','CCO.CCO','CCO.OCC','CCO.[H][H]',
 'CC(C)C','C1CC1','C1CCC1','C1CCCC1','C1CCCCC1','C1CCCCCC1',
 'C1CC2CCC1C2','C12C3C4C1C5C2C3C45','[Na+].[Cl-]',
 'F/C=C/F','F/C=C\\F','N[C@H](C)O','N[C@@H](C)O']
sources += ['O=C1OC[C@H]2OC(=O)C3=C(OC(=O)C4=CC(O)=C(O)C(O)=C4C4=C(O)C(O)=C(O)C=C4C(=O)O3)[C@@H]2OC(=O)C2=CC(O)=C(O)C(O)=C2C2=C(O)C(O)=C(O)C=C12']
sources=list(dict.fromkeys(sources));p=Chem.SmilesParserParams();p.removeHs=False
values=['# Original source; registered density radius1/2/3; removeHs=False']
keys=['# Original source; atom connectivity invariants; sparse keys at radius0/1/2/3']
contexts=['# source; atoms(symbol,H,charge,arom,noImplicit,isotope); bonds(a,b,RDKitType); symmetric bond rings']
for source in sources:
 m=Chem.MolFromSmiles(source,p);assert m is not None,source
 densities=[getattr(Descriptors,'FpDensityMorgan'+str(r))(m) for r in range(1,4)]
 values.append(source+'\t'+'\t'.join(format(x,'.17g') for x in densities))
 fingerprints=[rdFingerprintGenerator.GetMorganGenerator(r).GetSparseCountFingerprint(m).GetNonzeroElements() for r in range(4)]
 heavy=m.GetNumHeavyAtoms()
 assert densities == [len(x)/heavy if heavy else 0. for x in fingerprints[1:]]
 inv=rdMolDescriptors.GetConnectivityInvariants(m)
 keys.append(source+'\t'+','.join(map(str,inv))+'\t'+'\t'.join(','.join(map(str,sorted(x))) for x in fingerprints))
 atoms=';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())},{a.GetIsotope()}' for a in m.GetAtoms())
 bonds=';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{int(b.GetBondType())}' for b in m.GetBonds())
 rings=';'.join(','.join(map(str,r)) for r in m.GetRingInfo().BondRings())
 contexts.append(f'{source}\t{atoms}\t{bonds}\t{rings}')
for name,lines in [('descriptor_morgan_density.tsv',values),('morgan_sparse_keys.tsv',keys),('morgan_contexts.tsv',contexts)]:
 (fixture/name).write_text('\n'.join(lines)+'\n')
meta=dict(rdkit=rdkit.__version__,molecules=len(sources),scalar_values=3*len(sources),key_radii=[0,1,2,3],removeHs=False,includeChirality=False,includeRedundantEnvironments=False,fp='GetSparseCountFingerprint; unique keys, no bit-vector folding',native_acceptance_filter=False)
(root/'tasks/descriptors/morgan-density-controls.json').write_text(json.dumps(meta,indent=2)+'\n')
print(len(sources),'independent sources;',3*len(sources),'registered values and radius0..3 sparse key sets')
