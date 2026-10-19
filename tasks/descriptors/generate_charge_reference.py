"""Independent controls for the next Gasteiger/PEOE tranche (not shipped yet)."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors,rdPartialCharges,MolSurf
assert rdkit.__version__=='2026.03.6'
root=Path(__file__).resolve().parents[2]
smiles=list(dict.fromkeys(l.split('\t')[0] for l in (root/'crates/smiles/tests/fixtures/graph_contexts.tsv').read_text().splitlines() if not l.startswith('#')))
smiles+=['NC(=[NH2+])N','[O-]C(=O)C','N=C(N)N','[O-][N+](=O)c1ccccc1','C[S+](C)[O-]','[O-]S(=O)(=O)C','CS(C)=O','CS(C)(=O)=O','[SiH4]','B(O)(O)O','[B-](F)(F)(F)F','[Li+]','[Mg+2]','[Na+]Cl','C[Al](C)C','[BeH2]','[PH4+]','C[P+](C)(C)C','[N-]=[N+]=N']
smiles+=['[Fe]C.CCO','CCO.[Fe]C','[Og][Og].CCO','CCO.[Og][Og]']
smiles=list(dict.fromkeys(smiles));params=Chem.SmilesParserParams();params.removeHs=False
ids=['NumRadicalElectrons','MaxPartialCharge','MinPartialCharge','MaxAbsPartialCharge','MinAbsPartialCharge']+[f'PEOE_VSA{i}' for i in range(1,15)]
functions=dict(Descriptors._descList)
rows=[dict(id=i,title=i,group='Composition' if i=='NumRadicalElectrons' else 'Partial charge' if 'Charge' in i else 'Charge surface',unit='' if i=='NumRadicalElectrons' else 'e' if 'Charge' in i else 'Å²',digits=0 if i=='NumRadicalElectrons' else 4,origin='rdkit-standard',verification='pending-expanded-validation') for i in ids]
scalar=['# RDKit2026.03.6, removeHs=False, iterations=12']
atomic=['# SMILES, atom index, charge, attached-H total charge, hybridization, radical electrons']
conjugation=['# SMILES, bond index, conjugated']
contexts=['# RDKit 2026.03.6 sanitized contexts; same schema as graph_contexts.tsv']
for s in smiles:
 m=Chem.MolFromSmiles(s,params);assert m is not None,s
 rdPartialCharges.ComputeGasteigerCharges(m,nIter=12,throwOnParamFailure=False)
 atoms=';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())}' for a in m.GetAtoms())
 bonds=';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{b.GetBondTypeAsDouble()}' for b in m.GetBonds())
 contexts.append(f'{s}\t{atoms}\t{bonds}')
 for i in ids:scalar.append(f'{s}\t{i}\t{functions[i](m):.17g}')
 for a in m.GetAtoms():
  atomic.append(f'{s}\t{a.GetIdx()}\t{float(a.GetProp("_GasteigerCharge")):.17g}\t{float(a.GetProp("_GasteigerHCharge")):.17g}\t{a.GetHybridization().name}\t{a.GetNumRadicalElectrons()}')
 for b in m.GetBonds():conjugation.append(f'{s}\t{b.GetIdx()}\t{int(b.GetIsConjugated())}')
for name,lines in [('descriptor_charge.tsv',scalar),('charge_atomic.tsv',atomic),('charge_conjugation.tsv',conjugation),('charge_contexts.tsv',contexts)]:
 (root/'crates/smiles/tests/fixtures'/name).write_text('\n'.join(lines)+'\n')
positive={i for i in ids if any(l.split('\t')[1]==i and float(l.split('\t')[2])>0 for l in scalar[1:])}
report=dict(positive_ids=sorted(positive),status='native-model oracle-ready; integration validation in progress',rdkit=rdkit.__version__,descriptors=rows,molecules=len(smiles),scalar_values=len(scalar)-1,atomic_values=len(atomic)-1,conjugation_values=len(conjugation)-1,bins=list(MolSurf.chgBins))
(root/'tasks/descriptors/charge-definitions.json').write_text(json.dumps(report,indent=2)+'\n')
print(report['status'],report['molecules'],report['scalar_values'],report['atomic_values'])
