"""Pinned default legacy stereocounts, with original parser state retained."""
from pathlib import Path
import json
import hashlib
import csv
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, rdqueries

assert rdkit.__version__ == '2026.03.6'
assert Chem.GetUseLegacyStereoPerception() is True
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
sources = json.loads((root/'tasks/descriptors/extra50-controls.json').read_text())['source_smiles']
sources += ['CC(O)Cl', 'C[C@H](O)Cl', 'C[C@@H](O)Cl', 'C[C@H](C)O',
            'CC(F)(Cl)C', '[CH3:1]C(F)(Cl)[CH3:2]', '[CH3:0]C(F)(Cl)C',
            '[CH3:1]C(F)(Cl)[CH3:1025]', '[CH3:1][C@](F)(Cl)[CH3:2]',
            '[2H]C(F)(Cl)Br', 'F[C@]([2H])(Cl)Br', '[H][C@]([H])(F)Cl',
            'CN(CC)CCC', 'C[N+](CC)(CCC)CCCC', 'C[N@+](CC)(CCC)CCCC',
            'CP(CC)CCC', 'C[PH]CC', 'C[P@H]CC', 'C[AsH]CC',
            'C[S@](=O)CC', 'C[S+](CC)CCC', 'C[Se@](=O)CC',
            'CN1CC1', 'CN1C(C)C1', 'C1CCN2CCCC2C1', 'C1CCN2CCCCC2C1',
            'CC1CCC(O)CC1', 'C[C@H]1CC[C@@H](O)CC1', 'N[C@H]1CC[C@@H](O)CC1',
            'C[C@H](O)[C@H](O)C', 'C[C@H](O)[C@@H](O)C',
            'CC(O)C(O)C', 'CC(O)[C@@H](O)C',
            'C[C@](F)([C@H](O)Cl)[C@@H](O)Cl',
            'CC=CC', 'C/C=C/C', 'C/C=C\\C', 'CC=C(C)C',
            'C1CCCC=CC1', 'C1CCCC=CCC1', 'C1CCCC=CCCC1', 'CC=C=CC',
            '[H]C([H])([H])[H]',
            'CC(/C=C/C)(/C=C\\C)O', 'CC(/C=C/C)(/C=C/C)O',
            '[C@H](F)(Cl)Br', '[P@H](C)CC', 'CC[P@H]C', '[S@](=O)(C)CC']
sources += [r['smiles'] for r in json.loads((root/'tasks/descriptors/full100k-stereo-regressions.json').read_text())['inputs']]
# Predeclared104 original ring failures, retained in full after the SSSR port.
sources += [r['smiles'] for r in csv.DictReader((root/'tasks/evidence/catalog289-sps/ring-count104.csv').open())]
sources = list(dict.fromkeys(sources))
params = Chem.SmilesParserParams(); params.removeHs = False
functions = [('NumAtomStereoCenters', rdMolDescriptors.CalcNumAtomStereoCenters),
             ('NumUnspecifiedAtomStereoCenters', rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters)]
scalars = ['# source; registered id; pinned RDKit2026.03.6 default legacy, sanitize=True removeHs=False']
contexts = []
atom_vectors = ['# source; per-atom _ChiralityPossible; per-atom legacy CIP code (dash=absent)']
bridgehead_vectors = ['# source; per-atom rdqueries.IsBridgeheadQueryAtom (not CalcNumBridgeheadAtoms)']
bridgehead_query = rdqueries.IsBridgeheadQueryAtom()
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None and mol.GetNumAtoms(), source
    bridgehead_vectors.append(source+'\t'+','.join(str(int(bridgehead_query.Match(a))) for a in mol.GetAtoms()))
    for key, fn in functions:
        scalars.append(f'{source}\t{key}\t{fn(mol)}')
    contexts.append(dict(source=source,stereochem_done=mol.HasProp('_StereochemDone'),
                         atoms=[dict(z=a.GetAtomicNum(),isotope=a.GetIsotope(),
                                     map_present=a.HasProp('molAtomMapNumber'),map=a.GetAtomMapNum(),
                                     h=a.GetTotalNumHs(),aromatic=a.GetIsAromatic(),charge=a.GetFormalCharge(),
                                     tag=int(a.GetChiralTag()),possible=a.HasProp('_ChiralityPossible'),
                                     cip=a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None,
                                     rank=int(a.GetProp('_CIPRank')) if a.HasProp('_CIPRank') else None)
                                for a in mol.GetAtoms()]))
    atom_vectors.append(source+'\t'+','.join(str(a.HasProp('_ChiralityPossible')) for a in mol.GetAtoms())
                        +'\t'+','.join(a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else '-' for a in mol.GetAtoms()))
fixture = root/'crates/smiles/tests/fixtures/descriptor_stereo_counts.tsv'
fixture.write_text('\n'.join(scalars)+'\n')
(root/'crates/smiles/tests/fixtures/stereo_legacy_contexts.jsonl').write_text('\n'.join(map(json.dumps,contexts))+'\n')
(root/'crates/smiles/tests/fixtures/stereo_legacy_atoms.tsv').write_text('\n'.join(atom_vectors)+'\n')
(root/'crates/smiles/tests/fixtures/stereo_bridgehead_query.tsv').write_text('\n'.join(bridgehead_vectors)+'\n')
(root/'tasks/descriptors/stereo-controls.json').write_text(json.dumps(dict(
    rdkit=rdkit.__version__,legacy_default=True,source_count=len(sources),scalar_count=2*len(sources),
    source_smiles=sources,descriptor_ids=[k for k,f in functions],native_filter=False,
    fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest()),indent=2)+'\n')
print(len(sources),'sources;',2*len(sources),'independent legacy stereocounts')
