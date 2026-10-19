"""Independent pinned SPS pipeline; preserve source and every per-atom factor.

No native acceptance filter. No-heavy-atom division by zero is the sole named
undefined domain; all other oracle exceptions fail generation.
"""
from pathlib import Path
import hashlib
import json
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import SpacialScore

assert rdkit.__version__ == '2026.03.6'
assert Chem.GetUseLegacyStereoPerception() is True
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
sources = json.loads((root/'tasks/descriptors/stereo-controls.json').read_text())['source_smiles']
sources += ['C', '[H]', '[2H]', '[H][H]', '[H]C([H])([H])[H]',
            '[H]C([2H])([3H])F', '[H]C([H])(F)Cl', 'C[C@H]O',
            'C[C@H]C', '[C@H]', 'CO[C@H](O)O', 'CC(=O)C',
            'C1CC(C)C(C)C(C)C1', 'C1C[C@H](C)C(C)[C@H](C)C1',
            'C1C[C@H](C)[C@H](C)[C@H](C)C1',
            'CC1CC1C', 'CC1CC(C)C1', 'CC1CCC(C)CC1', 'CC1CCCCC(C)CCC1',
            'C1CCC2CCCCC2C1', 'C1CC2CCC1C2', 'C1CC2CC3CCC(C1)C23',
            'C1CCC(CC1)=C1CCCCC1', 'CC1CCC(=C(C)CC)CC1',
            'CC1CCC(=C(C)CC)C(C)C1', 'C1CCN2CCCC2C1',
            'C1=CCCCC1', 'C1=CCCCCC1', 'C1=CCCCCCC1', 'C1/CCCCCC=C/1',
            'C1/CCCCCC=C\\1', 'CC=C=C(C)Cl', '[H]C=C([H])C',
            '[2H]C=C([H])C', 'N->B', 'C[N+](C)(C)C', 'C[PH]CC',
            '[Be-2](F)(Cl)(Br)I', '[Si](F)(Cl)(Br)I', 'C[P+](F)(Cl)Br']
params = Chem.SmilesParserParams(); params.removeHs = False
# Symmetric substituted rings and atom-order variations stress dependent stereo.
for size in range(3, 13):
    for separation in range(1, size // 2 + 1):
        sources.append('CC1'+'C'*(separation-1)+'C(C)'+'C'*(size-separation-1)+'1')
for source in list(sources[-65:]):
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None, source
    sources += list(Chem.MolToRandomSmilesVect(mol, 3, randomSeed=1729))
# A separately retained regression selection from the difficult415 audit.
# Recalculate every oracle value; never use native values as reference.
regressions = json.loads((root/'tasks/descriptors/sps-regressions.json').read_text())['inputs']
sources += [r['smiles'] for r in regressions]
# First-ring caching is order-sensitive: validate macrocycle/ether/sulfur
# preparation after independent atom-order permutations, in Kekule notation.
for index in [1, 11, 12]:
    mol = Chem.MolFromSmiles(regressions[index]['smiles'], params)
    for seed in [1729, 2101]:
        sources += list(Chem.MolToRandomSmilesVect(mol, 3, randomSeed=seed, kekuleSmiles=True))
sources = list(dict.fromkeys(sources))
scalars = ['# source; SPS; RDKit2026.03.6; null=no-heavy-atom division by zero']
vectors = ['# source; modern tetrahedral mask; double-stereo non-NONE mask SORTED BY ENDPOINTS; hybrid factors; ring factors; degrees']
records = []
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None and mol.GetNumAtoms(), source
    if not mol.GetNumHeavyAtoms():
        try: SpacialScore.SPS(mol)
        except ZeroDivisionError: pass
        else: raise AssertionError('expected no-heavy SPS domain error')
        expected = None
        detail = SpacialScore._SpacialScore(mol, normalize=False)
    else:
        detail = SpacialScore._SpacialScore(mol)
        expected = SpacialScore.SPS(mol)
        assert expected == detail.score
    atom_mask = [int(i in detail.chiral_idxs) for i in range(mol.GetNumAtoms())]
    bond_mask = [int(b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() != Chem.BondStereo.STEREONONE)
                 for b in detail.mol.GetBonds()]
    bond_endpoints = [(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx())) for b in detail.mol.GetBonds()]
    sorted_bond_mask = [mask for edge,mask in sorted(zip(bond_endpoints,bond_mask))]
    hyb = list(detail.hyb_score.values()); ring = list(detail.ring_score.values())
    degree = [a.GetDegree() for a in detail.mol.GetAtoms()]
    scalars.append(source+'\t'+('null' if expected is None else repr(expected)))
    vectors.append(source+'\t'+'\t'.join(','.join(map(str, v)) for v in [atom_mask,sorted_bond_mask,hyb,ring,degree]))
    records.append(dict(source=source,expected=expected,modern_atom_mask=atom_mask,
                        double_stereo_mask=bond_mask,bond_endpoints=bond_endpoints,hybrid=hyb,ring=ring,degree=degree))
fixture = root/'crates/smiles/tests/fixtures/descriptor_sps.tsv'
fixture.write_text('\n'.join(scalars)+'\n')
(root/'crates/smiles/tests/fixtures/sps_contexts.tsv').write_text('\n'.join(vectors)+'\n')
(root/'tasks/descriptors/sps-controls.json').write_text(json.dumps(dict(
    rdkit=rdkit.__version__,legacy_default=True,sources=len(sources),native_filter=False,
    undefined_policy='no heavy atoms: explicit null; oracle ZeroDivisionError asserted',
    fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest(),records=records),indent=2)+'\n')
print(len(sources),'independent sources;',sum(r['expected'] is None for r in records),'named undefined SPS')
