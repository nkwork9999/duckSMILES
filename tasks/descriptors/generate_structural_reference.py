"""Pinned, independent controls for three registered structural counts."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
fixture = root / 'crates/smiles/tests/fixtures'
sources = [line.split('\t')[0] for line in (fixture/'bcut_original_contexts.tsv').read_text().splitlines() if not line.startswith('#')]
sources += [
    'NC=O', 'NC(=O)N', 'CN(C)C(=O)N(C)C', 'NC(=O)NC(=O)N',
    'O=C1NCCC1', 'O=C1NC(=O)CC1', 'NC(=S)N', 'N=C(N)N',
    'O=C(O)C', 'O=C(Cl)C', 'O=CNC(=O)C', 'O=C(N)O',
    'C[N+]1=CC=CC=C1', 'CC(=O)[NH3+]', 'CC(=O)[N-]C',
    'CC(=O)[15NH2]', 'CC(=O)N[2H]', 'CC(=O)N.CC(=O)N',
    'C[O+]=C(N)C', 'C1=[O+]CCC1', 'NC1=[O+]CCC1', 'O=C1N=CC=CC1',
    'C1CC1', 'C1CCC2(CC1)CCC2', 'C1CCC2(CC1)CCCC2',
    'C1CC2(C1)CC2', 'C1COC2(CC1)OCCO2', 'C1CC2(C1)CC1(CC2)CCC1',
    'C1CCC2(CC1)CC1(CC2)CCC1', 'C1CC2CCC1C2', 'C1CC2CCC1CC2',
    'C1CCC2CCCCC2C1', 'C1CC2CCC3CC(C1)C23', 'C1C2CC3CC1CC(C2)C3',
    'C12C3C4C1C5C2C3C45', 'C1CC2CCC1C2.C1CCC2(CC1)CCC2',
    'C1CCC2CCCCC2C1.C1CCC2CCCCC2C1', 'O=C1NCCC12CCCCC2',
    'C1CC2CC3CC1CC(C2)C3', 'C1C2C3C4C1C5C2C3C45',
]
sources = list(dict.fromkeys(sources))
ids = [('NumAmideBonds', 'Amide bonds', 'Functional groups'),
       ('NumSpiroAtoms', 'Spiro atoms', 'Rings'),
       ('NumBridgeheadAtoms', 'Bridgehead atoms', 'Rings')]
functions = dict(Descriptors._descList)
params = Chem.SmilesParserParams(); params.removeHs = False
scalars = ['# RDKit2026.03.6 registered counts; original SMILES/removeHs=False']
contexts = ['# Sanitized RDKit graph: source, atoms(symbol,H,charge,aromatic,noImplicit,isotope), bonds(a,b,order), symmetric bond rings']
positives = {id: 0 for id, _, _ in ids}
for smiles in sources:
    mol = Chem.MolFromSmiles(smiles, params)
    assert mol is not None, smiles
    for id, _, _ in ids:
        value = functions[id](mol)
        assert value >= 0 and int(value) == value
        scalars.append(f'{smiles}\t{id}\t{value}')
        positives[id] += value > 0
    atoms = ';'.join(f'{a.GetSymbol()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())},{int(a.GetNoImplicit())},{a.GetIsotope()}' for a in mol.GetAtoms())
    bonds = ';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{b.GetBondTypeAsDouble()}' for b in mol.GetBonds())
    rings = ';'.join(','.join(map(str, ring)) for ring in mol.GetRingInfo().BondRings())
    contexts.append(f'{smiles}\t{atoms}\t{bonds}\t{rings}')
(fixture/'descriptor_structural.tsv').write_text('\n'.join(scalars)+'\n')
(fixture/'structural_contexts.tsv').write_text('\n'.join(contexts)+'\n')
specs = [dict(id=id,title=title,group=group,unit='',digits=0,origin='rdkit-standard',verification='pending-expanded-validation') for id,title,group in ids]
(root/'tasks/descriptors/structural-definitions.json').write_text(json.dumps(dict(rdkit=rdkit.__version__,descriptors=specs,molecules=len(sources),scalar_values=len(scalars)-1,positive_controls=positives),indent=2)+'\n')
print(len(sources), 'independent contexts;', len(scalars)-1, 'scalar counts;', positives)
