"""Explicit-H and strict-query boundary controls; pinned registered oracle."""
from pathlib import Path
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')
root=Path(__file__).resolve().parents[2];fixture=root/'crates/smiles/tests/fixtures'
sources=[l.split('\t')[0] for l in (fixture/'descriptor_phase1.tsv').read_text().splitlines() if not l.startswith('#')]
sources += [
 '[H]OCC','[2H]OCC','[H]C(=O)NCC','[2H]C(=O)NCC','[H]OC(=O)C',
 '[H]N(C)CC','CC[H-]CC','CC[2H-]CC','CC(C)[H-]C(C)C','[2H]ON(=O)=O',
 '[H]C([H])([H])C([H])([H])[H]','CC([H])([H])C',
 'CC([H])([H])C([H])([H])C','[H]C([H])([H])C(C)(C)C',
 '[H]C([H])([H])C(C)(C([H])([H])[H])CC',
 'CC(C)(C)CC','FC(F)(F)CC','ClC(Cl)(Cl)CC','BrC(Br)(Br)CC','FC(F)(Cl)CC',
 'CC(=O)NCC','CCC(=O)OCC','CCC(=O)SCC','CCC(=O)S[H]',
 'CC(=[NH2+])NCC','CC(=[NH2+])OC','C=NCC',
 'c1ccccc1:c1ccccc1','c1ccccc1-c1ccccc1','c1ccccc1-O-c1ccccc1',
 'c1ccccc1N(C)C','S=C(NCC)NCC','CC(=NCC)NCC',
]
sources=list(dict.fromkeys(sources));p=Chem.SmilesParserParams();p.removeHs=False
lines=['# RDKit2026.03.6 Strict NumRotatableBonds; explicit H vertices retained']
for smiles in sources:
 m=Chem.MolFromSmiles(smiles,p);assert m is not None,smiles
 lines.append(f'{smiles}\t{Descriptors.NumRotatableBonds(m)}')
(fixture/'descriptor_rotors.tsv').write_text('\n'.join(lines)+'\n')
print(len(sources),'independent strict rotor values')
