"""Pinned QED properties and alert controls, independent of native acceptance."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import QED

assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
fixture = root / 'crates/smiles/tests/fixtures'
sources = [line.split('\t')[0] for line in
           (fixture / 'descriptor_morgan_density.tsv').read_text().splitlines()
           if not line.startswith('#')]
sources += ['N1C2=CC=CC=C2SC3=C1C=CC4=C3C=CC=C4',
            'OC1=C(O)C=C2C(=C1)OC3=CC(=O)C(=CC3=C2C4=CC=CC=C4)O',
            'CC(C)C1=CC2=C(C)C=CC2=C(C)C=C1',
            '[H]O[H]', '[H]N([H])C', '[2H]N([H])C',
            '[H][B-]([H])([H])[H]', '[H][Cu][H]',
            'F.F.F.F', 'COC(=O)C.COC(=O)C.COC(=O)C',
            'N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O',
            'CNC(=NCCSCc1nc[nH]c1C)NC#N',
            'CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12']
sources += ['[Hg]', '[Hf]', '[Ho]', '[Hs]', '[He]', '[Ru]', '[Rh]', '[Re]',
            '[Rn]', '[Ra]', '[Rb]', '[Dy]', 'c1cc[se]c1', 'c1cc[asH]c1',
            '[H]C([H])([H])C([H])([H])[H]',
            '[H]OOC', '[H]C(=O)C', '[H]/C=C/F', '[H]N=C(C)C',
            '[H-][Cu+]', 'C[H-]C', '[H:7]OC']
sources += json.loads((root / 'tasks/descriptors/aromatic-candidate-controls.json').read_text())['source_smiles']
sources = list(dict.fromkeys(sources))
params = Chem.SmilesParserParams()
params.removeHs = False
lines = ['# source; default QED; eight properties; 11 acceptor counts; 116 alert flags; Chem.RemoveHs default']
hydrogen_contexts = ['# source; retained atoms(symbol,isotope,totalH,charge,aromatic); bonds(a,b,RDKitType) after default RemoveHs']
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None, source
    h_removed = Chem.RemoveHs(mol)
    props = QED.properties(mol)
    acceptors = [len(h_removed.GetSubstructMatches(p)) for p in QED.Acceptors]
    alerts = [int(h_removed.HasSubstructMatch(p)) for p in QED.StructuralAlerts]
    assert props.HBA == sum(acceptors) and props.ALERTS == sum(alerts)
    lines.append('\t'.join([source, format(QED.qed(mol), '.17g'),
                           ','.join(format(v, '.17g') for v in props),
                           ','.join(map(str, acceptors)),
                           ','.join(map(str, alerts))]))
    if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        atoms = ';'.join(f'{a.GetSymbol()},{a.GetIsotope()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())}' for a in h_removed.GetAtoms())
        bonds = ';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{int(b.GetBondType())}' for b in h_removed.GetBonds())
        hydrogen_contexts.append(f'{source}\t{atoms}\t{bonds}')
(fixture / 'descriptor_qed.tsv').write_text('\n'.join(lines) + '\n')
(fixture / 'qed_hydrogen_contexts.tsv').write_text('\n'.join(hydrogen_contexts) + '\n')
metadata = dict(rdkit=rdkit.__version__, molecules=len(sources),
                property_order=list(QED.QEDproperties._fields),
                acceptor_patterns=QED.AcceptorSmarts,
                alert_patterns=QED.StructuralAlertSmarts,
                removeHs='Chem.RemoveHs default within QED.properties',
                source_removeHs=False, native_acceptance_filter=False,
                hydrogen_contexts=len(hydrogen_contexts)-1,
                source='https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/rdkit/Chem/QED.py')
(root / 'tasks/descriptors/qed-controls.json').write_text(json.dumps(metadata, indent=2) + '\n')
print(len(sources), 'independent QED/property/acceptor/alert controls')

# Element-name semantics need their own checks: e.g. [Te] is atomic number52,
# while [S] also requires aliphatic sulfur. Aromatic [te] is a single token.
queries = ['[Se]', '[se]', '[Te]', '[te]', '[As]', '[as]', '[S]', '[s]',
           '[P]', '[p]', '[C]', '[c]', '[N]', '[n]', '[O]', '[o]',
           '[Hf]', '[Hg]', '[Ho]', '[He]', '[Hs]', '[Ru]', '[Rb]', '[Ra]',
           '[Rn]', '[Re]', '[Rh]', '[Dy]', '[Cl]']
token_sources = ['c1cc[se]c1', 'c1cc[te]c1', 'c1cc[asH]c1', 'c1cc[pH]c1',
                 'c1ccsc1', 'c1ccoc1', 'c1ccncc1', 'CS', 'CP', 'CCl',
                 '[Hf]', '[Hg]', '[Ho]', '[He]', '[Hs]', '[Ru]', '[Rb]',
                 '[Ra]', '[Rn]', '[Re]', '[Rh]', '[Dy]']
tokens = ['# pinned RDKit2026.03.6: query; original source; matching atom indices, removeHs=False']
for query in queries:
    pattern = Chem.MolFromSmarts(query)
    assert pattern is not None and pattern.GetNumAtoms() == 1
    for source in token_sources:
        mol = Chem.MolFromSmiles(source, params)
        assert mol is not None
        hits = sorted({match[0] for match in mol.GetSubstructMatches(pattern)})
        tokens.append(query + '\t' + source + '\t' + ','.join(map(str, hits)))
(fixture / 'smarts_element_tokens.tsv').write_text('\n'.join(tokens) + '\n')
print(len(queries), 'queries;', len(token_sources), 'sources;', len(tokens)-1, 'element-token cells')
