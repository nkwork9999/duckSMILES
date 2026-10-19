"""Independent input-order atom H and directed dative-edge expectations."""
from pathlib import Path
import json
import rdkit
from rdkit import Chem, RDLogger
assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
p = Chem.SmilesParserParams(); p.removeHs = False
sources = ['N->[Cu+2]', '[Cu+2]<-N', 'N(->[Cu+2])C',
           '[NH3]->[Cu+2]<-[NH3]', '[13NH3]->[Cu+2]',
           '[Cu+2](<-N)(<-O)<-S', 'N->B', 'B<-N',
           '[Cu+2]1<-NCCN->1', 'N->1CCN[Cu+2]1',
           '[Cu+2]<-1NCCN1', '[Cu+2]1NCCN->1', 'N1CCN[Cu+2]<-1']
lines = ['# source; formula; atom attached-H vector; sorted directed dative edges; ordinary ring count']
from rdkit.Chem import rdMolDescriptors
for s in sources:
    m = Chem.MolFromSmiles(s, p); assert m is not None, s
    h = ','.join(str(a.GetTotalNumHs()) for a in m.GetAtoms())
    edges = sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE)
    lines.append('\t'.join([s, rdMolDescriptors.CalcMolFormula(m), h, ';'.join(f'{a},{b}' for a,b in edges), str(rdMolDescriptors.CalcNumRings(m))]))
root=Path(__file__).resolve().parents[2]
(root/'crates/smiles/tests/fixtures/dative_input.tsv').write_text('\n'.join(lines)+'\n')
(root/'tasks/descriptors/dative-input-controls.json').write_text(json.dumps(dict(rdkit=rdkit.__version__, sources=sources, native_filter=False),indent=2)+'\n')
print(len(sources), 'directed input controls')
