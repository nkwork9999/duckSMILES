"""Pinned H/mass controls for partial aromatic bond valence; no native filter."""
from pathlib import Path
import csv
import gzip
import json
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from compare_catalog import profile_reference

assert rdkit.__version__ == '2026.03.6'
RDLogger.DisableLog('rdApp.*')
root = Path(__file__).resolve().parents[2]
fixture = root / 'crates/smiles/tests/fixtures'
old = root / 'tasks/evidence/catalog232-coordination/chebi100k-v20/differences.jsonl.gz'
residuals = [json.loads(line) for line in gzip.open(old, 'rt')
             if 'formula' in json.loads(line)['differences']]
assert len(residuals) == 5
sources = [r['smiles'] for r in residuals]
sources += ['C:C', 'C:N', 'N:C', 'C:O', 'O:C', 'C:S', 'S:C',
            'C:N:C', 'C:O:C', 'C:S:C', 'C:N(C):C', 'C:[N+](C):C',
            'CC(=O):N', 'N:C(C)=O', '[13CH3]:N', '[H]C:N',
            'C1:C:C:C:C:C#1', 'C1:C:C:C:C:C=1', 'C1=CC#CC=C1',
            'C1:C:C:C:C:C:1', 'c1ccccc1', 'C1=CC=CC=C1',
            'C1=CNC=C1', '[nH]1cccc1', '[n+]1ccccc1',
            'CCO', 'O', 'N', 'P', 'CP(C)(C)C',
            '[H]O[H]', '[2H]N([2H])[2H]', '[NH4+].[Cl-]',
            '[H-][Na+]', 'C[H-]C', 'CN(=O)=O', '[OH2][Mg+2]']
sources = list(dict.fromkeys(sources))
params = Chem.SmilesParserParams()
params.removeHs = False
fields = ['MolWt', 'ExactMolWt', 'ds_ExplicitHydrogens',
          'ds_ImplicitHydrogens', 'ds_TotalHydrogens']
full_sources = {'C:C', 'C:N', 'N:C', 'C:O', 'O:C', 'C:S', 'S:C',
                'CC(=O):N', 'N:C(C)=O', '[13CH3]:N', '[H]C:N',
                'C:N:C', 'C:O:C', 'C:S:C', 'C:N(C):C', 'C:[N+](C):C',
                'C1:C:C:C:C:C#1', 'C1=CC#CC=C1'}
catalog = json.loads((root/'tasks/descriptors/phase1-catalog.json').read_text())['descriptors']
functions = dict(Descriptors._descList)
lines = ['# source; descriptor/formula; pinned RDKit2026.03.6 sanitize=True removeHs=False']
contexts = ['# source; atom total H(excludes neighboring H vertices); atom aromatic flag; bonds(a,b,type,aromaticflag); full reference diagnostics']
for source in sources:
    mol = Chem.MolFromSmiles(source, params)
    assert mol is not None, source
    atoms = list(mol.GetAtoms())
    explicit = sum(a.GetNumExplicitHs() for a in atoms) + sum(a.GetAtomicNum() == 1 for a in atoms)
    implicit = sum(a.GetNumImplicitHs() for a in atoms)
    values = [Descriptors.MolWt(mol), Descriptors.ExactMolWt(mol), explicit,
              implicit, explicit + implicit]
    lines.append(f'{source}\t@formula\t{rdMolDescriptors.CalcMolFormula(mol)}')
    for key, value in zip(fields, values):
        lines.append(f'{source}\t{key}\t{value:.17g}')
    if source in full_sources:
        profile = profile_reference(mol, include_extra=True)
        for spec in catalog:
            key = spec['id']
            if key in fields:
                continue
            value = ((functions[key](mol, force=False) if 'EStateIndex' in key else functions[key](mol))
                     if spec['origin'] == 'rdkit-standard' else profile[key])
            lines.append(f'{source}\t{key}\t{float(value):.17g}')
    hydrogens = ','.join(str(a.GetTotalNumHs()) for a in atoms)
    aromatic = ','.join(str(int(a.GetIsAromatic())) for a in atoms)
    bonds = ';'.join(f'{b.GetBeginAtomIdx()},{b.GetEndAtomIdx()},{int(b.GetBondType())},{int(b.GetIsAromatic())}' for b in mol.GetBonds())
    contexts.append(f'{source}\t{hydrogens}\t{aromatic}\t{bonds}')
(fixture / 'descriptor_valence_h.tsv').write_text('\n'.join(lines) + '\n')
(fixture / 'valence_h_contexts.tsv').write_text('\n'.join(contexts) + '\n')
meta = dict(rdkit=rdkit.__version__, molecules=len(sources),
            source_smiles=sources, formula_residual_original_indices=[r['index'] for r in residuals],
            selected_descriptors=fields, formula_values=len(sources), scalar_values=len(sources)*len(fields),
            additional_full_catalog_sources=sorted(full_sources), catalog_count=len(catalog),
            total_scalar_values=len(sources)*len(fields)+len(full_sources)*(len(catalog)-len(fields)),
            native_acceptance_filter=False,
            scope='Five retained formula residuals plus partial-aromatic/explicit-H/charge/ordinary controls. Full atom/bond diagnostics are not asserted as a full aromaticity correction.')
(root / 'tasks/descriptors/valence-h-controls.json').write_text(json.dumps(meta, indent=2)+'\n')
print(len(sources), 'sources;', len(sources)*len(fields), 'H/mass values plus formulas')
