"""Pinned RDKit Gasteiger parameter facts; no per-molecule expected answers."""
import re
from pathlib import Path
from urllib.request import urlopen
from rdkit import rdBase
from rdkit.Chem import MolSurf
assert rdBase.rdkitVersion=='2026.03.6'
url='https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/PartialCharges/GasteigerParams.cpp'
source=urlopen(url).read().decode()
rows=re.findall(r'^\s*"?([A-Z][a-z]?)\s+(\*|sp3|sp2|sp|so2|so)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',source,re.M)
assert len(rows)==32 and len({(r[0],r[1]) for r in rows})==32,len(rows)
s='//! Pinned RDKit2026.03.6 GasteigerParams.cpp facts; BSD-3-Clause notice applies.\n'
s+='pub(crate) fn parameters(symbol:&str,mode:&str)->[f64;3] {match (symbol,mode) {\n'
for element,mode,a,b,c in rows:s+=f'    ("{element}","{mode}")=>[{float(a)},{float(b)},{float(c)}],\n'
s+='    _=>[0.;3],\n}}\n'
s+=f'pub(crate) const BINS:[f64;{len(MolSurf.chgBins)}]={list(map(float,MolSurf.chgBins))!r};\n'
Path('crates/smiles/src/charge_data.rs').write_text(s)
print('pinned parameter rows',len(rows))
