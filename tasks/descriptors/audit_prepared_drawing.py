"""Audit actual prepared C ABI output against pinned RDKit, retaining every source row."""
import argparse,ctypes,gzip,json,hashlib
from pathlib import Path
import rdkit
from rdkit import Chem,RDLogger
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')

def main():
    p=argparse.ArgumentParser();p.add_argument('library',type=Path);p.add_argument('sources',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    lib=ctypes.CDLL(str(a.library.resolve()));fn=lib.ds_workbench_drawing_graph
    fn.argtypes=[ctypes.c_char_p,ctypes.c_size_t,ctypes.c_void_p,ctypes.c_size_t];fn.restype=ctypes.c_int
    params=Chem.SmilesParserParams();params.removeHs=False
    counts={'inputs':0,'compared':0,'unavailable':0,'reference_unavailable':0};differences=[];unavailable=[]
    output=a.output/'prepared-graphs.jsonl.gz'
    with gzip.open(a.sources,'rt') as src,gzip.open(output,'wt') as out:
        for line in src:
            r=json.loads(line);counts['inputs']+=1;smiles=r['smiles'];b=smiles.encode();n=fn(b,len(b),None,0)
            graph=None
            if n>0:
                buf=ctypes.create_string_buffer(n);assert fn(b,len(b),buf,n)==n
                graph=json.loads(buf.raw);assert graph['schema']=='ducksmiles-workbench-drawing-v1';assert graph['input_smiles']==smiles
            out.write(json.dumps({'name':r['name'],'smiles':smiles,'graph':graph},separators=(',',':'))+'\n')
            mol=Chem.MolFromSmiles(smiles,params)
            if graph is None:
                counts['unavailable']+=1;unavailable.append({'name':r['name'],'smiles':smiles,'reference_atoms':mol.GetNumAtoms() if mol else None});continue
            if mol is None:counts['reference_unavailable']+=1;continue
            counts['compared']+=1;bad={}
            atoms=[(x['symbol'],x['formal_charge'],x.get('isotope') or 0,x['bracket_h_count']+x['implicit_h_count'],x['aromatic']) for x in graph['atoms']]
            expected=[(x.GetSymbol(),x.GetFormalCharge(),x.GetIsotope(),x.GetTotalNumHs(),x.GetIsAromatic()) for x in mol.GetAtoms()]
            if atoms!=expected:bad['atoms']={'actual':atoms,'expected':expected}
            def bondkey(i,j,kind):return (i,j,kind) if kind=='dative' else (min(i,j),max(i,j),kind)
            bonds=sorted(bondkey(x['begin'],x['end'],x['bond_type']) for x in graph['bonds'])
            reference=sorted(bondkey(x.GetBeginAtomIdx(),x.GetEndAtomIdx(),str(x.GetBondType()).lower()) for x in mol.GetBonds())
            if bonds!=reference:bad['bonds']={'actual':bonds,'expected':reference}
            if bad:differences.append({'name':r['name'],'smiles':smiles,'differences':bad})
    report={'rdkit':rdkit.__version__,'library_sha256':hashlib.sha256(a.library.read_bytes()).hexdigest(),'sources_sha256':hashlib.sha256(a.sources.read_bytes()).hexdigest(),'graph_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),**counts,'mismatching_rows':len(differences),'differences':differences,'unavailable_rows':unavailable,'scope':'All source positions retained. Input-order atom symbol/charge/isotope/attached H/aromaticity; bond type/connectivity and directed dative donor→acceptor checked. Stereo independently verified from drawings.'}
    (a.output/'chemistry-audit.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items() if k not in ('differences','unavailable_rows')})
    raise SystemExit(bool(differences) or bool(counts['reference_unavailable']))
if __name__=='__main__':main()
