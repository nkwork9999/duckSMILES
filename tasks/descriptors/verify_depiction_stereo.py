"""Recover tetrahedral tags from actual 2D coordinates/wedges using RDKit.

The input molecule supplies connectivity and an independent expected tag only.
All input atom tags/bond directions are erased from the drawing reconstruction.
Original bond insertion order is retained so tag parity is comparable.
"""
import argparse,json
from pathlib import Path
import rdkit
from rdkit import Chem,RDLogger
assert rdkit.__version__=='2026.03.6'
RDLogger.DisableLog('rdApp.*')

def verify(record,counts=None):
    params=Chem.SmilesParserParams();params.removeHs=False
    original=Chem.MolFromSmiles(record['smiles'],params)
    assert original is not None
    drawing=record['drawing']; points=record.get('screenPoints',drawing['points'])
    assert len(points)==original.GetNumAtoms()
    wedges={(min(w['center'],w['neighbor']),max(w['center'],w['neighbor'])):w for w in drawing['wedges']}
    assert len(wedges)==len(drawing['wedges'])
    rw=Chem.RWMol()
    for atom in original.GetAtoms():
        a=Chem.Atom(atom);a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED);rw.AddAtom(a)
    for b in original.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx();w=wedges.get((min(i,j),max(i,j)))
        if w:i,j=w['center'],w['neighbor']
        rw.AddBond(i,j,b.GetBondType())
        out=rw.GetBondBetweenAtoms(i,j)
        out.SetIsAromatic(b.GetIsAromatic())
        if w:out.SetBondDir(Chem.BondDir.BEGINDASH if w['hashed'] else Chem.BondDir.BEGINWEDGE)
    mol=rw.GetMol();mol.UpdatePropertyCache(strict=False)
    conf=Chem.Conformer(len(points));conf.Set3D(False)
    for i,p in enumerate(points):conf.SetAtomPosition(i,(p['x'],-p['y'] if 'screenPoints' in record else p['y'],0))
    mol.AddConformer(conf)
    Chem.AssignChiralTypesFromBondDirs(mol,replaceExistingTags=True)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    bad=[]
    for w in drawing['wedges']:
        i=w['center'];expected=original.GetAtomWithIdx(i).GetChiralTag();actual=mol.GetAtomWithIdx(i).GetChiralTag()
        if counts is not None:
            key='reference_tetrahedral_centers' if expected in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW) else 'non_stereogenic_wedge_annotations'
            counts[key]=counts.get(key,0)+1
        if actual!=expected:bad.append(dict(center=i,expected=str(expected),actual=str(actual),hashed=w['hashed']))
    # Independently assigned RDKit stereo atoms identify the substituents to
    # inspect. The layout instead consumes raw slash/backslash annotations.
    mathematical=[(p['x'],-p['y'] if 'screenPoints' in record else p['y']) for p in points]
    for alkene in drawing.get('alkenes',[]):
        bond=original.GetBondBetweenAtoms(alkene['a'],alkene['b'])
        assert bond is not None
        expected=bond.GetStereo();neighbors=list(bond.GetStereoAtoms())
        if expected==Chem.BondStereo.STEREONONE:
            if counts is not None:counts['non_stereogenic_alkene_annotations']=counts.get('non_stereogenic_alkene_annotations',0)+1
            continue  # non-stereogenic duplicate substituents
        if counts is not None:counts['reference_stereogenic_alkenes']=counts.get('reference_stereogenic_alkenes',0)+1
        if len(neighbors)!=2:
            bad.append(dict(alkene=alkene,expected=str(expected),error='missing reference stereo atoms'));continue
        a,b=mathematical[bond.GetBeginAtomIdx()],mathematical[bond.GetEndAtomIdx()]
        left,right=mathematical[neighbors[0]],mathematical[neighbors[1]]
        axis=(b[0]-a[0],b[1]-a[1])
        cross=lambda p,c:axis[0]*(p[1]-c[1])-axis[1]*(p[0]-c[0])
        product=cross(left,a)*cross(right,b)
        opposite=expected in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOTRANS)
        if abs(product)<1e-8 or ((product<0)!=opposite):bad.append(dict(alkene=alkene,expected=str(expected),cross_product=product))
    return bad

def main():
    parser=argparse.ArgumentParser();parser.add_argument('input',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
    rows=centers=alkenes=screen_rows=0;bad=[];counts={}
    for line in args.input.open():
        record=json.loads(line);rows+=1;centers+=len(record['drawing']['wedges'])
        alkenes+=len(record['drawing'].get('alkenes',[]));screen_rows+=int('screenPoints' in record)
        failures=verify(record,counts)
        if failures:bad.append(dict(smiles=record['smiles'],failures=failures))
    report=dict(rdkit=rdkit.__version__,drawing_rows=rows,assigned_centers=centers,depicted_alkenes=alkenes,reference_counts=counts,screen_coordinate_rows=screen_rows,mismatching_rows=len(bad),failures=bad,
                scope='Source tetrahedral tags erased; recover from screen coordinates and actual wedge directions. Declared E/Z geometry checked against independently assigned RDKit stereo atoms. Omitted configurations excluded.')
    args.output.write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items() if k!='failures'})
    raise SystemExit(1 if bad or not centers else 0)
if __name__=='__main__':main()
