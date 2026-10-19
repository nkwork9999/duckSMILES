"""Numerical stability diagnostic: exact path-graph polynomial vs RDKit/NumPy.

This is not a tolerance policy. IPC remains outside the shipped catalog until
its mathematical and numerical compatibility policy is resolved.
"""
import json,math
from pathlib import Path
import rdkit,numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
assert rdkit.__version__=='2026.03.6'
def entropy(c):
 total=sum(c)
 return -sum((v/total)*math.log(v/total) for v in c if v)/math.log(2)
def sparse_ipc(m):
 n=m.GetNumAtoms();a=np.asarray(Chem.GetAdjacencyMatrix(m),dtype=float)
 an=a.copy();adj=[sorted(x.GetIdx() for x in atom.GetNeighbors()) for atom in m.GetAtoms()]
 c=[1.]
 for k in range(1,n+1):
  coef=(1./k)*sum(an[i,i] for i in range(n));c.append(abs(coef))
  for i in range(n):an[i,i]-=coef
  nxt=np.zeros((n,n))
  for i in range(n):
   for j in adj[i]:nxt[i]+=an[j]
  an=nxt
 avg=entropy(c)
 return sum(c)*avg,avg
rows=[]
for n in [20,40,80,128,200,256]:
 m=Chem.MolFromSmiles('C'*n)
 # Characteristic polynomial of path P_n:
 # coefficient of x^(n-2k) = (-1)^k choose(n-k,k).
 coefficients=[float(math.comb(n-k,k)) for k in range(n//2+1)]
 exact_avg=entropy(coefficients);exact=sum(coefficients)*exact_avg
 actual=Descriptors.Ipc(m);avg=Descriptors.AvgIpc(m)
 native,native_avg=sparse_ipc(m)
 rows.append(dict(atoms=n,analytic_ipc=exact,analytic_avg_ipc=exact_avg,
  rdkit_ipc=actual,rdkit_avg_ipc=avg,sparse_recurrence_ipc=native,sparse_recurrence_avg_ipc=native_avg,
  rdkit_relative_error_vs_analytic=abs(actual-exact)/exact,
  sparse_relative_error_vs_analytic=abs(native-exact)/exact))
print(json.dumps(dict(rdkit=rdkit.__version__,numpy=np.__version__,blas='Apple Accelerate (numpy.show_config)',
 conclusion='The division/trace recurrence is ill-conditioned on large path graphs. Do not broaden tolerance or call this parity. Both IPC descriptors remain unimplemented in the catalog.',rows=rows),indent=2))
