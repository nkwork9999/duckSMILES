"""Independent dense-matrix oracle for the documented heavy-graph extra50.

RDKit only parses the input. These are published descriptors, absent from the
pinned RDKit descriptor API, not new physical measurements or predictive models.
NumPy dense powers/Floyd-Warshall intentionally differ from native sparse DP/BFS.
"""
import math
import numpy as np

MDE = [(element, a, b) for element, maximum in [('C', 4), ('N', 3), ('O', 2)]
       for a in range(1, maximum + 1) for b in range(a, maximum + 1)]
IDS = ([f'ds_MDE{e}{a}{b}' for e, a, b in MDE]
       + [f'ds_GGI{k}' for k in range(1, 11)]
       + [f'ds_MWC{k}' for k in range(3, 11)]
       + [f'ds_SRW{k}' for k in range(3, 11)]
       + ['ds_Wiener', 'ds_Harary', 'ds_EccentricConnectivity', 'ds_DegreeDistance', 'ds_Gutman'])
assert len(IDS) == len(set(IDS)) == 50


def reference(mol):
    vertices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    n = len(vertices)
    if not n:
        return dict.fromkeys(IDS)
    remap = {v: i for i, v in enumerate(vertices)}
    elements = [mol.GetAtomWithIdx(v).GetSymbol() for v in vertices]
    A = np.zeros((n, n), dtype=np.float64)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in remap and b in remap:
            i, j = remap[a], remap[b]
            A[i, j] = A[j, i] = 1
    degree = A.sum(axis=1)
    D = np.where(A != 0, 1.0, np.inf)
    np.fill_diagonal(D, 0)
    for k in range(n):
        D = np.minimum(D, D[:, k:k+1] + D[k:k+1, :])
    values = {}
    for element, a, b in MDE:
        pairs = [D[i, j] for i in range(n) for j in range(i+1, n)
                 if elements[i] == elements[j] == element and math.isfinite(D[i, j])
                 and sorted((degree[i], degree[j])) == [a, b]]
        values[f'ds_MDE{element}{a}{b}'] = (
            len(pairs) * math.exp(-math.fsum(map(math.log, pairs))/len(pairs)) if pairs else None)
    reciprocal = np.zeros_like(D)
    np.divide(1., D, out=reciprocal, where=(D > 0) & np.isfinite(D))
    M = A @ (reciprocal**2)
    CT = M - M.T
    upper = np.triu(np.ones((n, n), dtype=bool), 1)
    for k in range(1, 11):
        values[f'ds_GGI{k}'] = float(np.abs(CT[(D == k) & upper]).sum())
    for k in range(3, 11):
        power = np.linalg.matrix_power(A, k)
        values[f'ds_MWC{k}'] = math.log1p(float(power.sum()))
        values[f'ds_SRW{k}'] = math.log1p(float(power.trace()))
    finite = upper & np.isfinite(D)
    values['ds_Wiener'] = float(D[finite].sum())
    values['ds_Harary'] = float(reciprocal[finite].sum())
    values['ds_EccentricConnectivity'] = float((np.where(np.isfinite(D), D, 0).max(axis=1)*degree).sum())
    values['ds_DegreeDistance'] = float(((degree[:, None] + degree[None, :])[finite]*D[finite]).sum())
    values['ds_Gutman'] = float(((degree[:, None]*degree[None, :])[finite]*D[finite]).sum())
    assert set(values) == set(IDS)
    return values
