"""Independent integer graph-information oracle for the two stable definitions.

Python arbitrary-width integers/full matrix recurrence contrast with the native
signed-limb/symmetric-triangle implementation. The larger domain delegates to
FLINT exact integer matrices; verification can explicitly select FLINT for all
sizes. Determinant/interpolation and analytic graph oracles separately verify
every coefficient on58 declared cases and four paths up to640 vertices.
RDKit supplies only graph connectivity, with explicit H vertices retained.
"""
import json
import math
from pathlib import Path
import numpy as np


def coefficients(adjacency, *, backend='python'):
    if backend not in ('python', 'flint'):
        raise ValueError('Unknown exact coefficient backend')
    n = len(adjacency)
    if backend == 'flint' or n > 256 or sum(map(len, adjacency)) > 1024:
        from flint_ipc_reference import coefficients as flint_coefficients
        return flint_coefficients(adjacency)
    if not 0 < n <= 256:
        raise ValueError("exact graph oracle domain is1..256 vertices")
    if any(len(row) != len(set(row)) or any(j < 0 or j >= n or j == i for j in row)
           for i, row in enumerate(adjacency)):
        raise ValueError("simple adjacency required")
    if sum(map(len, adjacency)) > 1024:
        raise ValueError("at most512 undirected edges")
    if any(i not in adjacency[j] for i, row in enumerate(adjacency) for j in row):
        raise ValueError("undirected adjacency required")
    if not any(adjacency):
        return [1] + [0] * n
    # Object arrays retain Python arbitrary-width integers. Only graph indices
    # use fixed-width np.intp; no coefficient passes through a float/int64 array.
    sources = np.array([i for i, neighbors in enumerate(adjacency) for _ in neighbors], dtype=np.intp)
    targets = np.array([j for neighbors in adjacency for j in neighbors], dtype=np.intp)
    diagonal = np.arange(n)
    previous = np.eye(n, dtype=object)
    result = [1]
    for k in range(1, n + 1):
        current = np.zeros((n, n), dtype=object)
        # Unbuffered accumulation includes every neighbor, including repeated
        # source indices. Full matrices independently contrast with native half
        # matrices and custom signed limbs.
        np.add.at(current, sources, previous[targets])
        trace = sum(current[i, i] for i in range(n))
        quotient, remainder = divmod(-trace, k)
        if remainder or type(quotient) is not int:
            raise ArithmeticError("integer characteristic trace division was not exact")
        result.append(quotient)
        current[diagonal, diagonal] += quotient
        previous = current
    if np.any(previous):
        raise ArithmeticError("Cayley-Hamilton identity failed")
    return result


def information(signed):
    absolute = [abs(c) for c in signed if c]
    total = sum(absolute)
    average = -math.fsum((v / total) * math.log2(v / total) for v in absolute)
    weighted = total * average
    if not math.isfinite(weighted) or not math.isfinite(average):
        raise ArithmeticError("nonfinite graph information")
    return {"Ipc": weighted, "AvgIpc": average}


def reference(mol, *, backend='python'):
    # Share the two reductions only within one fresh reference molecule.
    # There is no cross-input or cross-trial molecule-result cache.
    key = '_chemob_exact_information' if backend == 'python' else '_chemob_exact_information_' + backend
    cached = getattr(mol, key, None)
    if cached is not None:
        return cached
    adjacency = [[] for _ in mol.GetAtoms()]
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adjacency[a].append(b)
        adjacency[b].append(a)
    result = information(coefficients(adjacency, backend=backend))
    setattr(mol, key, result)
    return result


def self_test():
    source = Path(__file__).with_name("exact-ipc-reference.json")
    cases = json.loads(source.read_text())["cases"]
    for case in cases:
        adjacency = [[] for _ in range(case["atoms"])]
        for a, b in case["edges"]:
            adjacency[a].append(b)
            adjacency[b].append(a)
        signed = coefficients(adjacency)
        assert list(map(str, signed)) == case["coefficients"], case["name"]
        values = information(signed)
        for id, key in [("Ipc", "exact_ipc"), ("AvgIpc", "exact_avg_ipc")]:
            expected = float(case[key])
            assert abs(values[id] - expected) <= max(1e-12, 1e-12 * abs(expected)), (case["name"], id)
    print(len(cases), "Python full-integer matrices match determinant/analytic coefficients and Decimal85 entropy")


if __name__ == "__main__":
    self_test()
