"""Independent integer graph-information oracle for the two stable definitions.

Python arbitrary-width integers/full matrix recurrence contrast with the native
signed-limb/symmetric-triangle implementation. Determinant/interpolation and
analytic graph oracles separately verify every coefficient on58 declared cases.
RDKit supplies only graph connectivity, with explicit H vertices retained.
"""
import json
import math
from pathlib import Path


def coefficients(adjacency):
    n = len(adjacency)
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
    previous = [[int(i == j) for j in range(n)] for i in range(n)]
    result = [1]
    for k in range(1, n + 1):
        current = [[sum(previous[u][j] for u in neighbors) for j in range(n)]
                   if neighbors else [0] * n for neighbors in adjacency]
        trace = sum(current[i][i] for i in range(n))
        quotient, remainder = divmod(-trace, k)
        if remainder:
            raise ArithmeticError("integer characteristic trace division was not exact")
        result.append(quotient)
        for i in range(n):
            current[i][i] += quotient
        previous = current
    if any(any(row) for row in previous):
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


def reference(mol):
    # Share the two reductions only within one fresh reference molecule.
    # There is no cross-input or cross-trial molecule-result cache.
    cached = getattr(mol, "_chemob_exact_information", None)
    if cached is not None:
        return cached
    adjacency = [[] for _ in mol.GetAtoms()]
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adjacency[a].append(b)
        adjacency[b].append(a)
    result = information(coefficients(adjacency))
    mol._chemob_exact_information = result
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
