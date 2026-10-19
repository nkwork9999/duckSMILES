"""Experimental exact IPC oracle, not a catalog/compatibility policy.

Small graphs: determinant evaluation (integer Bareiss) and rational Newton
interpolation, independent of the proposed native Faddeev recurrence.
Large graphs: closed-form path/cycle/complete/bipartite polynomials and products.
Decimal entropy avoids using native floating-point results as expected values.
No native acceptance selection and no new dependency.
"""
import hashlib
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def determinant(matrix):
    n = len(matrix)
    if not n:
        return 1
    a = [row[:] for row in matrix]
    previous, sign = 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * value - a[i][k] * a[k][j]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                a[i][j] = quotient
            a[i][k] = 0
        previous = value
    return sign * a[-1][-1]


def determinant_interpolation(n, edges):
    matrix = [[0] * n for _ in range(n)]
    for a, b in edges:
        matrix[a][b] = matrix[b][a] = -1
    values = []
    for x in range(n + 1):
        for i in range(n):
            matrix[i][i] = x
        values.append(determinant(matrix))
    differences = []
    while values:
        differences.append(values[0])
        values = [b - a for a, b in zip(values, values[1:])]
    result = [Fraction(0) for _ in range(n + 1)]
    basis, factorial = [1], 1
    for k, difference in enumerate(differences):
        for j, value in enumerate(basis):
            result[j] += Fraction(difference * value, factorial)
        next_basis = [0] * (len(basis) + 1)
        for j, value in enumerate(basis):
            next_basis[j] -= k * value
            next_basis[j + 1] += value
        basis = next_basis
        factorial *= k + 1
    assert all(v.denominator == 1 for v in result)
    coefficients = [int(v) for v in reversed(result)]
    assert coefficients[0] == 1 and coefficients[1] == 0
    return coefficients


def path(n):
    c = [0] * (n + 1)
    for k in range(n // 2 + 1):
        c[2 * k] = (-1) ** k * math.comb(n - k, k)
    return c


def cycle(n):
    c = path(n)
    for k, value in enumerate(path(n - 2)):
        c[k + 2] -= value
    c[-1] -= 2
    return c


def product(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def entropy(coefficients):
    with localcontext() as context:
        context.prec = 85
        values = [Decimal(abs(c)) for c in coefficients if c]
        total = sum(values)
        average = -sum((v / total) * (v / total).ln() for v in values) / Decimal(2).ln()
        return str(average * total), str(average)


def main():
    cases = []
    def add(name, n, edges, coefficients=None):
        edges = sorted({tuple(sorted(e)) for e in edges})
        assert n > 0 and all(0 <= a < b < n for a, b in edges)
        assert n <= 256 and len(edges) <= 512
        if coefficients is None:
            coefficients = determinant_interpolation(n, edges)
        elif n <= 20:
            assert determinant_interpolation(n, edges) == coefficients, name
        assert len(coefficients) == n + 1
        ipc, avg = entropy(coefficients)
        cases.append(dict(name=name, atoms=n, edges=edges,
                          coefficients=list(map(str, coefficients)), exact_ipc=ipc, exact_avg_ipc=avg))

    for n in [1, 2, 3, 4, 7, 12, 20, 40, 80, 128, 200, 256]:
        add(f"path-{n}", n, [(i, i + 1) for i in range(n - 1)], path(n))
    for n in [3, 4, 5, 6, 7, 8, 12, 20, 64, 128, 256]:
        add(f"cycle-{n}", n, [(i, (i + 1) % n) for i in range(n)], cycle(n))
    for n in [2, 8, 32, 256]:
        c = [1, 0, -(n - 1)] + [0] * (n - 2)
        add(f"star-{n}", n, [(0, j) for j in range(1, n)], c)
    for n in [3, 8, 20, 32]:
        def choose(k):
            return math.comb(n - 1, k) if 0 <= k < n else 0
        c = [choose(k) - (n - 1) * choose(k - 1) for k in range(n + 1)]
        add(f"complete-{n}", n, [(i, j) for i in range(n) for j in range(i + 1, n)], c)
    for p, q in [(2, 3), (4, 8), (16, 16)]:
        n = p + q
        add(f"bipartite-{p}-{q}", n, [(i, p + j) for i in range(p) for j in range(q)], [1, 0, -p * q] + [0] * (n - 2))
    for copies, width in [(4, 3), (32, 6), (42, 3), (64, 4)]:
        edges, coefficients = [], [1]
        for k in range(copies):
            edges += [(k * width + i, k * width + (i + 1) % width) for i in range(width)]
            coefficients = product(coefficients, cycle(width))
        add(f"cycles-{copies}x{width}", copies * width, edges, coefficients)
    add("two-paths-128", 256,
        [(i, i + 1) for i in range(127)] + [(i, i + 1) for i in range(128, 255)], product(path(128), path(128)))
    for n in [1, 3, 20, 256]:
        add(f"isolated-{n}", n, [], [1] + [0] * n)
    # Deterministic irregular graphs, with integer determinant interpolation.
    for n in [4, 7, 12, 18, 24]:
        for seed in [3, 17, 41]:
            edges = [(i, j) for i in range(n) for j in range(i + 1, n)
                     if (i * 131 + j * 17 + seed * 31) % 101 < 22]
            add(f"irregular-{n}-{seed}", n, edges)
    root = Path(__file__).resolve().parents[2]
    destination = root / "tasks/descriptors/exact-ipc-reference.json"
    data = dict(scope="Experimental mathematical oracle, no pinned RDKit compatibility policy", cases=cases)
    destination.write_text(json.dumps(data, indent=2) + "\n")
    lines = ["# name; n; undirected edges; signed characteristic coefficients highest-first; exact IPC; exact AvgIPC"]
    for case in cases:
        lines.append("\t".join([case["name"], str(case["atoms"]),
                                ";".join(f"{a},{b}" for a, b in case["edges"]),
                                ",".join(case["coefficients"]), case["exact_ipc"], case["exact_avg_ipc"]]))
    destination.with_suffix(".tsv").write_text("\n".join(lines) + "\n")
    print(len(cases), "exact graph cases", hashlib.sha256(destination.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
