"""Additional exact coefficient oracle for larger graphs; reference environment only."""
from flint import fmpz_mat, __version__


def coefficients(adjacency):
    if __version__ != '0.9.0':
        raise RuntimeError('Reference environment requires python-flint0.9.0')
    n = len(adjacency)
    if not 0 < n <= 640:
        raise ValueError('Reference graph domain is1..640 vertices')
    if sum(map(len, adjacency)) > 2560:
        raise ValueError('At most1280 undirected edges')
    if any(len(row) != len(set(row)) or any(j < 0 or j >= n or j == i for j in row)
           for i, row in enumerate(adjacency)):
        raise ValueError('Simple adjacency required')
    if any(i not in adjacency[j] for i, row in enumerate(adjacency) for j in row):
        raise ValueError('Undirected adjacency required')
    entries = [0] * (n*n)
    for i, neighbors in enumerate(adjacency):
        for j in neighbors:
            entries[i*n+j] = 1
    polynomial = fmpz_mat(n, n, entries).charpoly()
    # FLINT's coefficient indices ascend from the constant term. Native and
    # existing determinant fixtures use the leading monic coefficient first.
    result = [int(polynomial[i]) for i in range(n, -1, -1)]
    if polynomial.degree() != n or result[0] != 1:
        raise ArithmeticError('Characteristic polynomial is not monic of degree n')
    return result
