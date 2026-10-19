"""Exact cut-vertex polynomial experiment; never used as verification oracle.

For components Ci of G-v, let qi=P(Ci), pi=P(Ci+v). Then
P(G)=sum_i pi*product(j!=i,qj) -(k-1)*x*product_i qi.
Small/biconnected leaves use FLINT in this prototype only. A native port must
use its independent integer recurrence and reproduce every signed coefficient.
"""
from functools import lru_cache
from math import prod

from flint_ipc_reference import coefficients as direct


def multiply(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if left:
            for j, right in enumerate(b):
                result[i+j] += left * right
    return result


def factored(adjacency, cutoff=32):
    calls = {'leaf_calls': 0, 'largest_leaf': 0, 'cut_calls': 0}

    def components(vertices, excluded=None):
        unseen = set(vertices)
        unseen.discard(excluded)
        result = []
        while unseen:
            seed = min(unseen)
            unseen.remove(seed)
            queue = [seed]
            for vertex in queue:
                for neighbor in adjacency[vertex]:
                    if neighbor in unseen:
                        unseen.remove(neighbor)
                        queue.append(neighbor)
            result.append(tuple(sorted(queue)))
        return result

    def balanced_cut(vertices):
        allowed = set(vertices)
        discovery, low, sizes, candidates = {}, {}, {}, []

        def visit(v, parent=None):
            discovery[v] = low[v] = len(discovery)
            sizes[v] = 1
            separated = []
            for neighbor in adjacency[v]:
                if neighbor not in allowed or neighbor == parent:
                    continue
                if neighbor not in discovery:
                    visit(neighbor, v)
                    sizes[v] += sizes[neighbor]
                    low[v] = min(low[v], low[neighbor])
                    if low[neighbor] >= discovery[v]:
                        separated.append(sizes[neighbor])
                else:
                    low[v] = min(low[v], discovery[neighbor])
            remaining = len(vertices) - 1 - sum(separated)
            if remaining:
                separated.append(remaining)
            if len(separated) > 1:
                candidates.append((max(separated), v))

        visit(vertices[0])
        return min(candidates)[1] if candidates else None

    @lru_cache(None)
    def solve(vertices):
        n = len(vertices)
        if not n:
            return (1,)
        parts = components(vertices)
        if len(parts) > 1:
            result = [1]
            for part in parts:
                result = multiply(result, solve(part))
            return tuple(result)
        cut = balanced_cut(vertices) if n > cutoff else None
        if cut is None:
            calls['leaf_calls'] += 1
            calls['largest_leaf'] = max(calls['largest_leaf'], n)
            remap = {v: i for i, v in enumerate(vertices)}
            graph = [[remap[j] for j in adjacency[v] if j in remap] for v in vertices]
            return tuple(direct(graph))
        calls['cut_calls'] += 1
        parts = components(vertices, cut)
        deleted = [solve(part) for part in parts]
        arms = [solve(tuple(sorted((*part, cut)))) for part in parts]
        prefix = [[1]]
        for p in deleted:
            prefix.append(multiply(prefix[-1], p))
        suffix = [[1] for _ in range(len(parts) + 1)]
        for i in range(len(parts)-1, -1, -1):
            suffix[i] = multiply(deleted[i], suffix[i+1])
        result = [-(len(parts)-1) * c for c in prefix[-1]] + [0]
        assert len(result) == n + 1
        for i, arm in enumerate(arms):
            term = multiply(arm, multiply(prefix[i], suffix[i+1]))
            assert len(term) == len(result)
            result = [a+b for a,b in zip(result, term)]
        return tuple(result)

    answer = list(solve(tuple(range(len(adjacency)))))
    calls['memoized_subgraphs'] = solve.cache_info().currsize
    return answer, calls


if __name__ == '__main__':
    import json
    from pathlib import Path
    cases = json.loads(Path(__file__).with_name('exact-ipc-reference.json').read_text())['cases']
    for case in cases:
        adjacency = [[] for _ in range(case['atoms'])]
        for a,b in case['edges']:
            adjacency[a].append(b); adjacency[b].append(a)
        actual, _ = factored(adjacency, cutoff=3)
        assert list(map(str, actual)) == case['coefficients'], case['name']
    print('All58 independent coefficient controls match the cut-vertex identity')
