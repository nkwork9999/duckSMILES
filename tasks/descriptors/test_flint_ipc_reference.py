import json
import math
from pathlib import Path
import unittest
from decimal import Decimal, localcontext

from flint_ipc_reference import coefficients


class FlintCoefficientReferenceTests(unittest.TestCase):
    def test_explicit_reference_backend_is_selected_without_silent_fallback(self):
        from exact_ipc_reference import coefficients as default_coefficients
        path3 = [[1], [0, 2], [1]]
        self.assertEqual(default_coefficients(path3, backend='flint'), [1, 0, -2, 0])
        with self.assertRaises(ValueError):
            default_coefficients(path3, backend='floating')

    def test_default_oracle_routes_extended_vertex_and_edge_domains(self):
        from exact_ipc_reference import coefficients as default_coefficients
        path257 = [[j for j in (i-1, i+1) if 0 <= j < 257] for i in range(257)]
        complete33 = [[j for j in range(33) if j != i] for i in range(33)]
        for adjacency in [path257, complete33]:
            with self.subTest(atoms=len(adjacency)):
                self.assertEqual(default_coefficients(adjacency), coefficients(adjacency))

    def test_all_determinant_and_analytic_controls_match_every_signed_coefficient(self):
        cases = json.loads(Path(__file__).with_name('exact-ipc-reference.json').read_text())['cases']
        self.assertEqual(len(cases), 58)
        for case in cases:
            with self.subTest(name=case['name']):
                adjacency = [[] for _ in range(case['atoms'])]
                for a, b in case['edges']:
                    adjacency[a].append(b)
                    adjacency[b].append(a)
                self.assertEqual(list(map(str, coefficients(adjacency))), case['coefficients'])

    def test_large_paths_match_the_closed_form_without_a_native_dependency(self):
        for n in [257, 384, 512, 640]:
            with self.subTest(atoms=n):
                adjacency = [[j for j in (i-1, i+1) if 0 <= j < n] for i in range(n)]
                expected = [0] * (n+1)
                for j in range(n//2+1):
                    expected[2*j] = (-1)**j * math.comb(n-j, j)
                self.assertEqual(coefficients(adjacency), expected)
                from exact_ipc_reference import information
                actual = information(expected)
                with localcontext() as ctx:
                    ctx.prec = 85
                    absolute = [abs(c) for c in expected if c]
                    total = sum(absolute)
                    average = -sum((Decimal(c)/Decimal(total)) * (Decimal(c)/Decimal(total)).ln()/Decimal(2).ln()
                                   for c in absolute)
                    for id, value in [('Ipc', average*Decimal(total)), ('AvgIpc', average)]:
                        numeric = float(value)
                        self.assertLessEqual(abs(actual[id]-numeric), max(1e-12, abs(numeric)*1e-12))

    def test_domain_errors_do_not_become_zero_polynomials(self):
        complete52 = [[j for j in range(52) if j != i] for i in range(52)]
        for adjacency in [[], [[] for _ in range(641)], [[0]], [[1, 1], [0]], [[1], []], complete52]:
            with self.subTest(adjacency_size=len(adjacency)):
                with self.assertRaises(ValueError):
                    coefficients(adjacency)


if __name__ == '__main__':
    unittest.main()
