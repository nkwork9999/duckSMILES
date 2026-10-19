"""Preserve whole-population coverage and the two independent accuracy gates."""
import unittest
from parallel_catalog import partitions, merge_metric


class ParallelCatalogTests(unittest.TestCase):
    def test_partitions_cover_each_original_position_once(self):
        for size in [1, 15, 34, 100_000]:
            parts = partitions(size, 4)
            self.assertEqual([i for start, end in parts for i in range(start, end)], list(range(size)))
            self.assertLessEqual(max(b-a for a,b in parts)-min(b-a for a,b in parts), 1)

    def test_metric_merge_keeps_failures_undefined_and_independent_policy(self):
        policy = dict(reference_source='exact-integer-characteristic-polynomial', tolerance=1e-10, relative_tolerance=1e-12)
        first = dict(**policy, compared=3, mismatches=1, reference_errors=0, undefined_native=1, undefined_reference=1, matching_undefined=1, max_abs_error=2.0, examples=[dict(index=1)])
        second = dict(**policy, compared=4, mismatches=2, reference_errors=1, undefined_native=0, undefined_reference=0, matching_undefined=0, max_abs_error=5.0, examples=[dict(index=8)])
        actual = merge_metric([first, second])
        self.assertEqual(actual['compared'], 7)
        self.assertEqual(actual['mismatches'], 3)
        self.assertEqual(actual['reference_errors'], 1)
        self.assertEqual(actual['matching_undefined'], 1)
        self.assertEqual(actual['max_abs_error'], 5)
        self.assertEqual(actual['reference_source'], policy['reference_source'])
        self.assertEqual([r['index'] for r in actual['examples']], [1,8])
        second['relative_tolerance'] = 1e-6
        with self.assertRaises(AssertionError):
            merge_metric([first,second])


if __name__ == '__main__':
    unittest.main()
