import unittest

import duckdb
from ducksmiles_embedded import DaltonMass, MassWindowSearch


class MassWindowTests(unittest.TestCase):
    def test_inclusive_bounds_and_reference_denominator(self):
        with duckdb.connect() as db:
            search = MassWindowSearch(db)
            values = tuple(DaltonMass(x) for x in (999, 1000, 1001, 1001.01, 999))
            matches = search.match((DaltonMass(1000),), values, 1000)
            self.assertEqual([(m.reference_index, m.peak_index) for m in matches],
                             [(0, 0), (0, 1), (0, 2), (0, 4)])
            for row, expected in zip(matches, (-1000, 0, 1000, -1000), strict=True):
                self.assertAlmostEqual(row.error_ppm, expected)

    def test_invalid_mass_unit_tolerance_and_empty(self):
        with duckdb.connect() as db:
            search = MassWindowSearch(db)
            for bad in (0, -1, float("inf"), float("nan")):
                with self.assertRaises(ValueError):
                    DaltonMass(bad)
            with self.assertRaises(TypeError):
                search.match((1000,), (DaltonMass(1000),), 1)
            with self.assertRaises(ValueError):
                search.match((), (), -1)
            self.assertEqual(search.match((), (DaltonMass(1000),), 1), [])
