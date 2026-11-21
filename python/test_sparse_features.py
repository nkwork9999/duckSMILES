import unittest

import duckdb

from ducksmiles_embedded import SparseFeatures


class SparseFeaturesTests(unittest.TestCase):
    def test_coordinates_round_trip_zeros_missing_and_duplicate_terms(self):
        with duckdb.connect() as db:
            features = SparseFeatures(db, {"first": 10, "second": 20, "third": 30})
            count = features.store("v1", [1, 2, 3], [
                (1, "first", 5.0), (1, "first", -2.0), (1, "second", None),
                (2, "first", -1.0), (2, "first", 1.0), (3, "third", 1.25),
            ])
            self.assertEqual(count, 3)
            self.assertEqual(features.dense("v1", [1, 2, 3], ["second", "first", "third"]), [
                (1, (None, 3, 0)), (2, (0, 0, 0)), (3, (0, 0, 1.25)),
            ])
            self.assertEqual(features.column("v1", "second"), [(1, None), (2, 0), (3, 0)])
            self.assertEqual(features.column("v1", "second", include_zeros=False), [(1, None)])
            self.assertEqual(features.coordinates("v1", [3, 1]), [(1, 10, 3), (1, 20, None), (3, 30, 1.25)])
            with self.assertRaises(KeyError):
                features.dense("v1", [99], ["first"])
            with self.assertRaises(ValueError):
                SparseFeatures(db, {"first": 90, "second": 20, "third": 30})
            with self.assertRaises(KeyError):
                features.store("v2", [1], [(1, "unknown", 1)])
            self.assertEqual(features.dense("v1", [1], ["first"]), [(1, (3,))])
