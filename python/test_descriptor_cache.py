import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import duckdb

from ducksmiles_embedded import DescriptorCache, DescriptorRevision, MoleculeReader


class DescriptorCacheTests(unittest.TestCase):
    def test_persists_negative_results_and_separates_versions_and_policies(self):
        artifact = os.getenv("DUCKSMILES_TEST_EXTENSION")
        config = {"allow_unsigned_extensions": bool(artifact)}
        samples = ["O", "CCO", None, "invalid", "CCO"]
        first = DescriptorRevision("molecules-v1", "identity-v1")
        with TemporaryDirectory() as directory:
            database = Path(directory) / "descriptors.duckdb"
            with duckdb.connect(database, config=config) as db:
                db.load_extension(artifact or "ducksmiles")
                expected = tuple(MoleculeReader(db).describe(samples))
                response = DescriptorCache(db).describe(first, samples)
                self.assertEqual(response.rows, expected)
                self.assertEqual(response.new_entries, 3)
            with duckdb.connect(database, config=config) as db:
                db.load_extension(artifact or "ducksmiles")
                cache = DescriptorCache(db)
                response = cache.describe(first, samples)
                self.assertEqual(response.rows, expected)
                self.assertEqual(response.new_entries, 0)
                self.assertEqual(cache.describe(first, ["O", "CCCC"]).new_entries, 1)
                self.assertEqual(cache.describe(DescriptorRevision("molecules-v2", "identity-v1"), samples).new_entries, 3)
                self.assertEqual(cache.describe(DescriptorRevision("molecules-v1", "stereo-v2"), samples).new_entries, 3)
                self.assertEqual(cache.describe(first, []).rows, ())
