import os
import unittest

import duckdb
from ducksmiles_embedded import MoleculeReader


class MoleculeReaderTests(unittest.TestCase):
    def setUp(self):
        artifact = os.getenv("DUCKSMILES_TEST_EXTENSION")
        config = {"allow_unsigned_extensions": bool(artifact)}
        if directory := os.getenv("DUCKDB_EXTENSION_DIRECTORY"):
            config["extension_directory"] = directory
        self.db = duckdb.connect(config=config)
        self.addCleanup(self.db.close)
        self.db.load_extension(artifact or "ducksmiles")
        self.reader = MoleculeReader(self.db)

    def test_batch_preserves_positions_and_nulls(self):
        rows = self.reader.describe(["CCO", None, "O", "CCO", "invalid"])
        self.assertEqual([r.formula for r in rows], ["C2H6O", None, "H2O", "C2H6O", None])
        self.assertAlmostEqual(rows[0].weight, 46.069, places=2)
        self.assertEqual(rows[0], rows[3])
        self.assertIsNone(rows[4].weight)

    def test_apostrophe_is_molecular_input(self):
        self.assertIsNone(self.reader.describe(["chemist's sample"])[0].formula)
        self.assertEqual(self.reader.describe(["O"])[0].formula, "H2O")

    def test_empty_and_wrong_types(self):
        self.assertEqual(self.reader.describe([]), [])
        for invalid in ("CCO", [123], {"smiles": "CCO"}):
            with self.subTest(invalid=invalid), self.assertRaises(TypeError):
                self.reader.describe(invalid)
