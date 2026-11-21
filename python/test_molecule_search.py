import os
import unittest

import duckdb

from ducksmiles_embedded import FingerprintFilter, MassRange, MoleculeSearch, MoleculeSearchSpec, PatternId


class MoleculeSearchTests(unittest.TestCase):
    def test_bounded_typed_search_can_combine_mass_fingerprint_and_registered_pattern(self):
        artifact = os.getenv("DUCKSMILES_TEST_EXTENSION")
        with duckdb.connect(config={"allow_unsigned_extensions": bool(artifact)}) as db:
            db.load_extension(artifact or "ducksmiles")
            db.execute("""
                CREATE TABLE molecules(molecule_id BIGINT PRIMARY KEY, smiles VARCHAR, mass DOUBLE);
                INSERT INTO molecules VALUES (10, 'CCO', 46.07), (3, 'CCO', 46.07), (5, 'CC', 30),
                    (2, 'CCCO', 60.1), (7, 'O', NULL), (9, 'invalid', NULL);
            """)
            oxygen = PatternId("contains oxygen")
            service = MoleculeSearch(db, "molecules", {oxygen: "[#8]"})
            result = service.search(MoleculeSearchSpec(mass=MassRange(40, 55), pattern=oxygen,
                                                        fingerprint=FingerprintFilter("CCO", 0.8), limit=1))
            self.assertEqual(result.state, "ok")
            self.assertEqual([(r.molecule_id, r.score) for r in result.hits], [(3, 1.0)])
            result = service.search(MoleculeSearchSpec(pattern=oxygen))
            self.assertEqual([(r.molecule_id, r.mass is None) for r in result.hits],
                             [(2, False), (3, False), (7, True), (10, False)])
            self.assertEqual(service.search(MoleculeSearchSpec(mass=MassRange(110, 120))).hits, ())
            self.assertEqual(service.search(MoleculeSearchSpec(fingerprint=FingerprintFilter("invalid", 0.8))).state,
                             "invalid_reference")
            with self.assertRaises(KeyError):
                service.search(MoleculeSearchSpec(pattern=PatternId("unknown")))
            with self.assertRaises(TypeError):
                MoleculeSearchSpec(pattern="contains oxygen")
