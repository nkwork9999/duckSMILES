import unittest

import duckdb

from ducksmiles_embedded import ConcentrationBasis, ConcentrationTable


class ConcentrationTests(unittest.TestCase):
    def test_mass_and_molar_units_share_only_an_explicit_compatible_basis(self):
        with duckdb.connect() as db:
            db.execute("""
                CREATE TABLE measurements(id BIGINT, value DOUBLE, unit VARCHAR, grams_per_mole DOUBLE);
                INSERT INTO measurements VALUES
                   (1, 1000, 'mg/L', NULL), (2, 1, 'g/L', NULL), (3, 2, 'mmol/L', 100),
                   (4, 2, 'mmol/L', NULL), (5, 2, 'seconds', 100), (6, NULL, 'g/L', NULL),
                   (7, 200, 'mg/L', 100);
            """)
            table = ConcentrationTable(db, "measurements")
            mass = table.read(ConcentrationBasis.MASS)
            self.assertEqual([(row.value, row.reason) for row in mass], [
                (1, None), (1, None), (0.2, None), (None, "molar_mass_required"),
                (None, "unknown_unit"), (None, "missing_value"), (0.2, None),
            ])
            self.assertEqual([row.unit for row in mass], ["g/L"] * 7)
            self.assertEqual(mass[-1].source_unit, "mg/L")
            self.assertEqual(mass[-1].rule_version, "concentration-1")
            molar = table.read(ConcentrationBasis.MOLAR)
            self.assertEqual(molar[0].reason, "molar_mass_required")
            self.assertAlmostEqual(molar[3].value, 0.002)
            self.assertAlmostEqual(molar[-1].value, 0.002)
