"""Prevent stable-definition agreement from being reported as pinned parity."""
import unittest
from compare_catalog import (FUNCTIONS, PARAMS, Chem, reference_source,
                             reference_value, tolerance, relative_tolerance,
                             allowed_error)


class ReferencePolicyTests(unittest.TestCase):
    def test_coefficient_backend_selection_preserves_reference_and_tolerances(self):
        mol = Chem.MolFromSmiles('CCCC', PARAMS)
        spec = {'id': 'Ipc', 'origin': 'rdkit-standard'}
        expected = reference_value(spec, mol, {}, {})
        self.assertEqual(reference_value(spec, mol, {}, {}, information_backend='flint'), expected)
        with self.assertRaises(ValueError):
            reference_value(spec, mol, {}, {}, information_backend='floating')

    def test_stable_information_and_pinned_reference_are_distinct(self):
        mol = Chem.MolFromSmiles("C" * 128, PARAMS)
        spec = {"id": "Ipc", "origin": "rdkit-standard"}
        actual = reference_value(spec, mol, {}, {})
        self.assertAlmostEqual(actual / 1.55012977156e27, 1, places=10)
        self.assertGreater(FUNCTIONS["Ipc"](mol) / actual, 100)
        self.assertEqual(reference_source(spec), "exact-integer-characteristic-polynomial")

    def test_relative_tolerance_only_applies_to_two_new_definitions(self):
        for id in ["MolWt", "MolLogP", "BCUT2D_MWHI", "ds_Wiener", "SPS"]:
            self.assertEqual(relative_tolerance(id), 0)
            self.assertEqual(allowed_error(id, 1e30), tolerance(id))
        self.assertEqual(allowed_error("Ipc", 1e30), 1e18)
        self.assertEqual(allowed_error("AvgIpc", 0), 1e-12)

    def test_old_registered_and_additional_reference_sources_remain_explicit(self):
        mol = Chem.MolFromSmiles("CCO", PARAMS)
        spec = {"id": "MolWt", "origin": "rdkit-standard"}
        self.assertEqual(reference_source(spec), "RDKit2026.03.6")
        self.assertEqual(reference_value(spec, mol, {}, {}), FUNCTIONS["MolWt"](mol))
        self.assertEqual(reference_source({"id": "ds_Wiener", "origin": "native-extra"}),
                         "independent-extra50-formulas")


if __name__ == "__main__":
    unittest.main()
