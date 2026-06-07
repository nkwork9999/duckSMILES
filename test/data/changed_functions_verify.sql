-- Verify that every public SQL function added or materially changed in the
-- recent duckSmiles implementation can be loaded and executed from DuckDB.
--
-- Usage from the repository root:
--   ./build/release/duckdb -unsigned < test/data/changed_functions_verify.sql

LOAD './build/release/extension/ducksmiles/ducksmiles.duckdb_extension';

CREATE OR REPLACE TEMP MACRO mol_v2000() AS E'  Ethanol\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n> <ID>\nV2000-1\n\n';

CREATE OR REPLACE TEMP MACRO mol_v3000() AS (
    'ethanol_v3000' || chr(10) ||
    '  test' || chr(10) ||
    chr(10) ||
    '  0  0  0     0  0            999 V3000' || chr(10) ||
    'M  V30 BEGIN CTAB' || chr(10) ||
    'M  V30 COUNTS 3 2 0 0 0' || chr(10) ||
    'M  V30 BEGIN ATOM' || chr(10) ||
    'M  V30 1 C 0.0000 0.0000 0.0000 0' || chr(10) ||
    'M  V30 2 C 1.5400 0.0000 0.0000 0 CFG=0' || chr(10) ||
    'M  V30 3 O 2.3100 1.3300 1.0000 0' || chr(10) ||
    'M  V30 END ATOM' || chr(10) ||
    'M  V30 BEGIN BOND' || chr(10) ||
    'M  V30 1 1 1 2 CFG=0' || chr(10) ||
    'M  V30 2 2 2 3' || chr(10) ||
    'M  V30 END BOND' || chr(10) ||
    'M  V30 END CTAB' || chr(10) ||
    'M  END' || chr(10) ||
    '> <ID>' || chr(10) ||
    'V3000-1' || chr(10) ||
    chr(10)
);

CREATE OR REPLACE TEMP MACRO sdf_text() AS (
    mol_v2000() || repeat(chr(36), 4) || chr(10) ||
    mol_v3000() || repeat(chr(36), 4) || chr(10)
);

CREATE OR REPLACE TEMP MACRO pdb_text() AS E'ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\nATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C\nATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C\nATOM      4  N   GLY A   2       4.000   5.000   6.000  1.00  0.00           N\nATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00  0.00           C\nATOM      6  N   ALA B   1       6.000   7.000   8.000  1.00  0.00           N\nEND\n';

CREATE OR REPLACE TEMP TABLE ducksmiles_changed_checks AS
SELECT *
FROM (
    VALUES
        -- Canonical SMILES and aromaticity perception.
        ('canonical_smiles_kekule_aromatizes', canonical_smiles('C1=CC=CC=C1') = 'c1ccccc1'),
        ('canonical_smiles_aromatic_stable', canonical_smiles('c1ccccc1') = 'c1ccccc1'),

        -- Murcko scaffolds and ring systems.
        ('murcko_scaffold_toluene', murcko_scaffold('Cc1ccccc1') = 'c1ccccc1'),
        ('murcko_scaffold_aspirin', murcko_scaffold('CC(=O)Oc1ccccc1C(=O)O') = 'c1ccccc1'),
        ('generic_scaffold_pyridine', generic_scaffold('c1ccccn1') = 'C1CCCCC1'),
        ('ring_systems_json_benzene', ring_systems_json('c1ccccc1') LIKE '%"num_atoms":6%"aromatic":true%'),

        -- MolHash full method family.
        ('mol_hash_methods_complete', mol_hash_methods() = '["canonical_smiles","formula","net_charge","degree_vector","atom_bond_counts","element_graph","bond_order_graph","anonymous_graph","murcko_scaffold","generic_scaffold"]'),
        ('mol_hash_canonical_smiles', mol_hash('CC(=O)[O-].[Na+]', 'canonical_smiles') = 'C(C)(=O)[O-].[Na+]'),
        ('mol_hash_formula', mol_hash('CC(=O)[O-].[Na+]', 'formula') = 'C2H3NaO2'),
        ('mol_hash_net_charge', mol_hash('CC(=O)[O-].[Na+]', 'net_charge') = '0'),
        ('mol_hash_degree_vector', mol_hash('CC(=O)[O-].[Na+]', 'degree_vector') = '0,1,1,1,3'),
        ('mol_hash_atom_bond_counts', mol_hash('CC(=O)[O-].[Na+]', 'atom_bond_counts') = '5:3'),
        ('mol_hash_element_graph', mol_hash('CC(=O)[O-].[Na+]', 'element_graph') = 'C(C)(O)O.[Na]'),
        ('mol_hash_bond_order_graph', mol_hash('CC(=O)[O-].[Na+]', 'bond_order_graph') = 'C(C)(=O)O.[Na]'),
        ('mol_hash_anonymous_graph', mol_hash('CC(=O)[O-].[Na+]', 'anonymous_graph') = 'C.C(C)(C)C'),
        ('mol_hash_murcko_scaffold', mol_hash('Cc1ccccc1', 'murcko_scaffold') = 'c1ccccc1'),
        ('mol_hash_generic_scaffold', mol_hash('c1ccccn1', 'generic_scaffold') = 'C1CCCCC1'),
        ('mol_hash_unknown_method_null', mol_hash('CCO', 'does_not_exist') IS NULL),

        -- Standardization and salt handling.
        ('largest_fragment', largest_fragment('CC(=O)[O-].[Na+]') = 'C(C)(=O)[O-]'),
        ('strip_salts', strip_salts('CCO.CN.[Cl-]') = 'C(C)O.CN'),
        ('neutralize_charges', neutralize_charges('CC(=O)[O-]') = 'C(C)(=O)O'),
        ('normalize_smiles', normalize_smiles('C1=CC=CC=C1') = 'c1ccccc1'),
        ('fragment_parent', fragment_parent('CC(=O)[O-].[Na+]') = 'C(C)(=O)O'),

        -- SMARTS matcher expansion and substructure JSON.
        ('smarts_isotope_predicate', mol_has_substructure('[13CH4]', '[13C]') = true),
        ('smarts_atom_map_predicate', mol_has_substructure('[CH4:7]', '[C:7]') = true),
        ('smarts_chirality_predicate', mol_has_substructure('F[C@H](Cl)Br', '[C@H]') = true),
        ('smarts_chirality_negative', mol_has_substructure('F[C@H](Cl)Br', '[C@@H]') = false),
        ('smarts_degree_D', mol_substructure_count('CCC', '[D2]') = 1),
        ('smarts_valence_v', mol_substructure_count('CCO', '[v4]') = 2),
        ('smarts_ring_connectivity_x', mol_substructure_count('c1ccccc1', '[x2]') = 6),
        ('smarts_ring_size_r', mol_substructure_count('c1ccccc1', '[r6]') = 6),
        ('smarts_combined_r_x', mol_substructure_count('c1ccccc1', '[r6;x2]') = 6),
        ('smarts_ring_non_ring_bond', mol_has_substructure('c1ccccc1-c1ccccc1', '*@*!@*@*') = true),
        ('smarts_ring_non_ring_bond_negative', mol_has_substructure('c1ccccc1', '*@*!@*@*') = false),
        ('smarts_not_in_ring_bond', mol_has_substructure('CC', '[#6]!@[#6]') = true),
        ('smarts_not_in_ring_bond_negative', mol_has_substructure('C1CC1', '[#6]!@[#6]') = false),
        ('smarts_not_aromatic_bond', mol_has_substructure('CC', 'C!:C') = true),
        ('smarts_recursive', mol_has_substructure('CCO', '[$([#6]-[#8])]') = true),
        ('smarts_recursive_negative', mol_has_substructure('CCN', '[$([#6]-[#8])]') = false),
        ('mol_substructure_matches_json', mol_substructure_matches_json('CCC', '[#6]~[#6]') LIKE '%"match":2%"target_atom":3%'),
        ('mol_substructure_invalid_null', mol_substructure_matches_json('not_a_molecule', 'C') IS NULL),

        -- MCS and scaffold network.
        ('mcs_smarts', mcs_smarts('CC(=O)O', 'CC(=O)N') = 'C(C)=O'),
        ('mcs_json', mcs_json('CC(=O)O', 'CC(=O)N') LIKE '%"num_atoms":3%"smarts":"C(C)=O"%'),
        ('scaffold_network_json_murcko', scaffold_network_json('Cc1ccccc1') LIKE '%"kind":"murcko_scaffold"%"smiles":"c1ccccc1"%'),
        ('scaffold_network_json_generic', scaffold_network_json('Cc1ccccc1') LIKE '%"kind":"generic_scaffold"%"smiles":"C1CCCCC1"%'),

        -- Descriptors, fingerprints, and similarity.
        ('logp_crippen', round(logp_crippen('CCO'), 4) = -0.0014),
        ('tpsa', round(tpsa('CC(=O)O'), 2) = 37.30),
        ('num_h_acceptors', num_h_acceptors('CC(=O)O') = 2),
        ('num_h_donors', num_h_donors('CC(=O)O') = 1),
        ('num_rotatable_bonds', num_rotatable_bonds('CCCO') = 1),
        ('ring_count', ring_count('c1ccccc1') = 1),
        ('num_aromatic_rings', num_aromatic_rings('c1ccccc1') = 1),
        ('num_heteroatoms', num_heteroatoms('CC(=O)O') = 2),
        ('fraction_csp3', round(fraction_csp3('CCO'), 3) = 1.000),
        ('add_hydrogens_valid', mol_is_valid(add_hydrogens('CCO')) = true),
        ('add_hydrogens_formula', mol_formula(add_hydrogens(add_hydrogens('CCO'))) = 'C2H6O'),
        ('morgan_fp_bits_default_len', octet_length(morgan_fp_bits('CCO')) = 256),
        ('morgan_fp_bits_custom_len', octet_length(morgan_fp_bits('CCO', 2, 128)) = 16),
        ('maccs_keys_len', octet_length(maccs_keys('CCO')) = 21),
        ('tanimoto_bit_self', round(tanimoto_bit(morgan_fp_bits('CCO'), morgan_fp_bits('CCO')), 3) = 1.000),
        ('tanimoto_bit_diff', round(tanimoto_bit(morgan_fp_bits('CCO'), morgan_fp_bits('c1ccccc1')), 3) = 0.000),

        -- MOL/SDF/V3000 structural extraction.
        ('mol_block_name', mol_block_name(mol_v2000()) = 'Ethanol'),
        ('mol_block_formula', mol_block_formula(mol_v2000()) = 'C2O'),
        ('mol_block_weight', round(mol_block_weight(mol_v2000()), 2) = 40.02),
        ('mol_block_num_atoms', mol_block_num_atoms(mol_v2000()) = 3),
        ('mol_block_num_bonds', mol_block_num_bonds(mol_v2000()) = 2),
        ('mol_block_property', mol_block_property(mol_v2000(), 'ID') = 'V2000-1'),
        ('mol_block_properties_json', mol_block_properties_json(mol_v2000()) = '[{"name":"ID","value":"V2000-1"}]'),
        ('mol_block_atoms_json', mol_block_atoms_json(mol_v3000()) LIKE '%"index":3%"symbol":"O"%'),
        ('mol_block_bonds_json', mol_block_bonds_json(mol_v3000()) LIKE '%"index":2%"bond_type":2%'),
        ('mol_block_json', mol_block_json(mol_v3000()) LIKE '%"name":"ethanol_v3000"%"has_3d":true%'),
        ('mol_block_has_3d', mol_block_has_3d(mol_v3000()) = true),
        ('mol_block_centroid_x', round(mol_block_centroid_x(mol_v3000()), 3) = 1.283),
        ('mol_block_centroid_y', round(mol_block_centroid_y(mol_v3000()), 3) = 0.443),
        ('mol_block_centroid_z', round(mol_block_centroid_z(mol_v3000()), 3) = 0.333),
        ('mol_block_radius_of_gyration', round(mol_block_radius_of_gyration(mol_v3000()), 3) = 1.240),
        ('mol_block_min_x', mol_block_min_x(mol_v3000()) = 0.0),
        ('mol_block_max_x', round(mol_block_max_x(mol_v3000()), 2) = 2.31),
        ('mol_block_min_y', mol_block_min_y(mol_v3000()) = 0.0),
        ('mol_block_max_y', round(mol_block_max_y(mol_v3000()), 2) = 1.33),
        ('mol_block_min_z', mol_block_min_z(mol_v3000()) = 0.0),
        ('mol_block_max_z', mol_block_max_z(mol_v3000()) = 1.0),
        ('sdf_count', sdf_count(sdf_text()) = 2),
        ('sdf_property_record_1', sdf_property(sdf_text(), 1, 'ID') = 'V2000-1'),
        ('sdf_property_record_2', sdf_property(sdf_text(), 2, 'ID') = 'V3000-1'),
        ('sdf_properties_json', sdf_properties_json(sdf_text()) LIKE '%"record":2%"value":"V3000-1"%'),

        -- PDB/CIF/XYZ-style structure metrics through auto-detect PDB input.
        ('structure_atom_count', structure_atom_count(pdb_text()) = 6),
        ('structure_chain_count', structure_chain_count(pdb_text()) = 2),
        ('structure_residue_count', structure_residue_count(pdb_text()) = 3),
        ('structure_model_count', structure_model_count(pdb_text()) = 1),
        ('structure_centroid_x', round(structure_centroid_x(pdb_text()), 3) = 3.500),
        ('structure_centroid_y', round(structure_centroid_y(pdb_text()), 3) = 4.500),
        ('structure_centroid_z', round(structure_centroid_z(pdb_text()), 3) = 5.500),
        ('structure_radius_of_gyration', round(structure_radius_of_gyration(pdb_text()), 3) = 2.958),
        ('structure_min_x', structure_min_x(pdb_text()) = 1.0),
        ('structure_max_x', structure_max_x(pdb_text()) = 6.0),
        ('structure_min_y', structure_min_y(pdb_text()) = 2.0),
        ('structure_max_y', structure_max_y(pdb_text()) = 7.0),
        ('structure_min_z', structure_min_z(pdb_text()) = 3.0),
        ('structure_max_z', structure_max_z(pdb_text()) = 8.0)
) AS v(name, ok);

SELECT
    count(*) AS total_checks,
    sum(CASE WHEN ok THEN 1 ELSE 0 END) AS passed_checks,
    sum(CASE WHEN ok IS DISTINCT FROM true THEN 1 ELSE 0 END) AS failed_checks
FROM ducksmiles_changed_checks;

SELECT name
FROM ducksmiles_changed_checks
WHERE ok IS DISTINCT FROM true
ORDER BY name;

SELECT
    CASE
        WHEN count(*) = 0 THEN 'all duckSmiles changed-function SQL checks passed'
        ELSE error('duckSmiles changed-function verification failed: ' || string_agg(name, ', '))
    END AS result
FROM ducksmiles_changed_checks
WHERE ok IS DISTINCT FROM true;
