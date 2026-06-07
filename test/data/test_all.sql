-- ==========================================================================
-- duckSmiles full function test via CSV
-- Usage: ./build/release/duckdb -unsigned -cmd "LOAD 'build/release/extension/ducksmiles/ducksmiles.duckdb_extension';" < test/data/test_all.sql
-- ==========================================================================

-- ============================
-- 1. SMILES functions
-- ============================
.print '=== SMILES: mol_is_valid, mol_formula, mol_num_atoms, mol_num_bonds, mol_weight, mol_exact_mass ==='

SELECT name,
       smiles,
       mol_is_valid(smiles) AS valid,
       mol_formula(smiles) AS formula,
       mol_num_atoms(smiles) AS atoms,
       mol_num_bonds(smiles) AS bonds,
       round(mol_weight(smiles), 2) AS weight,
       round(mol_exact_mass(smiles), 4) AS exact_mass
FROM read_csv('test/data/smiles.csv');

-- Verify invalid SMILES returns false/NULL
.print '=== SMILES: invalid input handling ==='
SELECT mol_is_valid('not_a_molecule') AS should_be_false,
       mol_formula('not_a_molecule') AS should_be_null,
       mol_num_atoms('not_a_molecule') AS should_be_null2;

.print '=== SMILES: Bemis-Murcko scaffold functions ==='

SELECT murcko_scaffold('CC(=O)Oc1ccccc1C(=O)O') AS aspirin_scaffold,
       generic_scaffold('c1ccccn1') AS pyridine_generic,
       ring_systems_json('c1ccccc1') AS benzene_ring_system;

.print '=== SMILES/SMARTS: substructure match JSON ==='

SELECT mol_has_substructure('CC(=O)O', 'C=O') AS has_carbonyl,
       mol_substructure_count('CCC', '[#6]~[#6]') AS cc_bonds,
       mol_substructure_matches_json('CC(=O)O', 'C=O') AS carbonyl_matches;

-- ============================
-- 2. InChI functions (11)
-- ============================
.print '=== InChI: inchi_is_valid, inchi_is_standard, inchi_version, inchi_formula ==='

SELECT name,
       inchi_is_valid(inchi) AS valid,
       inchi_is_standard(inchi) AS standard,
       inchi_version(inchi) AS version,
       inchi_formula(inchi) AS formula
FROM read_csv('test/data/inchi.csv');

.print '=== InChI: inchi_connections, inchi_hydrogens, inchi_charge ==='

SELECT name,
       inchi_connections(inchi) AS connections,
       inchi_hydrogens(inchi) AS hydrogens,
       inchi_charge(inchi) AS charge
FROM read_csv('test/data/inchi.csv');

.print '=== InChI: stereo functions ==='

SELECT name,
       inchi_has_stereo(inchi) AS has_stereo,
       inchi_num_stereo_centers(inchi) AS stereo_centers,
       inchi_stereo_bond(inchi) AS stereo_bond,
       inchi_stereo_tetrahedral(inchi) AS stereo_tet
FROM read_csv('test/data/inchi.csv');

-- Invalid InChI
.print '=== InChI: invalid input ==='
SELECT inchi_is_valid('not_an_inchi') AS should_be_false,
       inchi_formula('not_an_inchi') AS should_be_null;

-- ============================
-- 3. InChIKey functions (4)
-- ============================
.print '=== InChIKey: inchikey_is_valid, connectivity, stereo, protonation ==='

SELECT name,
       inchikey_is_valid(inchikey) AS valid,
       inchikey_connectivity(inchikey) AS connectivity,
       inchikey_stereo(inchikey) AS stereo,
       inchikey_protonation(inchikey) AS protonation
FROM read_csv('test/data/inchikey.csv');

-- Invalid InChIKey
.print '=== InChIKey: invalid input ==='
SELECT inchikey_is_valid('not-a-key') AS should_be_false;

-- ============================
-- 4. Comparison (1)
-- ============================
.print '=== inchi_skeleton_match ==='

SELECT inchi_skeleton_match(
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
) AS same_should_be_true,
inchi_skeleton_match(
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
    'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'
) AS diff_should_be_false;

-- ============================
-- 5. SELFIES functions (3)
-- ============================
.print '=== SELFIES: smiles_to_selfies, selfies_to_smiles, selfies_is_valid ==='

SELECT name,
       smiles,
       smiles_to_selfies(smiles) AS selfies,
       selfies_is_valid(smiles_to_selfies(smiles)) AS selfies_valid,
       selfies_to_smiles(smiles_to_selfies(smiles)) AS roundtrip
FROM read_csv('test/data/selfies.csv');

-- ============================
-- 6. MOL/SDF functions - inline test
-- ============================
.print '=== MOL block functions ==='

SELECT mol_block_name(E'  Ethanol\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n') AS name,
mol_block_formula(E'  Ethanol\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n') AS formula,
mol_block_num_atoms(E'  Ethanol\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n') AS atoms,
mol_block_num_bonds(E'  Ethanol\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n') AS bonds;

.print '=== MOL/SDF structural JSON ==='

WITH mol AS (
    SELECT 'ethanol_v3000' || chr(10) || '  test' || chr(10) || chr(10) ||
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
           '> <ID>' || chr(10) || 'V3000-1' || chr(10) || chr(10) AS block
)
SELECT mol_block_atoms_json(block) AS atoms,
       mol_block_bonds_json(block) AS bonds,
       mol_block_json(block) AS mol_json
FROM mol;

-- ============================
-- 7. Structure functions (4) - PDB inline test
-- ============================
.print '=== Structure (PDB) functions ==='

SELECT structure_atom_count('ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C
ATOM      4  N   GLY A   2       4.000   5.000   6.000  1.00  0.00           N
ATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00  0.00           C
ATOM      6  N   ALA B   1       6.000   7.000   8.000  1.00  0.00           N
END
') AS atom_count,
structure_chain_count('ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C
ATOM      4  N   GLY A   2       4.000   5.000   6.000  1.00  0.00           N
ATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00  0.00           C
ATOM      6  N   ALA B   1       6.000   7.000   8.000  1.00  0.00           N
END
') AS chain_count,
structure_residue_count('ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C
ATOM      4  N   GLY A   2       4.000   5.000   6.000  1.00  0.00           N
ATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00  0.00           C
ATOM      6  N   ALA B   1       6.000   7.000   8.000  1.00  0.00           N
END
') AS residue_count;

.print '=== ALL TESTS COMPLETE ==='
