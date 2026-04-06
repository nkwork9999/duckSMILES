# duckSmiles

Cheminformatics toolkit for DuckDB - SMILES, InChI, MOL/SDF, PDB, and SELFIES molecular analysis directly from SQL.

No external chemistry libraries required (pure Rust implementation, no RDKit dependency).

## Installation

```sql
INSTALL ducksmiles FROM community;
LOAD ducksmiles;
```

## Quick Start

```sql
-- Molecular formula from SMILES
SELECT mol_formula('CCO');
-- C2H6O

-- Molecular weight
SELECT round(mol_weight('c1ccccc1'), 2);
-- 78.11

-- Validate SMILES
SELECT mol_is_valid('CCO'), mol_is_valid('invalid');
-- true, false
```

## Functions

### SMILES (6 functions)

| Function | Return | Description |
|----------|--------|-------------|
| `mol_is_valid(smiles)` | BOOLEAN | Validate SMILES string |
| `mol_formula(smiles)` | VARCHAR | Molecular formula (Hill system) |
| `mol_num_atoms(smiles)` | INTEGER | Heavy atom count |
| `mol_num_bonds(smiles)` | INTEGER | Bond count |
| `mol_weight(smiles)` | DOUBLE | Average molecular weight |
| `mol_exact_mass(smiles)` | DOUBLE | Monoisotopic exact mass |

```sql
SELECT mol_formula('CC(=O)O'), mol_num_atoms('CC(=O)O'), round(mol_weight('CC(=O)O'), 2);
-- C2H4O2, 4, 60.05
```

### InChI (11 functions)

| Function | Return | Description |
|----------|--------|-------------|
| `inchi_is_valid(inchi)` | BOOLEAN | Validate InChI string |
| `inchi_is_standard(inchi)` | BOOLEAN | Check standard InChI (1S) |
| `inchi_version(inchi)` | VARCHAR | Version string |
| `inchi_formula(inchi)` | VARCHAR | Formula layer |
| `inchi_connections(inchi)` | VARCHAR | Connection layer (/c) |
| `inchi_hydrogens(inchi)` | VARCHAR | Hydrogen layer (/h) |
| `inchi_charge(inchi)` | VARCHAR | Charge layer (/q) |
| `inchi_stereo_bond(inchi)` | VARCHAR | Bond stereo (/b) |
| `inchi_stereo_tetrahedral(inchi)` | VARCHAR | Tetrahedral stereo (/t) |
| `inchi_has_stereo(inchi)` | BOOLEAN | Has stereochemistry |
| `inchi_num_stereo_centers(inchi)` | INTEGER | Number of stereocenters |

```sql
-- Stereochemistry detection (testosterone)
SELECT inchi_has_stereo('InChI=1S/C19H28O2/c1-18-9-7-13(20)11-12(18)3-4-14-15-5-6-17(21)19(15,2)10-8-16(14)18/h11,14-17,21H,3-10H2,1-2H3/t14-,15-,16-,17-,18-,19-/m0/s1');
-- true

SELECT inchi_num_stereo_centers('InChI=1S/C19H28O2/c1-18-9-7-13(20)11-12(18)3-4-14-15-5-6-17(21)19(15,2)10-8-16(14)18/h11,14-17,21H,3-10H2,1-2H3/t14-,15-,16-,17-,18-,19-/m0/s1');
-- 6
```

### InChIKey (4 functions)

| Function | Return | Description |
|----------|--------|-------------|
| `inchikey_is_valid(key)` | BOOLEAN | Validate InChIKey |
| `inchikey_connectivity(key)` | VARCHAR | First segment |
| `inchikey_stereo(key)` | VARCHAR | Second segment |
| `inchikey_protonation(key)` | VARCHAR | Third segment |

```sql
SELECT inchikey_connectivity('QTBSBXVTEAMEQO-UHFFFAOYSA-N');
-- QTBSBXVTEAMEQO
```

### Comparison (1 function)

| Function | Return | Description |
|----------|--------|-------------|
| `inchi_skeleton_match(a, b)` | BOOLEAN | Skeleton match (ignoring stereo) |

```sql
SELECT inchi_skeleton_match(
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
);
-- true
```

### MOL/SDF (6 functions)

| Function | Return | Description |
|----------|--------|-------------|
| `mol_block_formula(mol)` | VARCHAR | Formula from MOL block |
| `mol_block_weight(mol)` | DOUBLE | Weight from MOL block |
| `mol_block_num_atoms(mol)` | INTEGER | Atom count from MOL block |
| `mol_block_num_bonds(mol)` | INTEGER | Bond count from MOL block |
| `mol_block_name(mol)` | VARCHAR | Molecule name from header |
| `sdf_count(sdf)` | INTEGER | Count molecules in SDF |

### Structure - PDB/CIF/XYZ (4 functions)

| Function | Return | Description |
|----------|--------|-------------|
| `structure_atom_count(text)` | INTEGER | Atom count (auto-detect format) |
| `structure_chain_count(text)` | INTEGER | Chain count |
| `structure_residue_count(text)` | INTEGER | Residue count |
| `structure_model_count(text)` | INTEGER | Model count |

### SELFIES (3 functions)

| Function | Return | Description |
|----------|--------|-------------|
| `smiles_to_selfies(smiles)` | VARCHAR | Convert SMILES to SELFIES |
| `selfies_to_smiles(selfies)` | VARCHAR | Convert SELFIES to SMILES |
| `selfies_is_valid(selfies)` | BOOLEAN | Validate SELFIES string |

```sql
-- Roundtrip conversion
SELECT selfies_to_smiles(smiles_to_selfies('CCO'));
-- CCO
```

## Use Cases

- Molecular property extraction from chemical datasets
- Batch validation of SMILES/InChI in data pipelines
- Stereochemistry analysis and comparison
- Format conversion for ML pipelines (SMILES <-> SELFIES)
- Protein structure metadata extraction (PDB/CIF/XYZ)

## Architecture

- **Rust** (5 crates): Core molecular parsing and computation
- **C++**: DuckDB extension integration via FFI
- No external chemistry library dependencies

## License

MIT License
