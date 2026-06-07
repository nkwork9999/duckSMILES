# duckSmiles

**SQL to Molecules.** Cheminformatics toolkit for DuckDB.

Analyze SMILES, InChI, MOL/SDF, PDB, and SELFIES molecular structures directly from SQL — no Python, no RDKit, no setup.

**A lightweight, zero-dependency alternative to RDKit for common cheminformatics SQL queries.**

Instead of setting up Python + RDKit + pandas for molecular property extraction, just:

```sql
SELECT mol_formula(smiles), mol_weight(smiles) FROM molecules;
```

No external chemistry libraries, no Python environment, no compilation of RDKit from source.
Pure Rust implementation compiled into a single DuckDB extension.

### Who is this for?

- **Data scientists** working with chemical datasets (ChEMBL, PubChem, ZINC)
- **ML engineers** building molecular feature pipelines (SMILES <-> SELFIES conversion)
- **Cheminformatics teams** needing fast batch validation and property extraction
- **Anyone** who wants molecular analysis without the RDKit setup overhead

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

-- Normalize a common Kekule aromatic form
SELECT canonical_smiles('C1=CC=CC=C1');
-- c1ccccc1
```

---

## Function Reference

### SMILES Functions

SMILES (Simplified Molecular Input Line Entry System) is the most widely used text notation for molecules in cheminformatics. These functions parse SMILES strings and extract molecular properties.

#### `mol_is_valid(smiles) -> BOOLEAN`

Checks whether a SMILES string is syntactically valid and represents a parseable molecule. Returns `false` for malformed strings, empty input, or unsupported notation. Never throws an error — safe for batch processing over dirty data.

```sql
SELECT mol_is_valid('CCO');           -- true  (ethanol)
SELECT mol_is_valid('c1ccccc1');      -- true  (benzene, aromatic notation)
SELECT mol_is_valid('[Na+].[Cl-]');   -- true  (sodium chloride, disconnected)
SELECT mol_is_valid('not_a_molecule');-- false
SELECT mol_is_valid('');              -- false
```

**Supported SMILES features:** organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I), aromatic atoms (c, n, o, s, p, b), bracket atoms (`[NH3+]`, `[Fe+2]`, `[13C@@H]`), branches, ring closures, disconnected fragments (`.`), bond types (single `-`, double `=`, triple `#`), implicit hydrogens (valence-based), and conservative six-membered Kekule aromaticity perception for common benzene/pyridine-like rings. This is intentionally smaller than RDKit's full sanitization, valence model, stereochemistry handling, and canonicalization stack.

#### `mol_formula(smiles) -> VARCHAR`

Returns the molecular formula in Hill system order (C first, H second, then remaining elements alphabetically). Includes implicit hydrogens calculated from standard valence rules. Returns `NULL` for invalid SMILES.

```sql
SELECT mol_formula('O');              -- H2O
SELECT mol_formula('CCO');            -- C2H6O
SELECT mol_formula('c1ccccc1');       -- C6H6
SELECT mol_formula('CC(=O)O');        -- C2H4O2
SELECT mol_formula('[Na+].[Cl-]');    -- ClNa
SELECT mol_formula('CC(=O)Oc1ccccc1C(=O)O'); -- C9H8O4 (aspirin)
SELECT mol_formula('invalid');        -- NULL
```

#### `mol_num_atoms(smiles) -> INTEGER`

Returns the count of heavy atoms (non-hydrogen atoms). Implicit and explicit hydrogens are excluded. Returns `NULL` for invalid SMILES.

```sql
SELECT mol_num_atoms('CCO');          -- 3  (C, C, O)
SELECT mol_num_atoms('c1ccccc1');     -- 6  (6 carbons)
SELECT mol_num_atoms('[Na+].[Cl-]');  -- 2  (Na, Cl)
```

#### `mol_num_bonds(smiles) -> INTEGER`

Returns the total number of bonds in the molecule. Double bonds count as 1 bond, triple bonds count as 1 bond (bond count, not bond order sum). Ring closures are included. Returns `NULL` for invalid SMILES.

```sql
SELECT mol_num_bonds('CCO');          -- 2  (C-C, C-O)
SELECT mol_num_bonds('c1ccccc1');     -- 6  (benzene ring)
SELECT mol_num_bonds('C=C');          -- 1  (one double bond)
SELECT mol_num_bonds('C#C');          -- 1  (one triple bond)
```

#### `mol_weight(smiles) -> DOUBLE`

Returns the average molecular weight in Da (Daltons), using standard atomic weights (e.g., C=12.011, H=1.008, O=15.999). Includes the mass contribution of implicit hydrogens. Returns `NULL` (NaN internally) for invalid SMILES.

```sql
SELECT round(mol_weight('O'), 2);         -- 18.02  (water)
SELECT round(mol_weight('CCO'), 2);       -- 46.07  (ethanol)
SELECT round(mol_weight('c1ccccc1'), 2);  -- 78.11  (benzene)
SELECT round(mol_weight('CC(=O)Oc1ccccc1C(=O)O'), 2); -- 180.16 (aspirin)
```

#### `mol_exact_mass(smiles) -> DOUBLE`

Returns the monoisotopic exact mass, using the mass of the most abundant isotope of each element (e.g., C=12.0000, H=1.00783, O=15.99491). Used in mass spectrometry analysis. Returns `NULL` for invalid SMILES.

```sql
SELECT round(mol_exact_mass('O'), 4);     -- 18.0106
SELECT round(mol_exact_mass('CCO'), 4);   -- 46.0419
SELECT round(mol_exact_mass('C(Cl)(Cl)Cl'), 4); -- 117.9144 (chloroform)
```

#### `logp_crippen(smiles) -> DOUBLE`

Returns the **Wildman–Crippen LogP** (octanol/water partition coefficient prediction) using the atom-contribution method (S. A. Wildman, G. M. Crippen, *JCICS* 39, 868–873 (1999)). Implemented as a port of RDKit's `Crippen.txt` parameter table — 110 SMARTS patterns, 68 atom types, first-match-wins atom typing. Values match RDKit's `Crippen.MolLogP` exactly for small molecules and stay within the method's intrinsic ±0.68 log-unit accuracy on drug-like molecules. Returns `NULL` for invalid SMILES.

```sql
SELECT round(logp_crippen('C'), 4);            -- 0.6361   (methane)
SELECT round(logp_crippen('O'), 4);            -- -0.8247  (water)
SELECT round(logp_crippen('c1ccccc1'), 4);     -- 1.6866   (benzene)
SELECT round(logp_crippen('CCO'), 4);          -- -0.0014  (ethanol)
SELECT round(logp_crippen('CC(=O)Oc1ccccc1C(=O)O'), 4); -- ~1.31  (aspirin)
```

#### `tpsa(smiles) -> DOUBLE`

Returns the **Topological Polar Surface Area** using RDKit's default Ertl-style atom-contribution scope: nitrogen and oxygen atoms are included; sulfur and phosphorus are excluded. Useful for drug-likeness filters, permeability heuristics, and molecular ML features. Returns `NULL` for invalid SMILES.

```sql
SELECT round(tpsa('O'), 2);                            -- 31.50  (water)
SELECT round(tpsa('CCO'), 2);                          -- 20.23  (ethanol)
SELECT round(tpsa('CC(=O)O'), 2);                      -- 37.30  (acetic acid)
SELECT round(tpsa('CC(=O)Oc1ccccc1C(=O)O'), 2);        -- 63.60  (aspirin)
SELECT round(tpsa('c1ccncc1'), 2);                     -- 12.89  (pyridine)
```

#### `canonical_smiles(smiles) -> VARCHAR`

Returns a deterministic normalized SMILES string for the parser's supported graph subset. It also runs the conservative aromaticity perception pass, so common six-membered Kekule aromatic rings are emitted in lowercase aromatic form. Returns `NULL` for invalid SMILES.

This is useful for SQL-side deduplication and for normalizing common descriptor inputs, but it is not a replacement for RDKit's full canonical SMILES implementation: stereochemical canonicalization, full sanitization, and broad aromaticity models remain outside this lightweight subset.

```sql
SELECT canonical_smiles('C1=CC=CC=C1');  -- c1ccccc1
SELECT canonical_smiles('c1ccccc1');     -- c1ccccc1
SELECT canonical_smiles('C1CCCCC1');     -- C1CCCCC1
```

#### `add_hydrogens(smiles) -> VARCHAR`

Returns the input SMILES with every implicit H atom rewritten as an explicit `[H]` vertex (verbose bracket form). The output round-trips through `parse` and composes safely with other descriptors. Useful as a preprocessing primitive for SMARTS that match `[#1]`. Returns `NULL` for invalid SMILES.

```sql
SELECT add_hydrogens('CCO');
-- [C]([C]([O][H])([H])[H])([H])([H])[H]
```

#### `morgan_fp_bits(smiles [, radius, n_bits]) -> BLOB`

Returns a **Morgan / ECFP fingerprint** as a fixed-width bit vector (`ceil(n_bits/8)` bytes of BLOB). Defaults are radius=2 and n_bits=2048, i.e. **ECFP4 / 2048-bit** (the most common ML featurization choice). The 3-arg overload exposes full control: `radius` (ECFPn → radius=n/2) and `n_bits` (bit-vector width). Returns `NULL` for invalid SMILES.

Algorithm: layered BFS over atom neighborhoods + `hash_combine` + dead-atom dedup, ported from RDKit's `MorganGenerator`. **Not bit-exact RDKit-compatible** (uses an independent hash function), but algorithmically equivalent.

```sql
-- ECFP4 (default): inspect popcount via CAST to BIT
SELECT bit_count(CAST(morgan_fp_bits('CC(=O)Oc1ccccc1C(=O)O') AS BIT));  -- 26  (aspirin)

-- ECFP6 with 4096 bits
SELECT bit_count(CAST(morgan_fp_bits('CC(=O)Oc1ccccc1C(=O)O', 3, 4096) AS BIT));  -- 33

-- Tanimoto similarity via BIT operators (see also: tanimoto_bit below)
WITH x AS (
  SELECT CAST(morgan_fp_bits('CCO') AS BIT) AS a,
         CAST(morgan_fp_bits('CCN') AS BIT) AS b
)
SELECT bit_count(a & b)::DOUBLE / bit_count(a | b) AS tanimoto FROM x;
-- 0.3333
```

#### `tanimoto_bit(blob_a, blob_b) -> DOUBLE`

Computes **Tanimoto similarity** between two equal-length fingerprint BLOBs:
`popcount(a & b) / popcount(a | b)`. Operates directly on the raw BLOB bytes
(no `CAST AS BIT` round-trip), processing 8 bytes at a time via `u64::count_ones()`
so it lowers to POPCNT on x86_64 and CNT on aarch64. Algorithmically equivalent
to RDKit's `CalcBitmapTanimoto`.

- **Length mismatch** → `InvalidInputException` with both byte sizes in the message
  (better than a silent NULL — mixing different `n_bits` is a clear user error).
- **Both BLOBs all-zero** → `0.0` (matches RDKit's `union == 0 ? 0.0` convention).
- Bit-exact identical to the SQL-level `bit_count(CAST(a AS BIT) & CAST(b AS BIT))
  / bit_count(CAST(a AS BIT) | CAST(b AS BIT))` form, but skips two casts per row.

```sql
-- Pairwise similarity vs a reference (aspirin)
WITH ref AS (SELECT morgan_fp_bits('CC(=O)Oc1ccccc1C(=O)O') AS fp)
SELECT name,
       round(tanimoto_bit(morgan_fp_bits(smiles), (SELECT fp FROM ref)), 4) AS sim
FROM (VALUES
  ('aspirin',        'CC(=O)Oc1ccccc1C(=O)O'),
  ('salicylic acid', 'OC(=O)c1ccccc1O'),
  ('paracetamol',    'CC(=O)Nc1ccc(O)cc1'),
  ('caffeine',       'Cn1c(=O)c2c(ncn2C)n(C)c1=O'),
  ('methane',        'C')
) AS t(name, smiles)
ORDER BY sim DESC;
-- aspirin        1.0000
-- salicylic acid 0.2162
-- paracetamol    0.2162
-- caffeine       0.1020
-- methane        0.0000

-- Specialized vs SQL-native — bit-exact identical:
WITH q AS (SELECT morgan_fp_bits('CCO') AS a, morgan_fp_bits('CCN') AS b)
SELECT tanimoto_bit(a, b)
     = bit_count(CAST(a AS BIT) & CAST(b AS BIT))::DOUBLE
     / bit_count(CAST(a AS BIT) | CAST(b AS BIT)) AS match
FROM q;
-- true
```

#### `maccs_keys(smiles) -> BLOB`

Returns the **166 MACCS structural keys** as a fixed **21-byte (167-bit)** BLOB, ported from RDKit's `MACCS.cpp`. Bit `n` corresponds to RDKit MACCS key number `n` (bit 0 unused; bit 1, the isotope key, always 0). Returns `NULL` for invalid SMILES.

Unlike Morgan/ECFP (hashed local environments), each MACCS bit is a fixed yes/no structural question: most are SMARTS substructure matches, some are count thresholds (e.g. "≥2 oxygens"), and a few are special rules (per-atom element scan, "≥2 aromatic rings", multi-fragment). Lighter and more interpretable than Morgan, at coarser resolution.

**Bit-for-bit verified** against RDKit's `MACCSkeys.GenMACCSKeys` for aromatic SMILES written in lowercase form (`c1ccccc1`). Common six-membered Kekule aromatic forms (`C1=CC=CC=C1`) are now normalized by the conservative aromaticity perception pass, but broader RDKit aromaticity and sanitization behavior is intentionally not claimed.

```sql
-- 21-byte fixed width (vs Morgan's 256 bytes at 2048 bit)
SELECT octet_length(maccs_keys('CCO'));  -- 21

-- Drops straight into tanimoto_bit alongside Morgan
WITH ref AS (SELECT maccs_keys('CC(=O)Oc1ccccc1C(=O)O') AS fp)
SELECT name,
       round(tanimoto_bit(maccs_keys(smiles), (SELECT fp FROM ref)), 4) AS sim_maccs
FROM (VALUES
  ('aspirin',        'CC(=O)Oc1ccccc1C(=O)O'),
  ('salicylic acid', 'O=C(O)c1ccccc1O'),
  ('acetaminophen',  'CC(=O)Nc1ccc(O)cc1'),
  ('benzene',        'c1ccccc1'),
  ('ethanol',        'CCO')
) AS t(name, smiles)
ORDER BY sim_maccs DESC;
```

---

### InChI Functions (11)

InChI (International Chemical Identifier) is a layered textual representation developed by IUPAC. Each layer encodes different structural information (formula, connectivity, hydrogens, charge, stereochemistry). These functions parse InChI strings by layer extraction — no external InChI library required.

#### `inchi_is_valid(inchi) -> BOOLEAN`

Checks whether the string starts with `InChI=` and has a parseable structure.

```sql
SELECT inchi_is_valid('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- true
SELECT inchi_is_valid('not_an_inchi');                             -- false
```

#### `inchi_is_standard(inchi) -> BOOLEAN`

Checks whether the InChI uses the standard version prefix `InChI=1S/`. Standard InChI (`1S`) is the most widely used; non-standard (`1/`) may include additional layers.

```sql
SELECT inchi_is_standard('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- true
SELECT inchi_is_standard('InChI=1/C2H4O2/c1-2(3)4/h1H3,(H,3,4)');  -- false
```

#### `inchi_version(inchi) -> VARCHAR`

Returns the InChI version string (e.g., `1S` for standard, `1` for non-standard).

```sql
SELECT inchi_version('InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H'); -- 1S
```

#### `inchi_formula(inchi) -> VARCHAR`

Extracts the molecular formula layer — the first layer after the version prefix.

```sql
SELECT inchi_formula('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)');  -- C2H4O2
SELECT inchi_formula('InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2'); -- C6H12O6
SELECT inchi_formula('invalid');                                    -- NULL
```

#### `inchi_connections(inchi) -> VARCHAR`

Extracts the connection layer (`/c`), which describes the atom-to-atom connectivity (bond graph) of the molecule without hydrogen positions.

```sql
SELECT inchi_connections('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- 1-2(3)4
```

#### `inchi_hydrogens(inchi) -> VARCHAR`

Extracts the hydrogen layer (`/h`), which encodes fixed and mobile hydrogen positions. Includes both explicit H counts and mobile (tautomeric) H notation.

```sql
SELECT inchi_hydrogens('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- 1H3,(H,3,4)
```

#### `inchi_charge(inchi) -> VARCHAR`

Extracts the charge layer (`/q`). Returns the net formal charge of the molecule.

```sql
SELECT inchi_charge('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- NULL (neutral)
```

#### `inchi_stereo_bond(inchi) -> VARCHAR`

Extracts the bond stereochemistry layer (`/b`), encoding E/Z (cis/trans) double bond configuration.

```sql
SELECT inchi_stereo_bond('InChI=1S/C2H4/c1-2/h1-2H2/b2-1+'); -- 2-1+
```

#### `inchi_stereo_tetrahedral(inchi) -> VARCHAR`

Extracts the tetrahedral stereochemistry layer (`/t`), encoding R/S chirality at sp3 centers. Each stereocenter is listed with its parity (`+` or `-`).

```sql
-- Testosterone: 6 stereocenters
SELECT inchi_stereo_tetrahedral(
    'InChI=1S/C19H28O2/c1-18-9-7-13(20)11-12(18)3-4-14-15-5-6-17(21)19(15,2)10-8-16(14)18/h11,14-17,21H,3-10H2,1-2H3/t14-,15-,16-,17-,18-,19-/m0/s1'
); -- 14-,15-,16-,17-,18-,19-
```

#### `inchi_has_stereo(inchi) -> BOOLEAN`

Returns `true` if the InChI contains any stereochemistry layer (`/b`, `/t`, `/m`, or `/s`).

```sql
SELECT inchi_has_stereo('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- false (acetic acid)
SELECT inchi_has_stereo('InChI=1S/C19H28O2/c1-18-9-7-13(20)11-12(18)3-4-14-15-5-6-17(21)19(15,2)10-8-16(14)18/h11,14-17,21H,3-10H2,1-2H3/t14-,15-,16-,17-,18-,19-/m0/s1'); -- true (testosterone)
```

#### `inchi_num_stereo_centers(inchi) -> INTEGER`

Counts the number of tetrahedral stereocenters defined in the `/t` layer.

```sql
SELECT inchi_num_stereo_centers('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'); -- 0
SELECT inchi_num_stereo_centers('InChI=1S/C19H28O2/c1-18-9-7-13(20)11-12(18)3-4-14-15-5-6-17(21)19(15,2)10-8-16(14)18/h11,14-17,21H,3-10H2,1-2H3/t14-,15-,16-,17-,18-,19-/m0/s1'); -- 6
```

---

### InChIKey Functions (4)

InChIKey is a fixed-length (27-character) hash of an InChI string, formatted as `XXXXXXXXXXXXXX-XXXXXXXXXX-X` (three segments separated by hyphens). It enables fast exact-match lookups in databases.

#### `inchikey_is_valid(key) -> BOOLEAN`

Validates the InChIKey format: 14 uppercase letters, hyphen, 10 uppercase letters, hyphen, 1 uppercase letter.

```sql
SELECT inchikey_is_valid('QTBSBXVTEAMEQO-UHFFFAOYSA-N'); -- true
SELECT inchikey_is_valid('not-a-key');                    -- false
```

#### `inchikey_connectivity(key) -> VARCHAR`

Returns the first 14-character segment, which encodes the molecular connectivity (atom types and bonds, ignoring stereochemistry and charge). Two molecules with the same connectivity hash are constitutional isomers or identical.

```sql
SELECT inchikey_connectivity('QTBSBXVTEAMEQO-UHFFFAOYSA-N'); -- QTBSBXVTEAMEQO
```

#### `inchikey_stereo(key) -> VARCHAR`

Returns the second 10-character segment, encoding stereochemistry and isotope information. `UHFFFAOYSA` means "no stereochemistry defined."

```sql
SELECT inchikey_stereo('QTBSBXVTEAMEQO-UHFFFAOYSA-N');   -- UHFFFAOYSA (no stereo)
SELECT inchikey_stereo('WQZGKKKJIJFFOK-GASJEMHNSA-N');   -- GASJEMHNSA (glucose, has stereo)
```

#### `inchikey_protonation(key) -> VARCHAR`

Returns the third 1-character segment, encoding the protonation state. `N` = neutral, `M` = deprotonated (-1), `O` = protonated (+1), etc.

```sql
SELECT inchikey_protonation('QTBSBXVTEAMEQO-UHFFFAOYSA-N'); -- N (neutral)
```

---

### Comparison Function (1)

#### `inchi_skeleton_match(inchi_a, inchi_b) -> BOOLEAN`

Compares two InChI strings by their molecular skeleton: formula layer + connection layer + hydrogen layer. Stereochemistry, charge, and isotope differences are ignored. Useful for finding constitutional isomers or confirming structural identity regardless of stereochemistry.

```sql
-- Same molecule
SELECT inchi_skeleton_match(
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
); -- true

-- Different molecules (ethanol vs acetic acid)
SELECT inchi_skeleton_match(
    'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
    'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'
); -- false
```

---

### MOL/SDF Functions (6)

MOL blocks (V2000/V3000) are the standard file format for storing 2D/3D molecular structures with explicit atom coordinates. SDF (Structure Data Format) concatenates multiple MOL blocks with associated properties. These functions parse MOL/SDF text directly from VARCHAR columns.

#### `mol_block_name(mol) -> VARCHAR`

Extracts the molecule name from the first line of the MOL block header.

#### `mol_block_formula(mol) -> VARCHAR`

Computes the molecular formula from the atom block of the MOL file (counting each element symbol).

#### `mol_block_weight(mol) -> DOUBLE`

Computes the molecular weight from the atom types listed in the MOL block.

#### `mol_block_num_atoms(mol) -> INTEGER`

Returns the atom count parsed from the counts line or the atom block of the MOL file.

#### `mol_block_num_bonds(mol) -> INTEGER`

Returns the bond count parsed from the counts line or the bond block of the MOL file.

```sql
-- All MOL block functions (requires escaped newlines in SQL)
SELECT mol_block_name(mol_text) AS name,
       mol_block_formula(mol_text) AS formula,
       round(mol_block_weight(mol_text), 2) AS weight,
       mol_block_num_atoms(mol_text) AS atoms,
       mol_block_num_bonds(mol_text) AS bonds
FROM molecules;
```

#### `sdf_count(sdf) -> INTEGER`

Counts the number of molecules in an SDF file by counting `$$$$` delimiters. Useful for quickly checking dataset size without fully parsing each molecule.

```sql
SELECT sdf_count(sdf_text) AS num_molecules FROM sdf_files;
```

---

### Structure Functions — PDB/CIF/XYZ (4)

These functions parse protein/crystal structure files. The input format (PDB, mmCIF, or XYZ) is auto-detected from the content — no format flag needed.

**PDB** (Protein Data Bank): Standard format for macromolecular structures. Parsed via `ATOM`/`HETATM` records.

**mmCIF** (Macromolecular Crystallographic Information File): Newer format with `data_` prefix and `_atom_site.` schema. Used by the wwPDB as the primary archive format.

**XYZ**: Simple format listing atom symbol + x/y/z coordinates. Common in computational chemistry.

#### `structure_atom_count(text) -> INTEGER`

Counts atoms in the structure. For PDB: counts `ATOM` and `HETATM` records. For CIF: counts rows in `_atom_site` loop. For XYZ: reads atom count from header.

#### `structure_chain_count(text) -> INTEGER`

Counts unique chain identifiers (e.g., chain A, B, C in a protein complex). For XYZ format, returns 0 (no chain concept).

#### `structure_residue_count(text) -> INTEGER`

Counts unique residues (amino acids, nucleotides, ligands). Identified by unique (chain, residue number, residue name) combinations.

#### `structure_model_count(text) -> INTEGER`

Counts MODEL/ENDMDL blocks (NMR ensembles or multi-model files). Returns 1 for single-model structures.

```sql
-- Analyze a PDB file loaded as text
SELECT structure_atom_count(content) AS atoms,
       structure_chain_count(content) AS chains,
       structure_residue_count(content) AS residues,
       structure_model_count(content) AS models
FROM read_text('protein.pdb');
```

---

### SELFIES Functions (3)

SELFIES (Self-Referencing Embedded Strings) is a 100% valid molecular string representation designed for machine learning. Unlike SMILES, every SELFIES string decodes to a valid molecule — making it ideal for generative models (VAE, GPT, diffusion) where random token mutations must always produce valid chemistry.

#### `smiles_to_selfies(smiles) -> VARCHAR`

Converts a SMILES string to SELFIES notation. Each atom and structural feature is wrapped in brackets (`[C]`, `[=O]`, `[Branch1]`, `[Ring1]`).

```sql
SELECT smiles_to_selfies('CCO');            -- [C][C][O]
SELECT smiles_to_selfies('c1ccccc1');       -- [c][Ring1][C][c][c][c][c][c][Ring1][C]
SELECT smiles_to_selfies('CC(=O)O');        -- [C][C][Branch1][C][=O][O]
```

#### `selfies_to_smiles(selfies) -> VARCHAR`

Converts a SELFIES string back to SMILES notation. Enables roundtrip conversion for ML pipeline integration.

```sql
SELECT selfies_to_smiles('[C][C][O]');      -- CCO
SELECT selfies_to_smiles('[C][=O]');        -- C=O
```

#### `selfies_is_valid(selfies) -> BOOLEAN`

Checks whether a SELFIES string can be decoded to a valid molecule.

```sql
SELECT selfies_is_valid('[C][C][O]');       -- true
SELECT selfies_is_valid(smiles_to_selfies('CC(=O)O')); -- true
```

**Roundtrip example:**

```sql
-- SMILES -> SELFIES -> SMILES
SELECT smiles,
       smiles_to_selfies(smiles) AS selfies,
       selfies_to_smiles(smiles_to_selfies(smiles)) AS roundtrip
FROM (VALUES ('CCO'), ('CC(=O)O'), ('NCC(=O)O')) AS t(smiles);
```

---

## Real-World Examples: Kaggle Competition Data

### Molecular Data ML (Kaggle 2025) — 42/42 molecules parsed

```sql
-- Analyze photostability dataset directly from CSV
SELECT Batch_ID,
       Smiles,
       mol_formula(Smiles) AS formula,
       mol_num_atoms(Smiles) AS atoms,
       round(mol_weight(Smiles), 2) AS weight
FROM read_csv('train.csv');
```

### NeurIPS 2025 Open Polymer Prediction — 7,973/7,973 molecules parsed

```sql
-- Polymer SMILES use * for attachment points — replace with [H]
SELECT id,
       mol_formula(replace(SMILES, '*', '[H]')) AS formula,
       mol_num_atoms(replace(SMILES, '*', '[H]')) AS atoms,
       round(mol_weight(replace(SMILES, '*', '[H]')), 2) AS weight
FROM read_csv('train.csv');
```

Tested with real competition datasets containing complex aromatic systems, stereocenters, charged species, and polymer notation.

## Use Cases

- Molecular property extraction from chemical datasets (ChEMBL, PubChem, ZINC)
- Kaggle/NeurIPS competition feature engineering from SMILES
- Batch validation of SMILES/InChI in data pipelines
- Stereochemistry analysis and comparison
- Format conversion for ML pipelines (SMILES <-> SELFIES)
- Protein structure metadata extraction (PDB/CIF/XYZ)

## Architecture

- **Rust** (5 crates, ~3,200 lines): Core molecular parsing and computation
- **C++** (~330 lines): DuckDB extension integration via FFI
- No external chemistry library dependencies (no RDKit, no OpenBabel)
- Vectorized execution via DuckDB's `UnaryExecutor`
- Invalid input returns `NULL` — never crashes

## License

MIT License
