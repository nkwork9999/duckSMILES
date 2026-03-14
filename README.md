# duckSMILES

DuckDB extension for parsing [SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System) (Simplified Molecular Input Line Entry System) chemical notation directly from SQL.

Self-contained lightweight parser — no RDKit or external chemistry libraries required.

## Functions

| Function | Input → Output | Description |
|----------|---------------|-------------|
| `mol_is_valid(smiles)` | VARCHAR → BOOLEAN | Check if SMILES string is valid |
| `mol_formula(smiles)` | VARCHAR → VARCHAR | Molecular formula (Hill system, e.g. `C2H6O`) |
| `mol_num_atoms(smiles)` | VARCHAR → INTEGER | Heavy atom count (non-hydrogen) |
| `mol_num_bonds(smiles)` | VARCHAR → INTEGER | Bond count |

## Usage

```sql
LOAD 'smiles';

-- Validate SMILES
SELECT mol_is_valid('CCO');        -- true
SELECT mol_is_valid('not_valid');   -- false

-- Molecular formula
SELECT mol_formula('CCO');          -- C2H6O  (ethanol)
SELECT mol_formula('c1ccccc1');     -- C6H6   (benzene)
SELECT mol_formula('O');            -- H2O    (water)
SELECT mol_formula('CC(=O)O');      -- C2H4O2 (acetic acid)

-- Atom and bond counts
SELECT mol_num_atoms('c1ccccc1');   -- 6
SELECT mol_num_bonds('c1ccccc1');   -- 6

-- Batch processing
SELECT smiles, mol_formula(smiles) AS formula, mol_num_atoms(smiles) AS atoms
FROM (VALUES ('CCO'), ('c1ccccc1'), ('CC(=O)Oc1ccccc1C(=O)O')) AS t(smiles);
```

## Supported SMILES Features

- Organic subset atoms: `B C N O P S F Cl Br I`
- Aromatic atoms: `c n o s p b`
- Bracket atoms: `[NH3+]`, `[Fe+2]`, `[13C@@H]`
- Branches `()`, rings `1-9` / `%10-%99`, disconnected fragments `.`
- Implicit hydrogen (automatic, valence-based)
- Invalid SMILES → `false` / `NULL` (never crashes)

## Build

```bash
git clone --recurse-submodules https://github.com/nkwork9999/duckSMILES.git
cd duckSMILES
make release
```

Test:
```bash
./build/release/duckdb -c "
  LOAD 'build/release/extension/smiles/smiles.duckdb_extension';
  SELECT mol_formula('CCO');
"
```
