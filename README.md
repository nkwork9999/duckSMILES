# duckSmiles

**SQL to Molecules.** Cheminformatics toolkit for DuckDB.

[CPU analytics](#cpu-analytics-complete-sql-api-reference) adds statistical inference,
regression and GLMs, multivariate analysis, experimental design, quality/MSA,
distribution and survival methods, time-series forecasting and bounded predictive
models to molecular SQL. Every registered analytics function, parameter contract
and runnable example appears below; [JMP coverage](JMP_PARITY.md) records remaining work.

Formula-first reference indexing and explicit adduct mass conversions are
available through [`formula_exact_mass`, `adduct_neutral_mass`, and
`adduct_precursor_mz`](MASS_CONTEXT.md).

[`mol_graph_json`](NATIVE_GRAPH.md) exports a versioned native atom/bond graph
sidecar with input atom order, H, charge, isotope and raw stereo metadata. It is
not a pretrained GLACIER/ICEBERG tensor replacement.

Analyze SMILES, SMARTS, InChI, MOL/SDF, PDB, and SELFIES molecular structures directly from SQL — no Python, no RDKit, no setup.

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

The 0.6.0 source adds 254 CPU analytics functions (405 total SQL overloads).
Community packages depend on the version and source revision accepted upstream.
Before submission, the DuckDB 2 CI revision
`d54cbbe3df17e8e5ce68312f483b4155836112aa` passed all 13 SQL cases
(441 assertions) and all 405 catalog examples on macOS ARM64. The pinned Wasm
Rust compiler, `nightly-2025-02-01`, also passed the affected 58 analytics tests.
These local checks do not establish Wasm runtime SQL coverage or publication.
The submission fixes preserve argument names across DuckDB 1.5 and 2 and avoid
integer APIs unavailable in that pinned Rust compiler.

CI also builds and preloads JSON for the SQL regression suite. A required-extension
test loads the built DuckSMILES file and executes a JSON query, so missing
prerequisites cause a failure instead of silently skipping the tests. The DuckDB 2
preflight passes 14 SQL cases / 443 assertions. Community CI additionally loads
the built extension on test initialization and marks both extensions as required.

When CI supplies a prebuilt native DuckDB library, the build checks its embedded
source ID and version against DuckDB's version record for that commit before
setting the extension's CMake version. Missing or inconsistent metadata fails the
build; runtime version and ABI checks remain enabled. This handles prerelease
libraries supplied under a stable CI label.

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

-- SQL-native dedup/grouping hash
SELECT mol_hash('CC(=O)[O-].[Na+]', 'element_graph');
-- C(C)(O)O.[Na]

-- Salt handling and parent normalization
SELECT fragment_parent('CC(=O)[O-].[Na+]');
-- C(C)(=O)O
```

---

## Function Reference

Function descriptions, meaningful argument names and runnable examples are also
available from the SQL catalog after loading the extension:

```sql
SELECT function_name, parameters, parameter_types, description, examples
FROM duckdb_functions()
WHERE list_contains(categories, 'ducksmiles')
ORDER BY function_name, len(parameters);
```

To verify catalog coverage and execute every documented example against a locally
built extension, use a matching DuckDB CLI:

```sh
python3 scripts/check_function_catalog.py --duckdb /path/to/duckdb \
  --extension /path/to/ducksmiles.duckdb_extension
```

For experimental native Rust fingerprints, fragmentation features and streaming
command-line examples, see [CASMI_RESEARCH.md](CASMI_RESEARCH.md). These research
APIs are separate from the SQL functions documented below.

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

**Supported SMILES features:** organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I), aromatic atoms (c, n, o, s, p, b), bracket atoms (`[NH3+]`, `[Fe+2]`, `[13C@@H:7]`), isotope / atom-map / `@` / `@@` metadata retention, branches, ring closures, disconnected fragments (`.`), bond types (single `-`, double `=`, triple `#`), implicit hydrogens (valence-based), and conservative six-membered Kekule aromaticity perception for common benzene/pyridine-like rings. This is intentionally smaller than RDKit's full sanitization, valence model, stereochemical canonicalization, and broad aromaticity stack.

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

#### `mol_num_fragments(smiles) -> INTEGER`

Returns the number of disconnected molecular fragments (connected components).
This is useful for detecting salts and mixtures before standardization. Returns
`NULL` for invalid SMILES.

```sql
SELECT mol_num_fragments('CCO');          -- 1
SELECT mol_num_fragments('[Na+].[Cl-]');  -- 2
SELECT mol_num_fragments('C.C.O');        -- 3
```

#### `mol_num_bonds(smiles) -> INTEGER`

Returns the total number of bonds in the molecule. Double bonds count as 1 bond, triple bonds count as 1 bond (bond count, not bond order sum). Ring closures are included. Returns `NULL` for invalid SMILES.

```sql
SELECT mol_num_bonds('CCO');          -- 2  (C-C, C-O)
SELECT mol_num_bonds('c1ccccc1');     -- 6  (benzene ring)
SELECT mol_num_bonds('C=C');          -- 1  (one double bond)
SELECT mol_num_bonds('C#C');          -- 1  (one triple bond)
```

#### Fast molecule-profile descriptors (20)

These descriptors reuse the parsed molecular graph and are intended for wide,
vectorized profiling queries. Count functions return `INTEGER`; ratios and mass
return `DOUBLE`. Invalid SMILES become `NULL` at the DuckDB boundary.

| Group | SQL functions | Meaning |
|---|---|---|
| charge / H | `mol_formal_charge`, `mol_num_explicit_h`, `mol_num_implicit_h`, `mol_num_total_h` | Net formal charge and hydrogen counts; bracket H and sanitized pyrrolic N/P H are explicit |
| bond profile | `mol_num_single_bonds`, `mol_num_double_bonds`, `mol_num_triple_bonds`, `mol_num_aromatic_bonds` | Number of bonds by perceived order |
| ring / aromaticity | `mol_num_ring_atoms`, `mol_num_ring_bonds`, `mol_largest_ring_size`, `mol_num_aromatic_atoms` | Ring membership, largest SSSR ring, aromatic atoms |
| element profile | `mol_num_carbons`, `mol_num_nitrogens`, `mol_num_oxygens`, `mol_num_halogens` | C/N/O and F+Cl+Br+I counts |
| ratios / size | `mol_heteroatom_fraction`, `mol_aromatic_fraction`, `mol_heavy_atom_mass`, `mol_mean_degree` | Fractions over heavy atoms, isotope-aware heavy-atom mass, mean heavy-atom graph degree |

```sql
SELECT
  mol_formal_charge('[Cl-]') AS charge,        -- -1
  mol_num_total_h('CCO') AS hydrogens,         -- 6
  mol_num_aromatic_bonds('c1ccccc1') AS ar_bonds, -- 6
  mol_largest_ring_size('c1ccccc1') AS max_ring,  -- 6
  round(mol_heteroatom_fraction('CCO'), 3) AS hetero_fraction; -- 0.333
```

The 21 structural-profile functions above have a dedicated
[RDKit parity contract and reproducible test report](PROFILE_PARITY.md).
The checked-in corpus covers 4,510 inputs (including alternative SMILES and
element/isotope cases), with exact integer comparisons and absolute tolerance
`1e-9` for the four floating-point descriptors. This is a bounded compatibility
result, not a claim of full RDKit parser/sanitizer equivalence.

`mol_heavy_atom_mass` matches RDKit `HeavyAtomMolWt`: unlabelled atoms use
RDKit 2025.09.6 average atomic weights, isotope-labelled atoms use isotope
masses, and all hydrogen isotopes are excluded. It is not monoisotopic mass.
Legacy `mol_weight` and `mol_exact_mass` retain their existing tables.

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

#### `mol_mr(smiles) -> DOUBLE`

Returns the **Wildman–Crippen molar refractivity (MolMR)**, computed from the *same* atom-contribution table and first-match-wins typing as `logp_crippen` (the MR column of RDKit's `Crippen.txt`). MR captures molecular volume + polarizability and is one of the eight descriptors feeding drug-likeness models. Matches RDKit's `Crippen.MolMR` for small molecules. Returns `NULL` for invalid SMILES.

```sql
SELECT round(mol_mr('C'), 3);          -- 6.731   (methane)
SELECT round(mol_mr('c1ccccc1'), 3);   -- 26.442  (benzene)
SELECT round(mol_mr('CC(=O)Oc1ccccc1C(=O)O'), 2); -- ~44.6 (aspirin)
```

#### Ring-count descriptors `-> INTEGER`

RDKit's `rdMolDescriptors` ring breakdown, all over the SSSR ring set. A ring is **aromatic** when every ring bond is aromatic, **saturated** when every ring bond is single, **aliphatic** when it has ≥1 non-aromatic bond; **heterocycle**/**carbocycle** splits on whether any ring atom is non-carbon. (Aromaticity follows the parser's perception — see `num_aromatic_rings`.)

| Function | Counts rings that are… |
|----------|------------------------|
| `num_aliphatic_rings`        | not fully aromatic (≥1 non-aromatic bond) |
| `num_saturated_rings`        | all single bonds |
| `num_aromatic_heterocycles`  | aromatic **and** contain a heteroatom |
| `num_aromatic_carbocycles`   | aromatic, all-carbon |
| `num_saturated_heterocycles` | saturated **and** contain a heteroatom |
| `num_saturated_carbocycles`  | saturated, all-carbon |
| `num_aliphatic_heterocycles` | aliphatic **and** contain a heteroatom |
| `num_aliphatic_carbocycles`  | aliphatic, all-carbon |

```sql
SELECT num_aromatic_carbocycles('c1ccccc1');    -- 1  (benzene)
SELECT num_aromatic_heterocycles('c1ccncc1');   -- 1  (pyridine)
SELECT num_saturated_carbocycles('C1CCCCC1');   -- 1  (cyclohexane)
SELECT num_saturated_heterocycles('C1CCNCC1');  -- 1  (piperidine)
SELECT num_aliphatic_carbocycles('C1=CCCCC1');  -- 1  (cyclohexene, unsaturated but non-aromatic)
```

#### `qed(smiles) -> DOUBLE`

Returns the **QED drug-likeness score** (Bickerton et al., *Nat. Chem.* 2012) in `[0, 1]` — a weighted geometric mean of eight desirability-mapped properties (MW, ALOGP, HBA, HBD, PSA, ROTB, aromatic rings, structural alerts). Port of RDKit's `Chem.QED.default` (`weights_mean`). Returns `NULL` for invalid SMILES.

```sql
SELECT round(qed('CC(=O)Oc1ccccc1C(=O)O'), 3);   -- 0.550  (aspirin)
SELECT round(qed('CN1C=NC2=C1C(=O)N(C(=O)N2C)C'), 3); -- 0.529 (caffeine)
SELECT round(qed('N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O'), 3); -- 0.253 (QED paper ref)
```

Reproduces RDKit's reference values closely (paper molecules 0.253 / 0.234 land at 0.253 / 0.234). Two fidelity caveats, both minor: ROTB uses RDKit's *Strict* definition via a Rust bond pass (the engine doesn't parse bond-OR `-,:`), and 114 of 116 structural alerts are evaluated (the two skipped are multi-component `.` SMARTS). As RDKit's own `QED.py` notes, exact parity is bounded by the underlying LogP calculator anyway.

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

#### `murcko_scaffold(smiles) -> VARCHAR`

Returns the Bemis-Murcko scaffold: ring systems plus linker atoms that connect ring systems, with terminal side chains pruned. Acyclic molecules return an empty string. Returns `NULL` for invalid SMILES.

```sql
SELECT murcko_scaffold('Cc1ccccc1');                  -- c1ccccc1
SELECT murcko_scaffold('CC(=O)Oc1ccccc1C(=O)O');      -- c1ccccc1
SELECT murcko_scaffold('c1ccccc1CCc2ccccc2');         -- keeps both rings and the linker
```

#### `generic_scaffold(smiles) -> VARCHAR`

Returns the generic Bemis-Murcko scaffold by converting scaffold atoms to carbon and scaffold bonds to single bonds. This mirrors the common RDKit workflow of extracting a Murcko scaffold and making it generic for scaffold-class grouping.

```sql
SELECT generic_scaffold('c1ccccn1');  -- C1CCCCC1
```

#### `ring_systems_json(smiles) -> VARCHAR`

Returns ring systems as JSON. Each system includes 1-based atom indices, bond indices, ring count, aromatic flag, and the individual rings assigned to the system. Acyclic molecules return `[]`.

```sql
SELECT ring_systems_json('c1ccc2ccccc2c1');
```

#### `mol_has_substructure(smiles, smarts) -> BOOLEAN`

Returns `true` when the molecule contains at least one match for the SMARTS query. Returns `NULL` for invalid SMILES or unsupported SMARTS.

#### `mol_substructure_count(smiles, smarts) -> INTEGER`

Returns the number of unique SMARTS matches, deduplicated by the set of matched molecule atoms, mirroring RDKit's `SubstructMatch(..., uniquify=true)` behavior.

#### `mol_substructure_matches_json(smiles, smarts) -> VARCHAR`

Returns every unique SMARTS match as JSON. Each match includes `atom_indices` (1-based target atom indices in SMARTS pattern atom order) and an `atoms` mapping array with `query_atom`, `target_atom`, and target atom `symbol`.

```sql
SELECT mol_substructure_matches_json('CC(=O)O', 'C=O');
-- [{"match":1,"atom_indices":[2,3],"atoms":[{"query_atom":1,"target_atom":2,"symbol":"C"},{"query_atom":2,"target_atom":3,"symbol":"O"}]}]
```

SMARTS atom predicates include element / aromaticity, `#`, `H`, `X`, `D`, `v`, charge, `R/R0`, `r`, `x`, isotope, atom map (`:n`), bracket chirality (`@`/`@@`), ring bonds, recursive SMARTS, and boolean operators (`!`, `;`, `,`, implicit AND).

#### `mol_hash(smiles, method) -> VARCHAR`

Returns stable SQL grouping keys for molecular deduplication and graph bucketing. Supported methods are exposed by `mol_hash_methods()` and include `canonical_smiles`, `formula`, `net_charge`, `degree_vector`, `atom_bond_counts`, `element_graph`, `bond_order_graph`, `anonymous_graph`, `murcko_scaffold`, and `generic_scaffold`. Unknown methods and invalid SMILES return `NULL`.

```sql
SELECT mol_hash('CC(=O)[O-].[Na+]', 'formula');          -- C2H3NaO2
SELECT mol_hash('CC(=O)[O-].[Na+]', 'element_graph');    -- C(C)(O)O.[Na]
SELECT mol_hash('CC(=O)[O-].[Na+]', 'anonymous_graph');  -- C.C(C)(C)C
SELECT mol_hash_methods();
```

#### `largest_fragment(smiles) -> VARCHAR`
#### `strip_salts(smiles) -> VARCHAR`
#### `neutralize_charges(smiles) -> VARCHAR`
#### `normalize_smiles(smiles) -> VARCHAR`
#### `fragment_parent(smiles) -> VARCHAR`

SQL-native standardization helpers for dirty supplier or registry data. `largest_fragment` returns the largest connected component. `strip_salts` drops inorganic counterions when an organic component exists. `neutralize_charges` protonates common anions and removes removable cationic hydrogens while leaving metal ions charged. `normalize_smiles` neutralizes and canonicalizes with the local aromaticity model. `fragment_parent` strips salts, picks the main organic fragment, neutralizes it, and returns a canonical parent.

```sql
SELECT largest_fragment('CC(=O)[O-].[Na+]');     -- C(C)(=O)[O-]
SELECT strip_salts('CCO.CN.[Cl-]');              -- C(C)O.CN
SELECT neutralize_charges('CC(=O)[O-].[Na+]');   -- C(C)(=O)O.[Na+]
SELECT fragment_parent('CC(=O)[O-].[Na+]');      -- C(C)(=O)O
```

#### `mcs_smarts(smiles_a, smiles_b) -> VARCHAR`
#### `mcs_json(smiles_a, smiles_b) -> VARCHAR`
#### `scaffold_network_json(smiles) -> VARCHAR`

`mcs_smarts` returns an exact connected maximum common subgraph for two parsed molecules under the local atom/bond compatibility rules. `mcs_json` also returns 1-based atom mappings. `scaffold_network_json` emits molecule, Murcko scaffold, generic scaffold, and ring-system nodes plus scaffold edges as JSON.

```sql
SELECT mcs_smarts('CC(=O)O', 'CC(=O)N');  -- C(C)=O
SELECT mcs_json('CC(=O)O', 'CC(=O)N');
SELECT scaffold_network_json('Cc1ccccc1');
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

#### Other similarity metrics — `dice_bit`, `cosine_bit`, `tversky_bit`, … `-> DOUBLE`

The same raw-BLOB, no-cast machinery powers the rest of RDKit's
`DataStructs` similarity family. All take two equal-length fingerprint BLOBs
and share `tanimoto_bit`'s contract (length mismatch → `InvalidInputException`,
both-empty → `0.0`). Formulas are faithful ports of
`RDKit/Code/DataStructs/BitOps.cpp` (`c = popcount(a & b)`, `|a|`, `|b|` are the
on-bit counts):

| Function | Formula | Notes |
|----------|---------|-------|
| `dice_bit(a, b)`           | `2c / (\|a\| + \|b\|)`              | default metric for some RDKit FPs |
| `cosine_bit(a, b)`         | `c / sqrt(\|a\|·\|b\|)`             | Ochiai coefficient |
| `kulczynski_bit(a, b)`     | `c·(\|a\|+\|b\|) / (2·\|a\|·\|b\|)` | |
| `sokal_bit(a, b)`          | `c / (2\|a\| + 2\|b\| − 3c)`        | |
| `mcconnaughey_bit(a, b)`   | `(c·(\|a\|+\|b\|) − \|a\|·\|b\|) / (\|a\|·\|b\|)` | range `[−1, 1]` |
| `asymmetric_bit(a, b)`     | `c / min(\|a\|, \|b\|)`             | |
| `braun_blanquet_bit(a, b)` | `c / max(\|a\|, \|b\|)`             | |
| `russel_bit(a, b)`         | `c / N`                             | `N` = total bit width |
| `tversky_bit(a, b, α, β)`  | `c / (α\|a\| + β\|b\| + (1−α−β)c)`  | asymmetric; see below |

`tversky_bit` takes two extra weights `α, β ∈ [0, 1]` (out-of-range → error).
`α=β=1` reduces to Tanimoto; `α=β=0.5` to Dice; `α=0, β=1` gives a
substructure-style score (how much of `b` is contained in `a`) — handy for
"find molecules that contain this query fragment" searches.

```sql
-- Dice and Tversky vs a reference, same shape as tanimoto_bit
WITH ref AS (SELECT morgan_fp_bits('CC(=O)Oc1ccccc1C(=O)O') AS fp)
SELECT name,
       round(dice_bit(morgan_fp_bits(smiles), (SELECT fp FROM ref)), 4)            AS dice,
       round(tversky_bit(morgan_fp_bits(smiles), (SELECT fp FROM ref), 0.0, 1.0), 4) AS contains
FROM (VALUES
  ('aspirin',        'CC(=O)Oc1ccccc1C(=O)O'),
  ('salicylic acid', 'OC(=O)c1ccccc1O')
) AS t(name, smiles)
ORDER BY dice DESC;

-- α=β=1 is exactly Tanimoto:
WITH q AS (SELECT morgan_fp_bits('CCO') AS a, morgan_fp_bits('CCN') AS b)
SELECT tversky_bit(a, b, 1.0, 1.0) = tanimoto_bit(a, b) AS match FROM q;
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

### MOL/SDF Functions

MOL blocks (V2000/V3000) are the standard file format for storing 2D/3D molecular structures with explicit atom coordinates. SDF (Structure Data Format) concatenates multiple MOL blocks with associated properties. These functions parse MOL/SDF text directly from VARCHAR columns. V3000 support covers `COUNTS`, `ATOM`, and `BOND` CTAB blocks, including V3000 line continuation records, so the same MOL block functions work across V2000 and V3000 inputs. Structural JSON functions expose the complete parsed atom table, bond table, and SDF property block from SQL.

#### `mol_block_name(mol) -> VARCHAR`

Extracts the molecule name from the first line of the MOL block header.

#### `mol_block_property(mol, key) -> VARCHAR`

Extracts an SDF data field from a single MOL/SDF record. Property parsing follows the SDF data item structure: a `> <FIELD_NAME>` header, zero or more value lines, and a blank-line terminator. Multi-line values are preserved with newline characters; empty values are returned as an empty string; duplicate property names are preserved internally and this function returns the first matching value. Returns `NULL` if the molecule block is invalid or the property is absent.

```sql
SELECT mol_block_property(sdf_record, 'PUBCHEM_COMPOUND_CID') AS cid
FROM records;
```

#### `mol_block_properties_json(mol) -> VARCHAR`

Returns all SDF data fields from a single MOL/SDF record as an ordered JSON array of `{name, value}` entries. The array form intentionally preserves duplicate property names and original property order.

```sql
SELECT mol_block_properties_json(sdf_record) AS properties_json
FROM records;
-- [{"name":"ID","value":"123"},{"name":"NOTE","value":"line one\nline two"}]
```

#### `mol_block_atoms_json(mol) -> VARCHAR`

Returns the parsed atom block as an ordered JSON array. Each atom includes 1-based `index`, `symbol`, and numeric `x`, `y`, `z` coordinates. Works for both V2000 atom lines and V3000 `ATOM` records.

#### `mol_block_bonds_json(mol) -> VARCHAR`

Returns the parsed bond block as an ordered JSON array. Each bond includes 1-based `index`, 1-based `atom1`/`atom2` endpoints, and the numeric MOL bond type.

#### `mol_block_json(mol) -> VARCHAR`

Returns a complete parsed MOL/SDF record as a JSON object with `name`, `formula`, `weight`, `num_atoms`, `num_bonds`, `has_3d`, `atoms`, `bonds`, and `properties`.

```sql
SELECT mol_block_atoms_json(sdf_record) AS atoms,
       mol_block_bonds_json(sdf_record) AS bonds,
       mol_block_json(sdf_record) AS mol_json
FROM records;
```

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

The same functions accept V3000 CTAB records:

```sql
SELECT mol_block_formula(v3000_mol) AS formula,
       mol_block_num_atoms(v3000_mol) AS atoms,
       mol_block_num_bonds(v3000_mol) AS bonds,
       mol_block_has_3d(v3000_mol) AS has_3d
FROM molecules;
```

#### `sdf_count(sdf) -> INTEGER`

Counts the number of molecules in an SDF file by counting `$$$$` delimiters. Useful for quickly checking dataset size without fully parsing each molecule.

```sql
SELECT sdf_count(sdf_text) AS num_molecules FROM sdf_files;
```

#### `sdf_property(sdf, record_index, key) -> VARCHAR`

Extracts a property from a specific SDF record. `record_index` is 1-based, matching SQL row-numbering conventions. Returns the first matching value in that record, or `NULL` if the record or key is absent.

```sql
SELECT sdf_property(sdf_text, 2, 'PUBCHEM_IUPAC_NAME') AS second_name
FROM sdf_files;
```

#### `sdf_properties_json(sdf) -> VARCHAR`

Returns every parsed SDF record's property block as JSON. Each record includes its 1-based record number, MOL header name, and an ordered `properties` array. Multi-line values, empty values, duplicate field names, and field order are preserved.

```sql
SELECT sdf_properties_json(sdf_text) AS all_sdf_properties
FROM sdf_files;
-- [{"record":1,"name":"aspirin","properties":[{"name":"ID","value":"2244"}]}]
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

<!-- analytics-full-reference:start -->

## CPU analytics: complete SQL API reference

The analytics build adds numerical methods to the molecular SQL workflow. The
calculation core runs on the CPU, without a GUI, GPU runtime, Python or external
chemistry library. This section lists **every currently registered analytics SQL
function** (254 verified functions: the original 24 plus the method registry below), with its input contract and an
executable call. Existing community packages can predate this implementation;
use a build containing these registrations and check the loaded catalog:

```sql
SELECT function_name, description, examples
FROM duckdb_functions()
WHERE starts_with(function_name, 'chem_')
ORDER BY function_name;
```

Most functions use
`chem_function(data DOUBLE[][], parameters DOUBLE[]) -> VARCHAR` and return JSON
with `"version":1`. Each inner list is one observation or the explicitly documented
count table/subgroup. `chem_analyze(operation, data, parameters)` dispatches the
same numeric operations. `chem_scaffold_folds` instead takes a `VARCHAR[]` of
SMILES and is the sole molecular grouping wrapper in this API.

Parameters are positional; optional trailing parameters can be omitted. `[]`
means no numeric parameters; `[]::DOUBLE[]` makes the SQL type explicit. Inputs
must be finite rectangular lists without NULL elements. Unsupported shapes,
domains, numerical singularity/overflow or work limits return SQL `NULL`. Some
successful JSON fields are `null` when a documented quantity is undefined.
Standardize heterogeneous descriptors deliberately, encode categorical variables
explicitly, and handle missing observations before calling these numeric methods.

JSON indices are zero based; DuckDB list subscripts are one based. Preserve a
stable row order when materializing matrices, for example
`list([mw,logp,target] ORDER BY id)`. Keep feature order, units, preprocessing,
model JSON, seeds and input revisions together. Fit preprocessing inside training
folds; molecular scaffold grouping addresses one source of leakage and does not
replace assay, duplicate-measurement, temporal or target grouping.

### Original exploration, prediction, DOE and statistics

| SQL function | Input rows | Ordered parameters | Result / exact scope |
| --- | --- | --- | --- |
| `chem_analyze` | `operation VARCHAR`, numeric matrix | Same parameter vector as the selected operation | Generic dispatcher; scaffold grouping uses its separate wrapper |
| `chem_pca` | `[feature_0, ...]` | `[components=min(p,n-1), standardize=0]` | `mean`, `scale`, `components`, `scores`, `explained_variance`, `explained_variance_ratio` |
| `chem_kmeans` | `[feature_0, ...]` | `[clusters=min(2,n), max_iterations=300, seed=0]` | `labels`, `centers`, `inertia`, `iterations`, `converged` |
| `chem_ridge` | `[features..., target]` | `[lambda]`, nonnegative | `intercept`, `coefficients`, `stable_predict_parameters`, `feature_means`, `feature_scales`, `metrics` |
| `chem_pls` | `[features..., target]` | `[components]`, integer `1..p` | Same model fields, plus `components_used` |
| `chem_predict` | `[features...]` | Saved `stable_predict_parameters`, or `[intercept, coefficient_0, ...]` | `predictions` |
| `chem_cross_validate` | `[features..., target, fold_id]` | `[method, hyperparameter]`; method `0`=Ridge, `1`=PLS1 | `out_of_fold_predictions`, `rmse`, `r2`, `fold_metrics`, `preprocessing` |
| `chem_group_folds` | `[group_id]` | `[n_folds, seed]` | `folds`, `fold_sizes`, `n_groups` |
| `chem_scaffold_folds` | List of SMILES strings | `[n_folds, seed]` | Same fold result |
| `chem_factorial_design` | One `[low, high]` row per factor | `[center_points=0]` | `rows`, `coded_rows`, `run_count`, `factors` |
| `chem_central_composite_design` | One `[low, high]` row per factor | `[alpha=2^(p/4), center_points=1]` | Same fields, plus `alpha` |
| `chem_d_optimal_design` | Candidate feature rows | `[n_runs, max_iterations=50, seed=0]` | `indices`, `rows`, `log_determinant`, `means`, `scales`, `exchanges`, `converged` |
| `chem_response_surface` | `[features..., target]` | `[lambda=0]` | `terms`, `means`, `scales`, `coefficients`, `model_parameters`, `fitted`, `rmse`, `r2` |
| `chem_response_surface_predict` | Raw feature rows | Saved `model_parameters` | `predictions` |
| `chem_gp_predict` | Training prefix `[features..., observed_target]`, then candidates `[features..., dummy_target]` | `[n_train, length_scale, signal_variance, noise_variance, minimize, xi]` | `predictions`, `variances`, `expected_improvement`, `ranked_candidate_indices`, `feature_means`, `feature_scales`, `target_mean` |
| `chem_acquisition` | `[predictive_mean, predictive_sigma]` | `[best, minimize, xi]` | `expected_improvement` |
| `chem_diversity_select` | Existing feature prefix, then candidate feature rows | `[n_select, n_existing]` | `selected_indices` |
| `chem_summarize` | Numeric feature rows | `[]` | Column `means`, `sample_variances`, `standard_deviations`, `minima`, `maxima`; sample variance/SD are `null` for one row |
| `chem_welch_t_test` | `[value, group]`, group `0` or `1` | `[]` | Two-sided Welch test: group sizes/means, group-0 minus group-1 `mean_difference`, `t_statistic`, Satterthwaite `degrees_of_freedom`, `p_value` |
| `chem_one_way_anova` | `[value, group_id]` | `[]` | One-way between-groups F test: sums of squares, degrees of freedom, `f_statistic`, `p_value` |
| `chem_adjust_pvalues` | `[p_value]` in `[0,1]` | `[method]`: `0`=Holm, `1`=Benjamini-Hochberg | `adjusted_p_values` in original input order |
| `chem_process_capability` | `[value]` | `[lower_specification_limit, upper_specification_limit, within_process_sigma (optional)]` | Mean, overall sample SD, `pp`/`ppk`; `cp`/`cpk` only with supplied within-process sigma; descriptive mean ± 3 SD reference limits |
| `chem_kaplan_meier` | `[time, event]`; time nonnegative, event `1`=event, `0`=right-censored | `[]` | `times`, `at_risk`, `events`, `censored`, `survival`, `greenwood_variance`, `median_survival` |
| `chem_random_intercept` | `[value, group_id]` | `[]` | Balanced one-way random-intercept ANOVA variance estimator: `within_group_variance`, `between_group_variance`, `icc`, grand mean and group counts |

Every original function has a standalone call below. These examples use small
synthetic values; they are API demonstrations rather than scientific findings.

```sql
SELECT chem_analyze('acquisition', [[2.0,1.0]], [2.0,0.0,0.0]);
SELECT chem_pca([[1.0,2.0],[2.0,4.0],[3.0,6.0]], [1.0,0.0]);
SELECT chem_kmeans([[0.0],[1.0],[9.0],[10.0]], [2.0,100.0,42.0]);
SELECT chem_ridge([[0.0,1.0],[1.0,3.0],[2.0,5.0]], [1.0]);
SELECT chem_pls([[0.0,1.0],[1.0,3.0],[2.0,5.0]], [1.0]);
SELECT chem_predict([[1.0],[2.0]], [1.0,2.0]);
SELECT chem_cross_validate([[0.0,1.0,0.0],[1.0,3.0,0.0],[2.0,5.0,1.0],[3.0,7.0,1.0]], [0.0,1.0]);
SELECT chem_group_folds([[1.0],[1.0],[2.0],[3.0]], [2.0,42.0]);
SELECT chem_scaffold_folds(['c1ccccc1','Cc1ccccc1','C1CCCCC1','CCO'], [2.0,42.0]);
SELECT chem_factorial_design([[10.0,20.0],[1.0,3.0]], [1.0]);
SELECT chem_central_composite_design([[10.0,20.0],[1.0,3.0]], [1.0,1.0]);
SELECT chem_d_optimal_design([[-1.0],[-0.5],[0.5],[1.0]], [2.0,20.0,42.0]);
SELECT chem_response_surface([[-1.0,2.0],[0.0,1.0],[1.0,4.0]], [0.0]);
SELECT chem_response_surface_predict([[0.0],[1.0]], [1.0,0.0,1.0,1.0,1.0,2.0]);
SELECT chem_gp_predict([[0.0,2.0],[1.0,0.0]], [1.0,1.0,1.0,0.01,0.0,0.0]);
SELECT chem_acquisition([[2.0,1.0]], [2.0,0.0,0.0]);
SELECT chem_diversity_select([[0.0],[1.0],[5.0],[10.0]], [2.0,1.0]);
SELECT chem_summarize([[1.0],[2.0],[3.0]], []::DOUBLE[]);
SELECT chem_welch_t_test([[1.0,0.0],[2.0,0.0],[3.0,0.0],[2.0,1.0],[3.0,1.0],[4.0,1.0]], []::DOUBLE[]);
SELECT chem_one_way_anova([[1.0,0.0],[2.0,0.0],[3.0,0.0],[2.0,1.0],[3.0,1.0],[4.0,1.0]], []::DOUBLE[]);
SELECT chem_adjust_pvalues([[0.01],[0.04],[0.03]], [1.0]);
SELECT chem_process_capability([[1.0],[2.0],[3.0]], [-1.0,5.0]);
SELECT chem_kaplan_meier([[1.0,1.0],[2.0,0.0],[3.0,1.0]], []::DOUBLE[]);
SELECT chem_random_intercept([[1.0,0.0],[3.0,0.0],[5.0,1.0],[7.0,1.0]], []::DOUBLE[]);
```

### Additional statistical and modeling methods

The following table is derived from the extension's actual public method
registry. Each row includes the data shape, ordered parameter contract and
method limitations, followed by the same runnable call exposed in the SQL
catalog. A fitted-method family can still be partial even when its listed
limited kernel has passed tests.

<!-- analytics-registry-table:start -->

| SQL function / dispatcher operation | Input and ordered parameters / limited scope | Standalone SQL example |
| --- | --- | --- |
| `chem_one_sample_t` / `one_sample_t` | One-sample t inference: scalar rows; [null_mean,confidence=.95]. | `SELECT chem_one_sample_t([[1.0],[2.0],[4.0]], [0.0,0.95]);` |
| `chem_paired_t` / `paired_t` | Paired t inference: [a,b] rows; [null_difference=0,confidence=.95]. | `SELECT chem_paired_t([[1.0,0.0],[3.0,1.0],[6.0,3.0]], []);` |
| `chem_pooled_t` / `pooled_t` | Equal-variance two-sample t: [value,group0/1]; [null_difference=0,confidence=.95]. | `SELECT chem_pooled_t([[1.0,0.0],[2.0,0.0],[3.0,1.0],[5.0,1.0]], []);` |
| `chem_variance_test` / `variance_test` | Normal-model variance test and interval: scalar rows; [null_sd,confidence=.95]. | `SELECT chem_variance_test([[1.0],[2.0],[4.0]], [1.0]);` |
| `chem_variance_ratio_test` / `variance_ratio_test` | Normal-model F variance-ratio test: [value,group0/1]; [null_ratio=1,confidence=.95]. | `SELECT chem_variance_ratio_test([[1.0,0.0],[2.0,0.0],[3.0,1.0],[5.0,1.0]], []);` |
| `chem_proportion_test` / `proportion_test` | Probability-ordered exact binomial test and Wilson interval: binary rows; [null_p,confidence=.95]. | `SELECT chem_proportion_test([[1.0],[0.0],[1.0],[1.0]], [0.5]);` |
| `chem_contingency_test` / `contingency_test` | Pearson chi-square independence test, expected counts and Cramer V: observed count matrix; []. | `SELECT chem_contingency_test([[10.0,20.0],[20.0,10.0]], []);` |
| `chem_fisher_exact` / `fisher_exact` | Two-sided fixed-marginal Fisher exact test: 2x2 counts; []. | `SELECT chem_fisher_exact([[1.0,9.0],[11.0,3.0]], []);` |
| `chem_mann_whitney` / `mann_whitney` | Mann-Whitney asymptotic test with ties and continuity correction: [value,group0/1]; []. | `SELECT chem_mann_whitney([[1.0,0.0],[2.0,0.0],[3.0,1.0],[4.0,1.0]], []);` |
| `chem_wilcoxon_signed_rank` / `wilcoxon_signed_rank` | Paired Wilcoxon signed-rank normal approximation with zero/tie correction: [a,b]; []. | `SELECT chem_wilcoxon_signed_rank([[2.0,1.0],[4.0,2.0],[6.0,3.0]], []);` |
| `chem_kruskal_wallis` / `kruskal_wallis` | Tie-corrected Kruskal-Wallis test: [value,group_id]; []. | `SELECT chem_kruskal_wallis([[1.0,0.0],[2.0,0.0],[3.0,1.0],[5.0,1.0]], []);` |
| `chem_levene_test` / `levene_test` | Levene/Brown-Forsythe variance test: [value,group_id]; [center1median/0mean]. | `SELECT chem_levene_test([[1.0,0.0],[2.0,0.0],[4.0,0.0],[2.0,1.0],[4.0,1.0],[7.0,1.0]], [1.0]);` |
| `chem_equivalence_t` / `equivalence_t` | One-sample normal-model TOST equivalence: scalar rows; [lower_margin,upper_margin,alpha=.05]. | `SELECT chem_equivalence_t([[0.0],[0.1],[-0.1],[0.05]], [-1.0,1.0,0.05]);` |
| `chem_linear_model` / `linear_model` | OLS with coefficient, residual and overall-model inference: [features...,target]; [confidence=.95]. | `SELECT chem_linear_model([[0.0,1.0],[1.0,2.0],[2.0,4.0],[3.0,4.5]], [0.95]);` |
| `chem_logistic` / `logistic` | Binomial logit GLM, standardized penalty: [features...,binary_target]; [lambda=0,max_iterations=100,tolerance=1e-8]. | `SELECT chem_logistic([[0.0,0.0],[1.0,1.0],[2.0,0.0],[3.0,1.0]], [1.0]);` |
| `chem_poisson` / `poisson` | Poisson log-link GLM: [features...,count_target]; [lambda=0,max_iterations=100,tolerance=1e-8]. | `SELECT chem_poisson([[0.0,1.0],[1.0,1.0],[2.0,3.0],[3.0,4.0]], [0.1]);` |
| `chem_glm_predict` / `glm_predict` | GLM prediction from saved stable_predict_parameters or [link0identity/1logit/2log,intercept,coefficients...]. | `SELECT chem_glm_predict([[0.0],[1.0]], [1.0,0.0,1.0]);` |
| `chem_elastic_net` / `elastic_net` | Standardized coordinate-descent elastic net: [features...,target]; [lambda,alpha,max_iterations=10000,tolerance=1e-8]. | `SELECT chem_elastic_net([[0.0,1.0],[1.0,3.0],[2.0,5.0]], [0.1,0.5]);` |
| `chem_polynomial_features` / `polynomial_features` | Total-degree polynomial terms with exponent metadata: features; [degree1..3,include_bias0/1]. | `SELECT chem_polynomial_features([[1.0,2.0]], [2.0,1.0]);` |
| `chem_ir_chart` / `ir_chart` | Individuals/moving-range chart; ordered scalar rows; optional fixed [center,sigma]. | `SELECT chem_ir_chart([[1.0],[2.0],[4.0]], []);` |
| `chem_xbar_r_chart` / `xbar_r_chart` | Xbar-R chart using d2/d3 constants, equal subgroups of size2..10; []. | `SELECT chem_xbar_r_chart([[1.0,3.0],[2.0,4.0],[3.0,5.0]], []);` |
| `chem_xbar_s_chart` / `xbar_s_chart` | Xbar-S chart with c4 correction; equal subgroups of size2..1000; []. | `SELECT chem_xbar_s_chart([[1.0,3.0],[2.0,4.0],[3.0,5.0]], []);` |
| `chem_p_chart` / `p_chart` | Binomial fraction-nonconforming control limits: [defectives,sample_size]; []. | `SELECT chem_p_chart([[2.0,100.0],[4.0,100.0]], []);` |
| `chem_np_chart` / `np_chart` | Binomial count-nonconforming control limits, constant sample size: [defectives,n]; []. | `SELECT chem_np_chart([[2.0,100.0],[4.0,100.0]], []);` |
| `chem_c_chart` / `c_chart` | Poisson count control chart: count rows; []. | `SELECT chem_c_chart([[4.0],[6.0]], []);` |
| `chem_u_chart` / `u_chart` | Poisson rate control chart: [count,exposure]; []. | `SELECT chem_u_chart([[2.0,2.0],[6.0,2.0]], []);` |
| `chem_ewma_chart` / `ewma_chart` | EWMA with transient limits: ordered scalar rows; [lambda,L,known_center,known_sigma]. | `SELECT chem_ewma_chart([[1.0],[2.0],[3.0]], [0.2,3.0,0.0,1.0]);` |
| `chem_cusum_chart` / `cusum_chart` | Standardized two-sided CUSUM: ordered scalar rows; [k,h,known_center,known_sigma]. | `SELECT chem_cusum_chart([[1.0],[2.0],[3.0]], [0.5,5.0,0.0,1.0]);` |
| `chem_moving_average_chart` / `moving_average_chart` | Unweighted moving-average control chart, partial startup windows: scalar rows; [window,L,known_center,known_sigma]. | `SELECT chem_moving_average_chart([[1.0],[2.0],[3.0]], [2.0,3.0,0.0,1.0]);` |
| `chem_western_electric` / `western_electric` | Four Western Electric rules on fixed reference limits: scalar rows; [center,sigma]. | `SELECT chem_western_electric([[2.1],[2.2],[0.0]], [0.0,1.0]);` |
| `chem_gage_rr` / `gage_rr` | Balanced crossed random-effects ANOVA Gage R&R without pooling: [value,part,operator]; []. | `SELECT chem_gage_rr([[1.0,0.0,0.0],[3.0,0.0,0.0],[2.0,0.0,1.0],[4.0,0.0,1.0],[5.0,1.0,0.0],[7.0,1.0,0.0],[6.0,1.0,1.0],[8.0,1.0,1.0]], []);` |
| `chem_type1_gage` / `type1_gage` | Type1 Gage Cg/Cgk with20% tolerance fraction: scalar measurements; [reference,tolerance_width]. | `SELECT chem_type1_gage([[9.0],[10.0],[11.0]], [10.0,6.0]);` |
| `chem_agreement_kappa` / `agreement_kappa` | Unweighted two-rater Cohen kappa from a square count table; []. | `SELECT chem_agreement_kappa([[20.0,5.0],[5.0,20.0]], []);` |
| `chem_acceptance_oc` / `acceptance_oc` | Single binomial acceptance-sampling OC probabilities: defect probability rows; [n,c]. | `SELECT chem_acceptance_oc([[0.01],[0.1]], [10.0,0.0]);` |
| `chem_levey_jennings` / `levey_jennings` | Levey-Jennings fixed three-sigma limits: scalar rows; [reference_mean,reference_sigma]. | `SELECT chem_levey_jennings([[0.0],[3.1]], [0.0,1.0]);` |
| `chem_laney_p_chart` / `laney_p_chart` | Laney P-prime overdispersion limits estimated from standardized moving ranges: [defectives,n]; []. | `SELECT chem_laney_p_chart([[1.0,10.0],[9.0,10.0]], []);` |
| `chem_laney_u_chart` / `laney_u_chart` | Laney U-prime overdispersion limits estimated from standardized moving ranges: [count,exposure]; []. | `SELECT chem_laney_u_chart([[1.0,1.0],[9.0,1.0]], []);` |
| `chem_bootstrap_mean` / `bootstrap_mean` | IID percentile bootstrap mean: >=2 scalar rows; [replicates2..100000,confidence=.95,seed=0]; <=20000000 resampled values. | `SELECT chem_bootstrap_mean([[1.0],[2.0],[3.0]], [100.0,0.95,42.0]);` |
| `chem_permutation_test` / `permutation_test` | Mean-difference label permutation: [value,group0/1]; [replicates0exact,seed=0,alternative0two/1greater/2less]; bounded exact or Monte Carlo work. | `SELECT chem_permutation_test([[0.0,0.0],[1.0,0.0],[2.0,1.0],[3.0,1.0]], [0.0]);` |
| `chem_sample_rows` / `sample_rows` | Seeded row sampling: finite matrix; [count1..100000,replacement0/1,seed]; <=1000000 output cells. | `SELECT chem_sample_rows([[1.0,2.0],[3.0,4.0],[5.0,6.0]], [2.0,0.0,42.0]);` |
| `chem_simulate_normal` / `simulate_normal` | Box-Muller normal simulation: [mean,sd>=0] rows; [count_per_row,seed]; <=1000000 output samples. | `SELECT chem_simulate_normal([[0.0,1.0],[2.0,0.5]], [3.0,42.0]);` |
| `chem_power_z` / `power_z` | One-sample known-sigma normal power: [effect,sigma>0]; [n1..1000000000,alpha=.05,alternative0two/1greater/2less]. | `SELECT chem_power_z([[0.5,1.0]], [20.0,0.05,0.0]);` |
| `chem_sample_size_mean` / `sample_size_mean` | Normal-approximation mean margin planning: positive known sigma rows; [margin>0,confidence=.95]; sample size <=1000000000000. | `SELECT chem_sample_size_mean([[1.0],[2.0]], [0.5,0.95]);` |
| `chem_sample_size_proportion` / `sample_size_proportion` | Normal-approximation proportion margin planning: p in [0,1] rows; [margin>0,confidence=.95]; sample size <=1000000000000. | `SELECT chem_sample_size_proportion([[0.5],[0.2]], [0.05,0.95]);` |
| `chem_reliability_demonstration` / `reliability_demonstration` | Independent zero-failure binomial reliability demonstration: target reliability in (0,1) rows; [confidence=.95]; n<=1000000000000. | `SELECT chem_reliability_demonstration([[0.95],[0.99]], [0.95]);` |
| `chem_odds_ratio` / `odds_ratio` | Log-Wald odds ratio: positive 2x2 integer count table; [confidence=.95,continuity_correction=0]; correction must be explicit for zeros. | `SELECT chem_odds_ratio([[20.0,10.0],[10.0,20.0]], []);` |
| `chem_relative_risk` / `relative_risk` | Log-Wald relative risk: positive 2x2 integer count table; [confidence=.95,continuity_correction=0]; correction must be explicit for zeros. | `SELECT chem_relative_risk([[20.0,10.0],[10.0,20.0]], []);` |
| `chem_cohen_d` / `cohen_d` | Pooled-SD Cohen d and exact gamma-ratio Hedges g: [value,group0/1], >=2 observations per group; []. | `SELECT chem_cohen_d([[1.0,0.0],[2.0,0.0],[3.0,1.0],[5.0,1.0]], []);` |
| `chem_fit_normal` / `fit_normal` | Unweighted normal MLE: >=2 nonconstant scalar rows; []; <=100000 observations, absolute values<=1e150. | `SELECT chem_fit_normal([[1.0],[2.0],[4.0],[5.0]], []);` |
| `chem_fit_lognormal` / `fit_lognormal` | Zero-location lognormal MLE: >=2 varying strictly positive scalar rows; []; <=100000 observations. | `SELECT chem_fit_lognormal([[1.0],[2.0],[4.0],[8.0]], []);` |
| `chem_fit_exponential` / `fit_exponential` | Zero-location exponential MLE: nonnegative scalar rows with positive total exposure; []; <=100000 observations. | `SELECT chem_fit_exponential([[1.0],[2.0],[4.0]], []);` |
| `chem_fit_weibull` / `fit_weibull` | Zero-location Weibull MLE with bracketed shape solve: varying positive scalar rows; []; <=100000 observations. | `SELECT chem_fit_weibull([[1.0],[2.0],[3.0],[5.0]], []);` |
| `chem_fit_gamma` / `fit_gamma` | Zero-location Gamma MLE with bracketed shape solve: varying positive scalar rows; []; <=100000 observations. | `SELECT chem_fit_gamma([[1.0],[2.0],[3.0],[5.0]], []);` |
| `chem_fit_poisson` / `fit_poisson` | Poisson MLE with full count likelihood: exact nonnegative integer scalar rows<=9007199254740991; []; <=100000 observations. | `SELECT chem_fit_poisson([[0.0],[1.0],[2.0],[4.0]], []);` |
| `chem_normality_test` / `normality_test` | Jarque-Bera skew/kurtosis normality test with asymptotic chi-square p: >=3 nonconstant scalar rows; []; large samples recommended. | `SELECT chem_normality_test([[1.0],[2.0],[2.0],[3.0],[4.0],[4.0],[5.0],[7.0]], []);` |
| `chem_kde_gaussian` / `kde_gaussian` | Gaussian kernel density: scalar samples; [bandwidth>0,evaluation_points...], 1..4096 points and <=20000000 kernel terms. | `SELECT chem_kde_gaussian([[1.0],[2.0],[3.0]], [0.5,1.0,2.0,3.0]);` |
| `chem_box_cox` / `box_cox` | Box-Cox transform and profile likelihood: >=2 varying positive scalar rows; [lambda] or [] to estimate a bracketed finite optimum. | `SELECT chem_box_cox([[1.0],[2.0],[4.0]], [0.0]);` |
| `chem_survival_exponential` / `survival_exponential` | Zero-location exponential lifetime MLE with right censoring: [time>=0,event0/1], positive exposure and >=1 event; []. | `SELECT chem_survival_exponential([[1.0,1.0],[2.0,0.0],[3.0,1.0],[5.0,1.0]], []);` |
| `chem_survival_weibull` / `survival_weibull` | Zero-location Weibull right-censored lifetime MLE: [time>0,event0/1], identifiable positive-event shape; []. | `SELECT chem_survival_weibull([[1.0,1.0],[2.0,0.0],[3.0,1.0],[5.0,1.0]], []);` |
| `chem_logrank_test` / `logrank_test` | Two-group right-censored log-rank chi-square test: [time>=0,event0/1,group0/1]; []; same-time censors remain at risk. | `SELECT chem_logrank_test([[1.0,1.0,0.0],[2.0,0.0,0.0],[3.0,1.0,0.0],[2.0,1.0,1.0],[4.0,0.0,1.0],[5.0,1.0,1.0]], []);` |
| `chem_acf` / `acf` | Biased sample ACF: ordered nonconstant scalar rows; [max_lag0..min(n-1,10000)]; <=50000000 work units. | `SELECT chem_acf([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [3.0]);` |
| `chem_pacf` / `pacf` | Durbin-Levinson PACF from biased autocovariance: ordered nonconstant scalar rows; [max_lag0..min(n-1,2048)]. | `SELECT chem_pacf([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [3.0]);` |
| `chem_difference` / `difference` | Repeated lag differencing: ordered scalar rows; [lag=1,order1..32=1], must leave >=1 observation. | `SELECT chem_difference([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [1.0,1.0]);` |
| `chem_moving_average` / `moving_average` | Trailing simple moving average: ordered scalar rows; [window1..n]; no future observations are used. | `SELECT chem_moving_average([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [3.0]);` |
| `chem_exponential_smoothing` / `exponential_smoothing` | Simple exponential smoothing: ordered scalar rows; [alpha0..1,horizon1..10000], first-observation initialization. | `SELECT chem_exponential_smoothing([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [0.4,2.0]);` |
| `chem_holt` / `holt` | Holt linear/damped trend smoothing: >=2 ordered scalar rows; [alpha0..1,beta0..1,horizon1..10000,phi in (0,1]=1]. | `SELECT chem_holt([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [0.4,0.3,2.0,0.9]);` |
| `chem_ar_model` / `ar_model` | Conditional Gaussian autoregression via QR with intercept: ordered scalar rows; [order1..32,horizon1..10000]; requires identifiable residual variance. | `SELECT chem_ar_model([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [2.0,2.0]);` |
| `chem_periodogram` / `periodogram` | One-sided unwindowed periodogram density after constant detrending: 2..4096 ordered scalar rows; [sample_interval>0=1]. | `SELECT chem_periodogram([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [1.0]);` |
| `chem_seasonal_decompose` / `seasonal_decompose` | Classical additive centered-moving-average seasonal decomposition: ordered scalar rows; [period2..floor(n/2)]; unidentified edge terms are null. | `SELECT chem_seasonal_decompose([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [2.0]);` |
| `chem_ljung_box` / `ljung_box` | Ljung-Box autocorrelation diagnostic with chi-square approximation: ordered nonconstant scalar rows; [max_lag1..min(n-1,10000),model_df=0<lag]. | `SELECT chem_ljung_box([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [3.0,0.0]);` |
| `chem_covariance_matrix` / `covariance_matrix` | Centered covariance matrix: <=256 finite features; [sample1/0population=1], sample mode requires >=2 rows; <=100000000 numerical work units. | `SELECT chem_covariance_matrix([[1.0,2.0],[2.0,1.0],[3.0,6.0]], []);` |
| `chem_correlation_matrix` / `correlation_matrix` | Pearson correlation matrix: >=2 rows, <=256 finite nonconstant features; []; centered products and bounded work. | `SELECT chem_correlation_matrix([[1.0,2.0],[2.0,1.0],[3.0,6.0]], []);` |
| `chem_hierarchical_cluster` / `hierarchical_cluster` | Euclidean agglomeration with deterministic tie breaks: <=500 feature rows and bounded cubic work; [clusters,linkage0single/1complete/2average/3Ward=2]. | `SELECT chem_hierarchical_cluster([[0.0],[1.0],[4.0],[5.0]], [2.0,3.0]);` |
| `chem_mahalanobis` / `mahalanobis` | Squared Mahalanobis distance to training-prefix moments: <=256 features; [n_reference>=2,diagonal_ridge=0]; nonempty query suffix, singular covariance rejected. | `SELECT chem_mahalanobis([[-1.0,-1.0],[-1.0,1.0],[1.0,-1.0],[1.0,1.0],[2.0,0.0]], [4.0,0.0]);` |
| `chem_discriminant_analysis` / `discriminant_analysis` | Prefix-trained LDA/QDA: [features...,integer_label], query labels ignored; [n_train,method0LDA/1QDA=0,diagonal_ridge=0]; suffix posterior probabilities. | `SELECT chem_discriminant_analysis([[-2.0,10.0],[0.0,10.0],[2.0,20.0],[4.0,20.0],[-1.0,99.0],[3.0,99.0]], [4.0,0.0,0.0]);` |
| `chem_gaussian_naive_bayes` / `gaussian_naive_bayes` | Prefix-trained Gaussian naive Bayes with population class variances: [features...,integer_label]; [n_train,absolute_variance_floor>0=1e-9]; query labels ignored. | `SELECT chem_gaussian_naive_bayes([[-2.0,0.0],[0.0,0.0],[2.0,1.0],[4.0,1.0],[-1.0,99.0],[3.0,99.0]], [4.0,1e-09]);` |
| `chem_knn_predict` / `knn_predict` | K nearest neighbors with training-prefix-only population scaling: [features...,target]; [n_train,k1..n_train,task0regression/1classification]; suffix predictions. | `SELECT chem_knn_predict([[0.0,0.0],[1.0,0.0],[3.0,1.0],[4.0,1.0],[2.0,99.0]], [4.0,3.0,1.0]);` |
| `chem_dbscan` / `dbscan` | Euclidean DBSCAN in original feature units: [eps>=0,min_samples>=1 including self]; <=5000 rows and100000000 work units; noise label -1. | `SELECT chem_dbscan([[0.0],[0.1],[0.2],[2.0],[2.1],[2.2],[9.0]], [0.21,3.0]);` |
| `chem_mixed_model` / `mixed_model` | Unbalanced Gaussian random-intercept ML/REML with <=64 fixed features: [features...,y,integer_group]; [method0ML/1REML=1,maxiter1..10000=200,tolerance=1e-8]; bounded variance-ratio search. | `SELECT chem_mixed_model([[1.0,0.0],[2.0,0.0],[1.0,0.0],[2.0,0.0],[4.0,1.0],[5.0,1.0],[4.0,1.0],[5.0,1.0],[7.0,2.0],[8.0,2.0],[7.0,2.0],[8.0,2.0]], [1.0]);` |
| `chem_mixed_predict` / `mixed_predict` | Random-intercept prediction: [fixed_features...,integer_group]; saved mixed_model predict_parameters; known-group BLUP and zero-effect new-group predictions. | `SELECT chem_mixed_predict([[0.0],[1.0],[99.0]], [0.0,0.0,4.5,0.0,2.0,0.0,-1.0,1.0,1.0]);` |
| `chem_cox_ph` / `cox_ph` | Cox proportional hazards with <=32 features and bounded event work: [features...,time>=0,event0/1]; [ties0Breslow/1Efron=1,ridge=0,maxiter1..1000=100,tolerance=1e-8]. | `SELECT chem_cox_ph([[0.0,1.0,1.0,1.0],[1.0,0.0,1.0,1.0],[0.0,0.0,2.0,0.0],[2.0,1.0,3.0,1.0],[1.0,2.0,3.0,0.0],[2.0,0.0,4.0,1.0],[0.0,2.0,4.0,0.0],[1.0,1.0,5.0,1.0],[2.0,2.0,6.0,1.0],[0.0,1.0,6.0,0.0],[2.0,0.0,7.0,1.0],[1.0,2.0,8.0,0.0]], []);` |
| `chem_cox_predict` / `cox_predict` | Centered Cox survival prediction, <=32 features: feature rows; [baseline_survival0..1,coefficients(p),feature_origins(p),mean_offsets(p)]. | `SELECT chem_cox_predict([[0.0],[1.0],[2.0]], [0.8,0.6931471805599453,0.0,0.0]);` |
| `chem_relative_risk_predict` / `relative_risk_predict` | Centered Cox relative-risk prediction, <=32 features: feature rows; [coefficients(p),feature_origins(p),mean_offsets(p)]. | `SELECT chem_relative_risk_predict([[0.0],[1.0],[2.0]], [0.6931471805599453,0.0,0.0]);` |
| `chem_arima` / `arima` | Conditional-sum-of-squares seasonal ARIMA: ordered scalar rows; [p0..3,d0..2,q0..3,horizon0..1000,constant=1,maxiter=1000,tolerance=1e-8,P0..1=0,D0..2=0,Q0..1=0,season0or2..52=0]. | `SELECT chem_arima([[2.0],[4.0],[3.0],[5.0],[7.0],[6.0],[8.0],[9.0]], [0.0,1.0,0.0,2.0]);` |
| `chem_decision_tree` / `decision_tree` | CART squared-loss/Gini tree with <=64 features: [features...,target], training prefix; [n_train,depth1..16,min_leaf>=1,task0regression/1classification]. | `SELECT chem_decision_tree([[0.0,1.0],[1.0,3.0],[2.0,5.0],[3.0,7.0],[1.5,0.0]], [4.0,3.0,1.0,0.0]);` |
| `chem_random_forest` / `random_forest` | Seeded bootstrap CART forest with feature subsampling: <=64 features; [n_train,trees1..200,depth1..16,min_leaf,task0regression/1classification,seed]; bounded work/model size. | `SELECT chem_random_forest([[0.0,1.0],[1.0,3.0],[2.0,5.0],[3.0,7.0],[1.5,0.0]], [4.0,10.0,3.0,1.0,0.0,42.0]);` |
| `chem_gradient_boosting` / `gradient_boosting` | CPU first-order tree boosting for squared-error/binary log loss: <=64 features; [n_train,iterations1..500,rate in (0,1],depth1..8,min_leaf,task0squared/1binary,seed]. | `SELECT chem_gradient_boosting([[0.0,1.0],[1.0,3.0],[2.0,5.0],[3.0,7.0],[1.5,0.0]], [4.0,20.0,0.1,2.0,1.0,0.0,42.0]);` |
| `chem_mlp` / `mlp` | CPU single-hidden-layer tanh Adam network with training-prefix-only scaling: [n_train,hidden1..128,epochs1..10000,rate in (0,1],task0regression/1binary,seed,tolerance=1e-8]; bounded work. | `SELECT chem_mlp([[0.0,1.0],[1.0,3.0],[2.0,5.0],[3.0,7.0],[1.5,0.0]], [4.0,4.0,100.0,0.02,0.0,42.0]);` |
| `chem_box_behnken` / `box_behnken` | Box-Behnken two-factor-edge design: 3..64 [low,high] factor bounds; [center_points0..100000=1]; <=100000 runs and1000000 output cells. | `SELECT chem_box_behnken([[-1.0,1.0],[-1.0,1.0],[-1.0,1.0]], []);` |
| `chem_latin_hypercube` / `latin_hypercube` | Seeded jittered Latin hypercube: 1..64 [low,high] factor bounds; [runs1..100000,seed]; <=1000000 output cells. | `SELECT chem_latin_hypercube([[-1.0,1.0],[0.0,10.0]], [4.0,42.0]);` |
| `chem_mixture_lattice` / `mixture_lattice` | Lower-bound-offset simplex lattice filtered by upper bounds: 1..64 [nonnegative_low,high] rows with low<=high; [degree1..1000,total>0=1]; <=1000000 enumerated compositions. | `SELECT chem_mixture_lattice([[0.0,1.0],[0.0,1.0],[0.0,1.0]], [2.0]);` |
| `chem_fractional_factorial` / `fractional_factorial` | Positive-product fractional factorial: <=64 [low,high] bounds; [base_factors1..16,center_points,generator_masks...]; distinct masks select >=2 base columns; alias resolution not inferred. | `SELECT chem_fractional_factorial([[-1.0,1.0],[-1.0,1.0],[-1.0,1.0],[-1.0,1.0]], [3.0,1.0,7.0]);` |
| `chem_orthogonal_screening` / `orthogonal_screening` | Sylvester-Hadamard main-effect screening: 1..64 [low,high] factor bounds; []; interactions may be aliased, not a general definitive or Plackett-Burman design. | `SELECT chem_orthogonal_screening([[-1.0,1.0],[-1.0,1.0],[-1.0,1.0],[-1.0,1.0],[-1.0,1.0]], []);` |
| `chem_randomize_design` / `randomize_design` | Independent seeded shuffle per replicate then equal-size blocks: finite design rows; [replicates1..1000,seed,blocks=1]; rows divisible by blocks and <=100000 output runs. | `SELECT chem_randomize_design([[0.0],[1.0],[2.0],[3.0],[4.0],[5.0]], [2.0,42.0,2.0]);` |
| `chem_design_evaluate` / `design_evaluate` | Full-rank linear-model design diagnostics with added intercept: <=64 model-term columns, n>=terms+1; []; original/coded information matrices and coded spectral condition number. | `SELECT chem_design_evaluate([[-1.0,0.0],[0.0,1.0],[1.0,0.0],[2.0,1.0],[-2.0,-1.0]], []);` |
| `chem_mds` / `mds` | Classical metric Gower MDS, not optimized/nonmetric Stress1: 2..128 feature or distance rows; [dimensions=min(2,n-1),input0features/1distances=0,standardize0/1=0]. | `SELECT chem_mds([[0.0],[1.0],[2.0]], [1.0]);` |
| `chem_correspondence_analysis` / `correspondence_analysis` | Simple contingency-table SVD: nonnegative 2..64 by2..64 counts with positive margins; [dimensions=min(2,r-1,c-1),scaling0principal/1symmetric=0]. | `SELECT chem_correspondence_analysis([[20.0,10.0],[10.0,20.0]], [1.0]);` |
| `chem_factor_analysis` / `factor_analysis` | Iterated principal-axis factoring, not ML: <=64 variables; [factors,rotation0none/1varimax=1,maxiter1..2000=200,tol=1e-8,prior0ONE/1SMC=1,input0rows/1correlation=0]; bounded work. | `SELECT chem_factor_analysis([[1.0,0.64,0.64],[0.64,1.0,0.64],[0.64,0.64,1.0]], [1.0,0.0,500.0,1e-10,1.0,1.0]);` |
| `chem_canonical_correlation` / `canonical_correlation` | Two-block sample-covariance whitening with <=32 variables per block: numeric [X...,Y...] rows; [x_columns,dimensions=min(p,q),standardized_ridge=0]; no MANOVA inference. | `SELECT chem_canonical_correlation([[-1.0,-2.0],[0.0,0.0],[1.0,2.0],[2.0,4.0]], [1.0]);` |
| `chem_gaussian_mixture` / `gaussian_mixture` | Full-covariance constrained EM with <=32 features: [components1..32,maxiter1..2000=200,tol=1e-6,standardized_eigenvalue_floor=1e-6,seed=0,n_init1..16=1]; bounded work/model size. | `SELECT chem_gaussian_mixture([[-1.0],[0.0],[1.0]], [1.0]);` |
| `chem_gaussian_mixture_predict` / `gaussian_mixture_predict` | Saved Gaussian mixture probabilities/densities on same-width feature rows: unchanged gaussian_mixture predict_parameters, <=16384 doubles; reuses fit-time preprocessing. | `SELECT chem_gaussian_mixture_predict([[-1.0],[0.0],[1.0]], [1.0,1.0,0.0,0.0,1.0,1.0,0.0,1.0]);` |
| `chem_self_organizing_map` / `self_organizing_map` | PCA-initialized rectangular batch local-linear Gaussian SOM, <=64 features and256 nodes: [grid_rows,grid_columns,maxiter=100,bandwidth=0.5,standardize=1,tol=1e-6]; bounded work/model size. | `SELECT chem_self_organizing_map([[-1.0],[0.0],[1.0]], [1.0,1.0]);` |
| `chem_som_predict` / `som_predict` | Saved SOM nearest-node predictions on same-width feature rows: unchanged self_organizing_map predict_parameters, <=16384 doubles; distances use fit-preprocessed feature units. | `SELECT chem_som_predict([[-2.0],[0.0],[2.0]], [1.0,2.0,0.0,0.0,1.0,-1.0,1.0]);` |
| `chem_regression_metrics` / `regression_metrics` | [observed,prediction,optional_weight>=0]; []; <=100000 rows, numeric outcome/prediction absolute values<=1e150. Weighted RMSE/MAE/bias/R2 and residuals; constant-target R2 is null. | `SELECT chem_regression_metrics([[1.0,2.0],[2.0,2.0],[3.0,1.0]], []);` |
| `chem_confusion_matrix` / `confusion_matrix` | [actual_integer_label,predicted_integer_label,optional_weight>=0]; []; <=256 observed/predicted levels; actual rows, predicted columns. Undefined class precision/recall null; macro F1 includes union of positive-weight levels. | `SELECT chem_confusion_matrix([[10.0,10.0],[10.0,20.0],[20.0,20.0],[30.0,20.0]], []);` |
| `chem_binary_metrics` / `binary_metrics` | [outcome0/1,probability0..1,optional_weight>=0]; [threshold=.5,probability_clip=1e-15]. Inclusive threshold, clipped natural-log loss, Brier and confusion metrics. | `SELECT chem_binary_metrics([[0.0,0.1],[1.0,0.4],[0.0,0.4],[1.0,0.8]], [0.4]);` |
| `chem_roc_curve` / `roc_curve` | [outcome0/1,finite_score,optional_weight>=0]; []; positive weight in both classes. Descending inclusive thresholds grouped by equal score, trapezoid AUC; first threshold null denotes above maximum. | `SELECT chem_roc_curve([[0.0,0.1],[1.0,0.4],[0.0,0.4],[1.0,0.8]], []);` |
| `chem_precision_recall_curve` / `precision_recall_curve` | [outcome0/1,finite_score,optional_weight>=0]; []; both classes required. Descending grouped thresholds, initial precision1/recall0; step average_precision and trapezoid_area separately. | `SELECT chem_precision_recall_curve([[0.0,0.1],[1.0,0.4],[0.0,0.4],[1.0,0.8]], []);` |
| `chem_lift_curve` / `lift_curve` | [outcome0/1,finite_score,optional_weight>=0]; []; both classes required. Cumulative selection/gain/response/lift with equal-score batches; lift at zero selection null. | `SELECT chem_lift_curve([[0.0,0.1],[1.0,0.4],[0.0,0.4],[1.0,0.8]], []);` |
| `chem_calibration_curve` / `calibration_curve` | [outcome0/1,probability0..1,optional_weight>=0]; [bins1..1000=10]. Equal-width bins include probability1 in last bin; empty-bin means null; ECE and maximum error. | `SELECT chem_calibration_curve([[0.0,0.0],[1.0,0.5],[1.0,1.0]], [4.0]);` |
| `chem_profit_decision` / `profit_decision` | [outcome0/1,probability0..1,optional_weight>=0]; [profit00,profit01,profit10,profit11], actual rows/decision columns. Max expected-profit binary decision, exact ties choose0; totals outside f64 range null. | `SELECT chem_profit_decision([[0.0,0.1],[1.0,0.6],[1.0,0.2]], [2.0,-1.0,-5.0,10.0]);` |
| `chem_ensemble_predict` / `ensemble_predict` | Rows of member predictions; [] equal weights or [nonnegative_weight_per_member...]. Up to1024 members, stable normalized weighted average; no fitting or automatic weight selection. | `SELECT chem_ensemble_predict([[1.0,3.0],[3.0,7.0]], [1.0,3.0]);` |
| `chem_model_compare` / `model_compare` | [observed,prediction1,...] <=64 models; [task0regression/1binary=0,metric0RMSEorLogLoss/1MAEorBrier=0]. Rank caller-supplied predictions, no refitting; binary clipped log loss. | `SELECT chem_model_compare([[0.0,0.0,1.0],[1.0,1.0,0.4],[2.0,2.0,1.0]], []);` |
| `chem_holdout_split` / `holdout_split` | Finite rows, scalar integer class labels if stratified; [train_fraction in(0,1),seed,stratified0/1=0]. Per-group floor clamped to keep train/test nonempty, every stratum>=2; zero-based indices. | `SELECT chem_holdout_split([[0.0],[0.0],[1.0],[1.0]], [0.5,42.0,1.0]);` |
| `chem_rolling_origin` / `rolling_origin` | Chronologically ordered finite rows; [initial_train>=1,horizon>=1,step>=1=1,window0expanding=0]. Full-horizon folds only; training precedes validation, <=1000000 output indices. | `SELECT chem_rolling_origin([[0.0],[1.0],[2.0],[3.0],[4.0],[5.0],[6.0],[7.0]], [3.0,2.0,2.0,3.0]);` |
| `chem_curve_fit` / `curve_fit` | Nonlinear residual least squares, 10 curve families: [x,y] rows; [family0..9,max_iterations=300,tolerance=1e-9,optional_initial_coefficients...]. | `SELECT chem_curve_fit([[0.0,2.03],[0.5,2.64],[1.0,3.7],[1.5,4.88],[2.0,6.68],[2.5,8.94]], [0.0,500.0,1e-10,1.8,0.4]);` |
| `chem_curve_predict` / `curve_predict` | Saved nonlinear curve prediction: [x] rows; fit.predict_parameters=[family,coefficients...]. | `SELECT chem_curve_predict([[0.0],[1.0]], [0.0,2.0,0.6]);` |
| `chem_bspline_basis` / `bspline_basis` | Clamped Cox-de Boor B-spline basis: [x] rows; [degree1..3,full_clamped_knots...]. | `SELECT chem_bspline_basis([[0.0],[0.5],[1.0],[2.0]], [3.0,0.0,0.0,0.0,0.0,1.0,2.0,2.0,2.0,2.0]);` |
| `chem_spline_fit` / `spline_fit` | B-spline regression with second-difference coefficient penalty: [x,y] rows; [degree1..3,lambda,full_clamped_knots...]. | `SELECT chem_spline_fit([[0.0,1.0],[0.25,1.5],[0.5,2.0],[0.75,2.5],[1.0,3.0]], [1.0,0.0,0.0,0.0,0.5,1.0,1.0]);` |
| `chem_spline_predict` / `spline_predict` | Saved penalized B-spline prediction: [x] rows; fit.predict_parameters=[degree,knot_count,response_origin,knots...,offset_coefficients...]. | `SELECT chem_spline_predict([[0.0],[0.5],[1.0]], [1.0,4.0,0.0,0.0,0.0,1.0,1.0,1.0,3.0]);` |
| `chem_loess` / `loess` | Tricube local polynomial smoothing with optional robust bisquare weights: [x,y] rows; [span,degree0..2=1,robust_iterations0..10=0]. | `SELECT chem_loess([[0.0,1.0],[1.0,3.0],[2.0,5.0],[3.0,7.0],[4.0,9.0]], [1.0,1.0,0.0]);` |
| `chem_loess_predict` / `loess_predict` | Saved LOESS neighborhood prediction: [x] rows; fit.predict_parameters=[span,degree,n_train,x_train...,y_train...,robust_weights...]. | `SELECT chem_loess_predict([[0.5],[2.5]], [1.0,1.0,5.0,0.0,1.0,2.0,3.0,4.0,1.0,3.0,5.0,7.0,9.0,1.0,1.0,1.0,1.0,1.0]);` |
| `chem_spectral_snv` / `spectral_snv` | Per-spectrum standard normal variate: curve rows on a common grid; [variance_ddof=1]. | `SELECT chem_spectral_snv([[1.0,2.0,3.0],[10.0,12.0,14.0]], []);` |
| `chem_spectral_msc` / `spectral_msc` | MSC reference regression: curve rows; [] fits their mean reference, or args are the saved reference spectrum. | `SELECT chem_spectral_msc([[5.0,7.0,11.0,19.0]], [1.0,2.0,4.0,8.0]);` |
| `chem_savitzky_golay` / `savitzky_golay` | Fixed Savitzky-Golay filter or derivative: equally spaced curve rows; [odd_window,degree=2,derivative=0,spacing=1]. | `SELECT chem_savitzky_golay([[0.0,1.0,4.0,9.0,16.0,25.0,36.0]], [5.0,2.0,1.0,0.5]);` |
| `chem_curve_derivative` / `curve_derivative` | Finite curve derivatives: equally spaced curve rows; [order=1,spacing=1]. | `SELECT chem_curve_derivative([[0.0,1.0,4.0,9.0,16.0]], [2.0,0.5]);` |
| `chem_spectral_detrend` / `spectral_detrend` | Global polynomial baseline removal: curve rows; [degree=1] on a normalized equally spaced grid. | `SELECT chem_spectral_detrend([[2.0,3.0,6.0,11.0,18.0]], [2.0]);` |
| `chem_haar_transform` / `haar_transform` | Orthonormal Haar forward/inverse: power-of-two curve rows; [inverse=0,levels=full]. | `SELECT chem_haar_transform([[1.0,2.0,3.0,4.0]], []);` |
| `chem_fourier_transform` / `fourier_transform` | Direct discrete Fourier forward/inverse: real curve rows or packed complex rows; [inverse=0]. | `SELECT chem_fourier_transform([[1.0,2.0,3.0,4.0]], []);` |
| `chem_dynamic_time_warping` / `dynamic_time_warping` | Pairwise DTW distance and path: [value,series_id0or1] rows; [window_radius=unbounded]. | `SELECT chem_dynamic_time_warping([[1.0,0.0],[2.0,0.0],[3.0,0.0],[1.0,1.0],[1.0,1.0],[2.0,1.0],[3.0,1.0]], []);` |
| `chem_functional_pca` / `functional_pca` | Common-grid functional PCA: curve rows; [components,spacing=1], trapezoid L2 metric. | `SELECT chem_functional_pca([[1.0,2.0,3.0],[2.0,4.0,6.0],[3.0,6.0,9.0]], [1.0]);` |
| `chem_peak_find` / `peak_find` | Sampled local peaks with global prominence: curve rows; [minimum_prominence=0,minimum_distance=1]. | `SELECT chem_peak_find([[0.0,3.0,1.0,4.0,0.0,2.0,0.0]], []);` |
| `chem_manova` / `manova` | Whole-model numeric MANOVA with intercept: [X...,Y...] rows, <=32predictors/responses; [x_columns]; H/E and4criteria, exact rank-one F otherwise JMP approximate F with Roy p lower bound. | `SELECT chem_manova([[0.0,-1.0,-1.0],[0.0,-1.0,1.0],[0.0,1.0,-1.0],[0.0,1.0,1.0],[1.0,1.0,0.0],[1.0,1.0,2.0],[1.0,3.0,0.0],[1.0,3.0,2.0]], [1.0]);` |
| `chem_canonical_correlation_inference` / `canonical_correlation_inference` | Sequential Wilks/Bartlett asymptotic chi-square: sorted unregularized correlations as one-column rows; [n_samples>p+q,p1..64,q1..64]; all min(p,q) correlations required; not raw-data fitting. | `SELECT chem_canonical_correlation_inference([[0.8],[0.6]], [100.0,2.0,2.0]);` |
| `chem_multiple_correspondence_analysis` / `multiple_correspondence_analysis` | Complete integer-coded MCA: 2..32categorical variables, <=64observedcategories; [dimensions=min(2,n-1,categories-variables)]; indicator/Burt inertias, principal scores, Benzecri/Greenacre adjustments. | `SELECT chem_multiple_correspondence_analysis([[0.0,0.0],[0.0,0.0],[1.0,1.0],[1.0,1.0]], [1.0]);` |
| `chem_factor_analysis_ml` / `factor_analysis_ml` | Gaussian ML correlation-factor model with2..24variables: [factors,rotation0none/1varimax=1,maxiter1..1000=200,tol=1e-8,uniqueness_floor=1e-6,input0rows/1correlation=0,n_correlation_samples=0]; bounded local optimizer. | `SELECT chem_factor_analysis_ml([[1.0,0.49,0.49,0.49],[0.49,1.0,0.49,0.49],[0.49,0.49,1.0,0.49],[0.49,0.49,0.49,1.0]], [1.0,0.0,300.0,1e-9,1e-6,1.0,200.0]);` |
| `chem_svm_fit` / `svm_fit` | Binary soft-margin C-SVC with dual KKT diagnostics: [features...,label0/1]; [n_train,kernel0linear/1RBF=1,C=1,gamma=1/p,tolerance=1e-6,max_iterations=100000]. | `SELECT chem_svm_fit([[-2.0,0.0],[-1.0,0.0],[1.0,1.0],[2.0,1.0]], [4.0,0.0,10.0]);` |
| `chem_svm_predict` / `svm_predict` | Saved binary C-SVC prediction with fixed training preprocessing: feature rows; svm_fit.predict_parameters. | `SELECT chem_svm_predict([[-2.0],[2.0]], [1.0,0.0,1.0,0.0,2.0,0.0,0.0,1.0,-0.5,0.5,-1.0,1.0]);` |
| `chem_quantile_regression` / `quantile_regression` | Convex check-loss regression with ADMM residual convergence: [features...,target]; [quantile=.5,rho=1,max_iterations=10000,tolerance=1e-7]. | `SELECT chem_quantile_regression([[0.0,0.0],[1.0,1.0],[2.0,4.0],[3.0,3.0],[4.0,4.0]], [0.5,1.0,20000.0,1e-7]);` |
| `chem_huber_regression` / `huber_regression` | Fixed-scale Huber M-regression by IRLS: [features...,target]; [delta=1.345,scale=0OLS_MAD,max_iterations=200,tolerance=1e-8]. | `SELECT chem_huber_regression([[0.0,1.0],[1.0,3.0],[2.0,5.0],[3.0,27.0],[4.0,9.0],[5.0,11.0],[6.0,13.0]], [1.345,1.0]);` |
| `chem_robust_predict` / `robust_predict` | Saved Huber/quantile linear prediction with stable origin offsets: feature rows; fit.predict_parameters. | `SELECT chem_robust_predict([[2.0],[3.0]], [1.0,0.0,0.0,1.0,1.0,0.0,0.0,1.0,2.0]);` |
| `chem_marker_statistics` / `marker_statistics` | Diploid genotype call rows (individuals x markers); 0/1/2,missing=-1; []. Call rates, allele/MAF/heterozygosity/PIC, asymptotic Pearson HWE df1 and BH p-values; monomorphic/all-missing HWE null. | `SELECT chem_marker_statistics([[0.0,0.0],[1.0,1.0],[1.0,2.0],[2.0,-1.0]], []);` |
| `chem_marker_impute` / `marker_impute` | Diploid 0/1/2/-1 matrix; [method0mean/1mode/2LDknn/3HWE/4observed/5uniform/6specified=0,nearest_markers=10,nearest_samples=5,seed=0,specified=0]. All-missing loci rejected; mean/mode save values. LD-kNN squared dosage correlation, mean taxicab, majority call; no valid neighbors falls back mode. | `SELECT chem_marker_impute([[0.0,0.0],[1.0,1.0],[2.0,2.0],[-1.0,2.0]], [2.0,1.0,1.0]);` |
| `chem_marker_impute_predict` / `marker_impute_predict` | New dosage rows in [0,2] or missing=-1; saved mean/mode predict_parameters (one per marker). No query-driven fitting. | `SELECT chem_marker_impute_predict([[-1.0,-1.0],[2.0,0.0]], [0.5,1.5]);` |
| `chem_marker_relatedness` / `marker_relatedness` | Dosages in [0,2] or missing=-1; [method0VanRaden1/1VanRaden2/2diploidDominance=0]. <=1000individuals and100Mpairfeaturework. Additive missing mean-imputed; dominance requires discrete calls and missing centered effect0; monomorphic markers omitted. | `SELECT chem_marker_relatedness([[0.0,0.0],[1.0,1.0],[2.0,2.0]], [0.0]);` |
| `chem_distance_matrix` / `distance_matrix` | Numeric row vectors,abs<=1e150; [method0Euclidean/1booleanJaccard/2nonnegativeBrayCurtis/3numericGower/4Hamming=0,optional_missing_code]. <=1000rows,100Mpairfeaturework; pairwise available coordinates, no-common errors; Gower constant columns contribute zero. | `SELECT chem_distance_matrix([[0.0,1.0],[2.0,1.0]], [0.0]);` |
| `chem_marker_simulate` / `marker_simulate` | Scalar rows of allele frequencies in [0,1]; [samples,seed,inbreeding0..1=0]. Independent diploid loci with P(0)=q^2+Fpq,P(1)=2pq(1-F),P(2)=p^2+Fpq; <=1Moutputcells and1024markers. | `SELECT chem_marker_simulate([[0.25],[0.5]], [3.0,42.0,0.0]);` |
| `chem_marker_cross` / `marker_cross` | Complete parental genotype rows0/1/2; [parent0_index,parent1_index,offspring_count,seed]. Single-generation unlinked diploid meiosis, independent loci; <=1Moutputcells. | `SELECT chem_marker_cross([[0.0,2.0,1.0],[2.0,0.0,1.0]], [0.0,1.0,3.0,42.0]);` |
| `chem_genomic_blup` / `genomic_blup` | Rows[dosages...,phenotype], query suffix phenotype ignored; [n_train>=2,lambda1e-8..1e12=1,method0VanRaden1/1VanRaden2=0]. <=64markers,10000rows. Fixed residual/genetic variance ratio, intercept + ridge-equivalent marker effects; training-prefix means, no REML/GWAS inference. | `SELECT chem_genomic_blup([[0.0,1.0],[1.0,2.0],[2.0,5.0],[1.0,999.0]], [3.0,1.0]);` |
| `chem_genomic_blup_predict` / `genomic_blup_predict` | New dosage rows[0,2] or missing=-1; saved genomic_blup predict_parameters=[p,y_origin,y_mean_offset,training_marker_means(p),marker_effects(p)]. Fixed training preprocessing and missing mean substitution. | `SELECT chem_genomic_blup_predict([[0.0],[1.0],[2.0]], [1.0,1.0,1.6666666666666667,1.0,1.6]);` |
| `chem_aft_exponential` / `aft_exponential` | Right-censored exponential AFT: [features...,positive_time,event0/1]; [ridge=0,maxiter=200,tolerance=1e-8]. | `SELECT chem_aft_exponential([[1.0,1.0],[2.0,1.0],[3.0,0.0],[4.0,1.0],[5.0,0.0]], []);` |
| `chem_aft_weibull` / `aft_weibull` | Right-censored weibull AFT: [features...,positive_time,event0/1]; [ridge=0,maxiter=200,tolerance=1e-8]. | `SELECT chem_aft_weibull([[1.0,1.0],[2.0,1.0],[3.0,0.0],[4.0,1.0],[5.0,0.0]], []);` |
| `chem_aft_lognormal` / `aft_lognormal` | Right-censored lognormal AFT: [features...,positive_time,event0/1]; [ridge=0,maxiter=200,tolerance=1e-8]. | `SELECT chem_aft_lognormal([[1.0,1.0],[2.0,1.0],[3.0,0.0],[4.0,1.0],[5.0,0.0]], []);` |
| `chem_aft_predict` / `aft_predict` | Saved AFT survival/quantile prediction: [features...,query] rows; [mode0survival/1quantile,predict_parameters...]. | `SELECT chem_aft_predict([[0.0],[5.0],[10.0]], [0.0,0.0,0.0,0.0,0.0,1.6094379124341003]);` |
| `chem_nhpp_powerlaw` / `nhpp_powerlaw` | Single-system time-truncated power NHPP MLE: scalar exact recurrence times; [observation_end]. | `SELECT chem_nhpp_powerlaw([[1.0],[2.0],[4.0],[7.0]], [10.0]);` |
| `chem_nhpp_predict` / `nhpp_predict` | Power NHPP cumulative/intensity prediction: nonnegative scalar times; [shape,log_eta]. | `SELECT chem_nhpp_predict([[0.0],[5.0],[10.0]], [2.0,1.6094379124341003]);` |
| `chem_reliability_blocks` / `reliability_blocks` | Independent reliability blocks: component probability rows; [structure0series/1parallel/2kofn,k_if2]. | `SELECT chem_reliability_blocks([[0.8,0.7,0.9]], [2.0,2.0]);` |

| `chem_repeated_covariance` / `repeated_covariance` | Nine repeated covariance matrices: scalar increasing times; [structure_code, physical parameters...]. | `SELECT chem_repeated_covariance([[0.0],[1.0],[2.0]], [4.0,1.0,0.5]);` |
| `chem_repeated_model` / `repeated_model` | Gaussian repeated measures: [features...,y,subject,time]; [structure,method=1,maxiter=300,tolerance=1e-8,known covariance if method2]. | `SELECT chem_repeated_model([[0.0,1.1,0.0,0.0],[1.0,1.8,0.0,1.0],[2.0,2.1,0.0,2.0],[0.0,0.8,1.0,0.0],[2.0,2.4,1.0,2.0],[0.0,1.4,2.0,0.0],[1.0,1.9,2.0,1.0],[2.0,2.8,2.0,2.0]], [4.0,2.0,300.0,1e-08,0.7,0.4]);` |
| `chem_random_slope_model` / `random_slope_model` | Correlated Gaussian random slopes: [other_features...,slope,y,subject]; [method=1,maxiter=300,tolerance=1e-8]. | `SELECT chem_random_slope_model([[-1.5,-0.3999999999999999,0.0],[0.0,2.447801055611894,0.0],[1.5,4.800852965387551,0.0],[-1.43,0.9940080946479155,1.0],[-0.43,2.5060530045266307,1.0],[0.5700000000000001,3.8449574698911446,1.0],[1.57,5.7784240293359295,1.0],[-1.3599999999999999,1.6975969779199587,2.0],[-0.61,2.2441774366691205,2.0],[0.14,2.9939728231014513,2.0],[0.89,4.078268384613996,2.0],[1.6400000000000001,4.508643726327658,2.0],[-1.29,1.3353075168384074,3.0],[0.21000000000000002,1.9584611958404785,3.0],[1.71,3.0969290109247303,3.0],[-1.22,0.7311967164783243,4.0],[-0.21999999999999997,1.053909607842823,4.0],[0.78,2.041094518889131,4.0],[1.78,2.641025660787466,4.0],[-1.15,0.7112645756966788,5.0],[-0.39999999999999997,1.3007179446757398,5.0],[0.35000000000000003,2.4936964602489375,5.0],[1.1,3.0659463709565546,5.0],[1.85,3.871040404872158,5.0],[-1.08,1.113206127673739,6.0],[0.42000000000000004,3.4749909802800594,6.0],[1.92,6.188520632109679,6.0],[-1.01,1.183904059437987,7.0],[-0.25999999999999995,2.6856177920234843,7.0],[0.49000000000000005,4.176074956963551,7.0],[1.24,5.179436461116708,7.0],[1.99,6.850366871258961,7.0]], []);` |
| `chem_random_slope_predict` / `random_slope_predict` | Saved random slope prediction: [other_features...,slope,subject]; exact returned predict_parameters. | `SELECT chem_random_slope_predict([[0.0,0.0],[0.5,999.0]], [1,-0.3999999999999999,3.2327332447796517,1.71320244537735,-0.021116948503145182,-1.5,1.7559375000000002,1.134468100518366,0.8414294635617765,8,0,-0.028735991416393714,0.286619733520204,1,0.4445314172359469,0.20056681303699997,2,0.23636025116999793,-0.16587080212282412,3,-0.3765353663119541,-0.43652515986998724,4,-0.6876694303430354,-0.3976689798479615,5,-0.348895922670813,-0.1327743040990701,6,0.3024791285476626,0.27195226638001857,7,0.4584659137885829,0.37370043300261346]);` |
| `chem_glmm_binomial` / `glmm_binomial` | Bernoulli logit Gaussian random-intercept integrated ML: [fixed features...,y,group]; [maxiter=200,tol=1e-6,odd quadrature order=15,variance_upper=100]. | `SELECT chem_glmm_binomial([[0.0,0.0],[1.0,0.0],[0.0,1.0],[1.0,1.0]], []);` |
| `chem_glmm_poisson` / `glmm_poisson` | Poisson log Gaussian random-intercept integrated ML: [fixed features...,y,group]; [maxiter=200,tol=1e-6,odd quadrature order=15,variance_upper=100]. | `SELECT chem_glmm_poisson([[1.0,0.0],[3.0,0.0],[1.0,1.0],[3.0,1.0]], []);` |
| `chem_glmm_predict` / `glmm_predict` | Saved Gaussian random-intercept GLMM prediction: [features...,group] and unchanged predict_parameters; separate conditional,posterior and population means. | `SELECT chem_glmm_predict([[7.0],[99.0]], [1.0,1.0,0.0,1.0,3.0,0.5,0.6931471805599453,7.0,0.2,0.2,0.1,-0.34772255750516606,0.2,0.747722557505166,0.16666666666666666,0.6666666666666666,0.16666666666666666]);` |
| `chem_propensity_fit` / `propensity_fit` | Training-prefix logistic propensity model: [features...,treatment0/1]; [n_train,ridge=0,max_iterations=100,tolerance=1e-8]. Suffix labels ignored; saved preprocessing. | `SELECT chem_propensity_fit([[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,1.0],[1.0,0.0],[1.0,1.0],[1.0,1.0],[1.0,1.0]], [8.0]);` |
| `chem_propensity_predict` / `propensity_predict` | Saved training-only logistic propensity scoring: feature rows; propensity_fit.predict_parameters. | `SELECT chem_propensity_predict([[0.0],[1.0]], [1.0,0.0,0.0,0.5,0.5,1.0986122886681098]);` |
| `chem_causal_balance` / `causal_balance` | Raw-score overlap and weighted covariate balance: [features...,treatment,propensity]; [estimand0ATE/1ATT=0,explicit_clip=0]. Fixed baseline pooled SD, ESS; both groups>=2. | `SELECT chem_causal_balance([[0.0,1.0,0.5],[2.0,1.0,0.75],[1.0,0.0,0.25],[3.0,0.0,0.5]], []);` |
| `chem_propensity_match` / `propensity_match` | Greedy1:1 ATT matching with explicit absolute caliper: [outcome,treatment,propensity]; [caliper=1,replacement1/0=1,metric0probability/1logit=0]. Matched-subset effect; SE null. | `SELECT chem_propensity_match([[5.0,1.0,0.25],[8.0,1.0,0.75],[9.0,1.0,0.875],[1.0,0.0,0.125],[2.0,0.0,0.375],[4.0,0.0,0.875]], [0.2,1.0,0.0]);` |
| `chem_ipw_effect` / `ipw_effect` | Hajek IPW ATE/ATT with fixed-nuisance delta-method SE: [outcome,treatment,propensity]; [estimand0ATE/1ATT=0,explicit_clip=0]. | `SELECT chem_ipw_effect([[4.0,1.0,0.5],[6.0,1.0,0.75],[1.0,0.0,0.25],[3.0,0.0,0.5]], []);` |
| `chem_aipw_effect` / `aipw_effect` | Canonical AIPW ATE/ATT scores with influence SE: [outcome,treatment,propensity,m0,m1]; [estimand0ATE/1ATT=0,explicit_clip=0]. External nuisance fitting required. | `SELECT chem_aipw_effect([[4.0,1.0,0.5,2.0,4.0],[6.0,1.0,0.75,3.0,6.0],[1.0,0.0,0.25,1.0,3.0],[3.0,0.0,0.5,3.0,5.0]], []);` |
| `chem_uplift_curve` / `uplift_curve` | Supplied-score uplift/Qini evaluation: [factual_outcome,treatment,predicted_effect_score,assignment_propensity]; []. Grouped-score HT gain; no uplift-tree fitting. | `SELECT chem_uplift_curve([[3.0,1.0,0.9,0.5],[1.0,0.0,0.8,0.5],[2.0,1.0,0.8,0.5],[2.0,0.0,0.1,0.5]], []);` |
| `chem_assay_normalize` / `assay_normalize` | Seven assay row normalizations: [method0standardize/1CPM/2RPKM/3TPM/4CLR/5percent/6ratio,method-specific args]. Samples rows/features columns; lengths in base pairs. | `SELECT chem_assay_normalize([[1.0,3.0],[2.0,6.0]], [3.0,1000.0,2000.0]);` |
| `chem_quantile_normalize` / `quantile_normalize` | Equal-rank across-sample mean quantiles; [] args. Equal input ties receive the mean over tied reference ranks. Saves reference_quantiles. | `SELECT chem_quantile_normalize([[5.0,1.0,3.0],[4.0,2.0,6.0]], []);` |
| `chem_quantile_normalize_predict` / `quantile_normalize_predict` | New samples scored against saved sorted reference_quantiles matching row width, with tie averaging; no refitting. | `SELECT chem_quantile_normalize_predict([[4.0,4.0,0.0]], [1.5,2.5,5.5]);` |
| `chem_marker_admixture` / `marker_admixture` | Diploid independent-locus binomial admixture EM: [K,maxiter=1000,tolerance=1e-8,seed=0,n_init=3,frequency_floor=1e-8]. Genotypes0/1/2,-1missing. Joint Q/F local likelihood, convergence diagnostic, saved fixed-F inference. | `SELECT chem_marker_admixture([[0.0,0.0],[1.0,2.0],[2.0,-1.0]], [1.0,1000.0,1e-9,42.0,1.0]);` |
| `chem_admixture_predict` / `admixture_predict` | Infer ancestry of new diploid rows with saved fixed allele frequencies; exact marker_admixture predict_parameters; each starting ancestry weight strictly positive. | `SELECT chem_admixture_predict([[0.0,2.0],[1.0,1.0]], [2.0,2.0,10000.0,1e-10,1e-8,0.5,0.5,0.1,0.9,0.9,0.1]);` |
| `chem_rcbd_design` / `rcbd_design` | Randomized complete blocks: distinct treatment rows; [blocks,seed=0]. | `SELECT chem_rcbd_design([[0.0],[1.0],[2.0],[3.0]], [3.0,42.0]);` |
| `chem_balanced_incomplete_block` / `balanced_incomplete_block` | Exact balanced incomplete blocks: distinct treatments<=9; [block_size,lambda=1,seed=0]. | `SELECT chem_balanced_incomplete_block([[0.0],[1.0],[2.0],[3.0],[4.0],[5.0],[6.0]], [3.0,1.0,7.0]);` |
| `chem_resolvable_incomplete_block` / `resolvable_incomplete_block` | Resolvable affine incomplete blocks: q^2 treatment rows,q=2/3/5/7; [seed=0]. | `SELECT chem_resolvable_incomplete_block([[0.0],[1.0],[2.0],[3.0],[4.0],[5.0],[6.0],[7.0],[8.0]], [17.0]);` |
| `chem_split_plot_design` / `split_plot_design` | Classical two-level split plot: factor rows[low,high,role0whole/1subplot]; [replicates=1,seed=0]. | `SELECT chem_split_plot_design([[0.0,1.0,0.0],[2.0,3.0,1.0],[4.0,5.0,1.0]], [2.0,31.0]);` |
| `chem_strip_plot_design` / `strip_plot_design` | Classical crossed strip plot: factor rows[low,high,role0row/1column]; [replicates=1,seed=0]. | `SELECT chem_strip_plot_design([[0.0,1.0,0.0],[2.0,3.0,1.0],[4.0,5.0,1.0]], [2.0,31.0]);` |
| `chem_taguchi_orthogonal_array` / `taguchi_orthogonal_array` | Catalog orthogonal arrays: scalar factor levels2/3; [array4/8/9/18,seed=0]. | `SELECT chem_taguchi_orthogonal_array([[2.0],[3.0],[3.0],[3.0],[3.0],[3.0],[3.0],[3.0]], [18.0,11.0]);` |
| `chem_association_rules` / `association_rules` | Binary transaction rows/items columns,including empty baskets; [minimum_support=.1,minimum_confidence=.4,minimum_lift=1.2,max_antecedents=3,max_rule_size=4,max_rules=10000]. Complete Apriori within explicit limits,arbitrary disjoint consequents. | `SELECT chem_association_rules([[1.0,1.0,0.0],[1.0,1.0,1.0],[1.0,0.0,1.0],[0.0,1.0,0.0],[0.0,0.0,0.0]], [0.4,0.6,1.5,1.0,2.0]);` |
| `chem_basket_affinity` / `basket_affinity` | Binary transaction rows,<=32items; []. All distinct pairs: support,directional confidence,lift,Jaccard,leverage; undefined divisions JSONnull. | `SELECT chem_basket_affinity([[1.0,1.0,0.0],[1.0,1.0,1.0],[1.0,0.0,1.0],[0.0,1.0,0.0],[0.0,0.0,0.0]], []);` |
| `chem_item_reliability` / `item_reliability` | Complete finite respondents/items; [standardized0raw/1correlation=0]. Cronbach alpha,deleted-item alpha,corrected item-total correlations; negative alpha retained,undefined or nonrepresentable quantities null. | `SELECT chem_item_reliability([[0.0,0.0,0.0],[1.0,2.0,1.0],[2.0,4.0,0.0],[3.0,6.0,1.0]], []);` |

| `chem_kinship_model` / `kinship_model` | Gaussian y=Xbeta+g+e with supplied PSD kinship K; rows[y,fixed covariates...,K row]; [p,method0ML/1REML=1,h_upper=.99999999,grid=128,tol=1e-10]. Bounded nonnegative variance components; profile likelihood and BLUP. | `SELECT chem_kinship_model([[0.0,1.0,1.0,0.0,0.0,0.0,0.0],[2.0,1.0,1.0,0.0,0.0,0.0,0.0],[4.0,0.0,0.0,1.0,1.0,0.0,0.0],[6.0,0.0,0.0,1.0,1.0,0.0,0.0],[8.0,0.0,0.0,0.0,0.0,1.0,1.0],[10.0,0.0,0.0,0.0,0.0,1.0,1.0]], [0.0,1.0]);` |
| `chem_kinship_predict` / `kinship_predict` | Saved Gaussian genomic model scoring: rows[fixed covariates(p),cross-kinship to n training individuals]; unchanged kinship_model.predict_parameters. Reuses fixed effects and training kernel weights. | `SELECT chem_kinship_predict([[1,1,0,0,0,0],[0,0,0,0,0,0]], [1,6,0,0,0.5,10,0,0.9999999999999999,-0.18749999984346133,-0.18749999984346133,0,0,0.18749999984346133,0.18749999984346133]);` |
| `chem_irt_fit` / `irt_fit` | Bounded1PL/2PL/3PL marginal item likelihood under fixed N(0,1): item rows,optional last integer frequency; [model=2,freq=0,missing=0,maxiter=300,tol=1e-6,oddq=41,cmax=.5,guess_penalty=0]. | `SELECT chem_irt_fit([[0.0,0.0,30.0],[0.0,1.0,20.0],[1.0,0.0,20.0],[1.0,1.0,30.0]], [1.0,1.0,0.0,300.0,1e-08,61.0,0.5,0.0]);` |
| `chem_irt_predict` / `irt_predict` | Saved IRT characteristic/information curves: scalar ability rows and unchanged model_parameters; D1 logistic probabilities,individual/test Fisher information. | `SELECT chem_irt_predict([[0.5],[-1e+300],[1e+300]], [1.0,3.0,1.0,0.0,61.0,2.0,0.5,0.2]);` |
| `chem_irt_score` / `irt_score` | Posterior EAP ability scoring: binary item response rows,optional -1missing under saved policy,and unchanged model_parameters; finite extreme scores with quadrature diagnostics. | `SELECT chem_irt_score([[0.0,0.0],[0.0,1.0],[1.0,1.0]], [1.0,1.0,2.0,0.0,61.0,1.0,0.0,0.0,1.0,0.0,0.0]);` |
| `chem_cfa_model` / `cfa_model` | Gaussian simple-structure CFA: rows[S_MLrow(p),mean]; [N,q,correlated,maxiter,tolerance,assignments(p)]. | `SELECT chem_cfa_model([[1.0,0.7200000000000001,0.5599999999999999,0.23799999999999996,0.20999999999999996,0.242,10.0],[0.7200000000000001,1.0,0.63,0.26775,0.23625000000000002,0.20475000000000002,-2.0],[0.5599999999999999,0.63,1.0,0.20824999999999996,0.18374999999999997,0.15924999999999997,0.0],[0.238,0.26775,0.20825,1.0,0.6375,0.5525,5.0],[0.20999999999999996,0.23624999999999996,0.18374999999999997,0.6375,1.0,0.48750000000000004,1.0],[0.242,0.20475,0.15924999999999997,0.5525,0.48750000000000004,1.0,-3.0]], [300.0,2.0,1.0,500.0,1e-08,0.0,0.0,0.0,1.0,1.0,1.0]);` |
| `chem_cfa_scores` / `cfa_scores` | Saved CFA factor scores: observation rows; [method0regression/1Bartlett,score_parameters...]. | `SELECT chem_cfa_scores([[2.6,3.4,4.8]], [1.0,3.0,1.0,1.0,2.0,3.0,0.8,0.7,0.9,1.0,0.36,0.51,0.19]);` |
| `chem_sem_fit_indices` / `sem_fit_indices` | Gaussian SEM fit indices: rows[S_MLrow(p),Sigma row(p),observedmean,modelmean]; [N,kfree]. | `SELECT chem_sem_fit_indices([[1.0,0.3,1.2,0.2,1.0,1.1],[0.3,1.0,0.2,0.8,-2.0,-2.2]], [120.0,2.0]);` |
| `chem_sem_ram` / `sem_ram` | Supplied RAM mean/covariance: rows[Brow(p),Omega row(p),intercept]; no parameters. | `SELECT chem_sem_ram([[0.0,0.0,0.0,1.0,0.0,0.0,1.0],[0.5,0.0,0.0,0.0,0.8,0.0,2.0],[0.2,0.7,0.0,0.0,0.0,0.6,3.0]], []);` |
| `chem_degradation_path` / `degradation_path` | Transformed fixed-unit linear paths: [time,response,unit_id]; [intercept_mode=0,slope_mode=0,response_transform=0,time_transform=0]. | `SELECT chem_degradation_path([[0.0,10.1,1.0],[1.0,7.9,1.0],[2.0,5.9,1.0],[3.0,4.1,1.0],[0.0,12.1,2.0],[1.0,9.9,2.0],[2.0,7.9,2.0],[3.0,6.1,2.0]], [0.0,1.0,0.0,0.0]);` |
| `chem_degradation_predict` / `degradation_predict` | Saved path predictions and threshold inverse confidence sets: [time,unit_id,threshold]; [confidence,side,model...]. | `SELECT chem_degradation_predict([[2.0,7.0,6.0]], [0.95,0.0,1.0,0.0,0.0,1.0,7.0,10.0,-2.0,0.04,-0.01,0.01,0.16,8.0]);` |
| `chem_stability_analysis` / `stability_analysis` | Batch linear stability model selection and shelf life: [time,response,batch_id]; [lower_on,lower,upper_on,upper,pooling_alpha,confidence,pool_separate_MSE,horizon]. | `SELECT chem_stability_analysis([[0.0,10.1,1.0],[1.0,7.9,1.0],[2.0,5.9,1.0],[3.0,4.1,1.0],[0.0,12.1,2.0],[1.0,9.9,2.0],[2.0,7.9,2.0],[3.0,6.1,2.0]], [1.0,3.0,0.0,0.0,0.25,0.95,0.0,20.0]);` |
| `chem_destructive_degradation` / `destructive_degradation` | Independent-unit Gaussian/lognormal censored degradation ML: [time,stress,lower,upper,type]; [family=0,time_transform=0,stress_mode=0,maxiter=300,tolerance=1e-8]. | `SELECT chem_destructive_degradation([[0.0,0.0,10.1,10.1,0.0],[1.0,0.0,7.9,7.9,0.0],[2.0,0.0,5.9,5.9,0.0],[3.0,0.0,4.1,4.1,0.0]], []);` |
| `chem_destructive_predict` / `destructive_predict` | Saved independent-unit destructive distribution: [time,stress,threshold]; unchanged destructive_degradation.predict_parameters. | `SELECT chem_destructive_predict([[5.0,0.0,2.718281828459045]], [1.0,1.0,0.0,0.0,0.5,2.0,1.0,0.0,0.0,-0.2]);` |
| `chem_cumulative_damage_fit` / `cumulative_damage_fit` | Step-stress cumulative-exposure ML: [lo_time,hi_time,type,pattern_id]; [distribution,relation,continuation,maxiter,tolerance,pattern_count,patterns...]. | `SELECT chem_cumulative_damage_fit([[1.0,1.0,0.0,1.0],[3.0,3.0,0.0,1.0],[2.0,2.0,0.0,2.0],[6.0,6.0,0.0,2.0]], [0.0,0.0,1.0,300.0,1e-09,2.0,1.0,1.0,10.0,0.0,2.0,1.0,10.0,1.0]);` |
| `chem_cumulative_damage_predict` / `cumulative_damage_predict` | Saved cumulative-exposure reliability: [time,pattern_id,probability]; unchanged cumulative_damage_fit.predict_parameters. | `SELECT chem_cumulative_damage_predict([[4.0,7.0,0.5]], [1.0,0.0,0.0,1.0,0.0,0.6931471805599453,1.0,1.0,7.0,2.0,2.0,0.0,3.0,1.0]);` |

| `chem_kalman_smoother` / `kalman_smoother` | Known scalar-observation Gaussian Kalman filtering/RTS smoothing: [y,observed0/1] rows; [k,horizon,confidence,m0,P0,F,Q,H,R] row-major; k<=8,n<=5000,h<=500. | `SELECT chem_kalman_smoother([[1.0,1.0],[2.0,1.0],[1.5,1.0]], [1.0,2.0,0.95,0.0,1.0,1.0,0.2,1.0,0.5]);` |
| `chem_state_space_smoothing` / `state_space_smoothing` | Fit 6 additive innovations ETS families: scalar rows; [family0..5,period,horizon,confidence,maxiter,tol]. Estimated initial states,shared innovations ML and forecast intervals; bounded partial ETS. | `SELECT chem_state_space_smoothing([[1.0],[2.0],[1.0],[2.0],[1.0],[2.0],[1.0],[2.0]], [0.0,0.0,3.0,0.95,1000.0,1e-08]);` |
| `chem_state_space_predict` / `state_space_predict` | Reuse tagged Gaussian or additive-innovations ETS prediction_model: one ignored scalar placeholder row; [horizon,confidence,model...],h<=500. Models preserve final state and uncertainty convention. | `SELECT chem_state_space_predict([[0.0]], [2.0,0.95,1.0,0.0,1.0,1.0,0.25,1.0,0.2,1.0,0.5]);` |
| `chem_arima_forecast_intervals` / `arima_forecast_intervals` | Conditional normal ARIMA forecast intervals from supplied means and A/C/difference polynomials: [confidence,sigma2,degreeA,degreeC,nD,A...,C...,difference_lags...],h<=1000. | `SELECT chem_arima_forecast_intervals([[10.0],[10.0],[10.0]], [0.95,2.0,0.0,0.0,1.0,1.0,1.0,1.0]);` |
| `chem_cross_correlation` / `cross_correlation` | CCF for complete input/output rows: [maxlag<=500,normalization0global/1overlapPearson]; positive lag means input leads output; +/-2 white-noise SE bands. | `SELECT chem_cross_correlation([[1.0,0.0],[0.0,1.0],[-1.0,0.0],[0.0,-1.0]], [2.0,0.0]);` |
| `chem_arima_prewhitening` / `arima_prewhitening` | Input-only exact stationary ARMA prewhitening: [X,Y] rows; [n_train,p<=2,q<=2,maxlag,maxiter,tol]. Apply the same innovations operator to both series; suffix ignored. | `SELECT chem_arima_prewhitening([[-1.0,0.0],[1.0,1.0],[-1.0,2.0],[1.0,3.0],[-1.0,4.0],[1.0,5.0],[-1.0,6.0],[1.0,7.0],[-1.0,8.0],[1.0,9.0],[-1.0,10.0],[1.0,11.0]], [12.0,0.0,0.0,2.0,100.0,1e-08]);` |
| `chem_transfer_function_model` / `transfer_function_model` | Single-input rational omega/delta transfer with exact stationary ARMA noise: [X,Y] rows inclfutureX; [n_train,delay,numdegree<=2,dendegree<=2,p<=2,q<=2,h,confidence,maxiter,tol]. | `SELECT chem_transfer_function_model([[-1.0,1.7],[1.0,2.7],[-1.0,1.3],[1.0,2.3],[-1.0,1.7],[1.0,2.7],[-1.0,1.3],[1.0,2.3],[-1.0,1.7],[1.0,2.7],[-1.0,1.3],[1.0,2.3],[-1.0,1.7],[1.0,2.7],[-1.0,1.3],[1.0,2.3],[-1.0,1.7],[1.0,2.7],[-1.0,1.3],[1.0,2.3],[2.0,0.0],[-2.0,0.0]], [20.0,0.0,0.0,0.0,0.0,0.0,2.0,0.95,100.0,1e-08]);` |
| `chem_remove_linear_trend` / `remove_linear_trend` | Remove least-squares intercept+slope from equally spaced scalar series: [n_train=all],training-only stable centering/scaling and extrapolated trend. | `SELECT chem_remove_linear_trend([[4.0],[4.0],[6.0],[10.0]], []);` |
| `chem_remove_cycle` / `remove_cycle` | Remove one fixed-period cosine/sine cycle: scalar series; [period>2,include_constant=1,n_train=all]. Return amplitude/phase; phase null for numerical zero cycle. | `SELECT chem_remove_cycle([[5.0],[6.0],[-1.0],[-2.0],[5.0],[6.0],[-1.0],[-2.0]], [4.0,1.0]);` |
| `chem_functional_grid` / `functional_grid` | Per-function irregular common-grid interpolation: [x,y,id,optional role0train/1validation]; [],[size],or[size,lower,upper,outside_policy]. | `SELECT chem_functional_grid([[4.0,8.0,1.0],[0.0,0.0,1.0],[2.0,4.0,1.0],[1.0,3.0,2.0],[2.0,6.0,2.0],[3.0,9.0,2.0],[4.0,12.0,2.0]], [4.0]);` |
| `chem_functional_bin` / `functional_bin` | Explicit equal-width function bins: [x,y,id,optional role]; [number,lower,upper,outside0reject/1drop=1]. | `SELECT chem_functional_bin([[0.0,1.0,1.0],[0.25,3.0,1.0],[1.0,5.0,1.0],[2.0,9.0,1.0],[0.0,2.0,2.0],[2.0,10.0,2.0]], [3.0,0.0,2.0]);` |
| `chem_functional_thin` / `functional_thin` | Sorted per-function stride selection: [x,y,id,optional role]; [stride=2,offset=0,keep_last=0]. | `SELECT chem_functional_thin([[5.0,5.0,7.0],[4.0,4.0,7.0],[3.0,3.0,7.0],[2.0,2.0,7.0],[1.0,1.0,7.0],[0.0,0.0,7.0]], [2.0,1.0,0.0]);` |
| `chem_functional_snip` / `functional_snip` | SNIP peak clipping on uniform grids: [x,y,id,optional role]; [max_halfwindow,transform0raw/1Ryan-loglog=1,direction0increasing/1decreasing=0]. | `SELECT chem_functional_snip([[0.0,1.0,7.0],[1.0,1.0,7.0],[2.0,9.0,7.0],[3.0,1.0,7.0],[4.0,1.0,7.0]], [1.0,1.0,0.0]);` |
| `chem_functional_arwls` / `functional_arwls` | Baek2015 arPLS/arWLS baseline: [x,y,id,optional role]; [lambda=1e5,weight_tolerance=1e-6,maxiter=100,negative_std_ddof=0]. | `SELECT chem_functional_arwls([[0.0,0.0,7.0],[1.0,1.0,7.0],[2.0,0.0,7.0]], [1.0,1e-12,1.0]);` |
| `chem_functional_baseline` / `functional_baseline` | Subtract an explicit baseline function: [x,y,id,optional role]; [baseline_id,outside0null/1constant/2reject=2]. | `SELECT chem_functional_baseline([[0.0,1.0,9.0,0.0],[2.0,3.0,9.0,0.0],[0.0,4.0,1.0,0.0],[1.0,7.0,1.0,0.0],[2.0,9.0,1.0,0.0],[0.0,8.0,2.0,1.0],[2.0,12.0,2.0,1.0]], [9.0,2.0]);` |
| `chem_functional_dtw_matrix` / `functional_dtw_matrix` | All-pairs function DTW matrix: [x,y,id,optional role]; [window=-1,normalize0none/1population=0,cost0squared/1absolute=0]. | `SELECT chem_functional_dtw_matrix([[0.0,1.0,1.0],[1.0,2.0,1.0],[2.0,3.0,1.0],[0.0,1.0,2.0],[1.0,1.0,2.0],[2.0,2.0,2.0],[3.0,3.0,2.0],[0.0,2.0,3.0],[1.0,2.0,3.0],[2.0,2.0,3.0]], [2.0,0.0,0.0]);` |
| `chem_functional_golden_match` / `functional_golden_match` | Supplied golden/reference curve distance ranking: [x,y,id,optional role]; [target_id,metric0RMSE/1MAE/2DTW=0,window=-1,normalize=0]. | `SELECT chem_functional_golden_match([[0.0,1.0,9.0],[1.0,2.0,9.0],[2.0,3.0,9.0],[0.0,2.0,1.0],[1.0,3.0,1.0],[2.0,4.0,1.0],[0.0,1.0,2.0],[1.0,2.0,2.0],[2.0,3.0,2.0]], [9.0,0.0,-1.0,0.0]);` |
| `chem_text_dtm` / `text_dtm` | Aggregate supplied numeric tokens into a document-term count matrix; rows[document_id0based,term_id0based,count]; args[n_documents,n_terms]. Caller tokenizes and defines immutable vocabulary. | `SELECT chem_text_dtm([[0,2,1],[0,2,2],[2,0,4]], [4,3]);` |
| `chem_text_weight` / `text_weight` | Fit corpus dictionary/filter and exact JMP five DTM weightings; nonnegative integer count matrix; args[mode0binary/1ternary/2frequency/3log10/4TFIDF=4,min_corpus_count=1,max_terms=all]. | `SELECT chem_text_weight([[2,1,0],[0,1,3],[0,0,0]], [4,1,3]);` |
| `chem_text_weight_predict` / `text_weight_predict` | Apply saved dictionary and training IDF to new count rows; args unchanged text_weight.weight_parameters; training vocabulary width/order required. | `SELECT chem_text_weight_predict([[0,4,2]], [1,3,3,4,3,2,0,1,0.47712125471966244,0.47712125471966244,0.17609125905568124]);` |
| `chem_text_lsa` / `text_lsa` | Partial SVD of supplied nonnegative weighted DTM; args[k,preprocess0uncentered/1centered/2centered_scaled=2]; returns U,S,V and scores US with saved transform. | `SELECT chem_text_lsa([[3,0],[0,2],[0,0]], [2,0]);` |
| `chem_text_projection` / `text_projection` | Score supplied weighted query DTM with saved LSA/topic projection; args unchanged projection_parameters. Training mean,scale,divisor and loadings reused. | `SELECT chem_text_projection([[1.5,1.5]], [1,2,2,3,0,0,1,1,1,0,0,1]);` |
| `chem_text_topic_svd` / `text_topic_svd` | True orthogonal varimax of supplied-DTM SVD loadings VS; document topic vectors UR preserve reconstruction; args[k,preprocess=2,max_sweeps=200,tol=1e-8,kaiser0/1=1]. | `SELECT chem_text_topic_svd([[3,0],[0,2],[0,0]], [2,0,200,1e-10,1]);` |
| `chem_text_latent_class` / `text_latent_class` | Conditional-independence Bernoulli mixture EM on supplied count presence; args[k,seed=0,max_iter=200,tol=1e-8,probability_floor=1e-9]. Reports genuine likelihood,mixture/term probabilities and saved posterior model. | `SELECT chem_text_latent_class([[2,1,0],[0,1,3],[0,0,0]], [1,0,200,1e-10,1e-08]);` |
| `chem_text_latent_class_predict` / `text_latent_class_predict` | Score saved Bernoulli mixture using presence and absence of each fitted term; args unchanged text_latent_class.predict_parameters; returns posterior membership and log probability. | `SELECT chem_text_latent_class_predict([[1,0],[0,1],[0,0]], [1,2,2,1e-08,0.5,0.5,0.9,0.1,0.1,0.9]);` |
| `chem_spatial_covariance` / `spatial_covariance` | Four isotropic JMP spatial kernels: locations in1..3 dimensions; [kind0spherical/1exponential/2Gaussian/3power,partial_sill,nugget,rho]. Independent record nugget; covariance and semivariance matrices. | `SELECT chem_spatial_covariance([[0],[1],[0]], [0,4,1,2]);` |
| `chem_spatial_variogram` / `spatial_variogram` | Empirical isotropic Matheron semivariogram: rows[coordinates1..3,response]; [bins=10,lower=min_pair_distance,upper=max_pair_distance]. Pair counts,mean distance and half mean squared differences; empty classes null. | `SELECT chem_spatial_variogram([[0,0],[1,2],[2,4]], [2]);` |
| `chem_spatial_gls` / `spatial_gls` | Gaussian GLS with supplied absolute isotropic spatial covariance; rows[coordinates,fixed features,response]; [dimensions1..3,features0..12,kind,sill,nugget,rho]. Fixed effects,covariance,likelihood and reusable universal-kriging parameters; no covariance estimation. | `SELECT chem_spatial_gls([[0,1],[1,2],[2,6]], [1,0,0,2,1,1]);` |
| `chem_spatial_predict` / `spatial_predict` | Saved spatial GLS universal-kriging prediction: rows[coordinates,fixed features]; unchanged spatial_gls.predict_parameters. Mean,latent and independent-observation variances include fixed-effect uncertainty given supplied covariance. | `SELECT chem_spatial_predict([[0],[10]], [1,3,1,0,0,2,1,1,3,0,1,2,1,1,1,-2,-1.0000000000000002,3,1,0,1,0,0,1,0.33333333333333337]);` |

<!-- analytics-registry-table:end -->

<!-- analytics-continuation-contracts:start -->

### Continuation contracts and exact method choices

The verified continuation batch adds **62 operations across eight modules**,
bringing the catalog from 119 to **181 analytics functions** (332 total extension
overloads). The matching loaded-extension gate passed on macOS arm64 CPU.
Every registration has its own call in the table above. Numeric outputs and
saved scoring support these workflows; the implemented contracts still cover
subsets of their corresponding JMP method families.

**Latent methods and multivariate inference.** `chem_mds` uses classical
Gower double-centered squared distances. JMP's optimized Stress1 and nonmetric
methods are separate unfinished estimators. Simple correspondence analysis uses
contingency-table SVD with principal or symmetric scaling. MCA uses complete
integer-coded categories, indicator/Burt conventions, and separately reports
Benzécri and Greenacre adjustments. Principal-axis and maximum-likelihood factor
analysis are separate operations; each supports no rotation or Kaiser-normalized
varimax. ML factoring reports variance-boundary, convergence and fit diagnostics;
oblique rotation, confirmatory SEM and missing-data ML remain unfinished.
Canonical correlation whitens two numeric covariance blocks; sequential inference
consumes unregularized fitted correlations and uses Bartlett's asymptotic
chi-square approximation. MANOVA currently tests the whole centered numeric
predictor design against an intercept model, with four criteria and explicit
approximate F conventions; arbitrary effects/contrasts and repeated responses are
not supplied.

Gaussian mixture EM fits full covariances under a minimum standardized covariance
eigenvalue, with seeded starts, component ordering, likelihood history and portable
probability/density scoring. The constraint is an eigenvalue floor, not additive
ridge. Component-count selection and all covariance/noise options remain open.
SOM uses a rectangular PCA-initialized batch grid and fixed-bandwidth Gaussian
local linear smoothing; counted singular local fits fall back to weighted means.
Other topologies, schedules and guaranteed identical JMP codebooks remain open.
Stored mixture/SOM parameters preserve training origins, offsets and scales.

**Evaluation and splits.** Metrics accept caller-supplied predictions and
optional nonnegative weights. Confusion matrices put actual classes in rows and
predictions in columns; undefined per-class precision/recall is JSON null. ROC,
PR and lift group equal scores before integrating inclusive descending thresholds.
The initial ROC threshold is null (above every score). PR reports step average
precision and trapezoid area separately. Probability bins include one in the last
bin and preserve empty-bin nulls. Profit decisions use the supplied actual-by-
decision 2x2 profit matrix; ties choose decision zero. Ensemble prediction takes
fixed supplied weights, and model comparison ranks supplied predictions without
refitting or choosing weights. Holdout splits use `train_fraction` and optional
stratification; rolling-origin folds place all training indices before validation.
They generate indices and cannot establish that a caller fitted preprocessing
inside each split. Nested tuning, provenance and all model-selection methods
remain open.

**Nonlinear and smoothing.** `chem_curve_fit` minimizes ordinary response
residual squares by a scaled analytic-Jacobian Levenberg–Marquardt solver.
Parameters are `[family, max_iterations=300, tolerance=1e-9, optional_starts...]`.
When supplied, starts must follow all three leading values. Automatic starts are
heuristics and local optima are possible. Inspect `converged` and
`locally_identified`; covariance/SE are local Jacobian information estimates,
conditional on the chosen family, and are null when unidentified. These are the
complete currently supported family IDs and coefficient orders:

| Family ID | Coefficient order | Formula / domain |
| --- | --- | --- |
| 0 Exponential 2P | `a, b` | `a*exp(b*x)` |
| 1 Exponential 3P | `a, b, c` | `a+b*exp(c*x)` |
| 2 Logistic 4P | `rate, midpoint, lower, upper` | `lower+(upper-lower)/(1+exp(-rate*(x-midpoint)))` |
| 3 Gompertz 4P | `lower, upper, rate, midpoint` | `lower+(upper-lower)*exp(-exp(-rate*(x-midpoint)))` |
| 4 Power 3P | `a, b, c` | `a+b*x^c; x>0` |
| 5 Michaelis-Menten | `a, b` | `a*x/(b+x); x>0,b>0` |
| 6 Gaussian peak 3P | `height, center, width` | `height*exp(-0.5*((x-center)/width)^2); width>0` |
| 7 Logistic 2P | `rate, midpoint` | `1/(1+exp(-rate*(x-midpoint)))` |
| 8 Logistic 3P | `rate, midpoint, upper` | `upper/(1+exp(-rate*(x-midpoint)))` |
| 9 Gompertz 3P | `upper, rate, midpoint` | `upper*exp(-exp(-rate*(x-midpoint)))` |

`chem_curve_predict` takes the saved `[family, coefficients...]`. Arbitrary model
formulas, parameter bounds, weighted/user loss, profile likelihood, group
comparisons and the remaining pharmacokinetic/dissolution presets are unfinished.
B-spline basis/fit support degrees 1–3, full clamped nondecreasing knots with
`degree+1` endpoint repetitions, and at most 64 basis columns. Queries outside the
closed knot domain fail. `chem_spline_fit` minimizes SSE plus lambda times squared
second coefficient differences. This P-spline penalty is not JMP's integrated-
curvature smoothing-spline lambda. It reports `trace(S)` effective df,
`n-2trace(S)+trace(S²)` variance df and conditional sandwich covariance; automatic
knots/lambda are unfinished. LOESS uses nearest-neighbor tricube local polynomials
of degree 0–2, `max(ceil(span*n),degree+2)` neighbors, optional bisquare weights
based on six times median absolute residual, and at most 2,000 training rows.
Zero median residual stops reweighting; rank-deficient neighborhoods fail.
Saved LOESS state retains training rows and final robust weights. It is not an
assertion that span/edge conventions match JMP defaults.

**Functional and spectral numerics.** Spectrum rows use a common, equally
spaced grid unless specified otherwise. SNV centers and scales each row with the
chosen variance denominator. MSC uses the supplied reference or fits the mean
reference once; preserve that reference for future observations. Savitzky–Golay
uses fixed windows and polynomial boundary fits, with derivative and spacing
specified explicitly. Finite differences use declared second-order edge stencils.
Global polynomial detrending removes a fitted baseline on a normalized grid; it
does not implement SNIP or ARWLS. Haar is an orthonormal power-of-two Mallat
transform with explicit forward/inverse levels. Fourier is a direct DFT with
packed real/imaginary coefficients and a true complex inverse. DTW aligns two
possibly unequal scalar series under squared cost and reports a monotone path;
it does not construct a dataset-wide DTW matrix. Common-grid FPCA uses a sample
covariance decomposition under trapezoid L2 weights and returns eigenfunctions,
scores and reconstruction. Smoothed/sparse/irregular FPCA, automatic basis/BIC,
new-curve projection and functional DOE remain unfinished. Peak finding reports
sample locations, plateaus, global prominences and bases for each spectrum;
parametric peak decomposition is a separate unfinished method.

**SVM, quantile and robust regression.** Binary C-SVC minimizes
`0.5*||w||² + C*sum(max(0,1-y*f))` using feasible dual pair-coordinate SMO.
Linear and `exp(-gamma*||x-z||²)` RBF kernels use training-prefix population feature
scaling. The saved model retains support vectors, signed dual coefficients,
bias and preprocessing. Fits report primal/dual objectives, duality gap, box and
equality residuals, and the final recomputed up/low KKT gap; convergence requires
the requested gap, including after the last allowed update. Multiclass, SVR,
probability calibration, sample weights and automatic tuning remain unfinished.

Quantile regression is unpenalized convex check-loss minimization. ADMM requires
both primal and dual residual criteria; their reported units are standardized
response/design units, while objective/residuals are in original response units.
An iteration-limited result remains `converged:false`. Huber IRLS keeps scale
fixed: a supplied positive response-unit scale or an OLS-residual MAD estimate
computed once. Zero estimated MAD requires a supplied scale unless the fit is
numerically exact. Its score-sandwich covariance is fixed-scale HC0, with original
slope SE and standardized coefficient covariance; this is not simultaneous robust
scale or high-breakdown MM estimation. Both linear estimators share offset-safe
`chem_robust_predict` parameters. For an intercept-only model, supply one ignored
scalar dummy per prediction row. Quantile coefficient inference, arbitrary
categorical formulas and JMP robust-inference equivalence remain unfinished.

**Diploid marker numerics.** Marker calls use `0/1/2` with missing sentinel
`-1`. Statistics report allele/genotype counts, heterozygosity, PIC and Pearson
HWE approximation. Imputation supports mean, mode, LD-kNN, HWE/observed/uniform
random calls and a specified value. LD-kNN uses squared dosage correlation,
mean taxicab distance and majority calls, with explicit mode fallback. Mean/mode
models preserve fit values for future rows; other methods do not imply saved
query scoring. Relatedness includes additive VanRaden and dominance kernels.
Distance matrices have explicit Euclidean/Jaccard/Bray–Curtis/Gower/Hamming
numeric conventions. Seeded unlinked simulations and crosses obey diploid
Mendelian models, without phased recombination maps. GBLUP fixes the supplied
residual/genetic variance ratio and fits ridge-equivalent marker effects;
training-prefix marker means and phenotype centering are reused for suffix/saved
prediction. Certain-allele simulation endpoints are checked with the exact
seed `608688947055533`: probability zero always produces zero dosage and
probability one always produces two, including the formerly rounded RNG endpoint. This marker-ridge operation fixes lambda; the separate kinship contract below adds supplied-kernel ML/REML. GWAS testing, linked/polyploid admixture and genomic file parsers remain gaps; the separate omics contract adds independent-locus diploid admixture.

**Reliability continuation.** Exponential, Weibull and lognormal AFT fits
use full right-censored density/survival likelihoods with numeric location
covariates and a common scale. Exponential fixes SEV scale one; Weibull uses
SEV scale `1/shape`; lognormal uses normal log-time sigma. User-selected stress
transforms can be supplied as predictors; automatic stress-law specification and
scale covariates remain unfinished. Fits require positive times, bounded
`log_sigma`, identified observed information, and successful convergence.
Unpenalized fits report observed-information SE, 95% Wald intervals, AIC/BIC and
p-values. With ridge-penalized standardized slopes, curvature SE are descriptive
and calibrated p/CI/AIC/BIC fields are null. `chem_aft_predict` prepends mode zero
for survival at a nonnegative query time or one for a probability-to-time quantile
to the unchanged saved parameters. Nonrepresentable time quantiles are null while
finite log quantiles remain available. Left/interval censoring, delayed entry,
time-dependent/scale effects, frailty and competing risks remain unfinished.

Power-law NHPP fits one time-truncated system on `[0,T]` using
`shape=n/sum(log(T/t_i))` and `log_eta=log(T)-log(n)/shape`; predictions return
expected cumulative events and intensity. Failure-truncated modifications,
multiple systems/exposure, phases/reset and confidence bands remain unfinished.
Reliability blocks accept independent component success probabilities and compute
series, parallel or exact nonidentical-component k-of-n Poisson-binomial success
probabilities. They do not model arbitrary graph dependencies, common-cause
failures, component importance or repair simulation. These additions extend the
previous KM, zero-location lifetime and Cox methods under distinct contracts.

Exact contracts, output layouts, per-operation exclusions and source/test links:

- [Latent variables and clustering](tasks/jmp-latent-api.json): 8 operations.
- [Prediction evaluation and validation splits](tasks/jmp-evaluation-api.json): 12 operations.
- [Nonlinear curves and local smoothing](tasks/jmp-nonlinear-api.json): 7 operations.
- [Functional and spectral numerics](tasks/jmp-spectral-api.json): 10 operations.
- [Multivariate inference](tasks/jmp-latent-inference-api.json): 4 operations.
- [SVM, quantile and Huber regression](tasks/jmp-robust-learning-api.json): 5 operations.
- [Diploid marker analytics and GBLUP](tasks/jmp-genomics-api.json): 9 operations.
- [Reliability extensions](tasks/jmp-reliability-api.json): 7 operations.

<!-- analytics-continuation-contracts:end -->

<!-- analytics-models-contracts:start -->

### Mixed, causal, omics, classical DOE and market contracts

The earlier 28-operation checkpoint brought the previous 181 to **209 analytics functions and 360 total documented overloads**. The function table above contains every individual call. These are bounded method contracts; full JMP parity remains incomplete.

#### Mixed Advanced

Exact API and remaining options: [contract](tasks/jmp-mixed-advanced-api.json).

**matrix.** Finite nonempty rectangular numeric data<=100000rows,1024columns,1000000cells; fitting further capped256rows/16fixednumericfeatures. Exact signed subject IDs<=2^53-1.

**covariance codes.** 0: heterogeneous diagonal [variance_per_level...]; 1: unstructured [lower triangle covariance including diagonal,row-major]; 2: compound symmetry [variance,rho]; 3: heterogeneous compound symmetry [variances...,rho]; 4: AR1 [variance,rho]; 5: Toeplitz [variance,rho_lag1,...rho_lag(m-1)]; 6: heterogeneous Toeplitz [variances...,rho_lag1,...rho_lag(m-1)]; 7: antedependent [variances...,adjacent_rho0,...adjacent_rho(m-2)]; 8: equal-variance antedependent [variance,adjacent_rho0,...adjacent_rho(m-2)]

**time semantics.** Parameters refer to sorted global time levels. AR1 covariance=variance*rho^abs(actual_elapsed_time);negative rho requires integral gaps. Toeplitz requires uniform global grid and uses lag in grid increments;missing visits retain global lag. Antedependence multiplies all intervening global adjacent correlations,including omitted visits. Other heterogeneity parameters are categorical by ordered level.

**SPD.** Variances>0;correlations strictly(-1,1);CS rho>-1/(m-1);all matrices SPD at relative per-diagonal Cholesky pivot1e-12. No pseudoinverse. Full-rank whitened design relative QR threshold1e-11.

**likelihood.** Exact Gaussian conditional/profile likelihood;REML includes original-design determinant,so response and feature rescaling affects restricted likelihood by the mathematically required constants. AIC/BIC counts fixed coefficients plus covariance parameters;REML comparisons require identical fixed designs. Known-covariance method2 AIC/BIC=null.

**uncertainty.** coefficient_covariance/SE conditional on covariance estimate only;covariance parameter SE/CI/tests not supplied. Local profile Hessian is an identification diagnostic only. covariance_parameters_identified=false or boundary causes converged=false. optimizer_converged separately reports solver result.

**covariance boundary.** Search coordinate abs>24.9 (repeated coordinate cap25),relative variance<=1e-12,relative conditional pivot<=1e-8,or random G determinant/diagonal<=1e-8 criteria. These numerical boundary diagnostics do not provide boundary-mixture tests.

**fit limits.** Max30000likelihood evaluations and1000000000nominal blockdense work units;line search<=40trials;max1000iterations. Repeated covariance<=8levels,unstructured<=36covariance parameters. Numerical failure/invalid input/work exhaustion errors explicitly.

**saved model.** [p,response_origin,response_mean_offset,response_population_scale,standardized_intercept,feature_origins(p),feature_mean_offsets(p),feature_population_scales(p),standardized_fixed_slopes(p),n_subjects,(subject_id,normalized_intercept_BLUP,normalized_slope_BLUP)*n_subjects]. Here p includes the auto fixed slope;last fixed feature is the random slope covariate.

**outputs.** coefficients ordered intercept,other fixed slopes,auto slope;standard_errors same order. fitted_values are marginal fixed-effect means. conditional_fitted_values add BLUPs. blups columns are original-unit random intercept and raw-slope coefficient. covariance_parameters uses physical layouts above.

`chem_repeated_covariance`: **data:** 2..64 scalar strictly increasing global time levels. **args:** structure_code integer0..8; complete physical covariance parameters for the selected code **estimator:** Direct documented variance/correlation formulas; Cholesky verifies SPD with pivots relative to each original diagonal. **unsupported:** Covariance parameter estimation (use repeated_model); spatial covariance families; nonuniform Toeplitz grids

`chem_repeated_model`: **data:** 0..16 numeric fixed features followed by response, exact integer subject id, actual time. Intercept is automatic. <=256 rows,2..8 global levels,unique subject/time,>=2 subjects and within-subject replication. **args:** structure_code integer0..8; method0ML/1REML(default)/2known-covariance GLS; max_iterations integer1..1000,default300; tolerance1e-12..1e-3,default1e-8; gradient norm convergence sqrt(tolerance) clamped[1e-7,1e-5]; method2 only: complete physical covariance parameters; method0/1 prohibit extra values **estimator:** Independent-subject Gaussian likelihood. Subject-block Cholesky whitening and Householder QR compute GLS. ML profiles scale with n;REML with n-p and original-design determinant correction. Shape uses SPD Cholesky/PACF/partial-correlation coordinates and finite-difference BFGS. Method2 uses the supplied covariance without estimating a residual scale. Coefficient covariance conditions on estimated covariance parameters; no small-sample correction. Hessian of profile objective diagnoses local identification. Unconverged/boundary/unidentified solutions are explicit flags,never converged=true. **unsupported:** simultaneous random G-side plus structured R-side model; Kenward-Roger/Satterthwaite or covariance-parameter SE/CI/tests; automatic covariance selection; repeated-model portable prediction; likelihood comparison of REML models with different fixed designs

`chem_random_slope_model`: **data:** 0..15 other numeric features followed by numeric slope covariate,response,exact integer subject id. Intercept and slope are included as fixed and correlated random effects. <=256 rows,<=16 total fixed numeric features,>=3 subjects,within-subject replication. **args:** method0ML/1REML(default); max_iterations integer1..1000,default300; tolerance1e-12..1e-3,default1e-8; gradient norm convergence sqrt(tolerance) clamped[1e-7,1e-5] **estimator:** Gaussian marginal covariance sigma2 I + Z G Z^T with positive-semidefinite 2x2 G from a free lower factor. True profile ML/REML,origin-relative numerical standardization,conditional GLS covariance,group BLUPs G Z^T V^-1 residual. Random covariance/BLUPs reported for original raw slope units; portable model retains normalized coordinates. Boundary singular G returns variance_boundary=true and converged=false,even if optimizer converged. **unsupported:** crossed/nested multiple grouping factors; multiple random slope covariates; structured residual covariance jointly with G; GLMM; covariance-parameter SE/CI/tests or small-sample df corrections

`chem_random_slope_predict`: **data:** Same fixed features/slope order as training,followed by exact integer subject id. Finite rectangular matrix<=100000 rows/1000000 cells. **args:** exact predict_parameters from random_slope_model; layout in common_contract **estimator:** Reuses training feature/response references. Known subjects receive stored normalized random intercept+slope BLUPs;unknown subjects receive zero random effects. Returns conditional predictions,marginal predictions,and known_subjects flags. **unsupported:** refitting new subject effects; prediction intervals; parameter-uncertainty propagation

#### Glmm

Exact API and remaining options: [contract](tasks/jmp-glmm-api.json).

**finite.** Finite nonempty rectangular data; generic <=10000rows,128columns,1Mcells; fit additionally<=2000rows,<=16features,2..200groups. Exact signed group IDs abs<=2^53-1. Parameters<=16384finite doubles.

**identification.** Automatic unpenalized intercept; full-rank standardized fixed design; n>features+2; >=2groups with repeat observations. Both Bernoulli classes required; allzero Poisson rejected. A fitted-tolerance-independent bounded LP rejects Bernoulli complete/quasi separation and Poisson zero-count recession directions. Numerical LP reduced-cost/pivot tolerance1e-12,positive summed-improvement threshold1e-10,feasibility1e-8; ambiguous/work/iteration failures reject. Extremely near-degenerate finite-MLE existence remains a numerical assessment,not exact arithmetic. Guard-bound hits cannot report convergence.

**preprocessing.** Fit-only origin+mean-offset/population-scale feature normalization. Saved eta=eta_center+sum(original_slope*((x-origin)-mean_offset)); this avoids large-origin cancellation. Raw intercept is also reported and may lose digits; use saved model for stable prediction. Unrepresentable raw coefficients/intercepts reject explicitly.

**limits.** Nominal n*(features+2)*q*maxiter*4<=200000000; per-group mode root<=100steps with<=60bracket expansions per side, line search<=50trials/step, eigen<=64sweeps. This preflight is not a cumulative CPU instruction cap. Finite-MLE LP<=4096Bland pivots and<=200M nominal tableau-cell updates,otherwise reject. LP constraints use row-normalized fixed design and |direction_j|<=1.

**optimizer.** Constrained local projected BFGS; coefficients standardized bounded[-30,30] as numerical guard; variance>=0 and<=variance_upper. A separate v0 GLM and three positive starts are compared. Coefficient/upper-variance bound prevents overall convergence.

**quadrature.** Fit uses requested odd q3..31; final likelihood,EB moments,saved posterior useq+4<=35. Compare total and each group LL differences; tolerance=max(user_tol,1e-8)*(1+abs(total_LL)),group tolerance=total_tolerance/G. This is order consistency,not a rigorous error bound. q15 can be insufficient,as recorded in evidence.

**convergence.** converged=optimizer_converged && quadrature_converged && !coefficient_guard_bound && !variance_upper_bound. Variance zero itself is a legitimate boundary optimum. Finite maxiter or failure may return convergedfalse.

**uncertainty.** Interior covariance inverts finite-difference joint integrated-likelihood Hessian in fixed coefficients/log1p variance; fixed block transformed to original units. At v0, exact fixed GLM information. Non-SPD/unrepresentable information gives null covariance/SE. Estimates/SE from unconverged fits remain diagnostic. No variance SE/pvalue.

**portable layout.** [1,family0logit/1log,p,G,Q,variance,eta_center,origins(p),mean_offsets(p),slopes_original(p),group_ids(G),posterior_means(G),posterior_modes(G),posterior_variances(G),posterior_nodes(G*Q),posterior_weights(G*Q)]

**portable size.** 7+3p+4G+2GQ<=16384doubles; fit caps p16,G200,Q35 ->14855. Variance-zero saves Q1. Prediction validates uniqueIDs,nonnegative variance/weights,unit weight sums within1e-8,and finite consistent moments; normalizes accepted weights.

**prediction semantics.** conditional=known posterior-mean plug-in/new b0;posterior=known saved-posterior integral/new Gaussian-prior integral;population=Gaussian-prior integral for all. Logistic prior31/35consistency reported separately,not guaranteed for extreme variance.

**JMP differences.** JMP19.2 supports more distributions,summarized trials,probit,random slopes/multiple factors/repeated structures and permits noninteger Poisson responses/unbounded variance options; absent here. Do not claim identical JMP defaults or full GLMM parity.

**finite mle reference inference.** Bernoulli signed-X directions and Poisson directions unchanged on positive counts/nonincreasing on zero counts are standard GLM recession criteria. Their rejection in this random-intercept model follows because they improve conditional likelihood pointwise for every random intercept,thus also the positive Gaussian integral. This is a stated inference from the cited GLM theory,not a claim those papers establish full random-variance identification.

`chem_glmm_binomial`: **data:** 0..16 numeric fixed features then response and exact integer group ID; n>features+2,n<=2000,2..200 groups including >=2 repeated groups; full-rank fixed design. Binary0/1, or Poisson integer0..1e9. All rows fit. **args:** max_iterations integer1..1000,default200; tolerance1e-10..1e-3,default1e-6; odd quadrature_order3..31,default15; final validation uses q+4; variance_upper positive<=1000,default100 **estimator:** Actual Gaussian-random-intercept integrated maximum likelihood. Per-group safeguarded posterior mode/curvature and adaptive Gauss-Hermite quadrature; finite-difference projected BFGS over centered/scaled coefficients and log1p variance. An explicit variance-zero GLM candidate and variance starts0.25,1,4 are compared. Full Poisson factorial constants are included. AIC counts the variance component even at boundary. Fixed covariance is the fixed block of the inverse joint observed Hessian in the interior; exact fixed-only GLM information at variance zero. **unsupported:** multiple/nested/crossed random factors and random slopes; summarized binomial trials,probit,other response distributions; offset/exposure,weights,overdispersion,repeated residual covariance; REML/PQL/pseudo-likelihood/JMP default-estimator equivalence; variance-component Wald tests,likelihood-ratio boundary mixture tests,small-sample df corrections; global optimizer guarantee; local optima and finite quadrature remain diagnostic

`chem_glmm_poisson`: **data:** 0..16 numeric fixed features then response and exact integer group ID; n>features+2,n<=2000,2..200 groups including >=2 repeated groups; full-rank fixed design. Binary0/1, or Poisson integer0..1e9. All rows fit. **args:** max_iterations integer1..1000,default200; tolerance1e-10..1e-3,default1e-6; odd quadrature_order3..31,default15; final validation uses q+4; variance_upper positive<=1000,default100 **estimator:** Actual Gaussian-random-intercept integrated maximum likelihood. Per-group safeguarded posterior mode/curvature and adaptive Gauss-Hermite quadrature; finite-difference projected BFGS over centered/scaled coefficients and log1p variance. An explicit variance-zero GLM candidate and variance starts0.25,1,4 are compared. Full Poisson factorial constants are included. AIC counts the variance component even at boundary. Fixed covariance is the fixed block of the inverse joint observed Hessian in the interior; exact fixed-only GLM information at variance zero. **unsupported:** multiple/nested/crossed random factors and random slopes; summarized binomial trials,probit,other response distributions; offset/exposure,weights,overdispersion,repeated residual covariance; REML/PQL/pseudo-likelihood/JMP default-estimator equivalence; variance-component Wald tests,likelihood-ratio boundary mixture tests,small-sample df corrections; global optimizer guarantee; local optima and finite quadrature remain diagnostic

`chem_glmm_predict`: **data:** Same feature order as training then exact integer group ID; known or new IDs. <=10000 rows,<=1Mcells. **args:** Unchanged version1 predict_parameters from either GLMM fit; full layout in common_contract **estimator:** Known-group conditional prediction uses inverse link at saved posterior mean random effect; new-group conditional uses zero. Posterior prediction integrates saved discrete posterior for known groups and Gaussian prior for new groups. Population prediction integrates the Gaussian prior for every group. Poisson prior mean is exact exp(eta+variance/2); Bernoulli prior uses31/35-point Gaussian quadrature with order-difference diagnostic. Near-unit posterior weights are normalized; Poisson weighted means are accumulated in log space. **unsupported:** new-group posterior estimation; prediction intervals or parameter uncertainty propagation; posterior integration beyond the saved finite quadrature; automatic extrapolation/overflow repair; unrepresentable means return error

#### Causal

Exact API and remaining options: [contract](tasks/jmp-causal-api.json).

**numeric.** Finite rectangular nonempty <=100000rows,<=66columns,<=1000000cells,<=16384finite args. Per-method work and dimensions stricter.

**positivity.** All causal estimation/matching/evaluation input propensities must be strictly in(0,1). Optional clipping truncates only interior scores, reports count, does not cure failed identification and may bias inference.

**causal assumptions.** Consistency and well-defined interventions/no interference, conditional exchangeability given measured pre-treatment covariates, positivity and independent observations. These are not established by a balanced sample or good logistic fit.

**nuisance provenance.** Supply proper out-of-fold or independent nuisance predictions when needed; input cannot verify training provenance. Plain IPW conditional SE omits fitted nuisance uncertainty; AIPW needs orthogonality rate/regularity conditions.

**outcomes.** One actual observed factual outcome per observation; m0/m1 are predictions, never simultaneously observed potential outcomes.

**scope.** No actual JMP executable comparison, continuous treatment, missing outcomes, multiarm, clustered rows or finite-sample causal guarantee.

`chem_propensity_fit`: **data:** [features...,treatment]; prefix training labels0/1, suffix labels ignored; <=32 features, <=10000 training rows; n_train>parameter_count **args:** n_train: required integer3..min(n,10000); ridge: 0..1e12,default0; summed binomial negative log likelihood +ridge/2*squared standardized slopes; max_iterations: integer1..2000,default100; tolerance:1e-12..1e-2,default1e-8; final gradient max <= tolerance*n_train **estimator:** Analytic Newton/IRLS with backtracking, unpenalized intercept, population feature scales fit on training prefix, original-unit display coefficients and stable saved vector. Unpenalized separated designs reject finite-MLE claims. Return actual final gradient/objective/converged flag; no uncertainty or causal identification is inferred. **unsupported:** Multinomial/continuous treatments; Automatic covariate selection/tuning; Unmeasured confounding identification; Outcome model fitting; Automatic cross-fitting/provenance

`chem_propensity_predict`: **data:** [features...] with same ordered predictors; intercept-only model uses ignored scalar dummy **args:** unchanged saved vector[p,standardized_intercept,origins[p],mean_offsets[p],population_scales[p],standardized_slopes[p]] **estimator:** Stable origin-plus-offset feature normalization and logistic scoring. Query rows do not refit centering/scales. Returned probabilities can numerically reach0/1; downstream causal estimators reject endpoints rather than silently asserting positivity. **unsupported:** Model refitting; Probability calibration; Causal assumption verification

`chem_causal_balance`: **data:** >=1 numeric covariate, treatment0/1, propensity strictly(0,1); >=2rows per treatment group **args:** estimand:0ATE(default),1ATT; explicit_clip:0none(default),otherwise0<clip<.5; truncate interior input scores to[clip,1-clip] for weights only **estimator:** ATE weights T/e and(1-T)/(1-e); ATT weights T and(1-T)*e/(1-e). Weighted mean differences use fixed unweighted pooled sample SD=sqrt((sample_var_T+sample_var_C)/2), not a post-weighting denominator. Return unweighted/weighted SMD, weighted means, baseline SD, treated/control ESS=(sumw)^2/sumw^2 and raw input propensity common-support interval. Common-support=null if empirical ranges do not intersect; outside count uses raw scores even when weights are clipped. Undefined infinite SMD is null; exact constant-equal feature SMD0. Balance is a design diagnostic, not proof of exchangeability. **unsupported:** Distributional/full covariate overlap tests; Automatic covariate encoding; JMP all positivity diagnostics; Proof of conditional exchangeability

`chem_propensity_match`: **data:** [factual_outcome,treatment0/1,propensity strictly(0,1)], <=2000 rows, both groups present **args:** caliper: nonnegative absolute distance in selected score units,default1; replacement:1yes(default),0no; metric:0propensity(default),1logit propensity **estimator:** Process treated rows in input order; choose nearest admissible control, exact distance ties lowest original control index. Optional control reuse; unmatched treated rows excluded and reported. Estimate is mean factual treated-minus-matched-control outcome difference over matched treated subset. No generic independent-pair SE is returned; control reuse and estimated matching need valid matching inference. JMP caliper is multiplied by estimated score SD and its matched report fits a mixed model; neither convention is claimed here. **unsupported:** Optimal/global matching; Multiple controls per treated observation; Automatic score-SD caliper; Automatic min/max support trimming; JMP matched-sample GLMM inference; Abadie-Imbens matching variance

`chem_ipw_effect`: **data:** [factual_outcome,treatment0/1,propensity strictly(0,1)], both groups present **args:** estimand:0ATE(default),1ATT; explicit_clip:0none(default),otherwise(0,.5); changes effective weighting and can introduce bias **estimator:** ATE compares separately normalized weighted means with w1=T/e,w0=(1-T)/(1-e); ATT uses w1=T,w0=(1-T)*e/(1-e). IF_i=n*[normalized_w1_i*(Y_i-mu1)-normalized_w0_i*(Y_i-mu0)]. SE is sample SD(IF)/sqrt(n); fixed95%normal interval. This conditions on fixed supplied nuisance scores and does not include fitted-propensity first-order uncertainty. Cross-fitting alone does not make plain IPW orthogonal. Correct nuisance identification, consistency/no interference, conditional exchangeability and positivity are caller assumptions. **unsupported:** Estimated propensity-model uncertainty correction; Clustered/dependent rows; Sampling/external weights; Continuous/multivalued treatments; Automatic cross-fitting

`chem_aipw_effect`: **data:** [factual_outcome,treatment0/1,propensity strictly(0,1),predicted_m0,predicted_m1], both groups present **args:** estimand:0ATE(default),1ATT; explicit_clip:0none(default),otherwise(0,.5); explicit truncation can alter robustness/target **estimator:** ATE score psi=m1-m0+T*(Y-m1)/e-(1-T)*(Y-m0)/(1-e); estimate=meanpsi,IF=psi-estimate. ATT numerator g=T*(Y-m0)-(1-T)*e/(1-e)*(Y-m0); estimate=sumg/n_treated,IF=(g-T*estimate)/treatment_prevalence. ATT does not use m1. SE=sample SD(IF)/sqrt(n),fixed95%normal interval. Nuisance predictions are externally supplied; calibrated inference needs independent observations, valid identification, suitable cross-fitting/Donsker and nuisance-product rate conditions. Input cannot verify these conditions. **unsupported:** Automatic nuisance model fitting/cross-fitting; All JMP ratio-adjustment/inferential conventions; Clustered/dependent rows; Continuous/multivalued treatments; Finite-sample causal guarantee

`chem_uplift_curve`: **data:** [actual observed outcome,treatment0/1,predicted effect ranking score,known/valid assignment propensity strictly(0,1)], both groups present **args:** no parameters **estimator:** Sort supplied scores descending and group exact ties. z=T*Y/e-(1-T)*Y/(1-e). Gain at selection fraction k/n=sum_selected(z)/n. Qini=gain-fraction*finalgain; report trapezoid AUUC and Qini area. Initial threshold null denotes above all scores. Appropriate held-out evaluation needs randomized known assignment or observational exchangeability/valid propensity, plus consistency/positivity. Actual outcomes are factual observations, not individual treatment effects; observed treated outcomes cannot be used as counterfactual truth. This evaluates a supplied ranking and is not the JMP interaction-significance uplift partition learner. **unsupported:** JMP uplift partition fitting; True individual-effect labels; Automatic score training/cross-fitting; Normalized perfect-model Qini; Confidence bands/policy optimization

#### Omics

Exact API and remaining options: [contract](tasks/jmp-omics-api.json).

**inputs.** Finite rectangular nonempty sample rows/feature columns,<=10000rows,<=1024columns,<=1Mcells; args<=16384. SQLNULL invalid; genotype missing=-1 explicitly.

**normalization.** Method0 row sample/population z-score [0,ddof=1] with ddof0/1, abs<=1e150, nonconstant rows. Method1CPM [1]. Method2RPKM [2,lengthbp per feature...], Method3TPM [3,lengthbp...]. Method4CLR [4,explicit nonnegative pseudocount=0]. Method5relative abundance percent [5], Method6ratio [6]. Counts nonnegative, row totals positive; CLR x+pseudocount>0. No library-size offsets or joint batch correction.

**quantiles.** Values abs<=1e150. Fit all supplied sample rows: sort each row, average per rank across rows, apply mean over ranks spanned by an equal input value. Reference sorted vector is portable. Prediction uses only supplied saved reference, no query refit.

**likelihood.** Independent unrelated diploid genotypes, HWE within each population, binomial2 likelihood including heterozygote ln2. Missing calls skipped; at least one observation per row and marker during fit, per row during scoring. Joint ancestry Q and allele F via seeded multiple-start EM; K1 closed-form. Joint fit nonconvex, labels canonicalized lexicographically by F; possible nonunique Q/F and no SE/global guarantee.

**fit args.** [K1..min(n,16),maxiter1..10000=1000,tol1e-12..1e-3=1e-8,seed exact <=2^53-1=0,n_init1..8=3,frequency_floor1e-12...01=1e-8]. Maximum estimated n*p*K*(maxiter+K+1)*n_init work200M, K1 one iteration. Full parameter vector <=16384.

**saved model.** [K,p,maxiter,tolerance,frequency_floor,positive_initial_Q(K),row_major_F(K*p)]. Positive starting weights normalize to1. Scoring estimates only Q, fixed F; per-row iterations/converged and likelihood. No joint prediction batches used for refitting.

**convergence.** Relative change in observed log likelihood <=tol*(1+abs(previousLL)); not a stationarity or global optimality certificate. Honest false on exhausted maxiter. Monotone finite joint likelihood guarded; fixed-frequency query objective concave but flat/nonunique ancestry weights possible.

`chem_assay_normalize`: **unsupported:** TMM/KDMM and other assay normalization estimators; polyploid/phased/linkage-aware admixture; admixture parameter uncertainty and automatic K selection; global optimum guarantee or same optimizer as JMP/ADMIXTURE

`chem_quantile_normalize`: **unsupported:** TMM/KDMM and other assay normalization estimators; polyploid/phased/linkage-aware admixture; admixture parameter uncertainty and automatic K selection; global optimum guarantee or same optimizer as JMP/ADMIXTURE

`chem_quantile_normalize_predict`: **unsupported:** TMM/KDMM and other assay normalization estimators; polyploid/phased/linkage-aware admixture; admixture parameter uncertainty and automatic K selection; global optimum guarantee or same optimizer as JMP/ADMIXTURE

`chem_marker_admixture`: **unsupported:** TMM/KDMM and other assay normalization estimators; polyploid/phased/linkage-aware admixture; admixture parameter uncertainty and automatic K selection; global optimum guarantee or same optimizer as JMP/ADMIXTURE

`chem_admixture_predict`: **unsupported:** TMM/KDMM and other assay normalization estimators; polyploid/phased/linkage-aware admixture; admixture parameter uncertainty and automatic K selection; global optimum guarantee or same optimizer as JMP/ADMIXTURE

#### Design Classical

Exact API and remaining options: [contract](tasks/jmp-design-classical-api.json).

**input.** Finite rectangular nonempty data<=4096rows,64columns,1000000cells;operation-specific stricter bounds. Treatment plans require distinct treatment rows;IDs refer to input row positions.

**output.** Design<=16384runs and1000000setting cells. IDs zero-based except run_order1based. Physical row order is execution order;blocks/whole plots remain contiguous.

**randomization.** SplitMix64(seed) with rejection-sampled unbiased Fisher-Yates;seed exact0..2^53-1,including0. Repeat seed and inputs reproduce output. RNG is deterministic construction,not cryptographic.

**block information.** IncidenceN[treatment,block],pair_concurrence=N N^T (diagonal treatment replication),treatment_information=diag(replication)-N N^T/block_size. Connected complete/BIBD designs have rankv-1;contrast information is adjusted for block effects.

**BIBD limits.** v3..9,k2..4,k<v,lambda1..3,b<=128;backtracking nodes<=100000 andcandidate feasibility checks<=5000000. Errors distinguish unsupported parameter scope,no_solution (violated necessary equations or exhaustive search completed),search_limit (budget exhausted;existence undetermined). Repeated blocks allowed when lambda>1.

**split strip.** Hard-factor treatment combinations havea=2^(role0factor count),easy/column combinationsb=2^(role1factor count). Replicate labels are blocking strata. inference_ready only diagnoses positive classical replication error df;response inference is never performed.

**OA.** L4<=3binary,L8<=7binary,L9<=4ternary,L18<=1binary+7ternary. Only catalog column subsets;coded levels0..L-1. Main effects orthogonal,interactions may be aliased. No robust-design response analysis.

**limitations.** No general custom-optimal design search,no ANOVA/REML fitting in these constructors,no arbitrary mixed-level/IBD/resolvability claims.

`chem_rcbd_design`: **data:** 2..256distinct treatment-setting rows,1..16numeric columns. Each input row is one treatment. **args:** blocks integer1..1024(required); seed exact nonnegative integer<=2^53-1,default0 **estimator:** Each block contains every input treatment once. Fisher-Yates randomizes treatment order within each complete block;blocks remain contiguous. Returns treatment/block/run IDs,replications,pair concurrence and block-adjusted treatment information. **unsupported:** arbitrary optimal blocking; response-model fitting; BIBD incomplete blocks

`chem_balanced_incomplete_block`: **data:** 3..9distinct treatment-setting rows,1..16numeric columns. **args:** block_size integer2..4,strictly less than treatment count; lambda integer1..3,default1; seed exact nonnegative integer<=2^53-1,default0 **estimator:** Exact pair/replication-constrained BIBD. Integer equations determine r=lambda(v-1)/(k-1),b=vr/k. Search may use repeated blocks and fixes first block by treatment-label symmetry;certified binary affine construction handles(v,k,lambda)=(8,4,3),cyclic Z9 difference-family handles(9,4,3). Randomized treatment relabeling,block order and within-block order preserve exact incidence. No solver dependency in Rust. **unsupported:** arbitrary unbalanced incomplete block optimization; v>9,k>4,lambda>3; guaranteed completion outside bounded work budget

`chem_resolvable_incomplete_block`: **data:** q^2distinct treatment-setting rows,prime q in{2,3,5,7};1..16numeric columns. **args:** seed exact nonnegative integer<=2^53-1,default0 **estimator:** Affine-plane lines over prime finite field Fq:vertical plus q slope classes. v=q^2,k=q,r=q+1,b=q(q+1),lambda1. Each parallel class partitions all treatments once. Seed randomizes treatment labels,class/block order and order within block;parallel_class_ids explicitly retain class structure. This is a resolvable BIBD family,not an arbitrary resolvable-design search. **unsupported:** composite finite fields such as q4; arbitrary resolvable BIBD parameters; generic meaning of RIBD as any randomized incomplete plan

`chem_split_plot_design`: **data:** 2..12factor rows[finite low<high,role0whole-plot/1subplot];both roles required. **args:** replicates integer1..256,default1; seed exact nonnegative integer<=2^53-1,default0 **estimator:** Complete two-level factorial treatment structure crossed with replicate blocks. Whole treatment combinations are randomized per replicate and held constant over contiguous whole plots. All subplot combinations randomize independently inside each whole plot. Returns replicate/whole-plot/treatment IDs and aggregate classical replication error df[(r-1)(a-1),a(r-1)(b-1)]. Replicate1produces a plan with zero replication error df and inference_ready=false. No response inference performed. **unsupported:** arbitrary custom D/I-optimal split plots; more than2levels per factor; split-split plots; response-model inference

`chem_strip_plot_design`: **data:** 2..12factor rows[finite low<high,role0row-strip/1column-strip];both roles required. **args:** replicates integer1..256,default1; seed exact nonnegative integer<=2^53-1,default0 **estimator:** Complete crossed two-level row/column treatment combinations in each replicate. Row and column strips randomize independently. Membership crosses rather than nests:every row strip intersects every column strip once. Returned design rows follow stage1order;stage2_order is a zero-based permutation regrouping units by column strip,randomizing row subunits within each column strip. Replication error df[(r-1)(a-1),(r-1)(b-1),(r-1)(a-1)(b-1)]. Replicate1zero df explicitly;no response inference. **unsupported:** arbitrary custom optimal strip plots; nested split-split interpretation; additional third-stage easy factors; response-model inference

`chem_taguchi_orthogonal_array`: **data:** 1..8scalar factor level counts2or3. Output categorical level codes are zero-based. **args:** array required integer4/8/9/18; seed exact nonnegative integer<=2^53-1,default0 **estimator:** L4/L8 binary projective strength2arrays,L9 ternary finite-field array,NIST tabulated mixed L18(2^1 3^7). Requested factors map to first unused compatible catalog columns. Run order randomized. Runtime verifies pair balance. Normalized Helmert contrast main_effect_gram and main-vs-two-factor-interaction alias matrix expose confounding. Interaction terms use[factor,contrast,factor,contrast] zero-based;alias excludes intercept. **unsupported:** full Taguchi robust parameter-design platform; inner/outer array crossing; signal-to-noise response analysis; arbitrary mixed-level orthogonal array search; interaction-free claims

#### Market

Exact API and remaining options: [contract](tasks/jmp-market-api.json).

**inputs.** Finite nonempty rectangular <=100000rows,1024columns,1Mcells. Baskets strictly0/1 with <=32columns; no frequency weights,missing rows,automatic text extraction.

**rules.** Apriori includes all frequent sets up to max_rule_size2..6,max antecedents1..5,min support1e-12..1,confidence0..1,lift>=0. Defaults .1,.4,1.2,3,4. Rule output limit1..100000 default10000; exceeding it errors rather than truncates. <=1Mjoins,100Mtransaction count work,100kfrequent sets,1Mrule candidates; may reject large dense supported input by explicit work limit. No global rule ranking/top-k truncation. Every rule proper disjoint nonempty LHS/RHS. Counts exact; row number denominator includes empty baskets.

**rule metrics.** support=count(X union Y)/n;confidence=countXY/countX;lift=n*countXY/(countX*countY);leverage=support-pX*pY;conviction=(1-pY)/(1-confidence). Confidence1 yields convictionnull with conviction_infinite=true when pY<1; if pY1,indeterminate null/false. All thresholds inclusive.

**reliability.** Complete-case respondents x items,>=2rows,2..64items,abs<=1e150,n*p^2<=100M. Raw alpha=p/(p-1)*(1-sum var(item)/var(total)). Standardized uses correlation-scale items and requires every variance positive. Deleted alpha null when fewer than2items remain or variance zero. Corrected correlation excludes the item from total. Undefined/nonrepresentable JSONnull; negative alpha is retained. No automatic reversal,significance interval,IRT,unidimensionality/validity conclusion.

Independent checks use separate units: mixed covariance/GLS/ML/REML/BLUP 630 cases; GLMM 31 check groups/551 assertions; causal 51,062 scalar assertions; assay/ancestry 556 cases/22,948 fields; classical designs 586 cases; baskets/reliability 340 cases/14,842 fields. GLMM prior predictive quadrature can remain unconverged even when its likelihood fit converged; both diagnostics are retained. Joint ancestry EM is not a global optimum certificate.

<!-- analytics-models-contracts:end -->

<!-- analytics-specialized-contracts:start -->

### Genomic variance, item response, SEM and degradation contracts

The earlier checkpoint added **16 CPU operations**, bringing the then-registered API to **225 analytics functions / 376 total extension overloads**. All individual calls appear in the table above. The following exact numeric and inference contracts remain subsets of JMP19.2 methods.

#### KINSHIP

Complete machine-readable [API](tasks/jmp-kinship-api.json) and [independent reference](tasks/evidence/jmp-kinship-reference.json).

**input.** Finite nonempty rectangular inputs,abs<=1e150,<=1Mcells; ntraining4..128,pfixed0..16,n>p+2. Exact row order retained for K and queries. No observation weights/missing response. Intercept automatic.

**kernel.** Supplied exactly symmetric PSD K; diagonal mean positive and representable; any unequal mirrored entries reject (including roundoff asymmetry), so callers must explicitly symmetrize before fitting. max(|Kij|)/maxdiag<=1e4; eigenvalues below -1e-12*max(1,maxeigen) reject. Roundoff-small negative eigenvalues are clipped to zero; likelihood and training BLUP use that effective PSD kernel, saved normalized weights are restricted to its strictly positive eigenspace. Symmetric original training cross-rows reuse the same fitted BLUP to eigensolver rounding. K normalized by original mean diagonal for scale-invariant h.

**likelihood.** V=sigma_g^2 K+sigma_e^2 I;normalized R(h)=(1-h)I+hK/mean(diagK),sigma_g^2=profile_scale*h/meanDiagK,sigma_e^2=profile_scale*(1-h). ML denominatorn,REMLdenominatorn-p-1. REML includes logdet of original-design X^TR^-1X; response and fixed-feature unit changes carry required likelihood constants.

**optimizer.** Evaluate uniform h grid including0/upper,refine all resolved interior and endpoint-neighbor minima with golden sections,retain exact boundaries. <=30000profileevaluations. No global optimum theorem; grid can miss unresolved narrow interior extrema. Flat observed profile relative1e-9 returns components_identified=false/converged=false. Active upper bound returns converged=false; lower zero-genetic boundary can converge. Conditional fixed effects and BLUP use fitted covariance.

**numerics.** Response/features use origin-relative dimensionless differences; eigenspace weighting and Cholesky of standardized whitened normal Gram with relative diagonal-pivot1e-12. Numerically perfectly explained response relative residual norm<=1e-12 rejects; no positive-scale Gaussian MLE claimed. Component variance products use logarithms. Portable prediction expands each physical term including y_scale into binary mantissa/exponent products and ratios and uses compensated floating expansions before conversion to binary64; finite zero/nonzero effects survive intermediate overflow, underflow, and large cancellation, including subnormal kernel/feature units. Finite precision and eigensolver tolerance still apply. Positive fitted variance underflowing physical binary64 rejects with rescaling guidance.

**saved model.** [version1,n,p,y_origin,y_offset_scaled,y_scale,standardized_intercept,kernel_mean_diagonal,feature_origins(p),feature_mean_offsets_scaled(p),feature_difference_scales(p),standardized_slopes(p),positive_support_normalized_kernel_weights(n)]. Length8+4p+n. Algebraic prediction=y_origin+y_scale*(y_offset+fixed_normalized+(crossK/kernel_mean_diagonal)dotweights), evaluated in physical term scales. CrossK rows must share the training kernel units/order. With clipped near-PSD training K, weights remove clipped/null training directions; query crossK is interpreted against that retained support. No query/query covariance or full joint-PSD check is performed.

**uncertainty.** No variance-component inference or fixed-effect SE. Normalized h is not a general biological heritability estimate unless the caller kernel normalization/model assumptions justify it. No installed-JMP head-to-head comparison.

`chem_kinship_model`: **data:** n training rows[y,fixed numeric covariates(p),K_i0...K_i(n-1)]. n4..128,n>p+2,p0..16,one phenotype and one supplied kinship. **args:** p integer0..16 required; method0ML/1REML default1; h_upper .01..0.999999999999 default.99999999; grid16..512 default128; tolerance1e-12..1e-4 default1e-10 **unsupported:** GWAS tests and marker-specific variance refits; Multiple kernels/multitrait/nested effects/phenotype missing; Unbounded negative variance components; Variance-component SE,CI or boundary tests; Guaranteed global optimum; Asymmetric kernels, including unequal mirrored entries within former roundoff tolerance; explicitly symmetrize before fitting

`chem_kinship_predict`: **data:** Query rows[p original fixed covariates,n cross-kinship values in exact training-row order]. No query outcomes or variance fitting. **args:** unchanged fitted predict_parameters **unsupported:** Prediction intervals and covariance-parameter uncertainty; Refitting new genetic effects from observed query phenotypes; Checking the full joint query/training PSD kernel without query/query covariance; Inferring a corrected joint query/training PSD kernel from query cross-covariances alone

#### IRT

Complete machine-readable [API](tasks/jmp-irt-api.json) and [independent reference](tasks/evidence/jmp-irt-reference.json).

**finite data.** Finite nonempty rectangular numeric matrix<=10000rows,128columns,1000000cells;fit<=2000rows and item counts2..12for1PL,3..12for2PL,5..12for3PL. Parameters<=16384finite doubles.

**frequency.** Optional final frequency column integer0..1000000. Zero-frequency rows are validated then ignored. Identical response patterns collapse. Effective frequency excludes allmissing omitted patterns, must exceed kind*items+1and be<=10000000. No fractional case weights.

**missing.** NaN/SQLNULL are not response missing values; finite-1sentinel requires policy1incorrect or2omit. Default0reject. JMP missing->incorrect is matched only with policy1; omission is an explicit extension.

**assumptions.** One latent standard-normal ability,positive item discrimination,conditional local item independence. Ability location/scale fixed0/1;no learned person distribution. No automatic diagnostics prove substantive unidimensionality/local independence.

**model parameterization.** P_j(theta)=c_j+(1-c_j)*sigmoid(a_j*(theta-b_j));D=1.1PL fixes a1,c0;2PL fixes c0;3PL estimates c in[0,cmax].

**bounds.** Difficulty[-8,8],a[.05,5]numerical guards;contact prevents combined convergence. Guessing bound0..cmax(.01..0.8)defines constrained3PL model. Active zero/upper guessing coefficients have null SE/covariance rows;other uncertainty conditions on that face.

**penalty.** Explicit lambda/2*sum(c_j^2)on TOTAL log likelihood,not average. Returned penalized_log_likelihood=LL-penalty. It is not asserted to reproduce JMP proprietary guessing-penalty formula. Scaling all frequencies and lambda equally preserves estimates.

**identification.** Each item needs both observed classes.2PLrequires>=3items,3PL>=5;effectiveN>parameter_count+1. These are preconditions,not generic/global-identification proofs. Positive definite observed free information required for combined convergence;negative/nearzero relative eigenvalues<=1e-10yield information_identifiedfalse.

**quadrature.** Fit optimized at requested oddq11..81;reportedLL,posterior pattern moments,andinformation useq+8<=89. optimizer_gradient_quadrature_order names lowq;reported_likelihood_quadrature_order andinformation_quadrature_order name highq. ChecktotalLLdelta<=max(tol,1e-9)*(1+absLL)and eachpatternlogprobdelta<=max(tol,1e-9)*(1+abslogprob). Consistency is not rigorous error bound.

**convergence.** Combined converged=optimizer_converged&&quadrature_converged&&!numerical_guard_bound&&information_identified. Optimizer uses average penalized NLL gradient tolerance. Local multistart search,not global proof. Remaining finite iterates are diagnostic.

**inference.** parameter_covariance/standard_errors are null unless combined convergence;AIC/BIC also requirelambda0. information_identified/minimum_information_eigenvalue may still diagnose local curvature of a failed iterate. inference_valid reports validated converged curvature availability. Lambda>0uncertainty_kind explicitly denotes penalized-objective curvature,not MLE sampling covariance or a sandwich estimator. No p-values/CIs or boundary-mixture tests.

**covariance order.** Difficulties,then original-unit discriminations if2/3PL,then guessing if3PL. Hessian differentiates integrated scores in b,loga,c;inverse free block transformed by da/dloga=a. Active guessing rows/columns null;free stencil shrinks to25%of distance tobounds and rejects ifunresolvable.

**portable layout.** [1,kind1/2/3,items,missing_policy,saved_odd_q, a0,b0,c0,a1,b1,c1,...]. Savedq=fitq+8<=89;length5+3J<=41. Portable manual fixed-item models allow1..12items;fit still usesminimum2/3/5.

**score.** EAP under fixed N(0,1),not JMPpersonMLE. Score compares savedq andq+8<=97;maximum of logprob/mean/variance order differences must be<=1e-7*(1+abslogprob). Completely omitted response patterns exactlyprior. Conditions on supplied item parameters.

**sharp tail limit.** Measured12items,a5,b8,allcorrect,q81/89 returns mean8.39323 versus independentquad8.41561 and variance.04058 versus.03433;quadrature_convergedfalse withdelta.07225. Such returned values are diagnostic and must not be treated as converged precision. No adaptive integration fallback.

**resource.** Nominal collapsed_patterns*items*q*starts*maxiter<=500000000;maxiter<=2000,1or3starts,line search<=50trials,Jacobi<=64sweeps with separate boundedworkcheck. Nominal preflight is not cumulative CPU instruction count;observed-information evaluations are additional bounded work. Output matrices fit<=12items and<=10000curverows.

**JMP gaps.** Bounded local optimizer,explicit ridge penalty andposterior EAP scoring differ from JMPdefaults/abilityformula. NoGUI,model-selection/reportplots,multidimensional/polytomous response,latent regression,multiple-groupDIF,equating/adaptive testing,fullfitindices orWald/LR inference.

`chem_irt_fit`: **data:** Binary item columns, optional -1missing sentinel with explicit policy and optional final integer frequency. <=2000 rows;1PL2..12items,2PL3..12,3PL5..12. Each item needs correct/incorrect effective observations. No person/item covariates. **args:** model1Rasch/2logistic/3guessing;default2; frequency_column0no/1last;default0; missing0reject/1incorrect/2omit;default0;only finite-1sentinel; maxiter integer1..2000;default300; projected average-negative-objective gradient tolerance1e-10..1e-3;default1e-6; odd normal-Hermite quadrature q11..81;default41,validationq+8; guessing upper constraint cmax.01..0.8;default.5; total-loglikelihood penalty lambda/2*sum(c_j^2),lambda0..1e6;default0;3PLonly **estimator:** Fixed standard normal ability and independent Bernoulli items conditional on ability. Normal Gauss-Hermite MML with analytic integrated scores and bounded projected BFGS;1start1PL,3starts2/3PL. Discrimination exp(loga)>0,1PL a1/c0;2PL c0. Difficulty[-8,8],discrimination[.05,5]are numerical guard bounds and prevent combined convergence on contact. Guessing[0,cmax]is the explicitly constrained3PL model;active guessing faces may legitimately converge. Local optimum only. Frequencies collapse identical patterns; allcorrect/allincorrect people contribute to calibration. **unsupported:** unconstrained a,b or guessing1 bounds; 1/2PL guessing penalty,3PL with<5items; multidimensional or polytomous IRT,latent regression,group DIF/equating; automatic item reversal/selection/local-independence tests; conditional/joint ML and exact JMP optimizer or guessing-penalty equivalence; person ML ability estimates or full item parameter CI/Wald/LR inference

`chem_irt_predict`: **data:** One finite numeric ability per row,<=10000rows. **args:** unchanged model_parameters returned by irt_fit; portable layout in common_contract **estimator:** P=c+(1-c)*sigmoid(a*(theta-b)),D=1. Conditional Bernoulli Fisher information=(dP/dtheta)^2/(P(1-P)),summed across items. Stable boundary formulas give zero information at saturated tails. Also returns expected total score and conditional total-score variance. **unsupported:** ability uncertainty propagation or parameter uncertainty bands; adaptive testing/item selection; multidimensional curves

`chem_irt_score`: **data:** Same ordered item columns as fitted model,without frequency;<=10000rows. Portable fixed-item models may contain1..12items. **args:** unchanged model_parameters returned by irt_fit **estimator:** Posterior mean/variance under fixed N(0,1) prior,conditional on supplied item parameters. Uses savedq andq+8normal Hermite rules,returns higher-order moments/logpatternprobability. Both allcorrect/allincorrect patterns have proper finite EAP;allmissing under omission returns exactly mean0/variance1/logprobability0. This differs from JMP person-MLE ability scoring,which omits extreme patterns. Item parameters remain fixed. **unsupported:** JMP person-MLE/Newton ability formula; posterior intervals/plausible values or item-parameter uncertainty; refitting or scoring additional traits; automatic adaptive quadrature; sharp tails may remain unresolved and flagged

#### SEM

Complete machine-readable [API](tasks/jmp-sem-api.json) and [independent reference](tasks/evidence/jmp-sem-reference.json).

**input.** Finite nonempty rectangular numeric matrix<=100000rows,256columns,1000000cells;operation-specific tighter shape limits;<=4096finite args.

**covariance.** ML covariance uses divisorN. N-1normalization conversion is explicitly caller responsibility. Covariances must be symmetric and SPD at relative Cholesky pivot1e-12.

**CFA identification.** One active loading per indicator;>=3indicators/factor;latentvariance1;first indicator loading positive. Factor covariance is SPD and unit diagonal. Local observed-information SPD diagnoses identification;no general algebraic identification or global optimum certificate.

**CFA parameter order.** standard_errors/parameter_covariance order:active loading per indicator(p),uniqueness per indicator(p),free factor correlations lower triangle row-major if correlated. Manifest means are free parameters and count toward k;mean_covariance=implied covariance/N is reported separately. CFA k=3p+q(q-1)/2 when correlated,3p otherwise.

**CFA inference.** Optimizer convergence,identified,Heywood uniqueness boundary,latent boundary separate. SE/parameter covariance,overall fit p,RMSEA interval and information criteria null unless converged identified interior. Fitted means/covariance and discrepancy/indices remain descriptive candidate values. No parameterwise Wald p/CI.

**CFA limits.** <=20manifest variables,<=4factors,max1000iterations,100000likelihood/gradient evaluations or1e9nominal p^3work units. Line search<=50trials. Log anchor/uniqueness coordinates[-25,25],active normalized loadingsabs<=100,factor-correlation free coordinatesabs<=1000. Solver gradient threshold=sqrt(tolerance)clamped[1e-7,1e-5]. Heywood uniqueness<=1e-8 in standardized units;latent anchor<=1e-8,conditional factor covariance pivot<=1e-8,or coordinateabs>999boundary.

**saved model.** score_parameters=[p,q,means(p),Lambda row-major(pq),Phi row-major(q²),uniquenesses(p)]. Prefix method0/1when calling cfa_scores.

**fit indices.** df=p(p+3)/2-k includes means. Independence baseline mean=observedmean,cov=diag(S),df0=p(p-1)/2. Baseline discrepancy within64eps*p is treated as numerical zero. Undefined zero-df or baseline denominators returnnull. AICcnull whenN<=k+1.

**RMSEA convention.** Uses displayed official equation sqrt(max(chi2-df,0)/(N*df)),without arbitrary upper clipping. Help prose saysbounded0..1while its displayed formula lacks upper clipping;the unbounded formula choice is an explicit implementation interpretation,not a verified installed-JMP behavioral claim. CFI/TLI/NFI are clipped[0,1]as prose specifies;RNI is raw/unbounded.

**RMSEA CI.** Poisson mixture of central gamma CDFs,built around Poisson mode;normalizes retained weight mass. Max noncentrality1e7,max50000mixture offsets,gamma convergence<=200000iterations,bracket/bisection75steps. CentralCDF<=target gives zero noncentrality boundary;zero modeldf givesnull. Cap/no convergence givesnull interval plus explicit status.

**RMR SRMR.** Unique lower triangle covariance residuals plus all mean residuals,denominatorp(p+3)/2. SRMR covariance residuals divided bysqrt(Sii*Sjj),mean residuals bysqrt(Sii). No omission of nonzero mean terms.

**RAM.** p<=16;linear-system inverses use elimination with pivot rejection. Static invertibility and discrete dynamic stability distinct. Stability solvesP-BPB^T=I by p^2-dimensional linear system;SPD P=>strictly stable,nonSPD=>not strictly stable,singular/numeric unresolved=>null. This extra algebraic diagnostic is not an SEM identification/causality claim.

**unsupported.** No general SEM fitting,cross-loadings,residual correlations,FIML,multigroup/measurement invariance,robust correction,modification indices,or GUI.

`chem_cfa_model`: **data:** p3..20summary rows,each p ML covariance entries followed by sample mean. N-divisor covariance is required;N-1conversion is caller responsibility. **args:** N exact integer>p,<=1e9; q integer1..4; correlated flag0/1(required); max_iterations integer1..1000(required); tolerance1e-12..1e-3(required); p indicator assignments integer0..q-1;every factor>=3indicators **estimator:** True Gaussian summary-data ML with freely estimated manifest means. Each indicator loads on its specified one factor;latent variances fixed1 and first indicator loading positive to identify scale/sign. Positive uniquenesses and SPD correlated factor matrix. Analytic Gaussian score with correlation-coordinate chain rule;BFGS local search and Richardson derivatives of analytic score provide observed information. Physical-parameter delta covariance accounts for sign/scale constraints. Heywood/latent boundary,optimizer convergence and local identification are separate. No global optimum guarantee. **unsupported:** cross-loadings; correlated indicator residuals; general latent structural paths; FIML/missing data; multigroup/invariance; sandwich robust correction; parameterwise Wald tests/confidence intervals; globally certified ML optimum

`chem_cfa_scores`: **data:** Finite raw observed-variable rows,n×p;order matches stored means/loadings. <=100000rows/1000000cells. **args:** method0regression posterior mean/1Bartlett weighted least squares; exact returned score_parameters=[p,q,means(p),Lambda row-major(pq),Phi row-major(q²),uniquenesses(p)] **estimator:** Regression scores Phi Lambda^T Sigma^-1(x-mu),with conditional latent posterior covariance. Bartlett scores solve Lambda^T Theta^-1 Lambda and return measurement-error covariance. Reuses supplied saved parameters;no refit or propagation of fitted parameter uncertainty. **unsupported:** factor-score parameter uncertainty; missing observations; posterior sampling; general structural latent variable prediction

`chem_sem_fit_indices`: **data:** p1..20rows,each2p+2columns:observed ML covariance,model implied covariance,observed mean,model mean. Both covariances SPD. **args:** N exact integer>p,<=1e9; k_free exact integer0..p(p+3)/2,includes all free covariance and mean parameters **estimator:** Evaluates supplied Gaussian mean/covariance model. Chi-square=N*[logdetSigma-logdetS+tr(Sigma^-1S)+(mean difference)^T Sigma^-1(mean difference)-p]. Degrees of freedom=p(p+3)/2-k. Independence baseline has observed means and diagonal S. Reports central chi-square p,CFI/RNI/TLI/NFI,RMSEA with90%noncentral chi-square interval,RMR/SRMR with explicit mean residuals,Gaussian loglikelihood,AIC/AICc/BIC/BICu. Caller supplies valid model parameter count and is responsible for ML estimation/model assumptions. **unsupported:** automatic model identification check; checking supplied model was estimated by ML; robust/scaled chi-square corrections; modification indices; noncentrality>1e7confidence bounds

`chem_sem_ram`: **data:** p1..16rows,2p+1columns:B path coefficient row,SPD disturbance covariance Omega row,intercept. B diagonal must be zero. **args:**  **estimator:** Solves simultaneous linear system with A=I-B:mean=A^-1alpha,covariance=A^-1Omega A^-T,total=A^-1-I,indirect=total-B. Static invertibility is separate from strict discrete dynamic stability;the latter uses a discrete Lyapunov solution SPD diagnostic and may be null at a singular/numerically unresolved boundary. Positive Omega/invertibility do not establish model identification,fit or causality. **unsupported:** SEM parameter estimation; causal identification guarantee; dynamic stability from invertibility alone; nonsingular indefinite disturbance covariance; more than16variables

#### DEGRADATION

Complete machine-readable [API](tasks/jmp-degradation-api.json) and [independent reference](tasks/evidence/jmp-degradation-reference.json).

**data.** Finite nonempty rectangular <=10000rows,<=128columns,<=1000000cells; stricter per-op caps apply.

**args.** <=16384 finite doubles; IDs exact integers within absolute2^53-1. All fit rows are used, no implicit validation split.

**inference.** Normal independent errors or independent destructive/event units as specified. Serial dependence within a unit is not adjusted. Uncertainty from model choice is not included.

**path work.** n*design_rank^2<=200000000; <=128 unit IDs, QR rank threshold1e-11*sqrt(n) on column-scaled design.

**destructive limits.** <=4000rows;<=4 location coefficients, maxiter<=1000; standardized coefficient guards +/-1e6 and logσ[-30,30]. Actual likelihood counter capped at200000000 n*location_dimension units; exhausted budget rejects.

**cumulative limits.** <=4000rows;<=128patterns,<=128segments each,<=1024segments total;200000000 nominal likelihood units and per-evaluation work counter; guards centeredβ +/-1e4 and logσ[-10,10].

**numerical scope.** Nonconvergence flags are diagnostics; three-start BFGS is not a proof of global maximum. Null indicates documented nonrepresentability/unboundedness/unidentified curvature, never JSON NaN/Infinity.

**portable path df.** Saved path residual df must be integer1..10000, matching bounded fit resource domain; Student critical values solve the small survival tail directly.

**inverse numeric.** Cumulative time quantiles search log time, retaining relative accuracy far below a segment duration; stopped-pattern terminal rounding is clamped internally, query times remain strict.

**ML inference validity.** For both ML fits, optimizer_converged reports the projected-gradient criterion; parameter_guard_bound reports any active numerical coefficient/scale guard. information_identified reports a positive-definite local observed information matrix at the returned iterate, independently of convergence. inference_valid (also converged) requires optimizer_converged, no parameter guard, and identified local information. Covariance and AIC are null otherwise. Local curvature alone never validates inference from an unconverged estimate, and no global maximum proof is claimed. Coefficients, attained log likelihood, convergence diagnostics and saved predictions remain available.

`chem_degradation_path`: **data:** [nonnegative time,response,exact integer unit ID]; <=128 IDs; positive time for log time; positive response for log response; n>design rank, rank identified. **args:** intercept_mode:0=separate per ID,1=common at transformed time zero,2=fixed zero; default0; slope_mode:0=separate per ID,1=common; default0; response_transform:0=identity,1=natural log; default0; time_transform:0=identity,1=natural log,2=square root; default0 **estimator:** Twice reorthogonalized QR OLS for g(Y)=a_ID+b_ID*f(time) under the selected equality constraints. Pooled SSE/(n-rank) residual variance and marginal coefficient covariance. Common intercept/separate slopes retains a common intercept at f(time)=0. All rows fit; no train/test prefix inference. Transformed normal independent homoskedastic errors are required for inference. Correlation among repeated observations is not adjusted. **unsupported:** Bayesian hierarchical repeated-measures estimation; nonlinear paths; correlated residuals/random effects; arbitrary/custom transforms; X-level-within-group equality constraints

`chem_degradation_predict`: **data:** [time,known saved unit_id,response-scale threshold]; same transform domains as fitted model; no response or refitting. **args:** confidence:required in(0.5,1); side:required0=two-sided mean band,1=lower mean bound,2=upper mean bound; saved model unchanged: [1,response_transform,time_transform,G,(ID,a,b,covaa,covab,covbb,MSE,df)*G] **estimator:** Student-t mean and individual-response prediction bounds on the transformed scale, back-transformed to the original response scale. Threshold point solves a+b*u=g(threshold). Inverse confidence set contains u for which the threshold belongs to the pointwise mean confidence interval/bound. Stable quadratic/sign interval partition retains empty, singleton, disjoint and unbounded sets, intersected with physical time domain [0,infinity) or (0,infinity) for log time. Interval endpoints are mapped back to time; null upper with upper_unbounded=true denotes positive infinity. Nonrepresentable finite endpoints/point crossings may be null. These sets invert pointwise bands, not simultaneous bands or individual first-passage distributions. Back-transformed log response point is a median. **unsupported:** unknown ID prediction; joint simultaneous bands; first-passage lifetime distributions; serial-correlation-adjusted inference; uncertainty from model selection

`chem_stability_analysis`: **data:** Identity-scale [nonnegative time,response,batch_id], >=2 batches and each batch supports an identified separate line with residual degrees of freedom. **args:** lower_on:required0/1; lower:required finite lower specification (placeholder ignored when lower_on0); upper_on:required0/1; at least one specification enabled; upper:required finite upper specification; lower<upper when both enabled; pooling_alpha:required in[0.01,0.5); ICH-style conventional choice0.25; confidence:required in(0.5,1); conventional choice0.95; pool_separate_MSE:required0=individual batch MSE/df for full-model bands,1=pooled MSE/df; horizon:required positive finite upper time searched **estimator:** Fit full separate intercept/slopes, separate intercept/common slope, and common intercept/slope. First F-test the common-slope restriction against full with full-model MSE; if p<alpha choose full. Otherwise F-test common intercept against the common-slope model using that unrestricted model MSE; if p<alpha retain separate intercept/common slope, else pool all. Full model bands default to individual batch MSE/df; comparisons still use nested homoskedastic F tests. One enabled specification uses one-sided t(confidence) mean bound; both specifications use two-sided t((1+confidence)/2). Shelf life is the earliest mean-bound departure from specifications for any batch after time0, within horizon; initial noncompliance gives0, no crossing gives null+beyond_horizon. Model selection uncertainty and regulatory extrapolation restrictions are not incorporated. **unsupported:** nonlinear/time-transformed stability; regulatory approval/automatic ICH extrapolation limits; multiplicity-adjusted simultaneous bands; correlated observations; Bayesian models

`chem_destructive_degradation`: **data:** <=4000 independent units, one measurement per destroyed unit. [time,stress,lower_response,upper_response,censor_type]. type0 exact requires lower=upper; type1 left requires lower placeholder0; type2 right requires upper placeholder0; type3 finite interval lower<upper. Lognormal used bounds positive. More exact measurements than location coefficients, exact-design rank identified. **args:** family:0=normal on original Y,1=normal on log(Y), default0; time_transform:0=identity,1=ln,2=sqrt; default0; stress_mode:0=no stress (stress column unused),1=linear stress,2=Arrhenius Celsius11604.5181215503/(stress+273.15),3=ln(stress); default0; maxiter:integer1..1000,default300; tolerance:1e-12..1e-3,default1e-8; average negative-likelihood score **estimator:** Mean location μ=b0+b1*f(time) when no stress; otherwise μ=b0+b1*f(time)+b2*r(stress)+b3*f(time)*r(stress). This is a bilinear location model in transformed stress, not an exponential Arrhenius-rate kinetic path. Gaussian density includes -logσ-log(2π)/2; lognormal exact density additionally includes -logY Jacobian. Left/right/interval censor likelihood uses CDF/SF/CDF differences; no density Jacobian is added to probabilities. Analytic location/log-scale scores, bounded BFGS with backtracking, stable normal log tails, observed-information covariance in centered standardized coordinates. σ MLE uses n, not residual df. Guard/iteration/gradient and covariance diagnostics are explicit. Converged additionally requires identified observed information; an uncalibrated/null covariance cannot claim an identified fit. Covariance and AIC are available only when inference_valid=true; optimizer projected-gradient convergence, active guard status and local information identifiability are reported separately. An active numerical coefficient/log-scale guard always invalidates inference, even when its projected gradient vanishes. **unsupported:** arbitrary nonlinear/kinetic paths; heterogeneous scale; random unit effects/repeated measures; censor-only identifiable designs without sufficient exact responses; measurement error in time/stress; Bayesian inference; MLE covariance/information criteria for unconverged or numerical-guard-bound estimates

`chem_destructive_predict`: **data:** [time,stress,response threshold] with same ordered transform domains; marginal response probability under fitted location-scale model. **args:** saved model [1,family,time_transform,stress_mode,sigma,centered_intercept,p,feature_origins[p],mean_offsets[p],original_feature_slopes[p]] **estimator:** Return response median, mean,0.025/0.975 distribution quantiles, and P(response<=threshold). For lognormal, median expμ differs from mean exp(μ+σ²/2). Percentile bands describe the response distribution, not parameter confidence intervals. Probability below a response threshold is a marginal soft-failure probability at a fixed time, not the first time an individual path crosses. **unsupported:** parameter confidence bands; individual lifetime first passage; unit-level repeated prediction/random effects

`chem_cumulative_damage_fit`: **data:** <=4000 event units. type0 exact0<lo=hi; type1 left lo0<hi; type2 right lo>0,hi placeholder0; type3 interval0<=lo<hi. Pattern IDs must exist. n>parameter count and at least parameter-count exact failures. >=2 distinct transformed stress levels and identified observed information for converged fit. **args:** distribution:required0=exponential,1=Weibull,2=lognormal; relation:required0=linear,1=ln(stress),2=Arrhenius Celsius11604.5181215503/(stress+273.15),3=reciprocal stress; continuation:required0=stop at pattern end,1=extend final level,2=repeat full pattern; maxiter:required integer1..1000; tolerance:required1e-12..1e-3 for average score; pattern_count:required integer1..128; each pattern: [ID,segment_count,(positive_duration,stress)*segment_count];1..128 segments each,<=1024 total; no trailing arguments **estimator:** Under an explicit cumulative-exposure assumption, log scale μ(s)=b0+b1*r(s), A(t)=integral_0^t exp[-μ(stress(u))]du. Step integration is exact using log-sum-exp; response density is f_A(A(t))*exp[-μ(stress(t))]. Standardized exposure distribution is Exp(scale1), Weibull(shape1/σ,scale1), or Lognormal(log-location0,log-scaleσ); no extra shape is estimated for exponential. Censor likelihood uses F(A),S(A) or interval differences. Analytic scores include exposure and density-rate derivatives. Standardized stress optimization with three starts, bounded BFGS, covariance from observed information; no global optimum guarantee. Event density at an interior stress change is right-continuous; stop endpoint uses final level. Report failed convergence or unidentified curvature honestly. Covariance and AIC are available only when inference_valid=true; optimizer projected-gradient convergence, active guard status and local information identifiability are reported separately. An active numerical coefficient/log-scale guard always invalidates inference, even when its projected gradient vanishes. **unsupported:** ramp/sinusoidal/piecewise-ramp stress; Box-Cox/custom stress relationship; loglogistic/Frechet lifetime; history effects beyond cumulative exposure; delayed entry/competing risks/frailty; Bayesian repeated measures; MLE covariance/information criteria for unconverged or numerical-guard-bound estimates

`chem_cumulative_damage_predict`: **data:** [nonnegative time,known pattern ID,quantile probability strictly(0,1)]. Stop continuation rejects prediction times after final pattern end. **args:** saved model [1,distribution,relation,continuation,b0,b1,sigma,pattern_count,patterns...]; exponential requires sigma1 **estimator:** Return cumulative effective age, log age, survival/failure probability, event-time density, hazard/log hazard, and time quantile obtained by monotone inversion of exact cumulative exposure. Density includes dA/dt. Shape-specific density and hazard limits at time0 are handled; mathematically infinite/nonrepresentable values serialize as null, log quantities retained where finite. A quantile beyond a stopped pattern is null; it is not silently extrapolated. These contracts model accelerated event times under cumulative exposure, distinct from Miner fatigue damage. Quantile inversion is on log time; positive lifetimes below the minimum representable double become null. **unsupported:** unknown/new stress patterns without an explicit new carrier; parameter uncertainty/confidence bands; ramp/sinusoidal schedules; mechanical Miner fatigue model

IRT posterior scores can truthfully have `quadrature_converged:false`, including sharp high-ability tails; this does not certify their precision. Fit convergence, local information identification, active numerical bounds and valid inference are distinct diagnostics. Unconverged fits keep numeric estimates/likelihood and reusable model parameters, while unsupported inference fields are null. Kernel symmetry, PSD support, original response/feature units and subnormal arithmetic follow the explicit repaired kinship contract.

<!-- analytics-specialized-contracts:end -->

<!-- analytics-streams-contracts:start -->

### Time series, functional preprocessing, text and spatial contracts

This checkpoint adds **29 CPU operations**, for **254 registered analytics functions / 405 total extension overloads**. Every call appears in the table above. Implementation stopped at this verified checkpoint under the user instruction to push and stop. Full JMP19.2 parity remains incomplete.

#### State Space

Complete numeric [API](tasks/jmp-state-space-api.json) and [independent reference](tasks/evidence/jmp-state-space-reference.json).

**finite data.** Complete finite equally spaced numeric data;no SQL NULL/NaN missing sentinel except explicit Kalman observed0/1 flag. All matrix lengths and model versions checked;finite output failures reject.

**normal intervals.** Confidence .5..0.9999,normal quantile approximation;conditional process/measurement uncertainty only. No fitted parameter/variance-estimation or future-input uncertainty.

**gaussian portable layout.** prediction_model=[1,0,k,last_filtered_mean(k),last_filtered_cov(k*k),F(k*k),Q(k*k),H(k),R],row-major,length4+2k+3k*k,k1..8. P may be PSD;Q PSD,R>0.

**ets portable layout.** prediction_model=[1,1,k,y_origin,y_scale,standardized_sigma2,last_state(k),F(k*k),w(k),g(k)],row-major,length6+3k+k*k,k<=32. Positive scale/sigma2. Tag1 is shared additive innovations;tag0 independent Gaussian. F/w/g may represent a supplied linear additive innovations model;predict does not verify one of the six fitted ETS family structures.

**optimizer.** Deterministic bounded Nelder-Mead with3startsETS/prewhite,4transfer. Stop requires maximum simplex coordinate spread<max(tol^.75,1e-8) and objective spread<.01*tol*(1+abs objective),then finite difference projected gradient<max(sqrt(tol),1e-7) in raw u coordinates. h=1e-5 score check. Objective is profile log(sigma2)+logdet/n (ETS logdet0). maxiter perstart;returned iterations is winning start. Local optimizer,not proof of global optimum.

**identification.** Initial/GLS design reorthogonalized QR,relative column threshold1e-11. Local profile shape curvature central differencesh1e-4 on free face,Cholesky relative pivot threshold1e-14. Active legitimate faces condition the curvature;guards disqualify convergence. No SEs or global identification guarantee.

**resource.** Gaussian k<=8,n<=5000,h<=500. ETS period<=24,n*k*k<=2000000,n*k*k*maxiter<=250000000. Prewhite n<=2000,n*maxiter*(1+p+q)<=40000000. Transfer n<=2000,n*maxiter*(1+den+p+q)<=60000000. These are nominal preflight work limits;multi-start/line-search/curvature evaluations add bounded work. Params<=16384 at native layer.

**JMP gaps.** X11 is NOT implemented or renamed from moving-average decomposition. Multiplicative ETS/automatic search/full admissibility,seasonal/integrated/multiple-input transfer and prewhitening,automatic lag selection,parameter uncertainty,and full report diagnostics remain. Existing CSS ARIMA module is untouched. Exact JMP CCF normalization/default optimizer/transfer numerator parameterization equivalence not asserted.

**external covariance validation.** External initial/process/saved Gaussian covariances must be exactly symmetric in binary64. Scaled semidefinite Cholesky rejects any negative diagonal pivot or nonzero residual off-diagonal under a zero pivot; roundoff-ambiguous real-PSD matrices may reject, so callers must explicitly supply representable well-conditioned covariance. No near-PSD projection is performed. Initial covariance additionally needs scaled pivots>1e-14. Internal RTS rounding validation remains a separate tolerant numeric check.

`chem_kalman_smoother`: **data:** 1..5000 rows [finite response, observed flag0/1]. Flag0 ignores response placeholder but retains elapsed time; SQL NULL is invalid. **args:** k integer1..8; horizon integer0..500; confidence .5..0.9999; m0: k entries; P0: k*k proper positive-definite covariance row-major; F: k*k transition row-major; Q: k*k positive-semidefinite process covariance row-major; H: k observation loadings; R: positive measurement variance **estimator:** alpha0~N(m0,P0) is prior for observation0; alpha(t+1)=F alpha(t)+eta,Q; y(t)=H alpha(t)+epsilon,R. Independent process/measurement disturbances. Exact innovations Gaussian likelihood,Joseph filtering covariance,RTS smoothing via positive-definite solves. Singular predicted covariance needed by RTS rejects; no pseudoinverse. Future variance propagates final filtered covariance with Q; observation variance adds R. Output covariance matrices are flattened row-major inside one row per time. **unsupported:** parameter estimation/diffuse initialization; multivariate responses,time-varying matrices/non-Gaussian models; RTS with singular predicted covariance; parameter-estimation uncertainty; This general independent-noise Gaussian model is distinct from JMP ETS.

`chem_state_space_smoothing`: **data:** 8..5000 complete finite scalar rows; n>=2*k+4 with k=1+trend+seasonal_period. Nonconstant data and nondegenerate residual variance required. **args:** family0 ANN/1 AAN/2 AAdN/3 ANA/4 AAA/5 AAdA;default0; period2..24 if seasonal,else unused0..24;default4 seasonal/0 otherwise; horizon1..500;default10; confidence .5..0.9999;default.95; maxiter1..5000 per start;default1000; tolerance1e-10..1e-3;default1e-8 **estimator:** y_t=w*x_(t-1)+e_t,x_t=F*x_(t-1)+g*e_t with shared e iidN(0,sigma2). For fixed smoothing parameters,rank-safe QR profiles all initial states with sum-zero initial seasonals. Conditional innovations Gaussian ML sigma2=SSE/n,three-start bounded Nelder-Mead plus independent projected-score check. alpha=.0001+.9998u,beta=alpha*u,gamma=(1-alpha)*u,phi=.8+.2u for damped. Alpha guard contact prevents combined convergence;beta/gamma0 and phi.8/1 are legitimate constrained faces. Local free-face profile curvature must be PD for combined convergence;no parameter SE is returned. Nominal parameter count includes smoothing parameters,independent initial states and sigma2,even on an active face. AIC/BIC null if nonconverged/guarded/unidentified. Future observation variance sigma2*(1+sum_(j=0..h-2)(w*F^j*g)^2),conditional on estimated initial state and point parameters. **unsupported:** multiplicative error/trend/seasonality; automatic ETS model selection/full admissibility region or exact JMP optimizer defaults; parameter SE/Wald/CI; missing/diffuse estimates,parameter estimation uncertainty; bounds are an explicit partial ETS contract;no full JMP family equivalence

`chem_state_space_predict`: **data:** Exactly one finite scalar placeholder row [[0.]] for SQL/native use,ignored. Empty data accepted only at direct module layer. **args:** horizon1..500; confidence .5..0.9999; unchanged prediction_model from kalman_smoother or state_space_smoothing; portable layouts in common_contract **estimator:** Gaussian forecasts transition final filtered state/covariance and add process Q each step,measurement R to observation variance. ETS forecasts deterministic final state plus shared-innovation convolution. Gaussian latent variance is conditional state measurement variance;ETS latent variance excludes the new observation innovation. Scaling uses standard deviations before squaring to avoid overflowing intermediate scale^2. Positive physical variance below representable range rejects. **unsupported:** refitting/future-label learning; parameter uncertainty; saved transfer-function reuse is not provided by this operation

`chem_arima_forecast_intervals`: **data:** 1..1000 finite scalar point-forecast rows. **args:** confidence .5..0.9999; known innovation sigma2>=0; AR polynomial degree0..200; MA polynomial degree0..200; difference factors0..4; AR coefficients A0..degree,A0=1; MA coefficients C0..degree,C0=1; one integer lag1..100 per factor;combined denominator length<=500 **estimator:** D(B)=product(1-B^lag),psi(B)=C(B)/(A(B)D(B)),psi0=1. Conditional forecast error variance at h=sigma2*sum_(j=0..h-1)psi_j^2. Point means are supplied,not refitted. Existing arima.rs uses CSS;these extra intervals do not relabel its estimator exact ML. finite-horizon unstable polynomials allowed if representable. **unsupported:** parameter/innovation-variance estimation uncertainty; Box-Cox/backtransform/non-Gaussian intervals; future-exogenous input uncertainty; changing existing CSS ARIMA estimator

`chem_cross_correlation`: **data:** 3..10000 finite [input,output] rows;both nonconstant;overlap Pearson mode also rejects constant overlaps. **args:** maxlag0..min(500,n-2); normalization0centered global-energy(default)/1overlap Pearson **estimator:** Positive lag k pairs X_t with Y_(t+k). Global mode uses full-series centered energies in denominator. Pearson mode re-centers each overlap. Returns n-|k|overlaps and white-noise SE1/sqrt(n-|k|),two-SE marginal bands;these do not adjust multiple lag selection or prove input whiteness. Stable origin/scale normalization. **unsupported:** missing/unequal spacing; multiplicity-adjusted lag significance; exact JMP finite-sample normalization equivalence is unverified

`chem_arima_prewhitening`: **data:** Finite [input,output] rows<=10000,n_train12..min(2000,n). Both parameter estimation and output centering use training prefix only. Constant input rejects;constant output allowed with null CCF. **args:** n_train12..min(2000,n); AR order0..2; MA order0..2; CCF maxlag0..min(500,n_train-2);default10; maxiter1..5000;default1000; tolerance1e-10..1e-3;default1e-8 **estimator:** Fit input mean,ARMA coefficients and innovation variance by exact stationary Gaussian likelihood. Small companion-state covariance solves the stationary Lyapunov series by doubling;shared Kalman startup gains whiten all columns. Input mean profiles by GLS;output centers by ordinary training mean,and applies the SAME input-estimated whitening matrix,not a separate output fit. Returned innovations are physical-unit innovations with unit-noise-covariance whitening,not divided by fitted innovation sigma. CCF is nested under ccf;constant output returns null. PACF coordinates ensure stationary/invertible factors with guards abs(kappa)<=.98;guard contact prevents combined convergence. Same local optimizer/free-face curvature criteria as transfer. **unsupported:** integrated/seasonal/full-order ARIMA prewhitening; missing observations; output-specific filter estimation; automatic order selection,whiteness proofs or complete model diagnostics

`chem_transfer_function_model`: **data:** Exactly n_train+horizon finite [input,output] rows,n_train20..2000,horizon1..500. Training input/output must be nonconstant with nondegenerate fitted residual variance. Future input path required;future output values are ignored. No automatic lag/order selection. **args:** n_train20..2000; input delay0..24; numerator degree0..2; denominator degree0..2; stationary-noise AR order0..2; invertible-noise MA order0..2; horizon1..500; confidence .5..0.9999;default.95; maxiter1..5000;default1000; tolerance1e-10..1e-3;default1e-8 **estimator:** v_t=sum(delta_j*v_(t-j))+sum(omega_j*(X_(t-delay-j)-training_input_mean));Y_t=intercept+v_t+n_t;phi(B)n=theta(B)e. Zero centered-input presample transfer state is explicit;ARMA-noise initial covariance is stationary. For given stable delta and stationary/invertible ARMA,profile intercept/omega/sigma2 by exact Gaussian GLS QR;four-start bounded local optimization over shape PACFs. Returned delta and AR arrays use positive recurrence coefficients,MA array uses theta(B)=1+theta1B+...;numerator degree d has d+1 omega entries. Forecasts add conditional ARMA-noise mean to deterministic specified-input transfer. Exact stationary-noise Gaussian LL,not CSS/independent-error OLS. Parameter count=denorder+p+q+numerator_degree+3. AIC/BIC null when nonconverged/guarded/unidentified. Free profile-shape curvature is a local identification diagnostic,not full covariance or global optimality proof. **unsupported:** multiple inputs,seasonal/differenced noise or missing data; automatic lag/order selection or exact JMP numerator convention/default optimizer equivalence; parameter SE/Wald/LR/CI and nonidentification proofs; future input or parameter-estimation uncertainty; portable transfer prediction carrier is not yet exposed

`chem_remove_linear_trend`: **data:** 3..10000 finite scalar rows;zero-based implicit time index0..n-1. **args:** n_train3..n;defaultall **estimator:** Rank-safe QR on centered/scaled time;fit uses prefix only,extrapolates all row indices. Returns intercept,slope,trend,residuals and stable coefficient representation prediction_model=[1,y_origin,y_scale,b0,b1,time_mid,time_scale];this carrier is metadata,not accepted by state_space_predict. **unsupported:** automatic differencing/nonlinear trend; iid slope inference in autocorrelated errors; calendar/irregular time indices

`chem_remove_cycle`: **data:** 4..10000 finite scalar rows;zero-based implicit equally spaced time index. **args:** period>2 and<=1e6; include_constant0/1;default1; n_train4..n;defaultall **estimator:** Rank-safe QR on cos(2*pi*t/period),sin(...),optionalconstant. Component includes fitted constant when requested. Amplitude=hypot(coscoef,sincoef),phase=atan2(-sincoef,coscoef),so component=intercept+amplitude*cos(wt+phase). Phase null when amplitude<=1e-12*training_scale. Prefix fit predicts all indices;no automatic frequency selection. **unsupported:** automatic period discovery,full spectral selection; general seasonal/X11 decomposition; irregular time/calendar effects or iid inference

#### Functional Preprocess

Complete numeric [API](tasks/jmp-functional-preprocess-api.json) and [independent reference](tasks/evidence/jmp-functional-preprocess-reference.json).

**data.** Nonempty finite long rows[x,y,id] or[x,y,id,role], uniform width3/4,<=100000 observations,<=128functions,<=10000observations/function. IDs exact integers within absolute2^53-1. Role0=analysis/training(default),1=validation; a single ID cannot mix roles.

**args.** <=16384 finite doubles with exact operation-specific count/ranges.

**ordering.** Function IDs preserve first appearance; observations sorted by x, then original row index for ties. Grid/baseline methods require distinct x; binning/thinning/DTW retain repeated observations explicitly.

**baseline validation.** SNIP/arWLS/supplied-baseline transform role0 curves only. Role1 original values pass unchanged, baseline/weights zero placeholders, applied=false and arWLS convergence=null. This does not mean a zero baseline was estimated.

**uncertainty.** No inferential p-values, SEs or confidence bands are claimed for preprocessing or curve-distance rankings. Convergence concerns the specified numerical iteration only.

**resource limits.** SNIP sum(training lengths)*max_halfwindow<=200000000; arWLS sum(training lengths)*maxiter*16<=200000000; DTW pair n*m<=4000000 and total unordered or target/candidate cell pairs<=200000000.

**uniform grid tolerance.** 1e-8*spacing+8*machine_epsilon*max(abs(x)); no implicit resampling before baseline methods.

**numeric.** Reject unresolved/rank-singular penalty systems and nonrepresentable required arithmetic. Null output is reserved for declared empty bins, outside-support values, skipped diagnostics or nonconvergence diagnostics. No JSON NaN/Infinity tokens.

**inference boundary.** No GUI selection semantics, FPC model fitting, golden-curve optimization, Fourier/basis/FDA model expansion is added here.

**saved roundtrips.** name: saved_common_grid_new_validation_only; fit_operation: functional_grid; saved_fields: grid; replay_operation: functional_grid; replay_args: [len(grid),grid[0],grid[-1],2]; replay_data: 0; 1; 8; 1; 2; 5; 8; 1; example_fit_data: 0; 0; 7; 0; 1; 2; 7; 0; 2; 4; 7; 0; known_replay: grid.1: 1.0; values.0.1: 3.0; refits: False; name: saved_SNIP_baseline_as_supplied_reference; fit_operation: functional_snip; saved_fields: x.0; baselines.0; replay_operation: functional_baseline; replay_args: 900; 2; replay_data_rule: Zip saved x/baseline into [x,baseline,900,0] reference rows; append [0,4,800,0],[2,6,800,0],[4,8,800,0] for example constant baseline1.; known_replay: corrected.0.0: 3.0; corrected.0.1: 5.0; corrected.0.2: 7.0; refits: False; validation_policy: If new curve role1, correction deliberately skipped; original [4,6,8] retained.; name: saved_golden_reference; replay_operation: functional_golden_match; replay_args: 900; 0; -1; 0; saved_fields: Explicit supplied target x/y observations, not estimated model coefficients; replay_data: 0; 1; 900; 1; 2; 900; 2; 3; 900; 0; 2; 800; 1; 3; 800; 2; 4; 800; known_replay: candidate_ids.0: 800.0; scores.0: 1.0; target_excluded: True; refits: False

`chem_functional_grid`: **data:** Long function observations,>=2 distinct x per ID; internally sorted. Auto/size-only grid estimates use only role0 functions. Explicit bounds/size requires no training rows. **args:** []:size=min(4096,longest training function),domain=intersection of training supports,outside null; [size]:integer2..min(4096,longest training function); training-domain intersection; [size,lower,upper,outside_policy]:size2..4096,lower<upper; all quantities supplied, so validation-only replay and size larger than query sample count are supported; outside_policy:0=null beyond support,1=constant endpoint extension,2=reject; required in full layout **estimator:** Piecewise linear interpolation on a single evenly spaced grid. The full-argument layout replays a stored grid and estimates no new statistics. Auto bounds/length never use validation curves. Endpoints are exact, interior weights are convex; grids finer than representable x spacing reject. Null values explicitly mean unsupported query locations under policy0. Full explicit replay extends the JMP-style automatic longest-function size limit; it is an implementation contract, not an undocumented JMP default claim. **unsupported:** cubic/spline interpolation; duplicate-x aggregation; automatic extrapolation; implicit missing-data imputation; exact JMP domain/default interpolation equivalence

`chem_functional_bin`: **data:** Long observations; repeated x accepted and all observations counted, including multiplicities. Same explicit bin edges for training and validation. **args:** number:required integer1..4096; lower,upper:required finite lower<upper; outside_policy:0=reject any out-of-range value,1=drop/report count;default1 **estimator:** Average observed Y within left-closed/right-open bins, with final upper edge included. Empty bins preserve count0 and null mean. Report centers, edges, per-function counts and dropped observations. No interpolation or weighting by time interval. User-supplied bounds avoid fitting a domain to validation data. Edges/centers must be representably increasing. **unsupported:** quantile/variable-width bins; interval-integral means; automatic training-domain estimation; silent empty-bin imputation

`chem_functional_thin`: **data:** Long observations; each ID sorted by x and original row index for ties; repeated x allowed. **args:** stride:integer1..10000,default2; offset:integer0..stride-1,default0; keep_last:0/1,default0; union final point with stride selection without duplicate **estimator:** Retain sorted indices i with i mod stride=offset, optionally retaining the final point. Return selected x/Y and original zero-based input row indices. This explicit keep-every-stride convention removes ambiguity in thinning-rate wording; it does not estimate a smoother. An ID can remain with an empty selected array. **unsupported:** random/importance thinning; automatic anti-alias filtering; screen row selection state

`chem_functional_snip`: **data:** Training IDs require an evenly spaced strictly increasing x grid and n>=2*max_halfwindow+1. Validation IDs pass unchanged without fitting; their baseline arrays are zero placeholders with baseline_applied=false. Log-log input Y>=0; raw transform allows signed Y. **args:** max_halfwindow:required integer1..4999,also <=(n-1)/2 for every training ID; transform:0=identity,1=G(y)=ln(ln(1+y)+1);default1; original Ryan log-log transform has no square root; direction:0=widths1..max,1=max..1;default0; each pass reads the previous entire vector synchronously **estimator:** At each width w replace z_i by min(z_i,(z_(i-w)+z_(i+w))/2) where both neighbors exist; retain other points. Invert G after clipping and return baseline plus Y-baseline. This is an explicit span-scheduled SNIP variant based on peak clipping and Ryan double-log compression. It does not implement the full original energy-dependent FWHM/adaptive low-statistics filter, later LLS square-root variant, or unknown JMP settings. Clipping is distinct from a moving average. **unsupported:** original FWHM-dependent low-statistics filter; implicit irregular-grid interpolation; LLS square-root transform; automatic peak-width choice; JMP-default numerical equivalence; validation baseline estimation

`chem_functional_arwls`: **data:** Training IDs require >=3 uniformly spaced distinct x values. Validation rows remain unchanged with zero placeholder baseline/weights, convergence null, baseline_applied=false. **args:** lambda:1e-8..1e10,default1e5; penalty on index-based second differences; weight_tolerance:1e-12..1e-2,default1e-6; L2 relative proposed-versus-used weight change; maxiter:integer1..1000,default100; negative_std_ddof:0=population,1=sample;default0; explicit convention because exact software defaults are not assumed **estimator:** Solve (W+lambda D2^T D2)z=Wy with a true pentadiagonal Cholesky factorization. Start weights1; for d=y-z, estimate m,s from negative residuals and propose w_i=1/(1+exp(2*(d_i-(-m+2s))/s)), as Baek2015 arPLS (called ARWLS in JMP menu). Stop on relative weight change and small linear backward error. Returned weights are the weights that actually solved the returned baseline, not an unsolved next proposal. Zero/unidentified negative noise scale, iteration limits or numerical linear failures never claim convergence. No hidden ridge or weight floor. Affine signals within numeric tolerance have zero second-difference baseline correction and converge without iterations. **unsupported:** airPLS/AsLS substitution; moving-average substitution; implicit ridge/weight floor; nonuniform physical-spacing penalty; automatic lambda choice; baseline parameter confidence intervals; JMP-default numerical equivalence

`chem_functional_baseline`: **data:** Baseline ID must exist with role0,>=2 distinct x; its observations are a supplied reference and excluded from corrected IDs. Other training IDs may be irregular. Validation IDs pass unchanged. **args:** baseline_id:required exact ID of supplied reference function; outside_policy:0=null,1=constant endpoint extension,2=reject;default2 **estimator:** Linearly interpolate the supplied baseline to each training function x and subtract it. Return original data, interpolated baselines, corrected data and applied flags. Null outside support follows explicit policy0. The reference is neither estimated from other curves nor included in corrected candidates. A role1 validation function cannot be repurposed as a training baseline. **unsupported:** parametric baseline fitting/anchor regions; implicit validation correction; baseline extrapolation beyond selected policy; baseline uncertainty propagation

`chem_functional_dtw_matrix`: **data:** >=2 functions in within-ID sorted x order; unequal lengths permitted. Each ID has at least1 point. X determines sequence order, not a physical-time cost. **args:** window:-1=unconstrained,otherwise integer0..longest length; must reach every pair endpoint; normalize:0=original intensities,1=per-curve mean/population-standard-deviation standardization;default0; constant curves reject under1; cost:0=squared local amplitude difference and square root of summed path cost,1=absolute local difference and summed cost;default0 **estimator:** Sakoe-Chiba-band DTW with diagonal/up/left predecessor steps, two-row dynamic programming. Compute each unordered pair once and mirror it, diagonal0. Normalization is independently within each supplied curve, with scale-safe arithmetic; no training-set fitting. Window infeasibility or resource excess rejects rather than widening the window. DTW distances need not satisfy the triangle inequality. This contract supplies the numerical matrix for downstream clustering; it does not construct a JMP report script. **unsupported:** physical-x/time penalization; multivariate local costs; warping-path/warping-function export; soft DTW; automatic registration; triangle-inequality guarantees

`chem_functional_golden_match`: **data:** >=2 functions; explicit target ID excluded from candidate scores/ranking. RMSE/MAE require identical ordered x grids; DTW allows unequal lengths. Curve roles retained for candidate provenance. **args:** target_id:required exact ID of pre-specified reference/golden curve; metric:0=RMSE,1=MAE,2=squared-cost square-root DTW;default0; window:-1=unconstrained,otherwise valid DTW radius;default-1; normalize:0=original intensities,1=independent per-curve population standardization;default0 **estimator:** Compare each candidate with a supplied target, return scores in input function order and ascending ranking, deterministic ties by candidate order. Target is never included among candidates or used to estimate a common curve/model. Stable RMS avoids overflowing an intermediate unnormalized Euclidean norm. This is explicit distance matching, not FPC-score desirability profiling, model fitting, or process-setting optimization. **unsupported:** FPC profiler/desirability optimizer; automatic golden-target selection; affine amplitude calibration; counterfactual process outcomes; statistical match p-values/confidence intervals

#### Text Analysis

Complete numeric [API](tasks/jmp-text-analysis-api.json) and [independent reference](tasks/evidence/jmp-text-analysis-reference.json).

**input.** Nonempty finite rectangular numeric data,abs<=1e12,<=131072cells,columns1..4096;<=16384finite args. JSONversion1. Missing/NaN/Infinity reject. All term/row ordering caller-owned and stable; empty documents represented by zero count rows or explicit n_documents in text_dtm. No empty-column SQL inputs.

**tokenization.** No Unicode/string segmentation,normalization,language dictionaries,stemming,stopwords or synonym detection implemented. Caller supplies token IDs/counts and same fitted vocabulary for scoring; Japanese tokenization remains unsupported.

**weighting.** Binary 1[x>0];ternary min(x,2);frequency x;log frequency log10(1+x);TF-IDF x*log10(N/df) without smoothing or document-length scaling. Training N includes empty documents. Terms selected by corpus count descending with original-index ties. Saved IDF never refits query documents.

**svd normalization.** Follow literal JMP19.2 LSA help: divide uncentered matrix by n;centered/scaled by n-1. Centered-scaled uses sample standard deviation,constant columns scale1. This states the numeric contract; it is not sqrt(n) or sqrt(n-1) normalization. The help also describes PCA equivalence; no installed-JMP singular-value scale comparison has been performed,so that prose is not treated as independent validation.

**lsa scores.** For transformed A,approximation U*S*V^T. document_vectors=U,document_scores=US,term_vectors=V. Saved projection V scores US. Sign anchor is maximum absolute term-vector entry positive,source index tie. Repeated singular values allow alternative bases; compare subspace/reconstruction. Relative reconstruction error uses actual normalized residual,not subtracting close energies.

**topic scores.** L=V*S;orthogonal R from varimax;topic_loadings=L*R,document_scores=U*R,so scores*loadings^T=U*S*V^T. Projection=V*S^-1*R. Objective=sum_j[sum_i Z_ij^4-(sum_i Z_ij^2)^2/p],Z=L/max_abs(L),then optional Kaiser unit row norms. varimax_loading_scale records max_abs(L). rotation convergence concerns local sweeps only. Positive topic energy underflow or unrepresentable projection coefficients rejects with rescale guidance.

**latent class.** P(x)=sum_c pi_c product_t theta_ct^presence(1-theta_ct)^absence. Exact logsumexp E step;M step weighted means clamped to supplied probability bounds. Stop requires relative LL change<=tol and maximum parameter step<=sqrt(tol). Component mass<=1e-12 errors. Components canonically sorted by mixture probability then lexicographic term probability.

**diagnostics.** LCA separates optimizer_converged,components_distinct,necessary dimension bound(k*p+k-1 <= 2^p-1),and identification_verified=false. converged requires optimizer convergence and passes these necessary degeneracy checks; they do not prove statistical identification. AIC/BIC are descriptive conventional parameter penalties and null unless converged; no inference/calibration guarantee for singular mixtures. Varimax reports local rotation_converged separately.

**resources.** DTM output dimensions<=4096,cells<=131072. LSA n<=1024,p<=128,k<=32,n*p*p<=1M,100Jacobi sweeps. Topic k<=16,p*k*k*sweeps<=50M,max2000sweeps. LCA n2..2048,p<=128,k<=16,n*p*k*maxiter<=100M,maxiter2000. Scoring outputcells<=131072. Derived nonfinite results error.

**saved models.** operation: text_weight; result_path: weight_parameters; prediction_operation: text_weight_predict; layout: [version1,original_term_count p,selected_count k,weight_mode,training_n,selected_term_indices(k),training_IDF(k)],length5+2k; original ID universe retained; operation: text_lsa; result_path: projection_parameters; prediction_operation: text_projection; layout: [version1,input_terms p,output_dimensions k,training_divisor,training_means(p),training_scales(p),projection_matrix_row_major(p*k)],length4+2p+p*k; score=((query-mean)/scale/divisor)*projection; operation: text_topic_svd; result_path: projection_parameters; prediction_operation: text_projection; layout: [version1,input_terms p,output_dimensions k,training_divisor,training_means(p),training_scales(p),projection_matrix_row_major(p*k)],length4+2p+p*k; score=((query-mean)/scale/divisor)*projection; operation: text_latent_class; result_path: predict_parameters; prediction_operation: text_latent_class_predict; layout: [version1,nclasses k,nterms p,probability_floor,mixture_weights(k),Bernoulli_term_probabilities_class_major(k*p)],length4+k+kp

`chem_text_dtm`: **data:** Long token-count records,duplicates summed; IDs nonnegative integers within supplied dimensions,count integer0..1e9; aggregate cell<=1e9. Zero-count records supported; explicit empty documents preserved. **args:** n_documents integer1..4096; n_terms integer1..4096 **unsupported:** No string tokenization,Unicode normalization,Japanese dictionary,stemming,stopword or synonym fitting

`chem_text_weight`: **data:** Dense count matrix in immutable source vocabulary order,count integer0..1e9. Training frequency filters and top-term ordering are corpus count descending,source-index ascending ties; empty documents count in N. **args:** mode integer0..4 default4; min_corpus_count integer1..1e12 default1; max_terms integer1..p defaultp **unsupported:** No smoothed IDF,document-length normalization,automatic language preprocessing or dictionary from query data

`chem_text_weight_predict`: **data:** Dense query counts in entire original p-term universe; absent terms0; unselected training terms ignored. Unknown external words must be excluded/mapped by caller; out-of-range numeric IDs/width reject. **args:** unchanged weight_parameters layout [version1,original_term_count p,selected_count k,weight_mode,training_n,selected_term_indices(k),training_IDF(k)],length5+2k; original ID universe retained **unsupported:** No query IDF refit or inferred new vocabulary; OOV words are not automatically tokenized

`chem_text_lsa`: **data:** Supplied weighted DTM (use text_weight first if desired),n<=1024,p<=128,n*p*p<=1e6,k<=32 and numerical rank. Relative rank cutoff1e-10; one-sided Jacobi100sweeps; constant centered columns scale1. **args:** k integer1..min(n,p,32) required; preprocess integer0..2 default2 **unsupported:** No automatic TF-IDF fitting,model-rank selection,sparse unlimited corpus or installed-JMP equivalence of singular-value scale; repeated singular subspaces need not have unique vectors

`chem_text_projection`: **data:** Nonnegative weighted query matrix in exact fitted term order. Input columns<=128,output dimensions<=32,output<=131072cells; finite derived binary64 outputs required. TF-IDF must use original training model separately. **args:** unchanged projection_parameters layout [version1,input_terms p,output_dimensions k,training_divisor,training_means(p),training_scales(p),projection_matrix_row_major(p*k)],length4+2p+p*k; score=((query-mean)/scale/divisor)*projection **unsupported:** No query preprocessing refit,OOV inference,uncertainty or guaranteed finite result when derived binary64 terms overflow

`chem_text_topic_svd`: **data:** Same supplied-DTM SVD bounds; k<=16. Loadings initially VS; common max-absolute loading scaling before optional Kaiser row normalization,then pairwise orthogonal Jacobi varimax. Final maximum absolute term-loading anchor positive. **args:** k integer1..16 and numerical rank required; preprocess integer0..2 default2; max_sweeps integer1..2000 default200; tolerance1e-12..1e-3 default1e-8; kaiser integer0/1 default1 **unsupported:** No NMF/LDA/probabilistic-topic claim,global rotation maximum,automatic topic naming,GUI wordcloud or installed-JMP exact rotation reproduction

`chem_text_latent_class`: **data:** Count matrix binarized by presence,training n2..2048,p<=128,k<=16 and<=n,n*p*k*iterations<=100M. Deterministic SplitMix64 seeded shuffled distinct document indices initialize .75/.25 conditional probabilities. **args:** k integer1..min(n,16) required; seed integer0..2^53-1 default0; max_iterations integer1..2000 default200; tolerance1e-12..1e-3 default1e-8; probability_floor1e-12..0.1 default1e-9 **unsupported:** No global optimum or verified statistical identification,SE/CI,LDA/NMF,sparse JMP Text Explorer implementation identity,term-selection regression or sentiment dictionary

`chem_text_latent_class_predict`: **data:** Dense query count matrix in fixed p-term order,counts converted to presence; saved probabilities/mixing weights validated. Output<=131072cells. **args:** unchanged predict_parameters layout [version1,nclasses k,nterms p,probability_floor,mixture_weights(k),Bernoulli_term_probabilities_class_major(k*p)],length4+k+kp **unsupported:** No refitting,term-map inference,classification calibration or guarantee supplied arbitrary parameter vectors represent identified fitted models

#### Spatial

Complete numeric [API](tasks/jmp-spatial-api.json) and [independent reference](tasks/evidence/jmp-spatial-reference.json).

**input.** Finite nonempty rectangular data,abs<=1e100,at most1M cells; caller handles missing observations and coordinate units. Different methods have explicit row/dimension caps. All supplied parameter entries finite,<=16384.

**kernels.** kind0 spherical(max(1-d/rho,0)^2*(1+d/(2rho)) for d<rho),1 exp(-d/rho),2 exp(-(d/rho)^2),3 rho^d. Euclidean isotropy in1..3 coordinates; partial sill1e-100..1e100,nugget0..1e100,rho>0<=1e100 and power rho<1. Four19.2 equation images individually inspected. Nugget is an absolute variance, unlike JMP displayed nugget/residual scaling.

**measurement errors.** Nugget belongs to an observation record,not merely a coordinate. Distinct collocated records covariance=partial sill; diagonal=partial sill+nugget. Pair semivariance0 on self diagonal and nugget at distinct collocated observations. Saved prediction targets an independent new measurement at query coordinates; latent cross covariance never includes query nugget.

**GLS.** Fixed supplied absolute covariance only,not estimated from y. Exact fixed-response fits have finite likelihood under known V; no fictitious residual-variance estimate. Returned coefficient covariance includes intercept and known covariance only. Numerically singular covariance/design rejects at relative Cholesky diagonal1e-12; cov constructor may returnPSD/singular matrices with positive_definite_at_tolerance=false.

**prediction.** Universal-kriging variance includes estimation of the fixed trend with covariance known; no covariance-parameter estimation uncertainty. Tiny relative-negative rounding at1e-9 of contributing terms clamps to0; material negative variance rejects. Saved triangular factors and positive G inverse are syntax/PD validated,not cryptographically authenticated or reconstructed to prove consistency.

**saved model.** [1,n,d,p,kind,sill,nugget,rho,feature_origins(p),feature_mean_offsets_scaled(p),feature_difference_scales(p),bstd(p+1),training_coords(n*d rowmajor),normalized_design(n*(p+1) rowmajor),R_inverse_residual_original_y_units(n),lower_Cholesky_R(n(n+1)/2 rowmajor),lower_symmetric_inverse_Gram((p+1)(p+2)/2 rowmajor)]. R=V/(sill+nugget),bstd is coefficient vector in original response units against normalized centered features. Length8+3p+(p+1)+n*(d+p+2)+n(n+1)/2+(p+1)(p+2)/2,<=10580.

**variogram.** No response-trend removal,robust/Cressie estimator,minimum-pair threshold,spatial weights or variogram fitting. Default distance range uses observed minimum/maximum; exact collocated distances included; empty bins use JSON null. Width boundaries are the specified equal-width contract,not a claim to duplicate every JMP automatic bin adjustment.

**limits.** Kernel parameters and GLS coefficients must have finite physical binary64 outputs; otherwise return error/SQL NULL with rescaling guidance. No installed-JMP or cross-platform numerical validation.

**numerical repairs.** Direct semivariance uses expm1/analytic correlation complements and physical small-distance products,retaining nugget and representable tiny process increments even when covariance entries round equal. GLS coefficient covariance applies sqrt(total_sill) to transformation factors before quadratic products. Material-negative predictive variance uses a normalized finite comparison,never an overflowing absolute tolerance sum. Zero-nugget queries exactly matching a training coordinate and normalized trend row have structural zero latent/observation variance after material-negative validation.

`chem_spatial_covariance`: **data:** 1..256 rows of1..3 spatial coordinates **args:** Exactly[kind,sill>0,nugget>=0,rho>0]; power0<rho<1 **estimator:** Direct four published correlation functions; nugget added only on diagonal observation records **unsupported:** Covariance range/sill/nugget ML/REML or variogram-NLS estimation; Anisotropic/directional parameter fitting and random-effects G-side combinations; Covariance-parameter uncertainty,spatial F/df corrections and automated model selection; No installed-JMP numerical comparison

`chem_spatial_variogram`: **data:** 2..512 rows[coordinates1..3,response] **args:** 0..3 values: bins1..100=10,lower>=0=minimum observed distance,upper>=lower=maximum observed distance **estimator:** All unordered i<j pairs in inclusive bounds; equal-width classes left-closed/right-open except final inclusive right edge; identical bounds collapse to one bin **unsupported:** Covariance range/sill/nugget ML/REML or variogram-NLS estimation; Anisotropic/directional parameter fitting and random-effects G-side combinations; Covariance-parameter uncertainty,spatial F/df corrections and automated model selection; No installed-JMP numerical comparison

`chem_spatial_gls`: **data:** n>p+1,n<=128 rows[d coordinates,p features,response],intercept automatic **args:** Exactly[d,p,kind,sill,nugget,rho]; covariance parameters physically supplied and fixed **estimator:** Cholesky normalized by total sill; origin/difference-scaled feature and response GLS; standardized whitened Gram Cholesky(relative diagonal tolerance1e-12); known absolute V, no residual-scale fit **unsupported:** Covariance range/sill/nugget ML/REML or variogram-NLS estimation; Anisotropic/directional parameter fitting and random-effects G-side combinations; Covariance-parameter uncertainty,spatial F/df corrections and automated model selection; No installed-JMP numerical comparison

`chem_spatial_predict`: **data:** Query rows[d coordinates,p fixed features] in original units,exact width d+p **args:** Unchanged predict_parameters; never reads query outcomes or refits covariance **estimator:** x beta + c^T V^-1 residual; latent variance=sill-c^TV^-1c+delta^T(X^TV^-1X)^-1delta,delta=x-X^TV^-1c; observation variance adds independent query nugget **unsupported:** Covariance range/sill/nugget ML/REML or variogram-NLS estimation; Anisotropic/directional parameter fitting and random-effects G-side combinations; Covariance-parameter uncertainty,spatial F/df corrections and automated model selection; No installed-JMP numerical comparison

Numerical repairs are regression-tested: external Gaussian covariances reject asymmetric or negative-pivot inputs; spatial semivariance retains a small independent nugget and small-distance increments; physical GLS covariance avoids intermediate overflow; malformed saved spatial models reject material negative variance even at extreme scales; exact noiseless training records have zero prediction variance. These do not certify general hyperparameter estimation, uncertain-parameter intervals or whole-platform equivalence.

<!-- analytics-streams-contracts:end -->

### Worked workflows and saved models

```sql
-- Molecular descriptors become numeric features with an explicit row order.
WITH molecules AS (
    SELECT * FROM (VALUES (1,'CCO',1.0),(2,'CCCO',1.8),
                         (3,'CCCCO',2.1),(4,'CC(=O)O',1.3)) AS t(id,smiles,y)
), features AS (
    SELECT id, mol_weight(smiles) AS mw, logp_crippen(smiles) AS lp, y
    FROM molecules
)
SELECT chem_pca(list([mw,lp] ORDER BY id), [2.0,1.0]) AS chemical_space,
       chem_ridge(list([mw,lp,y] ORDER BY id), [1.0]) AS property_model
FROM features;

-- Store stable model parameters and score raw feature rows in original units.
WITH fit AS (
    SELECT chem_ridge([[0.0,1.0],[1.0,2.0],[2.0,4.0],[3.0,4.5]], [1.0]) AS model
)
SELECT chem_predict([[4.0],[5.0]],
    json_extract(model, '$.stable_predict_parameters')::DOUBLE[]) AS predictions
FROM fit;

-- A sample of right-censored lifetimes: event=0 means a surviving/censored unit.
SELECT chem_survival_weibull(
    [[0.5,1.0],[1.0,0.0],[1.5,1.0],[2.5,1.0],[4.0,0.0],[6.0,1.0]],
    []::DOUBLE[]);

-- ARIMA(0,1,0): no drift, three forecasts in the original response units.
SELECT chem_arima([[1.0],[3.0],[2.0],[5.0]], [0.0,1.0,0.0,3.0,0.0]);
```

Use the model-specific saved prediction parameter field rather than reconstructing
coefficients from a display. Ridge/PLS use `stable_predict_parameters` with
`chem_predict`; quadratic response surfaces use `model_parameters` with
`chem_response_surface_predict`; GLMs use their saved
`stable_predict_parameters` with `chem_glm_predict`. The complete molecular
feature-to-candidate workflow is also documented in [ANALYTICS.md](ANALYTICS.md).

### Interpretation, statistical assumptions and limits

- **Inference:** Parametric t/F/variance intervals assume the corresponding
  independent normal/common-variance model where required. Mann-Whitney,
  signed-rank and Kruskal-Wallis methods use their documented asymptotic
  calibrations and tie/zero corrections. Probability-ordered binomial and
  fixed-marginal 2x2 Fisher tests use explicit exact definitions. A one-sample
  TOST is not a complete two-sample, nonparametric or non-inferiority platform.
- **OLS, GLMs and penalties:** Inspect convergence and rank. Separate point
  prediction from coefficient and prediction uncertainty. Logistic/Poisson
  penalties act in standardized predictor units; an unpenalized separated model
  can fail. Elastic net is a fixed-parameter solver, not a complete tuning,
  categorical effect-selection or generalized-regression platform.
- **Control and measurement systems:** Limits are numeric outputs, not charts
  or evidence of process stability. IR estimates sigma from adjacent moving
  ranges; subgroup charts require the specified subgroup design. EWMA/CUSUM
  reference values must be supplied deliberately. Gage R&R is balanced crossed
  ANOVA with its declared interaction treatment, not arbitrary nested/unbalanced
  MSA. Cohen kappa is unweighted two-rater agreement.
- **Distribution and lifetime fits:** Gamma, Weibull, exponential and lognormal
  location is fixed at zero. Normal/lognormal scale uses the MLE denominator
  `n`. Uncensored distribution fitting and right-censored survival likelihoods
  are separate functions. The exponential lifetime rate includes all exposure;
  its standard error is asymptotic observed information. Weibull/Gamma shape
  roots must be identifiable and converge. No naive parameter-fitted KS p-value
  is returned. `chem_normality_test` is explicitly asymptotic Jarque-Bera, with a
  small-sample warning; it does not implement Shapiro-Wilk.
- **Transforms and densities:** Box-Cox requires positive nonconstant values and
  returns the profile likelihood including its Jacobian, with the Gaussian
  additive constant omitted. With no lambda it brackets and maximizes the
  profile; an unbracketed optimum fails explicitly. Gaussian KDE uses the
  caller's positive bandwidth and evaluation points, not an automatic selector.
- **ARIMA:** The method fits conditional sum of squared innovations, with zero
  presample deviations/innovations and a conditioning burn equal to the maximum
  composed AR/MA polynomial degree. It is not exact state-space likelihood.
  Reflection coefficients enforce stationarity/invertibility of the ordinary
  and seasonal factors. `include_constant=1` fits the mean of the differenced
  series; after one ordinary difference this acts as drift. Forecasts invert
  every ordinary/seasonal difference using stored histories. Inspect
  `converged`; an iteration-limited fit is not a certified global optimum.
  These AIC/BIC values use the conditional effective sample and should not be
  treated as interchangeable with exact-likelihood criteria or compared blindly
  across different differencing/conditioning samples. Parameter confidence
  intervals, estimated-parameter forecast uncertainty, exact likelihood and
  automatic order selection are separate unfinished work.
- **Additional DOE:** Box-Behnken enumerates two-factor edges; the jittered Latin
  hypercube stratifies each factor. The mixture lattice offsets lower bounds and
  filters upper bounds; this does not enumerate every constrained mixture design.
  Fractional factorial uses explicit positive-product generator masks without
  inferring alias resolution. Orthogonal screening uses Sylvester-Hadamard main
  effects and is not a general definitive-screening or Plackett-Burman designer.
  Randomization creates independent seeded replicates and equal-size blocks;
  design evaluation adds an intercept to supplied model-term columns. These
  operations do not supply split-plot error structures or automatic design power.
- **Original DOE and optimization:** Full factorial is two-level; CCD axial runs
  can extend outside factor bounds when alpha exceeds one. D-optimal exchange
  searches an explicitly supplied feasible candidate set and provides a local
  solution. GP inference uses fixed isotropic RBF hyperparameters and finite-pool
  EI; it does not perform continuous, constrained multiresponse or noisy-batch
  Bayesian optimization. The original random-intercept estimator remains a
  balanced one-way variance calculation, separate from broader mixed models.

Common SQL/C ABI bounds are 100,000 rows, 1,024 columns, 1,000,000 numeric cells,
16,384 numeric parameters and 64 MiB of JSON output. Method-specific limits are
stricter; excessively expensive or unrepresentable requests fail. Important
limits include:

| Methods | Additional bounded contract |
| --- | --- |
| PCA/k-means; Ridge/PLS | At most 256 features, plus per-operation work limits |
| DOE/response surfaces | At most 16 factors; candidate/design/term work is bounded; see the detailed API |
| GP/diversity | GP at most 512 observed rows and 10,000 candidates; diversity has an explicit work limit |
| Distribution/lifetime functions | At most 100,000 rows; finite values with absolute magnitude at most `1e150`; positive support and identifiability where required |
| Gaussian KDE | At most 4,096 evaluation points and 20,000,000 observation-point products |
| ARIMA | `p,q<=3`; `d,D<=2`; `P,Q<=1`; season `2..52` when seasonal terms are used; horizon `0..1000`; iterations `0..10000`; tolerance `1e-12..1e-3`; bounded work at most 500,000,000 |

These are in-memory numeric calls, not distributed or streaming statistical
solvers. An older published community package, another OS/architecture or Wasm
requires its own build and validation; successful local tests do not establish
cross-platform numerical parity.

### JMP/JMP Pro coverage and verification

[JMP_PARITY.md](JMP_PARITY.md) and the
[machine-readable feature matrix](tasks/jmp-feature-matrix.json) record official
JMP/JMP Pro 19.2 source URLs, edition evidence, data/assumption/test contracts and
missing features. All 502 leaf/method occurrences from 14 official capability
pages are retained, with selected reference suboptions expanded into 552 rows, including the separately sourced Item Reliability suboption.
These unequal-granularity rows cannot be converted into a percentage of JMP
coverage. This is a **19.2 reference baseline**, not a current-version parity claim.
The [20.0 help URL](https://www.jmp.com/support/help/en/20.0/jmp.html) responded,
but its retrieved index was byte-identical to the 19.2 index on 2026-10-21;
[the official documentation catalog](https://www.jmp.com/en/support/jmp-documentation)
still labeled 19.2 current. No distinct official 20.0 What's New was verified, so
20.0 additions remain an explicit unresolved follow-up rather than guessed rows.
See [the version audit](JMP_PARITY.md#200-urlの確認と版を混ぜないための記録).
A function name or a restricted implementation does not establish whole-platform
equivalence.

Large outstanding families include GLMM random slopes/multiple factors and other response families, spatial covariance-parameter fitting, general SEM/FIML/multigroup beyond simple CFA and supplied RAM, further causal estimators and uplift-tree fitting, irregular functional inference and automatic basis selection, full wavelets/X11/multiplicative ETS, seasonal or multivariate transfer, arbitrary DOE/mixture/restricted-randomization designs, Bayesian repeated-measures degradation and further reliability workflows, tokenization/sentiment/text regression and multidimensional IRT, choice/conjoint/MaxDiff, GWAS/multikernel inference and linked/polyploid admixture, further tuning/uncertainty/deployment/export options. The new numeric text/functional/time-series/spatial contracts are bounded subsets. DuckDB joins, aggregates,pivot/unpivot and SQL I/O remain reused host functionality. Pure GUI operations and GPU dependence are excluded; static figures/report formats remain separately tracked. No full-JMP coverage percentage is inferred from the unequal inventory rows. The user requested pushing this completed checkpoint and stopping; additional feature implementation is stopped.

The implementation workflow is test driven: observe an expected behavioral RED,
implement the bounded contract, then verify GREEN, numerical edge cases and
independent references before matching SQL/C ABI checks. Direct tests and a
loaded-extension test are distinct layers. No installed JMP was used for a
head-to-head numerical comparison.

<!-- analytics-validation-summary:start -->

- The current **254-function checkpoint** passed **1,302 Rust tests across 63 suites**, zero failures; two pre-existing external-data tests remain ignored. Strict analytics Clippy, scoped formatting and `git diff --check` passed.
- The matching macOS arm64 CPU Debug extension passed **13 native SQL cases / 439 assertions**, no skips, and every **405 documented catalog overloads / 405 executable examples**. C++ uses ASan/UBSan; Rust is not sanitizer-instrumented.
- All **29 new operations / 106 fields** passed the production exported C ABI. The general new native checker tested **29 examples / 118 fields**, **6 SQL saved-model round trips**, **39 invalid NULL cases**, selected/interleaved vectors; **45 strict JSON results / 11 queries**. The functional checker added **3 saved-reference round trips / 18 fields**, **43 invalid NULL cases**, explicit arWLS nonconvergence, selected/interleaved vectors; **22 JSON results / 7 queries**.
- All six prior numerical/model/design checkers passed against the matching build. [Current evidence](tasks/evidence/jmp-streams-validation.json) identifies the complete source graph, actual artifacts, logs and reproduction commands. [225-function](tasks/evidence/jmp-specialized-validation.json), [209-function](tasks/evidence/jmp-models-validation.json), [181-function](tasks/evidence/jmp-continuation-validation.json) and [119-function](tasks/evidence/jmp-sql-validation.json) records remain historical.

| Current added domain | Independent reference | Actual verification boundary |
| --- | --- | --- |
| Kalman/ETS/ARMA/transfer | 78 groups / 6,262 assertions | Actual production exported C ABI after covariance repair; [record](tasks/evidence/jmp-state-space-public-reference.json) |
| Functional preprocessing | 2,194 calls / 143,946 scalar assertions | Exact-source isolated C ABI, plus public Rust/native checks; arWLS177 converged/23 nonconverged retained; [record](tasks/evidence/jmp-functional-preprocess-reference.json) |
| Numeric token counts/weighting/LSA/varimax/LCA | 980 calls / 59,109 fields | Exact-source isolated wrapper plus public Rust/native checks; algebra error and EM fixed-point tolerance reported separately; [record](tasks/evidence/jmp-text-analysis-reference.json) |
| Spatial covariance/variogram/known-V GLS/kriging | 601 calls / 76,681 fields | Exact-source dense SciPy QR oracle, reviewed numerical boundaries, plus public/native checks; [record](tasks/evidence/jmp-spatial-reference.json) |

Earlier16-operation references (counts use different units):

| Earlier16-operation domain | Independent checks | Layer / evidence |
| --- | --- | --- |
| Kinship ML/REML / BLUP | 199 cases / 3,282 fields plus 3 valid-dimension validation rejects, including repaired Decimal/extreme cases | Isolated production-source numerical wrapper plus public/native gates; [record](tasks/evidence/jmp-kinship-reference.json) |
| 1PL/2PL/3PL IRT | 43 groups / 2,385 assertions | Production exported C ABI; [record](tasks/evidence/jmp-irt-public-reference.json) |
| CFA / factor scores / fit indices / RAM | 543 cases | Exact production-source numerical wrapper plus public/native gates; [record](tasks/evidence/jmp-sem-reference.json) |
| Degradation / stability / cumulative exposure | 820 scenarios / 66,812 assertions, 69 valid-fit Hessians | Exact production-source numerical wrapper plus public/native gates; [record](tasks/evidence/jmp-degradation-reference.json) |

Earlier28-operation independent references (their counts use different units):

| Added domain | Independent checks | Boundary / evidence |
| --- | --- | --- |
| Repeated covariance / random slopes | 630 cases | Production C ABI; [record](tasks/evidence/jmp-mixed-advanced-reference.json) |
| Bernoulli / Poisson GLMM | 31 groups / 551 assertions | Production C ABI and independent real-line integration/HiGHS; [record](tasks/evidence/jmp-glmm-reference.json) |
| Binary causal inference / uplift evaluation | 51,062 scalar assertions | Production C ABI; [record](tasks/evidence/jmp-causal-reference.json) |
| Assay normalization / admixture | 556 cases / 22,948 fields | Production exported C ABI; [record](tasks/evidence/jmp-omics-reference.json) |
| Classical blocked / split / strip / Taguchi designs | 586 cases | Isolated numerical wrapper plus public execute and native SQL; [record](tasks/evidence/jmp-design-classical-reference.json) |
| Market baskets / item reliability | 340 cases / 14,842 fields | Production exported C ABI, exhaustive enumeration and NumPy; [record](tasks/evidence/jmp-market-reference.json) |

The previous continuation references remain separate:

| Domain | Independent checks | Execution layer / evidence |
| --- | --- | --- |
| Latent/clustering | 4,369 scalar assertions | Actual production public-source snapshot compiled as a C ABI; [record](tasks/evidence/jmp-latent-reference.json) |
| Evaluation/splits | 1,140 cases / 22,209 fields | Exported production C ABI; [record](tasks/evidence/jmp-evaluation-reference.json) |
| Nonlinear/smoothing | 18,040 scalar assertions; 80 supplied-start fits, 80 automatic-start checks, 24 basis, 72 spline and 24 LOESS fits | Exported production C ABI; [record](tasks/evidence/jmp-nonlinear-reference.json) |
| Functional/spectral | 1,100 cases | Isolated Rust module wrapper plus separate public-crate tests; [record](tasks/evidence/jmp-spectral-reference.json) |
| Multivariate inference | 5,937 scalar assertions | Actual production public-source snapshot C ABI; [record](tasks/evidence/jmp-latent-inference-reference.json) |
| SVM/quantile/Huber | 2,062 scalar assertions / 60 fits; all oracle fits met their reported convergence criteria | Exported production C ABI; [record](tasks/evidence/jmp-robust-learning-reference.json) |
| Diploid markers | 853 cases / 137,076 fields | Exported production C ABI; [record](tasks/evidence/jmp-genomics-reference.json) |
| Reliability | 780 cases, including AFT information, survival/quantiles, NHPP and all-state block enumeration | Isolated Rust module wrapper plus separate public-crate tests; [record](tasks/evidence/jmp-reliability-reference.json) |

The nonlinear oracle deliberately preserves 2 of 80 supplied-start and 4 of 80
automatic-start `converged:false` results at tolerance `1e-12`; their admissible
numeric results are distinct from convergence certification. Earlier independent
references remain available for
[regression/multivariate/mixed](tasks/evidence/jmp-model-reference.json),
[inference/time-series/Cox](tasks/evidence/jmp-inference-timeseries-reference.json),
[distribution/ARIMA](tasks/evidence/jmp-distribution-arima-reference.json) and
[additional DOE](tasks/evidence/jmp-advanced-design-reference.json).
Check counts use different units and are not added into a JMP coverage percentage.
No installed JMP was used as a numerical oracle. Unfinished method families and
options remain explicit in the feature matrix.

<!-- analytics-validation-summary:end -->

For reproducible local checks, use the same task-owned build directory:

```sh
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-jmp-target cargo test --workspace
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-jmp-target cargo clippy -p ducksmiles_analytics --all-targets -- -D warnings
```

<!-- analytics-full-reference:end -->

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

Wasm toolchain, standard-library rebuild, and artifact checks:
[WebAssembly build contract](WASM_BUILD.md).

- **Rust** (6 crates): Molecular parsing/computation and a separate numerical analytics core
- **C++**: DuckDB extension registration, typed SQL wrappers and Rust FFI
- No external chemistry library dependencies (no RDKit, no OpenBabel)
- Vectorized execution via DuckDB's `UnaryExecutor`
- Invalid input returns `NULL` — never crashes

## License

MIT License
