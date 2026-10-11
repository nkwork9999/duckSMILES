# Formula and adduct mass context

Three additive SQL/Rust APIs implement `duck-formula-adduct-v1`. They do not
alter SMILES parsing, fingerprints, structure identity, or existing mass APIs.
No external package or runtime chemistry dependency is introduced.

```sql
SELECT formula_exact_mass('C8H10N4O2');              -- 194.0803755791 Da
SELECT adduct_precursor_mz(194.0803755791, '[M+H]+'); -- 195.08765203127 m/z
SELECT adduct_neutral_mass(195.08765203127, '[M+H]+'); -- 194.0803755791 Da
```

`formula_exact_mass(VARCHAR) -> DOUBLE` sums **neutral atomic** masses.
It accepts complete flat molecular formulas using
`C H N O S P Cl Br F I Si Se Na K B As Fe D`; `D` is deuterium.
An omitted count means one. Counts must be positive decimal integers without
leading zeroes. Repeated elements such as `CH3COOH` are summed. Inputs must be
ASCII, at most 4,096 bytes and at most 1,000,000 total atoms.

This versioned 18-symbol numerical table matches the public AFIX formula
routine on supported neutral flat formulas. It is not an average-weight table
or a claim of bitwise RDKit `ExactMolWt` parity. The existing
`mol_exact_mass(smiles)` table and behavior remain unchanged. Charged structures
and salts require an explicit, scientifically appropriate neutral formula;
the function cannot infer that choice.

`adduct_neutral_mass(DOUBLE, VARCHAR) -> DOUBLE` and
`adduct_precursor_mz(DOUBLE, VARCHAR) -> DOUBLE` use:

```text
delta = mass of added neutral atoms − mass of lost neutral atoms − z × electron_mass
precursor_mz = (n × neutral_mass + delta) / abs(z)
neutral_mass = (precursor_mz × abs(z) − delta) / n
electron_mass = 0.00054857990 Da
```

Here `n` is the number of neutral molecules and `z` is signed ion charge.
The complete supported set requires these exact spellings:

| Class | Labels |
|---|---|
| Positive, one M | `[M+H]+`, `[M+NH4]+`, `[M+Na]+`, `[M+K]+`, `[M-H2O+H]+`, `[M-2H2O+H]+`, `[M+2H]2+`, `[M]+`, `[M-H2O]+`, `[M+CH3OH+H]+`, `[M+CH3CN+H]+` |
| Negative, one M | `[M-H]-`, `[M-H2O-H]-`, `[M+CH2O2-H]-`, `[M+C2H4O2-H]-`, `[M+Cl]-`, `[M]-`, `[M-2H]2-`, `[M+Na-2H]-` |
| Two M | `[2M+H]+`, `[2M+Na]+`, `[2M+NH4]+`, `[2M+K]+`, `[2M-H]-`, `[2M+CH2O2-H]-`, `[2M+C2H4O2-H]-`, `[2M+Na-2H]-` |
| Three M | `[3M+H]+`, `[3M-H]-` |

These are arithmetic hypotheses, not validation that the molecule forms the ion
or has the atoms needed for a neutral loss. Do not silently relabel observed
adducts. Numeric inputs and outputs must be finite and positive.

SQL returns NULL for NULL, unsupported or invalid arguments. Rust functions in
`mass_context` return `Option<f64>`; the C ABI exports return NaN. FFI checks null
pointers, byte limits and UTF-8. Non-null pointers must still identify a readable
allocation of the supplied length.

| Rejection | Examples |
|---|---|
| Unknown symbol | `Xx`, `He` (not in this version's table) |
| Unsupported formula syntax | `(CH3)2`, `CuSO4.5H2O`, `[13C]H4`, `NH4+`, whitespace |
| Count or size | `C0`, `C01`, overflow, over 1,000,000 total atoms |
| Unknown/ambiguous adduct | `[M-2H]-`, `[M+H]`, `[M+Li]+`, aliases or spaces |
| Numeric value | NULL, NaN, infinity, zero, negative, overflow or nonpositive result |

## AFIX indexing and BIO candidate structures

The public **CASMI26: Metric Fusion Consensus Backfill** notebook's cell 9 embeds
`casmi_engine.py`. Its `formula_mass` computes neutral mass from
`molecular_formula`; `load_library` prioritizes it over precursor/adduct
arithmetic. This mitigates precursor rounding and adduct-label errors in the
**reference-library mass index**. It does not rewrite measured m/z, correct the
adduct label, or identify a structure. An explicit SQL view retains derivation:

```sql
CREATE VIEW reference_mass_context AS
WITH masses AS (
  SELECT *,
         formula_exact_mass(molecular_formula) AS formula_neutral_mass,
         adduct_neutral_mass(precursor_mz, adduct) AS precursor_neutral_mass
  FROM reference_metadata
)
SELECT *,
       coalesce(formula_neutral_mass, precursor_neutral_mass) AS index_neutral_mass,
       CASE WHEN formula_neutral_mass IS NOT NULL THEN 'formula'
            WHEN precursor_neutral_mass IS NOT NULL THEN 'precursor_adduct'
            ELSE 'unresolved' END AS mass_basis,
       'duck-formula-adduct-v1' AS mass_spec_id
FROM masses;
```

Keep raw formula/adduct/precursor m/z, source ID, version and original checksum.
Formula precedence is an indexing policy, not proof that the supplied formula
is correct. Unknown queries must not receive formulas from withheld answers.
Training/reference exclusions still apply before evaluation uses this view.

**BIO is a candidate source, not a descriptor.** It means ChEBI + LIPID MAPS,
added to COCONUT and then training structures after upstream key deduplication.
Retain source memberships even after canonical key deduplication. The supplied
6,930-bit vectors select RDKit Morgan r2/r3, RDKit path and MACCS bits; they are
not interchangeable with DuckSMILES Morgan vectors. This extension neither
downloads these databases nor imports their trained feature models. CASMI's
fixed RDKit scoring normalization remains unchanged.

Intentional differences: the public regex can accept partial formulas, whereas
this API rejects them. The public table gives `[M-2H]-` charge magnitude two;
this API requires `[M-2H]2-`. Sodium minus two hydrogens in a negative ion uses
one added electron, including `[2M+Na-2H]-`; the public `Na - 2*PROTON` expression
adds two. These cases are not called legacy parity.

The audited source snapshot was already present at
`docs/public_code_P001/sources/fusion/casmi26-metric-fusion-consensus-backfill.ipynb`.
Its SHA-256 and native verification are in `MASS_CONTEXT_VALIDATION.json`.
Notebook attribution identifies inherited prvsiyan components and candidate
dataset `prvsiyan/chebi-lipidmaps-casmi26`. No leaderboard gain is claimed.

## Validation

Six new Rust tests cover nitrogen-containing formulas, encoding/count errors,
all adducts, signed electron corrections, multimers, overflow, nonpositive results
and the C ABI. The workspace passes 716 tests, with one existing ignored test.

The actual SQL extension is tested with constants, 4,097-row vectors, filtered
rows and mixed NULLs. The independent verifier compares 1,022 formulas (including
1,000 deterministic nitrogen-containing compositions) to decimal arithmetic and
29 adducts at three masses to a grammar-based oracle in both directions. It also
executes all 150 catalog examples.

```sh
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-afix-target cargo test --workspace --offline
# Build using the existing DuckDB extension setup matching the host version.
/path/to/existing/python scripts/check_mass_context.py \
  --extension /path/to/ducksmiles.duckdb_extension --receipt /path/to/verification.json
```

Verified: macOS arm64, DuckDB 1.5.6, installed Rust/Clang and existing DuckDB
headers. No dependencies were installed or updated. Windows/Linux/Wasm were not
rebuilt. Repository-wide clippy reports existing C ABI lint errors; the added
module has no diagnostics. The new file is rustfmt-clean; existing repository
formatting differences remain.
