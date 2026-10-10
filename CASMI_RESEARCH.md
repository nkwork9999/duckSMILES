# CASMI molecular-feature research APIs

Experimental Rust APIs and streaming command-line examples, added against real
CASMI26 train/test. These are not replacements for the stable `morgan_fp_bits`
SQL function. No SQL/FFI registration, release, or model migration is implied.

## Native graph fingerprints

`morgan_experimental::fingerprints(smiles, width)` returns six sparse count
vectors: count2, count3, cycle2, shell2, cut2, hybrid2. They preserve isotope,
charge and aromatic labels, count repeated atom-layer environments, and vary
additional graph context. Raw @/@@ is deliberately not hashed. Stereo/CIP and
RDKit bit parity are not implemented. FNV64 and folding still permit collisions.

`fingerprints_v2(smiles, width)` returns four two-channel vectors: unchanged
count3 local features, followed by cycle, distance-shell, bridge-composition or
bridge-neutral-mass context. First channel indices are `[0,width)`; second channel
`[width,2*width)`. Apply square root, normalize each channel, then concatenate and
normalize. Default experiment uses width2048 and equal channel weights. The mass
channel is a structural hypothesis, not a predicted ion spectrum: it preserves
parent atom hydrogen counts, does not add caps or model rearrangements, and skips
unsupported or isotope-containing fragment masses. Input limit: 1..256 atoms.

The first globally mixed shell implementation performed poorly. Keeping local
features separate from bridge-mass context improved the fixed CASMI pilot; this
does not establish universal superiority over RDKit or competition performance.

## Paper-backed MinHash

`minhash::MinHash::new(dimensions, seed)` and `encode(&[u32])` implement the
affine MinHash equation in [Probst & Reymond 2018](https://doi.org/10.1186/s13321-018-0321-8),
also used by [MAP4](https://doi.org/10.1186/s13321-020-00445-4).
`similarity` counts equal coordinates. Never apply cosine to hash integer values.
Input tokens must already be distributed hashes (e.g. SHA-1 first four bytes),
not consecutive IDs. Empty sets yield the all-u32::MAX sentinel. Coefficient RNG
is versioned SplitMix64 and differs from the authors' NumPy/tmap implementations;
saved reference sketches are not bit-compatible.

The CASMI harness extracts MHFP6/MAP4 shingles with RDKit and feeds token hashes
to this native Rust compressor. MHECFP uses unfolded RDKit Morgan hashes. The
Rust compressor is implemented here; an independent native canonical rooted
SMILES/CIP implementation has NOT been claimed or shipped.

## Run

```sh
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-casmi-target cargo test -p ducksmiles_smiles --offline
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-casmi-target cargo build -p ducksmiles_smiles --examples --release --offline
```

`casmi_fingerprints` reads one SMILES per line and emits six sparse vectors as
JSON. `--v2` emits four two-channel vectors; `--legacy` emits the existing Morgan
r2/2048 bit vector in sparse form. Invalid inputs produce `null`, maintaining
alignment. `paper_minhash 256 42` reads comma-separated u32 tokens per line and
emits comma-separated MinHash coordinates. `--coefficients` as the third argument
prints the coefficient table for independent verification.

Experiment code, protocols and results are in the separate local repository
`casmi26-data-foundation/experiments/DM001/` through `DM004/` (a separate
research checkout; these experiment artifacts are not bundled here).
The harness keeps raw Parquet unchanged, excludes all400 public-test identities
from reference spectra, and records molecular splits and candidate coverage.
Only the small explicitly described pilot was fit/evaluated; existing Kaggle
neural models and official structure normalization were not changed.

## MS/MS features and additional50 descriptors

`ms_features::calculate(smiles, width)` (`duck-ms18-v1`) returns18 sparse count
channels, parent composition/properties, and connected fragments from one- and
two-bond cuts. Channels cover fragment masses/formulas, bond and cut environments,
neutral losses, acid/base proxies, ring cuts, nested paths, group distances,
heteroatom placement, graph topology and native properties. Two-cut fragments
must have both removed bonds on their actual boundary. Parent hydrogens are
retained; severed bonds are not capped. These are structural hypotheses, not
validated reaction products, BDE or gas-phase pKa estimates.

`paper_features::calculate(smiles, width)` (`duck-paper50-v2`) returns50 separately
defined sparse count vectors in `paper_features::NAMES` order. They add fragment
mass defects/Kendrick coordinates, element retention, oriented boundary patterns,
site distances, formula ambiguity, nested losses, nonlocal environments and
functional-group triangles. FIORA, FraGNNet, MIST, GLACIER, TDiMS and other primary
papers motivate these definitions; this does not reproduce50 published models.
The functions are custom combinations, not RDKit bit-compatible fingerprints.
RDKit graph primitives could be used to implement analogous functions.

Both APIs require one connected net-neutral molecule with supported elements
H/C/N/O/S/P/F/Cl/Br/I/B/Si/Se, no isotope labels, at most128 atoms and192 bonds.
Hash widths64..65536 are supported; experiments use2048. Unsupported inputs return
errors. No new dependencies, SQL bindings or FFI exports were added. Native
HBA/rotatable/LogP approximations differ from RDKit and are not used by the adopted
DM003 method. Version2 requires the acid OH, attached carbonyl C and its own
double-bonded O to survive together before emitting a CO2-loss hypothesis.

The `ms_features` CLI reads one SMILES per line and writes one JSON object with
`channels`, `properties`, `formula`, `mass`, and compact fragment records. The
`paper_features` CLI emits an array of50 sparse vectors; an optional0-based channel
number emits just that vector wrapped in a one-element array. Both emit a JSON
error object for unsupported rows and preserve row alignment.

DM003 improved a fresh480-molecule candidate-ranking pilot: Top1 37.3%→72.3%,
MRR .5008→.8309 against RDKitMorgan3+fragment mass. The selected pipeline combines
mass and cleavage-environment fingerprint prediction with direct exact-mass
fragment matching; this is not a Morgan-only algorithm improvement. Additional50
channels were individually compared on600 development molecules. Two selected
channels improved development but failed fresh confirmation on a different480:
MRR .84836→.84073, paired difference95% interval [-.01732,.00224]. They remain
available for research and are not adopted. The final method stays DM003.

See `experiments/DM003/RESULTS.md`, `experiments/DM004/RESULTS.md`,
`DM004/REGISTRY.md`, and `DM004/EFFECTS_50.csv` in the foundation repository for
definitions, primary sources and every channel's measured effect. Public-test
predictions exist for1,213 rows;959 protonated rows use the adopted method,254
other-adduct rows use explicitly marked baseline fallbacks. No Kaggle leaderboard
improvement or universal advantage over RDKit has been established.

## MS-FINDER hydrogen-rearrangement fingerprints (DM005)

`hr_features::calculate(smiles, width)` returns ten sparse channels in
`hr_features::NAMES` order: ion mass, ion formula, neutral-loss formula,
oriented cut environment, and nested ion pathway, for positive then negative mode.
Version `duck-hr10-v1` adapts the initial/subsequent hydrogen rules from
[Tsugawa et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC7063832/).
It is not the complete MS-FINDER ranker, a BDE model, or a probability model.

Use one/two-edge fragments with original parent H counts, single-bond CNOPS
boundaries, and zero formal charge on every parent atom. Existing element,
isotope, connected-component and size limits remain. Reject unsupported parents;
count unsupported fragment boundaries. Ion formulas conserve precursor elements,
with the H budget for [M+H]+/[M-H]-; subtract polarity times electron mass.
Nested paths retain the first charged endpoint and require a subsequent ±1H step.
No empirical Table1 frequency weights or external trained model are imported.
FNV64 folded counts are structural hypotheses and can collide; they are not
stereochemical/CIP fingerprints or certified reaction products.

`hr_features` reads one SMILES per line and emits a JSON object. `--positive`
returns only the five positive channels; `--audit` adds ion/fragment receipts.
Unsupported input emits an error object and preserves row alignment.
No SQL/FFI registration or new dependency was added.

All ten channels passed 100×5 random-SMILES invariance and independent RDKit
reconstruction on40 structures/12,500 ions. The positive-mode CASMI development
screen chose ion mass, yielding MRR .844366→.848015 and top1 444→449/600.
The development diagnostic95% interval for the difference includes0; independent
confirmation has not run. Keep DM003 adopted. Work stopped at the user's request;
see the foundation repository's `experiments/DM005/RESULTS.md` for that checkpoint.

## HR ion context and direct spectrum descriptors (DM006)

Work resumed on2026-10-10. `ion_context::calculate(smiles, width)` returns24
new fingerprints per polarity (`duck-ion-context48-v1`), inherited HR ion
receipts, and directed ion-state paths. See `ion_context::NAMES` for the order.
The channels encode cut count, H shift, formula ambiguity, conserved losses,
mass defects, oriented cut neighborhoods, retained functional groups and paths.
Positive channels occupy0..24; negative channels occupy24..48. Parent domain
and error handling remain those of DM005. Width2048 is used in the pilot.

`ion_context::match_smiles(smiles, mode, precursor, mz, intensity)` returns eight
direct descriptors (`duck-hr-match8-v1`), in `MATCH_NAMES` order. They measure
observed peak/intensity coverage, theoretical mass precision, peak F1, and
coverage from one-cut, two-cut, ring and linked-path hypotheses. Use mode1/-1
for protonated/deprotonated ions. Matching is nearest peak within max(.005Da,
20ppm), with lower-mz tie resolution. Exact duplicate observed masses are merged
before sqrt-intensity normalization; peaks at/above precursor−.02Da are excluded.
Theoretical precision merges masses at1e-5Da, while observed coverage is counted
once per peak. Empty spectra return zeros; invalid parameters return an error.

The `ion_context` streaming CLI takes one SMILES per line; `--positive` emits
24 channels and `--audit` also emits ion/path receipts. `--match` instead reads
five tab-separated fields: polarity, precursor, comma-separated mz,
comma-separated intensity, SMILES. Unsupported rows produce an error JSON object.
No dependencies, SQL functions or FFI exports were added.

These are32 custom definitions motivated by MS-FINDER, FIORA and FraGNNet, not
trained implementations of those models. RDKit graph primitives could implement
equivalent descriptors. Independent RDKit checks cover100×5 SMILES variants,
40 molecules/12,500 ions/12,218 paths and80 direct-match cases in both polarities.
Workspace release tests:638 passed,1 ignored. New files have no clippy diagnostics;
the repository's existing FFI errors and unrelated formatting differences remain.
The positive-mode600-molecule development screen is recorded in the foundation
repository's `experiments/DM006/RESULTS.md` and `EFFECTS_32.csv`. DM003 remains the
adopted method until a separately frozen policy passes fresh confirmation.
