# Next registered descriptor family: VSA / ASA

Implemented in v13: **all 44 entries**, bringing standard coverage to159/217.
Fixed references/full100k and selected iOS portrait verification are complete;
existing chemistry differences remain, so full parity is not achieved. Additional50
remains0/50. See `../evidence/catalog183-original-h/REPORT.md`.

## First tranche: 44 standards without a charge model

LabuteASA (1), SMR_VSA1..10 (10), SlogP_VSA1..12 (12),
EState_VSA1..11 (11), VSA_EState1..10 (10). After these, standard coverage would
be 159/217, not 217. Keep all current 139 IDs/order stable and append new IDs;
bump the native/cache revision, update coverage, dictionary groups and CSV.

Pinned primary definitions:
- [RDKit MolSurf.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Descriptors/MolSurf.cpp)
- [RDKit MolSurf.py](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/rdkit/Chem/MolSurf.py)
- [RDKit EState_VSA.py](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/rdkit/Chem/EState/EState_VSA.py)

Use the pinned PeriodicTable Rb0 radii and exact bond corrections; Labute's
include-H correction is its published algorithm, not an explicit-hydrogen
solvent simulation. EState_VSA bins surface contributions by EState, while
VSA_EState bins EState values by surface contribution. Preserve exact boundary
inclusion and explicit-H behavior from independent references.

## Shared dependencies and tests

1. Generate literal independent per-atom surface/Crippen fixtures and scalar
   references with the pinned wheel. Include isotopes, explicit H, ions,
   disconnected graphs, cages, heavy elements, fused heterocycles and values on
   both sides of every bin boundary. Observe missing-feature runtime failures.
2. Native Rust computes radii/ASA contributions once per graph. Reuse the
   existing EState vector. Crippen exposes original-atom contributions
   and expanded-H LogP/MR totals. Explicit-H inputs require an additional original-
   graph typing context because those definitions differ; independently verify these atom
   contributions before consuming them in surface bins.
3. Aggregate each family from shared arrays; undefined stays null. Avoid another
   parse, ring perception, hydrogen expansion or atom-typing pass for each bin.
4. Run full workspace, prior-error set and independent full-population gates.
   Repeat actual 10k Swift I/O/CSV and category UI on the enlarged catalog.
5. PEOE_VSA14 is a separate tranche requiring a validated Gasteiger charge model;
   it must not be supplied with arbitrary zero charge placeholders. BCUT charge
   variants can subsequently reuse the same validated charge array.

## Additional50 remains a separate requirement

Existing `paper_features.rs` has 50 research feature channels for MS. These are
hashed sparse vectors, not 50 validated scalar app properties. Do not silently
count channel names or summed hash buckets as completion of the requested50.
Reuse them only with explicit formulas, applicability, independent tests and a
useful app representation; do not invent experimental prediction claims.
