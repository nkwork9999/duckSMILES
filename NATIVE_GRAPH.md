# Native molecular graph sidecar

`mol_graph_json(smiles VARCHAR) -> VARCHAR` exports a deterministic JSON object
for a native DuckSMILES molecular input. The Rust API is
`native_graph::calculate(&str) -> Option<String>`; the C ABI is
`ds_mol_graph_json` in `src/include/ducksmiles.h`.

This is an input/descriptor sidecar for GLACIER, ICEBERG and other molecular
workflows. It is **not** their pretrained tensor format and does not replace
their RDKit featurizers, fragmentation graph, formula vocabulary, model
weights, adduct encoding, or collision-energy encoding. Arbitrary native
features cannot be appended to pretrained tensor inputs; a versioned adapter
and separately evaluated retraining are required for any such experiment.

```sql
SELECT mol_graph_json('N[C@@H](C)C(=O)O');
SELECT json_extract(mol_graph_json('[NH4+]'), '$.atoms[0].formal_charge');
-- 1
```

## JSON v1 contract

Top-level fixed strings are `schema="ducksmiles-native-graph-v1"`,
`backend="ducksmiles-native"`, `backend_version="0.3.0"` (Rust crate version),
`atom_order="input_parse_order"`, and
`hydrogen_policy="preserve_explicit_vertices_no_expansion"`.
`input_smiles` retains the exact accepted input string. `atom_count` and
`bond_count` count vertices and undirected bonds, including explicit H vertices.

`pretrained_tensor_compatible=false`. The fixed `computed` map has
`native_aromaticity=true` and `cip`, `ez`, `gasteiger`, `hybridization`,
`ring_basis`, `radical_electrons` all `false`. A raw token is not a computed
stereochemical assignment. Other existing native descriptors may be joined as
separate, versioned features; this API does not silently substitute them for
model features.

Each atom contains:

| Field | Meaning |
|---|---|
| `index`, `symbol`, `atomic_number` | 0-based parse index, element symbol and Z (1–118) |
| `isotope`, `atom_map` | Written isotope mass number and atom-map label, or JSON null |
| `formal_charge`, `aromatic`, `in_bracket` | Parsed charge, native-perceived aromaticity, bracket origin |
| `bracket_h_count` | Written H count in brackets, e.g. 4 for `[NH4+]`; not H vertices |
| `implicit_h_count` | Native valence-completed H count for unbracketed atoms only |
| `explicit_h_neighbor_count` | Number of adjacent H vertices, including labelled isotopes |
| `chirality_token` | Written `@`/`@@`, or null; no CIP R/S assignment |
| `stereo_neighbor_order` | Native parse-neighbor reference order for chiral atoms; `-1` denotes bracket H; empty without chirality |

Each bond contains `index`, `begin`, `end`, `bond_type` and `direction`.
`bond_type` is `single`, `double`, `triple` or `aromatic` after native aromatic
perception. One bond appears once; it is not duplicated as directed GNN edges.
`direction` is `1` for `/`, `-1` for `\`, and `0` when unspecified, read from
`begin` to `end`. Reverse traversal negates the sign. These are directional
input tokens, **not** resolved E/Z labels. Explicit single bonds between
aromatic atoms remain single where supported by the native parser.

## Atom identity, H and chemistry boundaries

No canonicalization, salt stripping, tautomer selection, charge neutralization
or H expansion/removal occurs. `O` has one vertex; `[H]O[H]` has three.
`[NH4+]` stores a bracket H count and one N vertex. Alternative SMILES can give
different atom orders and different JSON for the same molecule. Mapping a
SMILES through a canonicalizer or an H-removing toolkit requires a new atom
mapping; do not reuse previous graph indices. Native aromaticity perception
does change bond/atom aromatic flags from the written Kekule form.

Cache by exact input SMILES, exact graph JSON bytes, schema/backend version and
the extension/code SHA. The crate version alone does not distinguish builds.
An official CASMI IK14 may intentionally merge structures differing in H or
charge and is insufficient for model-input cache identity. Official scoring
continues to use its frozen RDKit implementation, independently of this API.

## Rejection and limitations

SQL NULL propagates. Invalid UTF-8, non-ASCII input, parse failure, unsupported
syntax or resource limits return SQL NULL / Rust None / C -1 consistently.
The sidecar is limited to 4,096 input bytes, 256 atom vertices and 512 bonds.
All 118 element symbols have a Z mapping; wildcard `*`, `[D]` (use `[2H]`), and
unknown elements are rejected. Recognition of an element is not a claim of a
validated valence, isotope mass or pretrained-model vocabulary for it.

The API-specific guard rejects dangling/repeated bonds, empty branches,
leading/trailing/repeated dots, dots inside branches, self/duplicate edges,
overflowing isotope/map annotations, and explicit bond/direction annotations
immediately before ring labels. The latter are a conservative supported-subset
restriction: this parser does not retain all opening-ring bond annotations.
Ordinary ring labels (including `%10`) and bond annotations between atoms are
supported. Other SQL functions keep their existing parser behavior unchanged.

This is not full RDKit sanitization or valence/radical/stereochemistry checking.
Native parser acceptance and aromaticity can differ from RDKit; record parity
failures, do not silently repair/relabel them. Unknown models' exact hydrogen,
atom-order, aromaticity, Gasteiger, ring basis, hybridization and stereo feature
requirements remain unverified. A caller must retain a failure receipt when
SQL returns NULL; there is no per-reason SQL diagnostic API in v1.

The C caller owns buffers. Non-null pointers must address valid, non-overlapping
ranges. A null or undersized output buffer returns the required byte count
without writes. A sufficiently sized buffer receives exactly that many bytes,
without a trailing NUL. No Rust allocation crosses the ABI; panics are caught.

## Verification and reproducible local build

No dependencies, model weights, external data or pretrained tensors are added.
Use the pinned DuckDB headers/runtime documented by `MASS_CONTEXT_VALIDATION.json`
and one task-owned `CARGO_TARGET_DIR` for release tests/build. Build recipe:

```sh
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-mspred-target \
  cargo test --workspace --release --offline
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-mspred-target \
  cargo build --workspace --release --offline
```

The C++ link and DuckDB metadata steps are the existing MASS_CONTEXT recipe,
with the target path above and extension version `native-graph-v1`. Run
`scripts/check_native_graph.py --extension /path/to/ducksmiles.duckdb_extension`
with an existing matching DuckDB Python. It loads the actual extension and
checks JSON, H/charge/isotope/stereo/input order, malformed inputs, vectorized
SQL, registered catalog examples, and all graph SQL regression queries. It
does not measure model accuracy or claim RDKit/pretrained parity.
