# Chemob / 139-entry checkpoint

## Scope

115/217 registered RDKit 2026.03.6 2D scalar descriptors, 24 native profile
fields. Additional 50 requested descriptors: 0/50. Not complete.

New: 85 registered functional-group SMARTS counts with capped atom-set-unique
matching, and four EState extrema with one per-atom calculation. Fixed SMARTS
bond Boolean precedence and explicit-H D degree. Independent goldens: 8,330
fragment counts + 1,230 nonfragment scalar values. BSD notices retained.

## Reproducible evidence

- `comparison.json`, `differences.jsonl.gz`: **v8** 139-entry native library,
  SHA256 `049e6ba0d039f0e3032ca30928578e2b75a90fc85064611929d8161ebef47c04`.
  100,000 distinct ChEBI concrete graphs; common accepted 99,755, all139+formula
  matching 98,784; accepted differences 971, native unsupported 245, formula
  differences six. Strict exit 1. `scope-delta.json`: same 1,216 differing-row
  set as old50, not proof each new descriptor matches globally.
- `bulk-io-summary.json`, `csv-validation.json`: **v8** actual Swift Workbench,
  10,000 inputs, 25 unsupported retained, 143 columns/139 numeric, 1,386,525
  numeric cells checked independently. Source text/order and CSV reimport passed.
  8,410,398 bytes. Cold 31.760s, warm 0.231s, no additional native calls.
  macOS host; excludes UI. Workload/load differs from prior checkpoints.
- `bulk-io-before-limit-fix.json`, `wide-csv-red.txt`: genuine 8 MiB reimport
  failure. Raised to 32 MiB; UTF-8 byte scanner preserves quoting/newlines/
  multilingual text. `swift-core.txt`: 11 release tests passed. Synthetic
  267-column test demonstrates format capacity, not chemical implementation.
- `rust-workspace-prepared.txt`: **v9** parser-ring-sharing refactor,
  735 passed / zero failed / one ignored.
- `parser-ring-reuse-timing.json`: **v8 vs v9**, 1,000 distinct source SMILES,
  999 supported, identical catalog and output; 139 fields, formula, JSON/FFI,
  one thread, no result cache. Three alternating trials, 2.940→1.764s median.
  After SHA256 `de4ae9a7bec812e35da8e251bce568279e51b9ef80048d5911958d82159c2bf9`.
  Mac observation with timing variability, no RDKit/device speed claim.
- `source-sha256.json` records the v8 feature checkpoint;
  `prepared-source-sha256.json` records v9 including parser/native_graph.
  Do not attribute the completed v8 full comparison to the later v9 build.

## Open work

102 standard descriptors, 50 meaningful additional descriptors, graph/aromaticity
and sanitize parity, property filters, proper 2D coordinates, full 10,000-row
cache persistence/UI, iPhone/device verification. New iPad portrait feature test
passed on iPad Air 11-inch M4 / iOS 26.5 Simulator (one test, zero failures),
with functional-group and EState screenshots visually inspected; older full test rotation failure remains unresolved.

The adopted ChEBI source selection and license documents are retained under
`../catalog50-chebi-concrete/`. Kaggle/CASMI foundation was inspected read-only
and not adopted for product use. App runtime remains native Rust with no Python
or RDKit runtime dependency. No 3D.
