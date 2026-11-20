# Molecular descriptor verification

DuckSMILES exposes the workbench's complete 291-column molecular profile to DuckDB.
The extension executes the existing Rust implementation. The offline Lean 4.19
executable checks its serialized results with the independently computed reference.
The catalog IDs, order, absolute and relative tolerances come from the repository's
[descriptor catalog](../../tasks/descriptors/phase1-catalog.json).
Per-column reference, evaluator, tolerance and test coverage appear in [coverage.csv](coverage.csv).

```sql
SELECT mol_descriptor_catalog_json();
SELECT mol_descriptor_values_json('CCO');
```

Element and graph summaries are calculated from sanitized observations in Lean
(28 columns). For 85 fragment descriptors Lean counts unique, capped SMARTS
embeddings supplied by RDKit. The remaining 178 numerical results come from
the project's independent oracle, including the exact-coefficient definitions
for Ipc and AvgIpc and dense graph formulas for 50 additional descriptors.
RDKit 2026.03.6 supplies molecular perception with explicit H retained.
The trusted observations, catalog and expected inputs belong to the verifier.
The candidate supplies only its results. Every column uses the same Lean check.

SQL NULL represents an invalid, unsupported or missing molecular input. A
JSON null *inside a profile* represents one undefined descriptor. For defined
numbers the absolute difference must be at most the larger of the absolute
tolerance and relative tolerance times the absolute reference. Lean compares
serialized decimal values by integer scaling. It rejects the wrong number
of rows or columns and changed input order. The kernel checks equivalence to
the declared integer comparison, presence agreement, full-coordinate coverage
and invariance of element counts under atom-list permutation. Molecular
perception, floating-point algorithms, serialization and the Rust binary have
finite differential-test coverage. Materials composition and viscosity APIs
are documented separately in [MATERIALS.md](../../MATERIALS.md).

## Reproduce

From this directory with the pinned Lean toolchain on PATH:

```sh
lake build
lake env lean Tests.lean
uv run --python 3.13 --with rdkit==2026.3.6 python reference.py CCO > /tmp/obs.json
.lake/build/bin/descriptor_check reference ../../tasks/descriptors/phase1-catalog.json /tmp/obs.json
(cd ../.. && cargo rustc --locked --release -p ducksmiles_smiles --lib --crate-type cdylib)
uv run --python 3.13 --with rdkit==2026.3.6 python integration.py ../../target/release/libducksmiles_smiles.so
```

On macOS use the .dylib artifact. If the bundled linker generates an executable
the system loader rejects, build with `LEAN_CC=/usr/bin/clang` and
`LIBRARY_PATH="$(lean --print-prefix)/lib"` and set MACOSX_DEPLOYMENT_TARGET
to the installed SDK's supported macOS version. For SQL integration add
`--with duckdb==1.5.5` to uv and pass `--extension /path/to/ducksmiles.duckdb_extension`
built for that exact DuckDB ABI. Add `--coverage coverage.csv --timings timings-m1.csv`
to write per-column results and the six timed trials. The inputs cover eight
valid molecules and two rejected strings. Five independent mutations exercise
error detection. These controls supplement the project's larger existing
[reference comparisons](../../tasks/evidence/catalog289-sps/REPORT.md).

## Timing

The [pipeline trials](timings-m1.csv) time a native Rust FFI batch, a SQL query
including result conversion, and a fresh Lean process reading prepared
observations and rendering JSON. Those paths include different work.
The [count trials](counts-m1.csv) use the same prepared atom labels and query
order in both languages, excluding input parsing, I/O and process startup.
Both programs check the same checksum. Run `lake build bench_counts` then
`uv run --python 3.13 --with rdkit==2026.3.6 python bench_counts.py` on the target host.
