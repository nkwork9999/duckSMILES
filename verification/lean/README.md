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

Lean computes all 291 output coordinates from prepared observations. An unknown
catalog ID fails the run. The observation file carries the molecular graph,
atom and bond perception, atomic contributions, initial Morgan invariants and
the eight QED component properties. RDKit 2026.03.6 prepares these inputs with
explicit hydrogen vertices retained, and separately prepares the hydrogen
representations required by Crippen, QED and BCUT. The independent numerical
oracle runs on fresh molecules during integration testing.

For the 85 functional groups Lean evaluates the compiled SMARTS, including
recursive queries, and counts distinct sorted vertex sets with a cap of 1000.
The kernel proves soundness and completeness of the enumerator for its
predicates, correctness of anchored recursive queries and membership and
uniqueness of the resulting images. A regeneration check pins the compiled
queries to the published SMARTS. A separate fixture covers all 85 with positive
examples across 98 molecules.

The remaining algorithms include the 50 heavy-graph indices, atomic VSA
binning, connectivity and shape factors, Morgan environment propagation,
characteristic-polynomial entropy, Burden matrix eigenvalue extrema and the
QED desirability function. The kernel relates each walk-count row to explicit
walk endpoint multiplicities, verifies exact trace divisions, proves the bin
index at the first strictly larger threshold and proves the partition of a
perceived ring class. The coverage table names the applicable result for
each column. General equivalence of the Rust binary with the Lean model,
chemical perception, component libraries and floating-point error bounds
remain open proof obligations. The executable comparisons test finite inputs.

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
uv run --python 3.13 --with rdkit==2026.3.6 python generate_fragments.py --check
uv run --python 3.13 --with rdkit==2026.3.6 python check_fragments.py
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
to write per-column results and the six timed trials. Pass `--extended` to
check the 228-input corpus used in CI and in the coverage table. Timings and
mutation controls use the smaller batch of eight valid and two rejected inputs. Five independent mutations exercise
error detection. These controls supplement the project's larger existing
[reference comparisons](../../tasks/evidence/catalog289-sps/REPORT.md).

## Timing

The historical 0.7.1 [pipeline trials](timings-m1.csv) time a native Rust FFI batch, a SQL query
including result conversion, and a fresh Lean process reading the earlier oracle-backed
observation format and rendering JSON. Those paths include different work.
The [count trials](counts-m1.csv) use the same prepared atom labels and query
order in both languages, excluding input parsing, I/O and process startup.
Both programs check the same checksum. Run `lake build bench_counts` then
`uv run --python 3.13 --with rdkit==2026.3.6 python bench_counts.py` on the target host.
