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
uniqueness of the resulting images. It also proves that every duplicate-free
enumeration of those images gives the same capped count. A regeneration check pins the compiled
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

## Comparison with registered RDKit functions

The [column inventory](rdkit-columns.csv) covers all 291 names. The
[observed differences](rdkit-differences.json) record the two input preparation
policies and the independent chain controls. On the 226 nonempty, parsable
inputs, Lean compared 49,042 outputs from all 217 registered default 2D
scalars using the catalog tolerances. Of these, 48,980 were finite, 60 were
nonfinite and two SPS calls raised `ZeroDivisionError` on hydrogen-only inputs.
All values agreed after mapping undefined results to descriptor null. The
other 24 graph summaries used public RDKit molecular properties; the extra 50
used independent formulas on RDKit graphs. Separately, the empty SMILES
parses to a zero-atom RDKit molecule, and DuckSMILES returns SQL NULL.

* **Hydrogen handling.** With `removeHs=False`, ordinary explicit hydrogen
  vertices remain in the graph. With RDKit's default SMILES parsing, four of
  six controls changed, affecting 136 values across 50 descriptor names.
  For `[H]O[H]`, Chi0 changed from 2.7071067811865475 to 0, and
  FpDensityMorgan1 from 4 to 1. Ethanol and isotope-labelled water were
  unchanged within tolerance.
* **Information content.** The registered Ipc and AvgIpc use the absolute
  characteristic-polynomial coefficients of the adjacency matrix. RDKit
  evaluates the recurrence in floating point. Lean and DuckSMILES obtain the
  coefficients using integer arithmetic before computing entropy. The two
  paths met tolerance for unbranched carbon chains of 20, 40 and 64 vertices.
  On 80 vertices, RDKit reported AvgIpc 3.466866198166465 and Lean
  3.4668641308764405. The exact path recurrence gives constant coefficient
  1, while the pinned RDKit numerical recurrence returned approximately
  −2.5746 billion on the recorded machine. The independent analytic check
  uses the path recurrence, Pₙ(x) = xPₙ₋₁(x) − Pₙ₋₂(x). Lean proves that
  its constant term on every path of size 2k is (−1)ᵏ, and tests the full
  order-80 molecular matrix independently.
  [RDKit Ipc](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/rdkit/Chem/GraphDescriptors.py)
  and [characteristic polynomial](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/rdkit/Chem/Graphs.py).
* **Undefined SPS.** The normalized spacial score divides by the heavy atom
  count. On `[H]` and `[H][H]` the pinned Python API raises an exception.
  DuckSMILES represents this coordinate as JSON null and preserves the rest
  of the vector. [RDKit SPS](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/rdkit/Chem/SpacialScore.py).

The audit runs each registered function on a fresh copy of the molecule. Lean
compares the returned numbers or nulls and records the ID of every mismatching
coordinate. The measured agreement is conditional on preprocessing, inputs,
registered version and tolerances.

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

To reproduce both RDKit input policies and the chain controls, run
`uv run --python 3.13 --with rdkit==2026.3.6 python compare_rdkit.py --output /tmp/rdkit-audit --library /path/to/libducksmiles_smiles.so`.
The output includes a per-descriptor CSV and the individual differences.

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
