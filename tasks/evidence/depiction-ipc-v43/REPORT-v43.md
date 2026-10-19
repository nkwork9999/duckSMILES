# Chemob v43 follow-up — local verification; physical device pending

2026-10-18. Core, chemistry, cache and stereo gates pass. Simulator Release
build-for-testing and unsigned iPhone arm64 Release build pass. Final UI
verification has passed the focused Simulator gates below. No commit, push, publication or release has occurred.
The immutable [v42 descriptor report](../catalog289-sps/REPORT-v42.md) remains valid
for its declared scope; v43 does not change descriptor arithmetic.

## Prepared drawing graph

`ds_workbench_drawing_graph` uses the same preparation pipeline as descriptors.
Legacy `ds_mol_graph_json` is unchanged. Input atom order and original SMILES are
preserved. Drawing schema `ducksmiles-workbench-drawing-v1`; graph cache has its
own namespace; layout namespace is `depiction-2d-v4`. Calculation revision is
still `ducksmiles-workbench-v42-cluster291-rdkit2026036`.

All 10,017 source positions (original first 10,000 plus 17 curated) were retained,
without selecting by native acceptance or drawing success. Compared 10,003 graphs
against RDKit 2026.03.6: atom symbol, charge, isotope, attached H, aromaticity,
bond type/connectivity and directed dative donor→acceptor all match (0 differing
rows). Remaining 14 inputs exceed the 256-atom drawing budget. Calculation uses
its separate 640-atom/1,280-bond budget. The former raw-graph 28 chemistry
differences are resolved by preparation, without altering the legacy raw API.

[Chemistry audit](prepared-final/chemistry-audit.json),
[prepared fixtures](prepared-final/prepared-graphs.jsonl.gz),
[audit runner](../../descriptors/audit_prepared_drawing.py).
Input/data rights remain those recorded in the v42 ChEBI release256 notice.

## 2D geometry and stereo

Original 10,000: 9,374 drawings, 612 quality rejections and 14 size rejections.
Of available originals, 261 explicitly flag some omitted stereo. All 17 curated
cases draw; 1 flags omitted stereo. Rejected drawings do not discard source or
calculated descriptors. Complex rings/bridges can remain unavailable; this is
not universal ChemDraw/RDKit drawing coverage.

Independent recovery from actual screen coordinates and wedge/hash direction
matches 39,089 reference tetrahedral centers and 3,593 reference stereogenic
alkenes, with 0 differing rows. One additional declared alkene annotation is
non-stereogenic in the reference and is counted separately. Omitted configurations
are excluded from these parity counts. Non-tetrahedral/CX limitations are flagged.
The earlier one-row E/Z fan regression is fixed and included in the full sweep.
Deliberately reversed wedge and E/Z controls are detected by the verifier.

[Layout summary](prepared-final/layout/summary.json),
[unavailable inputs](prepared-final/layout/unavailable.json),
[stereo round trip](prepared-final/stereo-roundtrip.json),
[E/Z negative control](ez-oracle-negative-control.json).

Timing is diagnostic only: a shared, memory-constrained Mac host, not controlled
hardware or physical iPhone measurement. Pure layout median 2.04 ms, p95 122.31 ms,
maximum 11.57 s (252 atoms); all 10,017 inputs plus fixture output took 587.61 s.
No owned compiler ran during this sweep; other host workload was uncontrolled.
Large first-time drawings still have a latency tail. Drawing work runs off the
main actor; both successful and unavailable results are reused. These figures
must not be presented as a release speed estimate or a new RDKit speed comparison.

## Calculation and cache regression

- Rust workspace: 811 tests passed, 0 failures, 0 warnings.
- Swift Core: 33 tests passed, 0 failures.
- Prepared API RED failures were observed before implementation; subsequent
  focused tests cover cleanup, direction, unsupported stereo and buffer/UTF-8
  contract, including unchanged short buffers and drawing budgets.
- Existing 71 Rust sources match frozen v42 exactly, except the new module
  declaration. Native descriptor/catalog/revision JSON bytes match v42 on
  1,041 bounded regression inputs. This is not a new full100k descriptor audit.
- Actual native→Swift cache probe passes attached-H labels, wedge/E/Z metadata,
  duplicate and disk reuse with no native/drawing work, and unavailable reuse.
- All 7 dictionaries have 224 matching keys/placeholders; translations are drafts.

[Workspace test summary](workspace-tests-summary.json),
[calculation source identity](calculation-source-identity-resumed.json),
[ABI byte parity](resumed-descriptor-abi-parity.json),
[cache proof](prepared-cache-proof.log).

## Ipc / AvgIpc

The 5,723 legacy numerical differences are retained explicitly. Exact-integer
stable coefficients match independent FLINT controls and remain invariant;
pinned RDKit floating Ipc changes approximately3.88x across equivalent atom
orders in source22. No tolerance widening or legacy compatibility mode is claimed.
CSV fields remain `Ipc_stable` / `AvgIpc_stable`.
[Diagnostic](information-order-diagnostic.json).

## App and physical-device status

Simulator Release build-for-testing and unsigned iPhone arm64 Release build passed.
The Simulator App executable is unchanged across the UI-harness refresh.
No new production layout or descriptor changes were made during final UI verification.

- Compact iPhone SE3 / iOS26.5 (375pt): chemical 2D examples and E/Z/bridge tests passed
  in the initial run. The property-column test failed because its diagnostic accessed
  an offscreen row before scrolling (run exit65, 2 passed/1 failed). The harness now
  scrolls before querying lower rows and preserves its geometry assertions.
- Final compact focused run: 2 passed/0 failed, exit0. Five chemical examples also
  check that property numbers and units do not overlap. Lower integer properties
  align with the initial right column within2pt; labels retain at least8pt separation.
  The separate E/Z/bridge pass used the identical App binary.
- iPad Air11 M4 / iOS26.5: 3 passed/0 failed, exit0. Chemical 2D, E/Z, bridged carbon,
  explicit partial-stereo notice and portrait/landscape inline results were exercised.
- Selected original-resolution exported screenshots were inspected. This is focused
  flow evidence, not a complete accessibility or all-language release audit.
- Normal Japanese iPad launch restored. All six saved input objects match the
  pre-test backup; normal UI shows six reused results, the 2D drawing and inline values.

Both paired phones remain unavailable (no tunnel/DDI). Physical execution and
physical speed remain unmeasured. No pairing/signing settings changed. Build time,
OS bootstrap time and XCTest runner installation time are not chemical calculation
benchmarks. Source changes remain uncommitted; no push, publication or release.

Build proof: [Simulator](app-simulator-build-proof.json), [iPhone arm64](app-iphone-build-proof.json).
UI evidence: [initial compact log](compact-resumed.log),
[final compact summary](compact-layout-check-summary.json),
[final compact screenshots](compact-layout-check-attachments/manifest.json),
[iPad summary](ipad-resumed-summary.json),
[iPad screenshots](ipad-resumed-attachments/manifest.json),
[normal Japanese screen](normal-ipad-final.png).

## Preserved artifacts and cleanup

Verified App ZIP CRCs and executable SHA256 against platform build proofs. Source
archive includes native/App source and tests, excluding build products.
[Artifact manifest](artifacts-v43.json), [cleanup record](final-build-cleanup.json),
[UI verification scope](ui-verification-scope.json).
Task-owned DerivedData, Apple Cargo outputs and Generated platform libraries were
removed after archive verification. Installed Japanese App and input data remain.
The existing AppIntents metadata extraction skipped warning remains; no Swift/Rust
compilation error was reported. Rebuild before additional automated UI tests.
