# Chemob / DuckSMILES v38 — validated CX wedge inputs

2026-10-18. Revision `ducksmiles-workbench-v38-cxwedge291-rdkit2026036`.
Frozen host library SHA256 `8cd3c5fcd76bab10c570f2f4b812a63640806b91c90ff0fe1201630057d9398a`.
Catalog remains **291 = 217 RDKit-derived definitions + 24 profiles + 50 additional**.

## Input support improvement

The prepared descriptor API now accepts a bounded, coordinate-free CX subset:
`wU`/`wD` wedge lists with valid unique bond references and incident atom indices.
Original strings remain App/import/cache/export identities. The current descriptors
do not measure axial handedness. This is not general CXSMILES interpretation,
axial stereochemical depiction, or support for a stereochemical graph sidecar.
Unknown, wiggly, radical, coordinate and other CX blocks are refused. Dative wedge
edges are excluded because reversing them could alter chemistry.

The raw native graph API keeps its previous grammar. No annotation is silently
discarded from the stored SMILES. Future configuration-sensitive descriptors need
an extended representation before using this projection.

The implementation was preceded by a genuine runtime RED: all 482 expected
supported controls were rejected. GREEN matches every declared descriptor and
formula against ORIGINAL annotated references for all 482 controls: 107 retained
population inputs plus 375 synthetic controls spanning endpoints, closures,
existing tetrahedral/double-bond stereo and atropisomer wedge combinations.
[Reference manifest](cx-wedge-reference.json), [runtime RED](cx-runtime-red.log),
[GREEN](cx-green.log). Full release Rust workspace: **794 passed, 0 failed,
2 ignored, 0 compiler warnings**. [Workspace summary](v38-workspace.json).

## Population evidence and its scope

The current whole-source **100000-input preparation/acceptance** run accepts
99859 inputs, leaving 141 unsupported. Exactly 107 previously unsupported
positions now pass; zero previously accepted positions regress.
[All positions and transitions](v38-acceptance100k.json).

A separate current C ABI comparison of all 107 new original annotated inputs
matches all 291 fields and formula within declared tolerances, with zero
reference errors. [Incremental comparison](v38-cx-positive/comparison.json).
Pinned IPC numerical compatibility differs on two of these rows, separately kept.

The earlier v37 whole100k numerical run remains a frozen baseline: all 99752
commonly accepted rows matched 291 references. **A whole100k v38 numerical rerun
has not been claimed.** The current change adds an input projection; it does not
change descriptor formulas or tolerances. [Frozen baseline](REPORT.md).

## App input/output and cache

SMI paste previously treated a CX annotation as a compound name. It now keeps
its full block before an optional name. Unknown, malformed and additional blocks
remain available to the engine instead of becoming accepted bare prefixes.
CSV/TSV comma lists and engine row-injection bounds are also tested. Swift Core:
**17 tests passed** after actual runtime RED and a mixed-space/tab regression.

Actual Swift Workbench processing of the original 10000 distinct compounds
passes import, source names/SMILES/order, calculation, CSV output/reimport and
warm reuse. Unsupported rows decrease from 25 to 16. Output has 295 columns,
2801641 finite numeric cells, and `Ipc_stable`/`AvgIpc_stable` headers.
Warm reuse restores all 10000 rows with zero new native calculations.
[Actual I/O summary](v38-bulk-io-summary.json),
[independent CSV parser/value checks](v38-swift10k-csv-values.json).

The actual cache test also loads a real frozen37 unsupported CX result and
recomputes it as valid under38, preserves duplicate reuse, and reopens the new
positive from disk. All earlier chemistry/revision/cache controls pass.
[Actual cache integration](v38-cache-integration.log).
The capacity-wide disk test is the separately retained v37 proof; do not relabel
that run as v38. No physical-device speed claim is made.

## Remaining work

141 unsupported source inputs remain, including resource-limited graphs and
non-tetrahedral/ambiguous metal coordination. A read-only non-tetrahedral
experiment shows why stripping @SP/@TB/@OH is unsafe: ten examples change SPS,
and six also change unspecified stereocenter counts. Correct geometry/stereo
semantics and canonical metal coordination require further implementation.
[Diagnostic experiment](v38-non-tetra-base-experiment.json).

The two stable information descriptors remain mathematically validated alternatives
to pinned RDKit numerical results; compatibility limits from37 remain visible.
The strict input-acceptance completion gate has not been achieved. Physical phones
remain unavailable, so actual device runtime and speed remain unmeasured.

## Simulator interface verification

iPad original-source paste/CX inline values and the original10000-row CSV flow
both pass. The actual GUI produces a29132263-byte CSV **byte-identical** to the
actual Mac output; an independent parser verifies all2801641 finite cells, source
strings/order, blanks and status. The visible system share view is captured, with
no sending/copy/save action.
[UI summary](app-v38/ipad-summary.json),
[whole-file parity](app-v38/simulator-export-byte-parity.json),
[value checks](app-v38/simulator-export-291-values.json).

Visual review exposed rounded negative zero (`-0.0014` displayed as `-0.00`).
An actual iOS runtime RED records four assertion failures before correction.
All four shared number-formatting tests, including German/Swiss/French decimals
and scientific values, and the corrected CX GUI test then pass. The displayed
rounded zero is normalized after locale-specific rounding; native/CSV values
are not modified. [Current5-test gate](app-v38/rounded-negative-zero-green-summary.json).

The iPhone17e Simulator sample/detail/input/reuse flow passes with the same final
App. A Release configuration mismatch initially requested an x86_64 slice despite
the aarch64 native build script; the project now declares arm64. Actual RED build
and subsequent green builds are retained. This is an Apple Silicon Simulator
build, not a claim of Intel Simulator support.
[Phone gate](app-v38/phone-summary.json).

The original normal Japanese iPad Simulator has been updated to38 without test
arguments. All six saved IDs/names/SMILES are semantically preserved. The revised
number display and inline lower panel are verified in the original image.
[Current normal Japanese screen](app-v38/original-ipad-japanese-six.png).
Simulator timings are uncontrolled observations, not physical-device benchmarks.

## Device build and retained artifacts

The final unsigned iPhoneOS arm64 App builds successfully with minimum iOS17,
iPhone/iPad families, all four required native exports, version0.1.0(1), and
three matching source notices. The ZIP also includes the root MIT license.
This proves build/link integration, not physical runtime. Both paired phones
still have unavailable tunnels/DDI services.
[Build proof](app-v38/iphoneos-build-proof.json),
[unsigned App ZIP](app-v38/Chemob-0.1.0-v38-unsigned-iphoneos.zip),
[device availability](app-v38/physical-device-availability.json).

No commit, push or release was made. Shared formatter tests are included in the
iOS test target as well as the standalone Swift Core package. Final source
archive includes all Rust workspace members, Core manifest/sources/tests, App
sources/resources/configuration, public10k test data/attribution/licenses, and
reference tools. It does not bundle the test data into the production App.
