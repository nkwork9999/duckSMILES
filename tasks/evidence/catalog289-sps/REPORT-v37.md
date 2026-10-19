> Latest input/App checkpoint: [v38 report](REPORT-v38.md). Preparation accepts
> 99859/100000;141 unsupported remain. The following is the frozen **v37 whole
> numerical** baseline, with its original source/library hashes and failed gates.

# Chemob / DuckSMILES v37 — descriptor verification

2026-10-18. **Current status: implementation and integration verified on the supported scope; the strict acceptance-parity gate remains FAILED (248 inputs).**

Native revision `ducksmiles-workbench-v37-valence291-rdkit2026036`.
Frozen macOS library SHA256 `3e2fa996604baba34d8a8b0b38f5c77a64434fed18a2b4da7457cf16224cb882`.

## Implemented scope

| Category | Count | Reference |
|---|---:|---|
| Registered RDKit-derived 2D definitions | 217 | 215 pinned numerical implementations + 2 stable mathematical definitions |
| Existing profile fields | 24 | Independent graph/property extraction |
| Additional published descriptors | 50 | Independent dense formulas; absent from the audited pinned descriptor API |
| Total exported native values | 291 | Catalog order retained; two new information fields appended |

The 50 additions comprise 19 molecular-distance edge descriptors, 10 topological
charge indices, 16 walk/closed-walk counts and 5 distance/connectivity indices.
[Definitions and limits](../../descriptors/EXTRA50-DESIGN.md),
[coverage](../../descriptors/coverage.json).

## Whole 100000-input result

[Machine-readable comparison](v37-full100k/comparison.json),
[compact audited summary](v37-full100k-summary.json).

- 100000 distinct original SMILES retained, with source positions checked exactly once.
- 99752 accepted: all 291 fields and molecular formula agree within declared tolerances.
- Zero descriptor mismatches, zero reference computation errors, zero formula mismatches.
- 248 RDKit-valid inputs remain unsupported. Their indices are unchanged from frozen35.
- The overall strict gate is **FAILED** for these 248 acceptance differences.
- The gate for values on commonly accepted inputs is **PASSED**.

Undefined values are checked as undefined and export as blank cells. They are not
silently replaced with zero. This verifies computational definitions against the
specified software/oracles; it does not establish experimental property accuracy.

### Information descriptors and model compatibility

Ipc/AvgIpc use exact integer characteristic-polynomial coefficients and binary64
entropy. Independent Python arbitrary-width integer full matrices verify the
native signed-limb/symmetric implementation. 58 determinant/analytic graph controls
and 85-digit Decimal entropies, 34 original molecular controls, 20000 signed i128
arithmetic controls and 270-digit carry/borrow/division boundaries also pass.

Both stable fields match their independent reference on all 99752 accepted inputs.
The declared absolute/relative limits are Ipc 1e-10/1e-12 and AvgIpc 1e-12/1e-12.
No existing descriptor tolerance was widened. The weighted Ipc can be enormous,
so its absolute error must be interpreted with its relative limit.

Pinned RDKit 2026.03.6/NumPy 2.5.4 gives different numerical results for these two
fields on 5617 accepted inputs. This is a separately recorded **FAILED numerical
compatibility gate**. For example, the 128-carbon path has exact Ipc≈1.55013e27;
the pinned result is about 112.4 times larger. Stable AvgIpc≈3.80581 versus≈1.14759.

The App labels these fields as stable and explains model compatibility. CSV column
names are **Ipc_stable** and **AvgIpc_stable**. Native catalog IDs remain Ipc/AvgIpc.
Models trained with the pinned numerical features need feature compatibility
validation. The 215 other registered standards use pinned RDKit values directly.
[Exact definitions](../../descriptors/ipc-definitions.json),
[pinned differences](v37-full100k/pinned-information-differences.jsonl.gz).

## Invalid-input correction

Frozen36 incorrectly accepted overvalent carbon/nitrogen/oxygen/fluorine inputs.
A prepared-only read-only explicit-valence guard now matches pinned Atom.cpp,
including charged P/S/As/Se and the two-coordinate hydride exception. The raw graph
API keeps its existing grammar/contract. The cache revision was advanced.

All 14603 unfiltered atom/ligand/aromatic/coordination controls match acceptance and
formula: 12290 accepted, 2313 rejected. These controls check acceptance/formula,
not every descriptor on all synthetic isolated-ion examples. Runtime RED and
GREEN evidence are retained. Whole 100k preparation still accepts 99752 and retains
the same 248 unsupported positions. A real frozen36 falsely accepted cached C5
entry is recomputed/rejected by the actual Swift37 store.

## Actual App input/output

Original 10000 distinct source compounds passed actual Swift Workbench import,
calculation, CSV output and reimport. Source names/SMILES/order and all 25 unsupported
rows are retained. The CSV has 295 columns and 2799118 finite numerical cells.
Warm reuse performs zero additional native calculations for all 10000 inputs.
The current v37 actual disk store also saves and reopens all 10000 distinct rows:
29042424 cache bytes, identical results/source order, and zero new native calls.
The 25 unsupported rows remain represented.
[Current disk-cache proof](v37-cachebulk10k.json).

The actual iPad Simulator UI imports the original CSV over an existing six-row
library, finds source row 9999 and shows its formula and Wiener 3360 inline. It
writes a 29109614-byte CSV with the two stable headers. Its entire data body is
**byte-identical** to the Mac Swift37 output. An independent Python CSV parser
checks every finite value, undefined blank, status, name and source position.
CSV SHA256 `37b11280a9939789d9ceba249490e276322d354ae3fde4816af18cf2febc2335`.
[Actual Simulator CSV checks](app-v37/simulator-export-291-values.json),
[byte comparison](app-v37/simulator-export-byte-parity.json).

The six-row system share/cancel flow checks the visible file title, cancellation,
and preservation of all six inputs across table pages. The system share UI may
follow the device language independently of the App language. Full-screen capture
is used for the remote activity view. Earlier locator/page assumptions and the
one discarded test launched after a failed duplicate-resource build are retained
explicitly; they are not accepted final proofs.

## Interface and platforms

- iPad six-flow suite passes: French decimals in table/details, actual portrait and
  landscape layout, cached numeric ranges/undefined-vs-zero, corrected rejection
  versus valid PF6, stable information/compatibility text, and original10k CSV.
- iPhone 17e/iOS 26.5 Simulator sample/search/detail/invalid-input/paste/reuse passes.
- Focused visible share/cancel flow passes. Original PNGs inspected.
- Original iPad normal no-argument Japanese launch preserves all six saved input
  records. JSON key serialization changed; all IDs/names/SMILES remain identical.
- Seven draft languages have 221 matching keys/placeholders:en,ja,id,hi,de,fr,it.
  Swiss decimal conventions use the interface locale. Native-speaker review remains.
- Version metadata 0.1.0(1) and three source notices are included in the App.
- The unsigned iPhoneOS arm64 build passes, targets iPhone/iPad with iOS 17 minimum,
  and contains the required four native exports. Physical runtime was not tested.
  [Build proof](app-v37/iphoneos-build-proof.json),
  [unsigned build archive](app-v37/Chemob-0.1.0-unsigned-iphoneos.zip).
  This archive needs signing before installation on a physical device.
- Full release Rust workspace 792 passed/0 failed/2 ignored; the external100k
  preparation test was invoked separately and passes. One older ignored test remains.
- Swift Core 14 tests pass. No Swift/compiler warnings attributable to current code.
  Xcode's no-AppIntents metadata messages are recorded separately.

[App evidence and original screenshots](app-v37/),
[latest normal iPad](app-v37/original-ipad-japanese-six.png).

## Dataset and conditions

ChEBI structure registry release 256, with source/selection hashes retained.
Source 215662 rows; 215636 have SMILES. SHA ordering scans 110202 candidates,
skips 4113 duplicate canonical isomeric graphs, 327 RDKit-invalid and 5762 generic
(dummy-atom) graphs, and retains 100000 definite distinct compounds. Original source
SMILES are not replaced by the identity key. Salts, charge, isotopes and stereo
remain. Selection uses no native acceptance filter. The existing Kaggle/CASMI
foundation was audited read-only; its reuse rights were not established, so it was
not used for this public fixture. Dataset provenance/rights were recorded earlier.
[Selection manifest](v37-source-selection.json).

Reference: Python 3.13.3/RDKit 2026.03.6/NumPy 2.5.4, macOS 26.2 arm64.
Four disjoint 25000-row ranges share frozen tools/library; merged counts, policies,
source positions and difference records are checked. This is a correctness run,
with no throughput claim. Archived earlier incomplete36 runs are explicitly discarded.

## Remaining limits

Supported input is bounded to 256 atoms, 512 bonds and 4096 UTF-8 bytes per structure.
There are 102 oversized-atom inputs in the unsupported cohort (eight also exceed
512 bonds) and 146 other unsupported grammar/preparation cases, chiefly CX wedge
annotations/non-tetrahedral coordination. Read-only diagnostics flag 108 CX-annotated
inputs and 36 non-tetrahedral tokens (flags overlap with resource limits); three
small metal-cluster inputs remain unclassified by these syntax flags. Maximum
size in the retained unsupported cohort is 598 atoms / 642 bonds. These flags are
a worklist, not proof of the only failing preparation stage.
[Diagnostic source positions](v37-unsupported-diagnostics.json). Their original records remain available.

The current picture is a labelled connectivity preview. Proper chemical 2D
coordinates, stereochemical depiction and SDF import are future features. 3D is
outside this descriptor milestone. Physical phones are paired but unavailable at
verification time; physical-device speed and runtime accuracy remain unmeasured.
Mac/Simulator timings taken during other verification work are observations only.
Earlier controlled same-output spectral optimization and pinned-workflow timing
have their own version/scope in the historical report.

[Historical checkpoints](REPORT-history-before-final-v37.md).

## Reproducible checkpoint

[Source archive](source-v37-final.tar.gz) contains all five Rust workspace members,
Swift Core package manifest/sources/tests, App sources/resources/benchmarks, project
configuration, and licenses. Reference tools and environment hashes are retained
with the full comparison. No commit, push or release was made.

Task-owned build products and raw UI bundles were removed after success/failure
summaries and selected original screenshots were retained. The two owned test
Simulators were shut down; the original Japanese iPad remains on the normal v37
App with its six saved records.
