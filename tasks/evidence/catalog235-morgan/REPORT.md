# Morgan densities and aromatic preparation — v22

2026-10-18. **211/217 registered RDKit standards +24 profile fields =235**.
Requested meaningful additional50: **0/50**. The overall goal is unfinished.
Revision: `ducksmiles-workbench-v22-macro235-rdkit2026036`.
Definitions and primary references: [NEXT-MORGAN](../../descriptors/NEXT-MORGAN.md).

## What changed

Added `FpDensityMorgan1/2/3`: unique unfolded sparse Morgan keys at each radius,
divided by heavy atoms. One radius3 pass supplies all three values. No bit-vector
folding or chirality is used. Initial invariants and complete radius0..3 sparse
key sets match 352 independently sanitized reference graphs exactly.

Workbench preparation now handles tested charged-radical aromatic cycles and
default-valence candidate rules. Standard HBA/HBD nitrogen query charge rules
were corrected. A held-out macrocyclic polyphenol exposed an existing ether-
bridge aromaticity error: neutral degree2 O/S in a first-encounter ring of nine
or more members cannot donate two aromatic electrons under the pinned rule.
The raw parser/native-graph contract remains separate. This is not a complete
RDKit sanitizer or an exact complex-ring/basis implementation.

## Independent chemical checks

- Pinned oracle: RDKit2026.03.6, sanitize=True, source removeHs=False.
- Runtime RED controls preceded fixes; failing logs are retained here.
- Full release Rust: **766passed,0failed,1existingignored** before the subsequent
  unregistered QED test was added. No new Rust warnings; scoped formatting passed.
- C ABI352-source set: all235 scalars plus formula pass strict comparison.
- Expanded26-source aromatic set: 132 differing scalar cells before the final
  macrocycle fix, zero afterward; all26 formulas also match.
- Reference tolerances unchanged; Morgan densities use1e-9 (observed max error0).

## Actual Swift App code: chemically distinct10k

Mac host, actual `App/Chemistry.swift` and store, **no UI rendering**.
Original names/SMILES/order,235 numeric columns,239 CSV columns and CSV reimport
pass independent validation: **2,342,433 numeric cells checked**. All25
unsupported rows remain with blank numeric values. CSV22,179,931bytes.

The new three Morgan columns match on **9,975/9,975 accepted inputs**, max error0.
The value gate passes, but the complete gate **fails** on25unsupported inputs.
This is a three-column comparison, not an all235 whole10k parity claim.

All10k warm results are reused with zero additional native calls. The last full
10k disk-capacity check wasv21:21,912,427bytes,10k exact disk restores/zero new
native calls. v22 verifies the newly corrected macrocycle by actual Swift disk
reload; the same literal fails withv21 and passes withv22. The actualv22 bulk
run checks10k in-memory reuse. These scopes are deliberately distinguished.

Only one actualSwift row changes fromv21→v22; the other9,999 outputs are exact.
All235+formula on the changed row improve51differences→0,51resolved/0introduced.
Reports retain source/output/library hashes before old binaries/raw runs were
removed. Latest frozen dylibSHA256:
`65d65e535b572ea7f8c2ce16d1f9e11a950e3a26bfc715d6e3cc65f2720cd2d8`.

Observed Mac import91.575s,warm0.164s,export/snapshot33.507s under concurrent
comparison/system load. These are not controlled before/after or RDKit speed
comparisons and are not physical-device timings. A one-second symbolicated Mac
CLI sample identifies BCUT spectral work and Bertz string formatting among
costs; it contains no iPad frame-time evidence.

## Full100k result identifies its older engine

The completed full100k run is **v20/catalog232**,208 standards+24 profiles:
99,585 all232+formula match;167 accepted rows differ;248 unsupported. Strict
gate fails. Five formula residuals remain. See
[saved results](../catalog232-coordination/chebi100k-v20/comparison.json).
On the common202-field/common-acceptance cohort versusv15,482→162 differing
rows,320resolved/0introduced value rows; **three newly unsupported rows** remain.
It is not a currentv22/all235 result.

## App and remaining work

Latest iPad Air11M4/iOS26.5 Simulator build and four targeted portrait UI tests
passed. Seven PNGs exported; the three new Morgan/radical/macrocycle screenshots
were inspected at original resolution. Normal Japanese launch restores the six
saved inputs. This does not resolve the older rotation or physical-device gates.
Seven language files have127 matching keys/placeholders; wording remains draft.
Physical iPhones are unavailable. Numeric filtering, proper chemical2D,
iPhone/current landscape/10k rendered UI and the older rotation failure remain
open. Six standards remain: QED,SPS,two stereo counts,two IPC forms. Additional
50 meaningful descriptors and whole-catalog current100k verification remain
required. No3D or new runtime dependencies were introduced.
