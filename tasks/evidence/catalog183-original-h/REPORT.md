# Chemob / v13 surface and original-H checkpoint

Catalog183 =159/217 registered RDKit standards +24 native profiles. Additional50:
0/50. Remaining58 standards and50 meaningful additional descriptors; not complete.

## Added44 standards

LabuteASA1; SMR_VSA10; SlogP_VSA12; EState_VSA11; VSA_EState10. Shared original
parse/rings, one Labute area array and one EState array. Crippen provides expanded-
H molecular totals and original-input atomic contributions. Explicit-H inputs
need two representation contexts; all surface bins reuse those arrays.
Pinned source definitions and BSD attribution are recorded in
`../../descriptors/NEXT-VSA.md`, `surface-definitions.json` and license notices.
No3D, Python or RDKit runtime dependency in the app.

## Genuine regression and correction

v12 full100k found **one new row**, index95722 `[2H]O`, on four surface bins.
Expanded-H typing correctly computes molecular LogP/MR but can change the type
of an original explicit H vertex. v13 types the original graph for those inputs,
while retaining expanded-H molecular totals. Atomic and combined-total regression
RED logs are retained. `[H]O` and five additional partial-H controls were added.

## Verified native / reference evidence

- Rust workspace **740 passed /0 failed /1 existing ignored**, no warnings.
- **138 molecules ×44 =6072 scalar references**, **4485 atomic references**;
  prior1230 nonfragment and8330 functional-group cells remain covered.
- Exact bin edges and either side tested; undefined per-field values remain
  null/dash/blank, not arbitrary zero. The independent verifier accepts matching
  RDKit nonfinite/native null, but rejects replacing defined ethanol ASA with
  null. v12 evidence records the verifier's original failing undefined test.
- SMR_VSA8 and SlogP_VSA9 cannot be positive with the pinned Crippen parameter
  table. Their registered zero bins are retained for compatibility; they are
  not counted as additional50. Other42 have positive fixture examples.
- `regression/comparison.json`: five explicit-H/zero-radius controls, all183
  values+formula match; strict passed.
- `nci/comparison.json`:4999 inputs/4900 raw unique;4991 common accepted,
  4986 all183+formula matching. Five value differences +eight acceptance
  differences remain; strict fails as before.

## Complete100k comparison

`chebi/comparison.json`:100,000 distinct source graphs,99,755 common accepted,
**99,223 all183+formula matching**,99,238 standard+formula matching.532 accepted
value/formula differences +245 unsupported remain; six mass/formula differences.
Strict exit1: full compatibility is **not achieved**.
`chebi/scope-delta.json`: no new/resolved differing rows versus v11's139 fields,
and every legacy field's mismatch count/max error is unchanged. v12's newly
introduced `[2H]O` row95722 is resolved. The added44 can still differ on graphs
already differing in the older engine; absence of new rows is not full parity.
Frozen Mac dylib SHA
`e0d985deff543a205f5bda55c724133dbbc7b53c4018da59931f0fb24e7b511d`.

## Actual Swift I/O/cache

`bulk-io-summary.json`, `csv-validation.json`: v13 Workbench10,000 distinct
inputs,25 unsupported preserved,187CSV columns/183numeric,1,825,425 numeric
cells independently checked. Original text/order and CSV reimport passed.
CSV14,122,923 bytes. Mac cold read/import/calculation18.028s, warm0.157s,
10,000 reused/no additional native calls. Timings overlapped correctness work;
not a fair version/RDKit comparison and not physical-device/UI speed.
`prefix-output-equality.json`: on this10k subset, prior139 fields and legacy
outputs match, and v12/v13 all183 values match. The changed case lies outside
this10k subset and is covered separately and in the full100k comparison.

`cache-bulk.json`: actual ChemistryStore save/reopen at10,000 distinct inputs,
14,048,974-byte cache, all10,000 reused, no additional native calculations,
all restored values and row order equal. Fits existing32MiB bound.

## iPad reflection

v13 arm64 iOS Simulator build passed on iPad Air11-inch M4/iOS26.5. Two selected
UI tests passed: native ASA/SMR/SlogP/EState surface values and page navigation;
actual pasted `[2H]O` input → molecular LogP−0.82 and SMR_VSA1=6.9097 /
SMR_VSA2=0.0000. Search remains hittable and no detail modal is needed.
Three PNGs were inspected; original-resolution inspection confirms full numeric
strings (resized previews were hard to read). `ipad-ui-summary.json` /txt/png
retain the result. Seven languages have116matching resource keys. Standard
scientific titles remain their canonical Latin names; translation wording is
still draft.

The older landscape rotation failure, physical-device speed, current iPhone,
large-input UI, numeric filters, proper chemical2D layout and full parity remain
open. The preview is a labeled connectivity diagram, not a chemical2D layout.
Goal remains active with58 standards +50 additional still required.
