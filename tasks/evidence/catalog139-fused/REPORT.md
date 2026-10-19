# Chemob / v11 fused-aromaticity checkpoint

Catalog unchanged: 139 entries = 115/217 RDKit standards + 24 native profiles.
Additional 50: 0/50. Goal remains incomplete; no 3D.

## Chemical interpretation

v10 fixed globally cyclic versus exocyclic donation and the exactly-one-shared-
bond / <=24-atom ring neighbor rule. Its complete 100k comparison reduced total
errors but introduced nine previously matching rows. **v10 is superseded**, not
an accepted final accuracy result. See `../catalog139-cyclic/chebi/`.

v11 computes global atom donation once, excludes ineligible rings from fused
candidate systems, checks connected subsystems of 1–6 rings, excludes atoms
buried in >=3 selected rings from the electron count, and marks only perimeter
bonds. It follows the pinned [RDKit aromaticity implementation](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Aromaticity.cpp).
Existing original aromatic inputs remain preserved; this is not full RDKit
sanitization, Kekulization or valence equivalence. BSD notices retained.

## Verified

- `aromatic-subsets-red.txt`: actual runtime failure on newly exposed graphs
  before subsystem implementation. Final fixture: 27 independent RDKit atom/
  bond aromatic-flag graphs, including all nine newly exposed v10 cases.
- `rust-workspace.txt`: 736 passed, zero failed, one existing ignored; no new
  warnings. Includes fragment 8,330 cells and nonfragment standard 1,230 cells.
- `nci/comparison.json`: 4,999 inputs (4,900 distinct raw SMILES), 4,991 common
  accepted, 4,986 all139+formula matching; five accepted differences + eight
  acceptance differences remain. Strict gate failed. v8 had 20 accepted
  differences; v11 has no new differing NCI rows.
- `regressions/`: deliberately error-enriched 1,223 rows, comprising all 1,216
  v8 differing rows and seven observed new v10 rows. **439 old differing rows
  now match**; all seven new v10 rows match. Remaining 777 includes 245
  unsupported. Six mass/formula differences remain. Not a random sample.
- `tail-regressions/`: the final two new v10 rows both match all139+formula.
  Together with the enriched set, all nine new v10 regressions are covered.
- `bulk-io-summary.json`, `csv-validation.json`: v11 actual Swift Workbench,
  10,000 distinct ChEBI source inputs; 25 unsupported preserved, 143 columns,
  139 numeric descriptors, **1,386,525 numeric cells** independently checked.
  CSV reimport, original names/SMILES/order and blank invalid values passed.
  Export 8,410,418 bytes. Cold import/calculation 16.096s, warm 0.186s,
  10,000 reused / no additional calculations. Mac host, no UI; overlapped
  correctness work, so not a fair cross-version or RDKit speed comparison.
- `csv-tamper-rejected.json`: the independent verifier rejects an actual
  exported functional-group value altered by +1 and emits no passing report.
- `../catalog139-cyclic/cache-bulk.json`: v10 actual ChemistryStore at 10,000
  distinct inputs; disk cache 8,669,292 bytes; reopening reuses all 10,000,
  zero native calculations, all values and row order equal. Separate revision
  from v11; not a v11 full-disk-cache measurement.

## Complete full-population comparison

`chebi/comparison.json`: 100,000 distinct source graphs, 99,755 both accepted,
**99,223 all139+formula matching**, 99,238 standard+formula matching, 532
accepted value/formula differences + 245 unsupported =777 total differing rows.
Six mass/formula differences remain. Strict gate exited **1**.
`chebi/scope-delta.json`: **439 old differing rows resolved, zero new differing
rows** against v8. This is not full RDKit compatibility.
Library SHA256 `cf14b860fd0030cc19f7f62391d1ea92bae7412291a39583ab86081ff5d91dc0`.

Formatting: affected files passed. Workspace-wide cargo fmt fails on 28
unchanged baseline files; see `format-summary.json` and `rust-fmt.txt`.
The aligned Crippen data table has an explicit rustfmt skip marker; no values
changed. Formatting-only edits were made after the scientific host build.

## App / remaining verification

The installed iPad Simulator app is now **v11**. iPad Air 11-inch M4 /
iOS 26.5 arm64 Simulator build passed. Two selected portrait UI tests passed:
functional-group paging / four EState values, and actual pasted NCI-1218 input
through analysis to the inline panel. Corrected aromatic bonds **11** and atoms
**10** (old 6/6), and displayed LogP **4.42** (old4.22) were asserted through
accessibility values. Search remains hittable and no detail modal appears.
Three screenshots were visually inspected; `ipad-portrait-summary.json` and
`ipad-portrait.txt` retain the result. The graph preview remains a circular
connectivity diagram, not a proper chemical 2D depiction.

Older full iPad rotation failure remains unresolved; it was not deleted or
weakened. Current iPhone, physical-device speed and large-input UI rendering
remain unverified. This portrait pass does not resolve those gates.

Native unsupported-source audit: 102 exceed the 256-atom bound, 107 have
CXSMILES extensions, 36 non-tetrahedral annotations. Original inputs were not
rewritten or stripped; categories are source-feature classifications, not a
proof that dropping annotations would preserve chemistry.

102 remaining standards, 50 meaningful additional descriptors, numerical
filters, proper 2D depiction and complete graph/sanitize parity remain work.
