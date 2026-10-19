# Chemob v18 /232 descriptor checkpoint — 2026-10-18

**208/217 standards +24 existing diagnostics =232. Additional50:0/50.**
The full goal is unfinished. Revision:
`ducksmiles-workbench-v18-structure232-rdkit2026036`.
Frozen Mac dylib SHA:
`ab80803b2f43de3e9fc4829f6fca1ea07b48caa613bb4b6d9848a61d39d68db3`.

## Added three registered counts

NumAmideBonds,NumSpiroAtoms,NumBridgeheadAtoms. Shared parsed atoms/adjacency
and symmetric-ring topology; no3D,new dependency or extra parse. Amide uses
the registered carbonyl/non-ring-O/N query. Spiro/bridgeheads use pairwise ring
intersections and union atom identities,not molecular degree alone.
Definitions,primary source and notices:
[`NEXT-STRUCTURAL.md`](../../descriptors/NEXT-STRUCTURAL.md).

Runtime RED:missing registered NumAmideBonds,then all **753 scalar counts across
251 source inputs pass**. Positive controls:amide23,spiro8,bridgehead8. No native
acceptance/agreement filtering. Charged/isotope/explicit H,urea,imides,ring O,
non-amides,disconnected copies,fused/bridged/spiro/cage cases included.

Full Rust workspace **757 passed /0 failed /1 existing ignored**,no new Rust
warnings;scoped rustfmt passed. Swift Core source is unchanged;thev17 eleven
Core tests remain its latest test run. v18 release Core module and actual App
bulk/cache executables built and ran. Seven languages123 matching keys,draft
translations. Native/App BSD notices include Lipinski attribution.

All251-source C ABI comparison:251 both accepted,249 all232+formula matching.
New3 counts all match;two existing aromaticity rows still differ in other
fields,so strict gate fails. Frozen-v17 old229 catalog metadata and all result
fields/values are exactly unchanged on these251 source controls.

## Actual App Swift 10k

10,000 distinct ChEBI source inputs. **232 numeric/236 CSV columns**;
**2,312,540 numeric cells** independently checked using Python CSV. Original
names,SMILES,row order,unsupported25 rows with blank values and successful
CSV reimport preserved. CSV21,720,528 bytes. Old229 values and all existing
result fields exactly unchanged on all10k rows;prefix receipt retained.

Actual disk cache21,456,002bytes:all10k restored,zero new native calls,identical
values/order. Fits32MiB. Host observation:cold35.814s,warm.147s,
export/snapshot13.499s. No other task build/corpus test overlapped the cold bulk
run;system load was not controlled. These are not an equivalent-workload
version/RDKit comparison,physical-device measurement or rendered UI timing.
In particular,thev17 run overlapped multiple jobs,so its138s must not be used
as a speedup baseline.

Independent RDKit check of **only the new3 in actual Swift rows**:
9,975 common accepted,9,960 all selected values match;**15 accepted value
differences +25 unsupported**. Per-field mismatches:amide1,spiro14,bridgehead12.
All15 are previously differing source rows in thev15 full100k evidence. Most
are coordination compounds;two large cage/macrocycle cases expose remaining
ring-set differences. New differing source rows:0. Value and overall gates
still fail. Scope-delta retains every source and field difference.

This is not an all232-field or whole100k parity claim. The completed full100k
remainsv15/catalog202/178 standards;final217+50 verification is still required.

## iPad App

iPad Air11M4,arm64,iOS26.5 Simulator build and **three selected portrait UI
tests pass**:new counts across inline pages,BCUT8,explicit/isotope-H original
surface and BCUT handling. Urea shows2 amide bonds,norbornane2 bridgeheads,
spiro compound1 spiro atom. Search remains hittable;no detail modal or screen
scrolling. Six PNGs exported,new-count three inspected at original resolution.
Installedv18 normally relaunched with original saved inputs.

Xcode emitted the same two AppIntents metadata-extraction-skipped warnings as
v17 (no AppIntents framework dependency);no new engine compiler warning.
Thev17 five-test portrait checkpoint remains historical evidence,not an extra
v18 five-test run. Current iPhone,physical-device speed,10k rendered UI,
landscape and the old rotation failure remain open. Preview is connectivity
only;numeric filters and proper chemical2D remain unimplemented.

- [Amide count](ipad-attachments/00A4A1CB-4B1D-495A-9B79-0511183A056C.png)
- [Bridgehead count](ipad-attachments/A52F5AD0-8D42-4B6F-A850-9F1CEDECFE89.png)
- [Spiro count](ipad-attachments/B362F12A-8D0B-4A0D-816D-94D4407F6059.png)

## Remaining

Nine standards:Morgan density3,stereo2,QED1,SPS1,IPC2. Additional50 remain0.
Chemical cleanup,aromaticity,coordination/dative,ring perception,mass/formula,
unsupported notation/large graphs and remaining App gates are still required.
Next precision work targets pinned N/P/oxo-halogen cleanup with independent
all-catalog controls. Existing Morgan is explicitly not bit-exact RDKit;do not
reuse its density as a supposedly compatible registered descriptor.
