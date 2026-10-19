# QED and element/preparation corrections — v24

2026-10-18. **212/217 RDKit standards +24 existing profiles =236**.
Requested additional50 **0/50**; overall goal unfinished. Revision
`ducksmiles-workbench-v24-qed236-rdkit2026036`.
[Definitions and primary references](../../descriptors/NEXT-QED.md).

## Native chemical evidence

QED uses the pinned eight-property/11acceptor/116alert definition, default H
removal and WEIGHT_MEAN. Existing properties and shared SMARTS context are
reused. Foundation fixes address two-letter element parsing, uppercase
non-organic SMARTS semantics, missing atomic numbers, bounded isolated As/N+
Kekulé cycles, aromatic candidate element scope and integral radical valence.
Raw parser/native graph consumer semantics remain tested and separate.

Runtime failing controls retained. Initial364 QED controls expose6scalarerrors;
expanded383 expose22. All235+formula aromatic35 controls expose128scalarerrors
before preparation/radical fixes, then pass. An additional C ABI35 test finds
one QED error on aromaticTe not present in the10k corpus. Final413 QED controls
and638element-token cells cover this, including five previously missing atomic
number entries. No numerical tolerance was relaxed to hide a discrepancy.

Final release Rust **770passed/0failed/1existingignored**, no new Rust warnings.
413score+3304property+4543acceptor+47908alert reference values match;51default-
H-removal retained graphs and638element-token truth-table cells match.
The first hydrogen graph test compared ring-closure endpoint order on an
undirected edge: a harness error, corrected without a chemistry change. Its
log is retained. Scoped formatting passed after formatting new test lines.

Frozenv24 C ABI: **413/413 and35/35 all236+formula strictpass**.
SHA256 `b88aec91c283940b285b1c94d590be3053775d958ecc22ffb1cbf31e5231a723`.
The earlier v23 C ABI383 pass and35/34fail are retained; v23was never built
into the App. This does not prove whole-corpus precision or full sanitization.

## Actual Swift10k and cache

Actual App code on Mac, noUIrendering. **236numeric/240CSVcolumns**,2,352,408
numeric cells independently checked. Names,originalSMILES,sourceorder,CSV
reimport pass; all25unsupported rows remain with blank values. QED matches
**9975/9975accepted** inputs; valuegate passes but completegate fails25unsupported.
This is a QED-only whole10k chemical comparison, not all236-field parity.

All235common fields/baseoutputs are exact on all10k v22→v23; all236/baseoutputs
are exact on all10k v23→v24. Thus these foundation changes do not change old
values in this particular10k set; separate held-out controls catch the affected
chemistries. Original data and library/output hashes are saved in receipts.

Currentv24 full10k disk store:22,105,571bytes;10kexactrestore,zero new native
calls, allfields/order identical,25unsupported retained. Focused actualSwift
QED/As/Se/Te cache and disk reload pass. Same final Te-specific harness linked
v23fails(-5),v24passes(0); prior As-specific harness failsv22 and passesv23.
All10k in-memory warm results also reuse with zero additional native calls.

Uncontrolled Macv24 import54.161s/warm.239s/export-and-snapshot17.819s are
observations, not controlled version comparisons or device timings. The prior
v23import47.651s occurred under different load and is not a regression claim.

## Same-workload RDKit timing

Fixed413raw-SMILES regression controls (not representative10k throughput),
fresh parse each molecule/trial,236fields+formula+baseoutputs+JSONencode/decode,
no result cache. Correctness/status gate passes on all inputs before timing,
no native-acceptance filtering. Five alternating paired trials after one warmup.
Mac median: **native0.725767s, RDKit2.848137s (~3.92x ratio)**. Raw times,
source/library hashes and load averages are in rdkit-paired-timing.json.
RDKit includes Python computation of the24profile fields; its successful BCUT
vector wrappers cache/reuse all eight outputs on each fresh molecule (verified
in pinned source). Native batches64 through C ABI, RDKit uses a Python rowloop.
Background hostload is uncontrolled. This is not physicaliPhone/iPad speed,
UI rendering speed or a general3.92x throughput claim.

## App and required remaining work

CurrentiPadAir11M4/iOS26.5Simulator v24build and3portraitUItests pass. Six
PNGs exported; two newQED/As originals visually inspected. Normal Japanese launch restores the six saved inputs; the original-resolution
screenshot was inspected. New QED,Hg,As,Te,H-removal values and
priorMorgan/macrocycle regressions reach the inline panel. Seven language files129matchingkeys/
placeholders, translationsdraft. No new runtime dependency or3D. Numericfilters,
properchemical2D, currentiPhone/landscape/10k renderedUI and olderrotationfailure
remain open. PhysicaliPhones unavailable. Five standards remain:stereo2,SPS,IPC2.
Meaningfulextra50 candidate research is not an implementation/count increase.

Completed full100k remainsv20/catalog232:99585all232+formula match,167accepted
value-differing rows,248unsupported,fiveformularesiduals; strictFAILED. Three
newunsupported inputs versusv15remain in the report. No currentv24whole100k
parity claim. Source/ring/notation/size accuracy gaps and final full-catalog
100k verification remain required.

## Error-enriched wider precision check

ALL415rows that differed or were unsupported in frozenv20whole100k were
rechecked withv24, without selecting by current native acceptance. Result:
19all236+formula matches,148acceptedvalue-differing rows,248unsupported retained;
strictFAILED. Fiveformula residuals remain. **QED differs on25accepted rows;
Morgan1/2/3 each differ on7**. These are largely complex-ring/aromatic/valence
foundation cases. This makes the limits of the413-controls/10k-QED passes
explicit. Do not infer whole100k v24 counts from this selected error cohort.
Next precision work must address those retained failures; do not hide them by
loosening tolerances or silently omitting sources.
