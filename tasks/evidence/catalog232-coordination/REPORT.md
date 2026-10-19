# Chemob v20 / directed coordination checkpoint — in progress

Catalog232:208/217 standard+24 existing profiles. Additional50:0/50. Goal unfinished.
Revision `ducksmiles-workbench-v20-coordination232-rdkit2026036`.
Frozen Mac SHA256 `523972a9d34245b792dd13cb14d53243783a2153e486c50fdcc008089c469b9e`.

## Verified host changes

Directed dative bonds separate numeric order1 from donor/acceptor valence0/1.
Default ring perception excludes these edges;degree/distance and descriptor-
specific semantics are retained. Workbench prepares the as-written graph before
aromatic perception. Raw native-graph contract/source text retained. Multiple-
metal canonical selection,arrow input grammar/full sanitization remain open.
[Definitions and primary sources](../../descriptors/NEXT-COORDINATION.md).

35 original valid sources:8120 scalar+35formula checks pass in Rust/API and
actual C ABI. Beforev19:30 differingrows/308fields→afterv20:0/0. Independent
atomic valence/H/charge/radicals/hybridization,conjugation,directed edge types
and endpoint-mapped ring sets also match. FullRust762pass/0fail/1existingignored;
no new Rust warnings. Scoped formatting passes. No dependency/3D added.
Genuine runtime failing controls:initial93→1→0,expanded3→0. A separate context
harness incorrectly assumed equal bond-array ordering;that harness was fixed
to compare edge identities,and is not counted as a production correction.

## Actual Swift10k I/O and disk cache

10000distinct ChEBI source strings,232numeric/236CSV columns. IndependentPython
CSV checks2312508numericcells,original names/SMILES/order,25unsupported with
blankvalues,and CSV reimport. CSV21719687bytes. Warmreuse10000/0newNative.
Cold71.268s/warm.209s/export28.994s were observed while the full100k comparison
was running:uncontrolled load,not a speedup/fairRDKit/physicaldevice/UI claim.

Actual10k disk cache21455273bytes;10000exactrestores,0extra nativecalculations,
allfields/undefinedvalues/roworder retained. FocusedactualSwift cache checks Ce
18HBA/12rotors/0rings/18single and independentBertzCT literal,Mg8BCUTnil. Same
Swiftlinkedoldv19fails(signal5);v20passeswith exactdiskreload. Catalog232metadata
unchanged. Beforehashes/projections saved before supersededv19binary/raw deletion.

## Changed-row independent numerical comparison

18actualSwift rows changed,9982exactunchanged. Selection is any actualresult
change before consultingRDKit,not acceptance/agreement. All208standard+24profile
+formula on18:before245differencecells/18rows→after8cells/1row;237resolvedcells,
0introducedfields. Ce regression fully corrected. Remainingoriginalindex3589,
ChEBI-structure-77631,cucurbituril+Pt:8ring-set fields. This is a changed-row
scope;it does not replace whole-corpus comparison. See changed-inputs/.

## Required evaluation still running

Whole100k independentall232+formulaCABI comparison finished on frozenv20.
99752bothaccepted;99585all232+formula match,167acceptedrowsdiffer and248are
nativeunsupported. 99586matchall208standard+formula (one profile-only difference).
Fiveformularesiduals remain. Strictgate fails. This is208standards/232fields in
v20,not the newer211standards/235fields. Inputselection preserves100000distinct
canonicalisomericgraphs and original source strings;no native acceptancefilter.
The248unsupported count includes the existingnotation/size limits and possible
newmulti-metal fallback;comparetheprior245source-set before claimingnewcountcause. Appv20 build and3selectedportraitUItests pass on iPadAir11M4/arm64/iOS26.5. FivePNGs exported;new2 inspected atoriginalresolution. NormalJapaneselaunch restored6savedinputs. Same2SDKAppIntentsmetadata warnings;test-runnerLLDBmetadata diagnostics logged despite3passingtests. No new engine/compiler warning. Sevenlanguages123matchingkeys/placeholders;the existing localization checker rejected blankseparators,and now skipsblanklines (runtimefail→pass);translationdraftstatusunchanged.
BothphysicaliPhones checked:unavailable. BeforebuildnormalJapaneseiPadv19saved
sixinputs were inspected. No currentphysicalspeed/iPhone/landscape/10krendered
UI claim. Connectivitypreview/proper2D,numericfilters and the oldrotationfailure
remain open. Nine standards,additional50 and otherchemical/notation/size gaps
remain required work. Do not mark the full goal complete.


## Current App evidence

- [Coordination ring count](ipad-attachments/D539A9B7-31D4-4316-92D2-2D66061DF8E4.png)
- [Undefined BCUT values](ipad-attachments/06C121DF-4335-471A-A9D6-BA216F1BE199.png)
- [Normal saved-input launch](normal-launch-saved-inputs.png)

Chemical uniqueness was separately audited read-only:all10000 inputs have
10000distinct canonicalisomericSMILES under pinnedRDKit/removeHs=False,with
isotope/charge/stereo retained and no salt/tautomer/neutralization merging.
No dataset filtering or rewriting was used for this audit.
