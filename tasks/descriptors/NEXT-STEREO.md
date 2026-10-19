# Implemented engine v30: legacy2 + modern SPS; full parity unfinished

Catalog289:215/217 standards+24profiles+50extras. Independent789 SPS values and per-atom/per-bond/factor contexts pass, including20 difficult-source regressions and18 macrocycle atom-order variations. The initial751 all289+formula strict pass is historical v29; expanded789 current v30 has783 full matches/6 remaining non-SPS chemical-preparation differences. Allthree new standards match all789 and all167 accepted difficult415 sources.

Current v30 corrected donor electron/degree eligibility, uppercase as-written atom flags, initial D2-root first-ring caching for macrocycle ether/sulfur donors, and incident-aromatic-bond valence in coordination cleanup. Raw perceived native graph remains unchanged. The D2 order correction retains the established relevant-cycle basis for later rings; it is not a general claim of canonical RDKit SSSR ordering.

FullRust779pass/0fail/1existingignored/0warnings. ActualSwiftv30 10k289numeric/293CSV,2779168cells checked,25unsupported; all289/base exactly equalv29. New3 independent9975 accepted checks carry forward. Cacheold29RED/new30GREEN includes corrected boron/macrocycle values and disk reload. Current App build succeeds; UI baseline/new-feature tests running, installed normal App verification still v26. Full100k current289 audit running; do not quote final counts until completion.

Difficult415:21→35 full matches,146→132 accepted value failures,248unsupported,0formula.736 old differing cells resolved,18 newly differing cells on6 already-failing rows remain and are recorded; no introduced failing rows in that selected cohort. No full precision claim. IPC2 remains unimplemented. Primary-source archive and notices retained.

--- Original research notes below ---

# Remaining standard stereocounts and SPS — research, not implemented

v26currently212/217standards+24profiles+50extras. No count increase from this
research. After current App gates/cleanup, generate independent broader scalar
and per-atom/per-bond controls, runtime missing-feature RED, then implement.

## Primary pinned sources and version state

Direct Release_2026_03_6 source SHA/line receipts:stereo-source-audit.json.
Cached web fetch of Lipinski.cpp returned530lines while direct pinned GitHub
returned553lines and the expected functions. Use the direct pinned file hashes
and actual2026.03.6 oracle, not stale cached search snippets or master docs.

Lipinski.cpp numAtomStereoCenters/numUnspecifiedAtomStereoCenters version1.0.1
count atoms carrying _ChiralityPossible; unspecified also requires CHI_UNSPECIFIED.
If _StereochemDone is absent, the function copies and runs assignStereochemistry
(clean=true,force=true,flagPossible=true). Parsed input normally already has
stereochemistry assigned. Actual pinned wheel legacy default is True. Counts
are not simply the number of @ tags; invalid/redundant tags must be cleaned.
Chirality.cpp assignAtomChiralCodes uses CIP-equivalence ranks and legal-center
rules, plus later ring/special-case processing. Do not substitute generic
hybridization or new FindPotentialStereo for the legacy count without testing.

SPS(version1.0.0) uses a molecule copy, FindPotentialStereoBonds, modern
FindMolChiralCenters(includeUnassigned=True,includeCIP=False,
useLegacyImplementation=False), then sums hybrid*stereo*ring*degree² divided by
heavy-atom count. hybrid=1SP,2SP2,3SP3,4otherwise; stereo=2for modern atomcenter
or a DOUBLE bond with stereo!=NONE(ANY included),else1; ring=2nonaromatic-ring
atom,1otherwise. Actual code GetDegree includes retained explicit H vertices,
despite a prose reference to heavy neighbors. Do not use extra50 heavy graph.

## Confirmed preliminary controls

stereo-preliminary-controls.json records actual pinned calculations:
C:SPS0; [H]C([H])([H])[H]:SPS64 with H vertices retained;
[H][H]:SPSraisesZeroDivisionError(no heavy atoms), stereocounts0.
CC(O)Cl:total1/unspecified1/SPS15.75;
C[C@H](O)Cl:total1/unspecified0/SPS15.75;
CC=CC:SPS9.5 and one STEREOANY double bond.
C1CCN2CCCC2C1:count1/unspecified1,SPS36.666666666666664,
modern chiral_idxs=[7] in this exact pinned copy pipeline. Do not infer a
historic renumbering issue's values for this version without an actual control.

A known no-heavy SPS domain error should be specified explicitly as undefined
in the FFI; do not hide arbitrary reference exceptions as success. Any exception
policy change needs named tests/source justification and clear comparator
reporting, not a broad tolerance/acceptance relaxation. Reference fixture must
also cover tetrahedral P/As,S/Se,N(quaternary,small-ring,bridgehead), explicit/
isotopic/bridgingH,dependent ring stereocenters, assigned-invalid tags, equivalent
ligands, enantiomers, atom permutations, doublebond ring-size thresholds,
stereo-free inputs and complex ChEBI examples. Keep original source/state.

Sources:
https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/Lipinski.cpp
https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Chirality.cpp
https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/FindStereo.cpp
https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/rdkit/Chem/SpacialScore.py

The two remaining IPC standards still have a separate numerical-instability/
policy issue documented in NEXT-GRAPH; user option was already asked previously,
do not repeat/block all independent work. Proper2D,numericfilters,library reset/
replace and wider acceptance/precision remain product tasks. Current imports
append only with a10000total cap; a normal saved6-input library cannot add a
full10000without a reset/replace flow. The10kUItest deliberately starts empty
under--ui-testing; it does not erase the user's normal saved library.
