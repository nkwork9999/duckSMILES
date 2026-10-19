# Three registered structural counts

Pinned oracle: RDKit2026.03.6 registered NumAmideBonds,NumSpiroAtoms and
NumBridgeheadAtoms. These are standard descriptors,not additional50.

The amide definition matches a nonaromatic C double-bonded to a nonaromatic,
non-ring O and connected to nonaromatic N by the default SMARTS bond. Urea
therefore counts two matches; a ring carbonyl oxygen is excluded. The direct
native graph calculation implements the registered query,not a generic count
of nitrogen near oxygen.

Spiro atoms are unique atoms shared by exactly one vertex in a pair of
symmetric-SSSR atom rings. Bridgeheads are deduplicated endpoints in a ring
pair's shared bond subgraph when more than one bond is shared. Molecular degree
alone is insufficient; fused rings sharing one bond do not count as bridgehead
pairs under this definition. The existing ring topology is reused.

Primary definition source:
[pinned RDKit Lipinski.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/Lipinski.cpp).
BSD notice is included in native/App notices. No new runtime dependency or3D.

## Controls

`generate_structural_reference.py`:251 independent source inputs and753
registered scalar expectations. Positive controls:amide23,spiro8,bridgehead8.
Includes charged/isotope/explicit-H amides,urea,imides,non-amides,ring O,
disconnected copies,fused/bridged/spiro systems,cages and multiple shared rings.
Independent sanitized graph and symmetric bond-ring fixtures are also retained.

The end-to-end production workbench test first failed at runtime on missing
NumAmideBonds,then all753 scalar expectations passed without filtering inputs.
The production implementation uses shared parsed molecular/ring data. Full
workspace757passed/0failed/1existing ignored. C ABI on all251source inputs:
249 all232+formula matches,with the same two known aromaticity rows differing.
New3 have zero differences. The old229 catalog prefix and all result fields
are exactly unchanged across all251source inputs versus frozenv17.

This does not establish full217+50 or whole-catalog100k compatibility. Current
expanded I/O/cache,selected App UI and independent actual-Swift subset checks
are recorded separately in `../evidence/catalog232-structural/` as they finish.
The existing native Morgan/QED algorithms are not assumed to match registered
functions; their precise defaults and numerical contracts require separate work.
