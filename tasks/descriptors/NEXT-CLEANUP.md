# Workbench chemical preparation and strict rotor correction

Pinned definitions:RDKit2026.03.6
[MolOps.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/MolOps.cpp),
[Lipinski.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/Lipinski.cpp).
Existing BSD notices cover both. No new runtime dependency or3D.

Preparation is scoped to `workbench::analyze`;source text and the raw native
graph contract remain unchanged. Neutral five-valence N has the pinned ordered
O-double/N-triple cleanup passes. Neutral five-valence,three-neighbor P uses
the selective O/C-or-N condition. Neutral Cl/Br/I oxo groups use the specified
valence/oxygen-neighbor condition. Update touched implicit H and rerun aromatic
perception with shared ring topology. Explicit bracket H is retained. This is
not the full sanitizer:coordination,validation and other stages remain work.

An unusual O=N#N input exercises both N passes,charge change and a newly added
implicit H;its expected output follows the pinned reference,not an assumption
that cleanup must conserve its original charge or H count.

Strict rotors now distinguish explicit vertex degree from heavy degree,count
all attached H for methyl exclusions,retain the full linkage exclusion when
the double-bonded N is charged,and accept non-ring aromatic bonds as the query
does. This is a registered query count,not a torsional barrier prediction.

## Independent evidence

49 source inputs,11,368 scalar expectations plus49 formula checks. Runtime
RED1,303 differences→47 (H refresh needed)→1 (rotor rule)→0. Separate73-input
strict-rotor reference exposed9 differences,all fixed. Full Rust760passed/
0failed/1existing ignored. Raw graph/idempotence contract also tested.
Full49-source C ABI strict comparison passes all232 fields and formula.

Actual10k Swift and supplemental15 changed-row comparison are in
`../evidence/catalog232-cleanup/`. The changed-row view checks all232 fields+
formula and does not replace whole-corpus evaluation. One Ce coordination
complex remains different and has a new single-bond-count difference due to
the still-missing dative representation;this is explicitly retained as work.

## Next: coordination semantics

A dative edge is directional. Pinned Bond.cpp gives numerical order1 but
valence contribution0 at donor and1 at acceptor. It is not a Single bond.
Default ring perception excludes it;adjacency/degrees and individual
descriptor definitions require separate handling. BCUT's registered wrapper
must remain undefined when its bond-type handling cannot represent it.

Primary sources:
[Bond.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Bond.cpp),
[QueryOps.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/QueryOps.cpp),
MolOps.cpp organometallic cleanup. QueryOps uses its pinned M-atom query,not a
hand-picked transition-metal list. Multi-metal ordering uses canonical ranks;
do not substitute input order and claim full equivalence. No coordination
production implementation has been added at this checkpoint.
