# Workbench coordination semantics — v20

Pinned RDKit2026.03.6. This implements a bounded part of chemical preparation,
not full sanitization or complete RDKit equivalence. Original input is retained.

## Representation and definitions

`BondOrder::Dative` stores donor=`a`,acceptor=`b`. Numeric order is1;endpoint
valence contribution is0/1. Graph adjacency/degrees/components retain the edge;
default ring perception excludes it. Do not globally erase it or equate it to
Single. SMARTS bond truth tables now include a fifth type;default single-or-
aromatic excludes dative,while wildcard/negation retain their meanings.

Workbench parses as-written with the existing bounded grammar,then applies
N/P/halogen cleanup,coordination inference,ring refresh if needed,and aromatic
perception. The raw native-graph API continues its prior contract. Neither the
raw grammar nor its schema has acquired explicit `->`/`<-` input support.

Donor eligibility uses effective atomic number and maximum allowed valence,
including the aromatic total-degree4 condition,and the pinned MDL M atom query
for metal membership. With one eligible metal neighbor the edge is directed
donor→metal. Multiple eligible metals currently return unsupported;matching
RDKit's canonical atom ordering and metal tie-break is still required. No
input-order substitute is claimed. No new dependencies or3D.

Electronic states use endpoint valence;orbital counts exclude donor dative
edges where the reference does. TPSA retains graph degree but does not count
dative as a single bond;Labute applies zero bond-radius correction. BCUT's
registered throwing-wrapper behavior is all eight values undefined,including
mass. BertzCT deliberately uses bond-type numeric17 in its connection counts
(the pinned Python definition),while weighted distances use numeric order1.
FractionCSP3 follows total-degree4 even for coordinated aromatic carbon.

## Independent checks

35 valid original inputs,8120 scalar values plus35 formulas. Includes13 prior
ChEBI coordination error rows and3 explicit/isotope-H carbon controls. Two
invalid aromatic source spellings are separately recorded as RDKit-rejected;
they are not presented as successful value comparisons. No native acceptance
or agreement filter. Beforev19:30 differing rows/308 fields;afterv20:0/0.

Initial API runtimeRED93→1→0,expandedRED3→0. Atomic H/charge/aromaticity/total
valence/radicals/hybridization,directed bond types/conjugation,and endpoint-
mapped ring sets also match. Parsers append ring closures in different bond
array positions;the context harness compares edge identities,not bond indices.
Full Rust workspace762passed/0failed/1existingignored;no new Rust warnings.
Actual C ABI35/35 strict all232+formula pass. Scalar tolerances unchanged.

Actual Swift cache test checks Ce18 acceptors/12 rotors/0 rings/18 single bonds,
its independent BertzCT literal,and Mg all8BCUT undefined. The same Swift test
linkedv19 traps;v20 passes and restores values exactly from disk. Actual10k I/O,
whole100k comparison and current App UI verification are ongoing;record their
final scope in the checkpoint report. No physical-device speed claim.

## Primary implementation references

- [Bond.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Bond.cpp)
- [MolOps.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/MolOps.cpp)
- [QueryOps.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/QueryOps.cpp)
- [ConjugHybrid.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/ConjugHybrid.cpp)
- [MolSurf.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/MolSurf.cpp)
- [Lipinski.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/Lipinski.cpp)
- [GraphDescriptors.py](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/rdkit/Chem/GraphDescriptors.py)

Upstream BSD attribution retained in crates/smiles/RDKIT-DATA-LICENSE.txt.
