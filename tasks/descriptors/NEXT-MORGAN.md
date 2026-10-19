# Registered Morgan density and prepared aromaticity

Pinned RDKit2026.03.6. Three standards added: `FpDensityMorgan1/2/3`. Catalog235
contains211standards+24existingprofiles;additional50 remains0. Definition:
unique keys in the unfolded sparse count fingerprint divided by heavy atoms,
with0 for zero heavy atoms. Chirality disabled. The legacy native fingerprint
module remains separate;its folded bit-vector density is a different operation.

## Native calculation

Radius0 connectivity vectors contain atomic number,total degree,total H
including explicit H vertices,formal charge,integer mass difference from average
atomic weight,and a trailing1 only for ring members. Hashes use portableu32
arithmetic. Each layer hashes its center and sorted bond/neighbor pairs. Bond
codes retain aromatic12 and dative17. Duplicate environments are bond sets;
the first hash/atom within a duplicate environment survives. Dead atoms have
zero next-layer invariants. The three cumulative key counts share one radius3
pass. Initial invariants and every radius0..3 key set match independent graph
fixtures (352sources after the macrocycle regression was added).

## Chemical preparation corrections

Workbench-only aromatic candidate checks follow adjusted default valence and
radical rules. An isolated even C/N ring with a two-coordinate unprotonatedN+
radical and exactly one double-bond requirement per atom gets an explicit
Kekulé representation before perception. Fused/odd/other-atom cases require
further work. Original input and the raw native graph consumer contract remain.

The standard nitrogen acceptor queryv3 has no charge constraint;donorv3 likewise
has no charge constraint,whilev4requires+1. Earlier extra native charge guards
were removed after independent failing controls. This reproduces the pinned
queries,including their unusual charged-radical and pyridinium results.

A held-out10k source revealed one Morgan mismatch from an already incorrect
macrocycle model. Native/RDKit shared the same seven ring atom sets. Pinned
Aromaticity.cpp disables neutral degree2 O/S without multiple bonds when its
first ring encounter has >=9members. Such atoms are ether bridges there.
This rule is now applied using shared native ring order. First-ring/basis order
and other complex aromaticity/SSSR cases still require whole-corpus evaluation.
Smallfuran/thiophene and chargedO controls retain their reference behavior.

## Evidence and scope

Before macrofix,v21:351whole235+formulaCABI sources allmatch;actualSwift10k
new3:9974/9975acceptedmatch,onevalue-differing macrocycle+25unsupported. The
macrocycle differed in51fields,including48older fields. Expanded26-source
all235+formula regression exposed132errors;afterv22fix all26match. v22also
matches all235+formula on352sources. FullRust766passed/0failed/1existingignored.
New3 tolerance1e-9;oldtolerances unchanged. Initial comparator refusal came
from missingnewID tolerance registration,not numeric differences. Sparse keys
are tested exactly. Independent tables never select by native acceptance.

ActualSwiftv21 10k235numeric/239CSV/2342433cells passed I/O checks,with25issue
rows retained. Allprevious232fields/baseoutputs were exact acrossall10k versus
v20. Cache capacity21912427bytes,10000exactdiskrestores/0newNative. v22focused
cache adds the corrected macrocycle literal;linkedv21fails,linkedv22passeswith
exactdiskreload. v22actual10kreplay passed235-column I/O and allnew3 on9975accepted inputs,
with25unsupported retained. Four portrait App tests pass. Finalscope is in
../evidence/catalog235-morgan/REPORT.md. Completedfull100k is frozenv20/catalog232,
not a currentv22/all235 parity claim.

A one-second symbolicated sample of actualMacBulkIO shows native BCUTspectral,
graph/Bertz formatting and fragment work. This CLI uses actualAppcode and has
no UI rendering. Its main-thread idle state does not measure iPad frame time.
Concurrent experiments invalidate simple before/after speed comparisons.

## Primary references

- [MorganGenerator.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Fingerprints/MorganGenerator.cpp)
- [FingerprintUtil.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Fingerprints/FingerprintUtil.cpp)
- [Portablehash](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/RDGeneral/hash/hash.hpp)
- [Aromaticity.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Aromaticity.cpp)
- [Lipinski.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/Lipinski.cpp)

Upstream notices retained in crates/smiles/RDKIT-DATA-LICENSE.txt. No newruntime
dependency or3D calculation. Remaining6standards,meaningful50extra descriptors,
fullchemical/notation/size parity and App gates remain required.
