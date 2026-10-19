# BCUT2D implementation and verification

Pinned oracle: RDKit 2026.03.6 registered `Descriptors.BCUT2D_*` functions.
This is eight standard descriptors, not part of the additional 50. No 3D,
third-party numerical library or new crate dependency is used at runtime.

## Definition

First remove **all hydrogen vertices**, including isotope and isolated H,
following `MolOps::removeAllHs`. Preserve heavy isotope masses. Recompute
attached H, electronic context, Gasteiger charges and atomic Crippen values on
that graph. Reuse existing charge/Crippen arrays when the source has no H vertex.
The original input graph and its VSA/molecular-total calculations are retained.

For each of mass, Gasteiger charge, atomic Crippen LogP and atomic Crippen MR:

- diagonal: that atomic property;
- bonded off-diagonal: single 1, double 0.7071067811865475,
  triple 0.5773502691896258, aromatic 0.8164965809277261;
- other off-diagonals: 0.001, including disconnected components;
- output: largest then smallest eigenvalue, in that property order.

These are the pinned upstream literals, not rounded user-facing constants.
An empty graph after H removal returns eight zeroes. Unknown Gasteiger modes
make all eight fields undefined, matching the registered Python wrapper's NaN
result when the underlying vector function throws. Do not retain a plausible
mass-only result or replace undefined charge with zero.

`spectral.rs` scales finite symmetric matrices, applies Householder reduction,
then uses an LDL/Sturm count and bisection for the two extreme eigenvalues.
It does not compute eigenvectors or the full spectrum. Invalid input or a
nonconverged result stays undefined. Descriptor absolute tolerance is 1e-9.

## Independent controls

`generate_bcut_reference.py` emits 216 independent sanitized graph contexts,
216 RemoveAllHs contexts and 1,728 scalar reference values. They include
explicit/isotope/isolated/bridging H, nondefault valence, charged components,
heavy isotopes and unsupported charge parameter modes. Pure graph-context
tests isolate the BCUT algorithm from known source-SMILES sanitization gaps.
End-to-end workbench tests retain the original 41 source inputs. A separate
216-source C ABI comparison and actual Swift 10k audit are required; graph
context parity alone is not source input compatibility.

The numerical solver has 26 independently generated NumPy eigvalsh references:
dense/random, diagonal, rank-one, path, star and clustered/degenerate cases,
up to 256 atoms. Additional shift/scale/permutation and invalid-input controls
exercise the production solver. The eigensolver, H removal and catalog API
each had a runtime failing test before implementation.

Current run evidence is saved in `../evidence/catalog229-bcut/`.
The full 217+50 descriptor goal and final 100k gate remain open.

## Primary sources

- [RDKit BCUT.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/BCUT.cpp)
- [RDKit AddHs.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/AddHs.cpp)
- [LAPACK DLAEBZ algorithm documentation](https://raw.githubusercontent.com/Reference-LAPACK/lapack/master/SRC/dlaebz.f)

RDKit BSD attribution is included in the native and App notices. The numerical
solver is a native implementation of the standard mathematical operations;
LAPACK is not linked or copied into the App.
