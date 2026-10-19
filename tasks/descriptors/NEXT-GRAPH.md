# Graph / topology tranche — v14

Graph19landed in v15 (178/217then,catalog202). Current v16 has197/217,
catalog221including24pre-existing
native diagnostics. Additional50 remains0/50. Full parity is not claimed.

## Pinned primary definitions

- [ConnectivityDescriptors.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Descriptors/ConnectivityDescriptors.cpp)
- [GraphDescriptors.py](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/rdkit/Chem/GraphDescriptors.py)
- [Subgraphs.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Subgraphs/Subgraphs.cpp)
- [ConjugHybrid.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ConjugHybrid.cpp)
- [Atom.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Atom.cpp)
- [Aromaticity.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Aromaticity.cpp)

These definitions are BSD-3-Clause; see the existing RDKIT-DATA-LICENSE.txt.

## Shared work and scientific meaning

Chi0/Chi1: degree-weighted connectivity. Chi0n..4n /Chi0v..4v: products of
inverse square-root valence deltas along paths; heavy-element correction for
v. These summarize branching and local connectivity rather than measured
experimental properties. Formal charge is not subtracted from the reference
connectivity delta definition. Chi0/1 include explicit-H graph vertices/bonds;
higher-order paths use useHs=False.

One traversal collects unique bond-set paths of1..4 bonds. Its final atom may
repeat any previous atom, closing a ring. Traversal stops there. This includes
ring-with-tail paths; restricting closure to the starting atom undercounts them.
The reference Chi implementation skips a repeated *starting* final atom only;
a final atom repeated elsewhere is multiplied again. Preserve this behavior.

HallKierAlpha: sum the pinned hybridization-dependent C/N/O/P/S/halogen
contributions, with Rb0 radius fallback for other elements. Its local context
includes conjugation, valence and radicals. The old conformer hybridization
helper cannot distinguish conjugated amide N/anisole O from simple SP3 atoms.
Kappa1/2/3 reuse alpha and path counts. Phi = Kappa1*Kappa2/heavyAtoms (zero
when no heavy atoms), a molecular flexibility index. No3D coordinates required.

BalabanJ and BertzCT share one weighted shortest-path matrix (bond length
1/order; disconnected sentinel1e8). BalabanJ retains RDKit's q-n+1 expression
for disconnected graphs. BertzCT classifies sorted distance rows using the
reference4-decimal symmetry keys/cutoff100, and weights pair/triple connections
by bond order. H, isotope and disconnected examples are independently covered.

## Oracle separation and open chemistry gaps

`generate_graph_reference.py` selects185reference inputs without consulting
native acceptance. It now emits3515scalar values,2290atomic contexts and740path
counts, and independent sanitized RDKit atom/bond graphs. Pure algorithm tests
use those graphs; integration tests use the pre-existing41 reference SMILES.
The complete177-source-SMILES C ABI comparison is retained separately and is
**strict-failed**, not filtered to native successes: v14 had172 complete matches,5
differing rows; v15 fixes the three FractionCSP3 radical-carbon gaps, leaving
175complete matches and two aromaticity/sanitization gaps (`[n+]1ccccc1`, hypervalent P ring). They affect
old descriptors as well as BalabanJ/BertzCT. Do not confuse the sanitized-context
algorithm test with full end-to-end SMILES parity. These gaps remain required
work, and the large-population harness retains them.

## IPC remains unimplemented, not counted

Ipc/AvgIpc use the characteristic-polynomial coefficients of the adjacency
matrix. The RDKit/NumPy trace/division recurrence becomes ill-conditioned on
large graphs. `investigate_ipc.py` compares it with the exact path-polynomial
coefficient |c_(n-2k)| = choose(n-k,k). At n128, the oracle Ipc is about112times
the analytic value; native sparse accumulation also diverges. At256 the error
is much larger. This is not ordinary last-digit tolerance. Experimental native
recurrence is retained outside the shipped crate/catalog in
`experimental_ipc.rs`; literals in `ipc-numerical-controls.tsv` and the numerical
report remain evidence. Need a stable characteristic-polynomial algorithm and
an explicit mathematical/numerical compatibility policy. Do not widen tolerance
or emit zero/null placeholders and call them implementations.

## Remaining at v15 (historical39)

Graph2 IPC; PEOE_VSA14; BCUT8; Lipinski5 (amide/stereo/bridgehead/spiro);
Descriptors8 (radicals1, charges4, Morgan density3); QED1; SPS1. Additional50 is
separate. Gasteiger should supply PEOE/charge extrema/BCUT charges from one
validated shared atom array. All additions must keep the first202 IDs stable,
bump cache revision, pass independent goldens and expose remaining full-corpus
differences honestly.

## v16 hydride correction

The completed v15full100k comparison found10new differing hydride rows.
RDKit's outer-electron minus GetTotalNumHs expression promotes to uint32 and
wraps for BH4−/AlH4− and related inputs. v16 uses explicit wrapping_sub rather
than producing NaN from a signed negative delta. Eight independent hydride
fixtures produced16runtime-RED Chi0n/v values before correction. All10actual
source regressions now match all221fields+formula (strict passed). This does
not replace a final current-catalog full100k gate. Additional50still0/50.

Current v16 implements19of those (radicals1,charge extrema4,PEOE14),leaving20standard fields plus the separate50additional requirement. See NEXT-CHARGE.md and README.md for current gates.
