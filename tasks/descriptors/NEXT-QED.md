# Registered QED — pinned definition and preparation

RDKit2026.03.6 `rdkit.Chem.QED.qed`, version1.1.0, default WEIGHT_MEAN.
v24/catalog236 adds one standard:212/217 standards+24existing profiles.
Requested meaningful additional50 remains0/50; overall goal unfinished.

## Calculation

QED applies eight ADS desirability curves and a weighted geometric mean to
MW,CrippenLogP,QED-specific HBA,HBD,N/O TPSA,Strictrotors,AROM,ALERTS. HBA sums
11 individual single-query atom match counts (default1000 cap per query), and
is intentionally different from the standard NumHAcceptors value. ALERTS counts
presence once for each of116patterns, including two disjoint multi-component
queries. AROM uses ordinary SSSR dimension after deleting atoms matching
`[$([A;R][!a])]`; it is not the NumAromaticRings field or symmetricSSSR count.
The residual cycle dimension is computed by union-find, excluding dative edges.

Default Chem.RemoveHs happens within the reference QED implementation. Native
calculation preserves isotopic,isolated,bridging,hydride,H-only-neighbor and
double-bond-defining directional Hvertices. Other explicit Hs transfer to
attached-H totals on a private graph. No changes to source/raw graph. Most
inputs need no H removal; these reuse existing MW,LogP,HBD,TPSA,Strictrotors and
ring topology. One PreparedMatcher shares adjacency/rings across the acceptor
and114single-component alert queries; only two multi-component alerts use the
separate disjoint-mapping matcher. No new runtime dependency.

## Failing controls and foundation corrections

The initial364-source independent score/property test exposed six scalar
differences on Hg and two Rucomplexes. Expanded383 controls exposed22differences,
including Hf/Ho/Ru/Rh/se/As and explicit-H peroxide rotors. Element prefixes
Hg/Hf/Ho and Ru/Rh were incorrectly interpreted as Hcounts or ring predicates.
Two-letter element recognition now precedes one-letter SMARTS predicates and
uses the full118-element table; bracket aromaticas/se are single element tokens.

Arsenic cannot be aromatic under the default pinned model (Z>18 exceptions are
Se/Te). Workbench preparation reifies tested isolated C/N/As cycles with zero/
one-double-bond requirements per atom, solving the cycle constraints rather
than assuming every cycle is even. It retains the prior N+radical regression
behavior. Fused/other-element/other-charge cases remain outside this bounded
rule; this is not a complete Kekulizer. Raw parser semantics are tested.

Expanded35-source all235+formula controls exposed128scalar differences before
that preparation fix. Four remaining NumRadicalElectrons differences on aromatic
Se/Te/P/S+ exposed a separate issue: assignRadicals sees integral Kekulé valence,
not the later1.5+1.5 restored aromatic bond sum. Aromatic electronic states now
use their adjusted integral valence; other states retain their prior rule.
The35-source all235+formula gate then passes. Definitions/tolerances were not
relaxed. QED catalog runtime test initially fails because the field is missing,
then passes all383 independent scores after registration/shared-property reuse.

The final expanded413-source internal fixture records all11acceptor counts and all116alert
flags, not just aggregate values. It records51default-H-removal retained graphs.
The first graph test compared undirected bond endpoint order and failed on a
ring closure: a harness error, not a chemistry RED. Its corrected comparison
canonicalizes ordinary edge endpoints/array order and preserves dative direction.
A separate35-source C ABI check then exposed one new QED mismatch on aromatic
Te, absent from the10k corpus. [Te]/[Se] constrain atomic number, while organic-
subset [S]/[P]/[C] also require aliphatic atoms. Bracket [te] is one aromatic
atom token. A new638-cell,29query×22source truth table exposed3such mismatches.
The broader atomic-number path then exposed five missing table entries (Ho,Hs,
Ra,Rn,Dy); its fallback now covers the shared118-element table while retaining
the common-element fast path. Runtime failing controls are saved. v23was never
built into the App. v24is the final candidate; final whole-suite/CABI/Swift/App
results belong in the evidence report.

## References

- [QED.py](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/rdkit/Chem/QED.py)
- [AddHs.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/AddHs.cpp)
- [RemoveHs defaults](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/MolOps.h)
- [Aromaticity.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Aromaticity.cpp)
- [Atom.cpp](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Atom.cpp)

Novartis/GregLandrum/contributor BSD notices retained in native and App resources.
Matching this empirical score's implementation does not establish experimental
property/predictive accuracy. It does not supply a new descriptor for the extra50.
