# Current v42 validation and practical-use metadata

All50 declared formulas now pass the complete original100000-input comparison,
actualSwift/iPad10000CSV and cache gates. Formula, domain, and an individual
practical_use are in extra50-definitions.json. They are structural-feature
candidates with the documented graph policy, not established predictive gains.
The native catalog finalization changes descriptive verification/practical_use
only: all70 calculation source files and the entire compiled SPECS prefix remain
byte-identical to the fully audited42 source. See v42-catalog-metadata-source-proof.json
and v42-catalog-metadata-binary-proof.json (1006 fresh rows, all291 exactly equal).
Final Simulator and unsigned device builds pass; the final Japanese normal iPad
launch preserves all six inputs and visibly reuses all six results. No numerical
formula, order, tolerance, graph policy, cache semantics, or original input changed.

[Complete results](../evidence/catalog289-sps/REPORT-v42.md).

--- Earlier design and controls ---

# Additional50 — published formulas, implementation and validation

v26 candidate implements50 scalar2D descriptors absent from the pinned RDKit
built-in descriptor API. These are published topological families, not novel
physical measurements or demonstrated predictive improvements. Existing24
profiles and optionalRDKitAUTOCORR2D are excluded. No new runtime dependency/3D.

Independent480 original-source inputs ×50 =24000cells pass, including no-heavy
undefined cases, salts,isotopic/explicit/bridgingH, rings, coordination, long
chains, seeded trees and atom permutations. All50 have positive independent
controls. Sparse native BFS/walk DP is compared against dense NumPy matrix
powers/Floyd-Warshall plus hand graph checks. Runtime missing-feature RED and
subsequent explicit-dative-input RED are retained. Thirteen separate arrow
input controls verify donor→acceptor direction, attached H, formulas and
covalent-ring counts while retaining legacy raw native-graph grammar.
Full catalog/100k,Swift/cache/CSV/App gates are still pending forv26.

## Built-in audit

extra50-builtin-audit.json records pinned2026.03.6 module names. RDKit already
exposes CalcAUTOCORR2D and192AUTOCORR2D channels outside the default217 list.
Moran/Geary/Moreau-Broto families are therefore EXCLUDED from this target.
Generic graph/matrix primitives can of course compose many formulas; absence
claims concern descriptor implementations, not mathematical expressibility.
The module inventory is supplemented by a pinned RDKit source/contrib term
audit (2006bounded text files, zero candidate-name hits; archiveSHA recorded).
This is an API/source term audit, not a proof against unnamed implementations. CalcEccentricity is3D and is not the proposed graph
EccentricConnectivityIndex.

## Implemented50 formulas

-19MolecularDistanceEdge: carbon degreepairs11,12,13,14,22,23,24,33,34,44 (10),
 nitrogen11,12,13,22,23,33 (6), oxygen11,12,22 (3). For n eligible pairs, define
 n/geometricmean(graphdistance). Compute via mean-log to avoid product overflow;
 no eligible pair means undefined, not zero. Elements and branching classes
 describe how separated comparable carbon/nitrogen/oxygen sites are.
-10raw Galvez topological-charge indices GGI1..GGI10. For A adjacency and R with
 R_ij=1/d_ij^2 (zero diagonal/disconnected), CT=A*R-(A*R)^T; GGI_k sums |CT_ij|
 over unordered pairs at distancek. This is topology asymmetry, not a Gasteiger
 or measured partial charge. Do not also count normalized JGI aliases.
-8molecular walk counts MWC3..MWC10, each ln(1+sum(A^k)). Repeated vertices/edges
 allowed, unweighted graph. Exclude first-order bond-count rescaling and raw/
 log duplicate registrations. Captures branching propagation over3..10bonds.
-8self-returning walk counts SRW3..SRW10, each ln(1+trace(A^k)), capturing cyclic
 and returning connectivity. One descriptor per order; no duplicate raw form.
-5published graph indices: Wiener=sum unordered finite distances; Harary=sum
 reciprocal finite distances; EccentricConnectivity=sum degree_i*eccentricity_i;
 DegreeDistance=sum(degree_i+degree_j)*distance_ij; Gutman=sum degree_i*degree_j*distance_ij.

These are published scalar descriptor families rather than arbitrary fingerprint
bits. The parameters have explicit bond-distance/walk-length meanings. Most
have redundant information on special graph families; do not claim50independent
physical measurements or predictive improvement without evaluation.

## Defined graph/domain policy

Hydrogen-suppressed heavy-atom graph (exclude all H vertices including isotopes),
original element identity, unweighted bonds including directed-coordination
connectivity viewed undirected for these topological functions. Ignore bond
order/charge/stereo in these50definitions; do not merge such source/cache keys.
No3D, tautomerization, desalting or silent graph repairs. Preserve salts; use
only finite within-component pairs. State this departure from sentinel-distance
implementations explicitly. No-heavy-atom inputs:undefined50values. No eligible
MDE pair:undefined; missing distanceorder in GGI/walk families:mathematicalzero.

MDE valence means graph degree in the source implementation (confirmed against
_graph_matrix.Valence), not chemical bond-order sum. Defined disconnected MDE
uses eligible finite pairs globally; do not silently adopt component summation.
Use stable logarithmic geometric means; do not reproduce NumPy product overflow.
Native sparse walk recurrence and independent Python/NumPy matrix powers can
validate one another. Exact small-graph integer/brute-force checks, relabeling,
reversed-SMILES, aromatic/Kekulé and disconnected graph controls are required.
No new runtime dependencies. Shared BFS/distances/adjacency plus sparse DP;
no new dense eigenproblem. Benchmark cost before making all50eager in the App.

## Primary implementation references (pinned v1.2.0)

-https://raw.githubusercontent.com/mordred-descriptor/mordred/v1.2.0/mordred/MolecularDistanceEdge.py
-https://raw.githubusercontent.com/mordred-descriptor/mordred/v1.2.0/mordred/TopologicalCharge.py
-https://raw.githubusercontent.com/mordred-descriptor/mordred/v1.2.0/mordred/WalkCount.py
-RDKit exclusion: https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/Code/GraphMol/Descriptors/AUTOCORR2D.cpp

MORDRED-NOTICE.txt and the bundled App notices retain BSD-3-Clause attribution.
The explicit DegreeDistance name avoids the ambiguous Schultz/MTI terminology.
Primary graph formula references:
https://www.pmf.kg.ac.rs/KJS/en/images/volumes/vpl40/kjs40azari33.pdf
https://link.springer.com/article/10.1007/s10878-013-9645-5
Full machine-readable formulas/applicability: extra50-definitions.json.

Confirmed source facts: Mordred _graph_matrix.Valence is adjacency row/column
sum, i.e. unweighted degree under the default useBO=False. WienerIndex.py uses
half the distance-matrix sum. EccentricConnectivityIndex.py uses sum(E*degree).
Pinned primary files:
https://raw.githubusercontent.com/mordred-descriptor/mordred/v1.2.0/mordred/_graph_matrix.py
https://raw.githubusercontent.com/mordred-descriptor/mordred/v1.2.0/mordred/WienerIndex.py
https://raw.githubusercontent.com/mordred-descriptor/mordred/v1.2.0/mordred/EccentricConnectivityIndex.py

Speed note for later Bertz work: the older handoff's proposed f64-bit keys are
NOT proven equivalent to4-decimal strings. Mathematically equal multiples1/6
can differ in their final bits after differently ordered floating additions.
Consider exact scaled/quantized distance-class keys with independent symmetry
and full-output controls; do not implement naive to_bits grouping or claim an
exact domain proof from the mathematical lattice alone.
