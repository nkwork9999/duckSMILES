# Current resolution — v37, 2026-10-18

Ipc and AvgIpc are now registered using exact integer coefficients and stable
entropy. Both match their independent reference on all 99752 commonly accepted
inputs; pinned numerical results differ on 5617. CSV uses `Ipc_stable` and
`AvgIpc_stable`, with App compatibility guidance. The remaining strict acceptance
gate has 248 unsupported inputs. [Current report](../evidence/catalog289-sps/REPORT.md).

---

The following is the historical design rationale before integration.

# IPC: exact coefficients before choosing numerical compatibility

The two missing standards are `Ipc` and `AvgIpc`. They remain unregistered.
The earlier optional policy question is unanswered; this experiment does not
select a default or relabel mathematical agreement as pinned numerical parity.

## Why a separate oracle is needed

The pinned RDKit/NumPy floating Faddeev recurrence can diverge greatly on long
path graphs. At128 vertices, the retained diagnostic is about112 times the
analytic result. Copying that numerical instability would be a poor basis for
a new scientific calculation. Conversely, silently replacing values in a
feature set used by an existing RDKit model can change its predictions.
Both effects need concrete measurements before deciding the shipped contract.

## Proposed exact experiment

Use the unweighted undirected adjacency of **all vertices**, including explicit
H vertices and dative connectivity. This matches the graph used by pinned IPC,
despite its hydrogen-suppressed wording. No3D or trained parameters.

1. Obtain characteristic coefficients with sparse integer Faddeev arithmetic.
   A small-integer fast path may restart with signed arbitrary-width integers
   on overflow. Trace division by the iteration number must be exact, with
   remainder asserted. Never continue after overflow or truncate coefficients.
2. Convert absolute exact coefficients to finite doubles only at the entropy
   step. Compute both variants together. Empty/invalid adjacency is an error;
   an edgeless nonempty graph has zero entropy and zero IPC.
3. Validate integer coefficients separately from final floating entropy.
4. Compare exact results and pinned results in separate columns. Keep original
   input rows, report relative as well as absolute differences, and do not
   widen existing descriptor tolerances to hide instability.

`generate_exact_ipc_reference.py` is an experimental independent oracle:
integer Bareiss determinants at integer points, rational Newton interpolation
for small irregular graphs; analytic path/cycle/complete/bipartite polynomials
and polynomial products for large graphs up to256 vertices/512edges. Entropy
uses85-digit Decimal arithmetic. It is not called by native production code.
The generator has run:58 exact graph cases. The first standalone integer
candidate passes all signed coefficient vectors and85-digit entropy references,
plus invalid adjacency checks. Runtime RED was retained before implementation.
This is experimental only:neither descriptor is registered or exposed by the App,
and no pinned-compatibility/default policy has been selected.

## Primary references

[Rehman and Ipsen, La Budde's Method](https://arxiv.org/abs/1104.3769) studies
coefficient computation through orthogonal reduction and polynomial recurrences,
including error analysis. It is a floating alternative worth comparing, rather
than a guarantee of exact graph coefficients.

[Johansson, fast nearly division-free characteristic polynomial algorithms](https://arxiv.org/abs/2011.12573)
reviews exact-small-integer-division methods and a faster dense variant. The
current bounded molecular graphs are sparse, so initial implementation should
measure the simpler sparse recurrence before adding algorithmic complexity.

No new runtime dependency is authorized or required for this experiment.
