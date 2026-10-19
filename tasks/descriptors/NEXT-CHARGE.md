# Gasteiger / charge tranche — v16

Implemented19standards:NumRadicalElectrons1,partial-charge extrema4,PEOE_VSA14.
Catalog221=197/217standards+24existing profiles. Additional50remains0/50.
Full chemical parity is not achieved; see `../evidence/catalog221-charges/REPORT.md`.

Primary definitions,pinned Release_2026_03_6:
- [GasteigerCharges.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/PartialCharges/GasteigerCharges.cpp)
- [GasteigerParams.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/PartialCharges/GasteigerParams.cpp)
- [MolSurf.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/Descriptors/MolSurf.cpp)

## Shared inputs / semantics

One electronic context (valence,radicals,conjugation,hybridization) feeds Graph
and charges. One12-iteration Gasteiger array feeds all18charge fields; radical
count comes from the same context. One original-atom Labute area array feeds
existing VSA and PEOE. No3D coordinates,no new crates,no runtime RDKit/Python.

DAMP=.5,DAMP_SCALE=.5,IONXH=20.02,EPS_DOUBLE=f64::EPSILON. Initial charge split
follows conjugated two-bond same-element pairs,not an entire-domain average.
Preserve repeated markers/input neighbour order. Compute all electronegativities
before each charge-update round; attached-H totals and separate-H atoms differ.
Parameter modes use complete hybridization,sulfur fallback and pinned unknown-
mode zero-parameter behavior. Preserve its NaN propagation,not arbitrary neutral
charges.32parameter rows are independently pinned facts with BSD attribution.

Python extrema update uses the new charge as the first min/max argument,so
NaN behavior depends on atom order; Rust f64::min/max would silently differ.
MinAbsPartialCharge is min(abs(minCharge),abs(maxCharge)),not the smallest
absolute atom charge. PEOE bisect_right uses13charge thresholds;NaN goes to the
last bin as in the pinned reference. Undefined extrema remain null/dash/blank.

## Independent gates and remaining gaps

207molecules/3933scalars/2372atom q/qH/hybrid/radical controls,conjugated bonds,
all bin edges and mixed finite/unknown components in both input orders. All14
bins have positive independent examples. Atomic zero-stub runtime RED retained;
then full model passes. A real native ethanol control passes the independent
subset verifier; changing a defined charge to null produces one mismatch.

Full207-source C ABI:205all221+formula matches,two pre-existing aromaticity
gaps. Actual10k Swift output,new19only:9975both accepted,9964complete matches,
11value-differing rows plus25unsupported. All11are already differing structural/
sanitization rows in v15; no new differing row in this subset. Gates still fail.
Keep formula/charge/bond cleanup,aromaticity and coordination-bond differences
in scope. Do not call the small sanitized-context test full SMILES compatibility.

Current10k CSV/cache and selected iPad portrait tests pass. v15full100k is an
older202-field/178-standard snapshot; it is not a current197-standard100k pass.
Final whole-catalog100k validation remains required. Next20standards:BCUT8,
Morgan density3,Lipinski5,QED1,SPS1,IPC2. BCUT can reuse the validated charge and
Crippen arrays; it still needs a validated native symmetric eigensolver. IPC
numerical-policy preference and additional50-use-case questions remain optional
and pending while independent work continues.
