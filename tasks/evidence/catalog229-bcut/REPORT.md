# Chemob v17 /229 descriptor checkpoint — 2026-10-18

**205/217 standard descriptors +24 existing diagnostics =229**. Additional
50: **0/50**. The full goal is unfinished. No 3D or new crate dependency.
Revision `ducksmiles-workbench-v17-bcut229-rdkit2026036`; frozen Mac library
SHA `aea81ade6c80e80d17c324215126a96cb07a2ce7481ba2fc74bd2a19c79eadc9`.

## Added eight BCUT descriptors

Mass, Gasteiger charge, atomic Crippen LogP and MR each supply a Burden matrix
diagonal. Return largest and smallest eigenvalues for each matrix. No 3D.
One existing charge/Crippen array is reused when there are no H vertices.
Otherwise remove all H vertices, including isotopes, and recompute that context
without changing the original input's surface or molecular-total calculations.
Unknown charge parameters make all eight fields undefined, matching the pinned
registered wrapper. Empty H-suppressed graphs have eight defined zeroes.

Native scaled Householder reduction plus Sturm bisection computes only extreme
eigenvalues. Definitions, source links and hydrogen policy:
[`NEXT-BCUT.md`](../../descriptors/NEXT-BCUT.md). RDKit BSD notices accompany
the native and App sources. No LAPACK dependency or copied implementation.

## Independent numerical and source input checks

Runtime RED before catalog integration, eigensolver and hydrogen-removal
implementation. Independently sanitized contexts:216 molecules,1,728 BCUT
values and exact RemoveAllHs graph/H/isotope controls pass. Solver:26 independent
NumPy eigvalsh matrices up to256, shift/scale/permutation and invalid-input
controls pass. Context tests do not stand in for source-SMILES compatibility.

Rust workspace **756 passed /0 failed /1 existing ignored**, no new warnings.
Swift Core **11 passed**; scoped rustfmt passed. Original41 end-to-end inputs
include BCUT references. Separate full216-source C ABI gate:216 both accepted,
214 all229+formula matches. Known aromaticity differences remain for
`[n+]1ccccc1` and `CP1(C)=CC=CN=C1C`; all eight BCUT fields are affected too.
Strict gate fails. Unknown charge modes have24 matching undefined examples per
BCUT field in this source set.

NCI4,999:4,991 both accepted,4,986 all229+formula matches;5 accepted value
differences and8 acceptance differences. Strict fails. Versus v16, the same13
source rows differ and existing field mismatch counts are unchanged. BCUT has
185 matching undefined cases per field. Retained:source-controls,nci and
nci/scope-delta.json. No acceptability filtering or tolerance relaxation.

## Actual Swift 10k I/O, cache and new-field oracle

10,000 distinct ChEBI source inputs. Actual App Workbench generated229 numeric
and233 CSV columns. Python CSV independently checked **2,282,615 numeric cells**,
original names/SMILES/order, unsupported25 rows with blank values and successful
Swift CSV reimport. CSV21,600,601 bytes. Old221 values and all existing JSON
fields are exactly unchanged in all10k rows; prefix receipt retained.

Actual disk cache:21,396,042 bytes, all10k restored with zero further native
calculations and identical values/order. Fits32MiB. Mac cold138.009s,warm.847s;
export/snapshot47.586s. These overlapped NCI verification, another bulk run and
compilation: **not a controlled version/RDKit benchmark or device/UI speed**.
BulkIO and CacheBulk had both exited0 by the user-requested pause; post-resume
independent checks completed without recomputing their native results.

RDKit audit of **only the new BCUT8 in actual Swift results**:9,975 common
accepted,9,962 all selected fields match;13 accepted value differences and25
unsupported. All13 are existing structural/sanitization differences in the
v15 full100k evidence; no new differing source row in this subset. They include
coordination compounds, nitro/amine oxide/oxo-halogen cleanup and aromaticity.
Value and overall gates fail. Actual-swift-bcut-10k retains all differences and
the source-row comparison. This is not an all229-field or full100k parity claim.

The completed full100k comparison remains **v15/catalog202/178 standards**.
Final217+50 whole-catalog100k verification is still required.

## iPad App verification

iPad Air11M4,arm64,iOS26.5 Simulator build succeeded. **Five selected portrait
UI tests pass**:all8 BCUT values visible inline,explicit/isotope H with original
surface preserved,charge pages,radical count and Graph pages. Search remains
hittable and no detail modal is required. Nine PNGs exported; the ethanol and
isotope-water BCUT captures inspected at original resolution. Seven languages
have120 matching keys; translations remain drafts. Installed v17 App normally
relaunched with original saved inputs. Build logs and XCTest summary retained.

No physical-device speed measurement:both registered iPhones were unavailable.
Landscape/old rotation failure,current iPhone and10k rendered UI remain open.
Preview is explicitly labelled connectivity-only, not a chemical2D depiction.
Small negative rounded values can still display as `-0.00`; native values are
preserved, and this is separate from scientific accuracy.

- [Eight BCUT values](ipad-attachments/3F7437FE-1E23-409A-ACD6-847D2A2CB1C5.png)
- [Isotope H removal for BCUT](ipad-attachments/5ABBBD92-4BE2-46D6-9DC5-00DD9C5DE89E.png)

## Remaining work

12 standards:IPC2,Morgan density3,Lipinski5,QED1,SPS1;meaningful additional50;
known chemistry precision and unsupported notation/size;numeric filters,
proper chemical2D and remaining device/UI gates. IPC numerical-policy and
additional50-use-case optional questions remain unanswered; other work proceeds.
