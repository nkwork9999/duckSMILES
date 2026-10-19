# Chemob v16 /221 descriptor checkpoint — 2026-10-17/18

**197/217 standards +24existing profile diagnostics =221**. Additional50 is
**0/50**. The goal is active and not complete. No3D or new crate dependency.
Revision `ducksmiles-workbench-v16-charges221-rdkit2026036`; frozen Mac dylib SHA
`495dd3c001b91dcf85944bab85acc9e735b1fd9d6450ebf64533540458422119`.

## Added19 and shared model

NumRadicalElectrons1, charge extrema4, PEOE_VSA14. One electronic context feeds
Graph, radical counts and Gasteiger; one12-iteration charge array feeds all18
charge fields; one original-atom Labute area array feeds surface/charge bins.
Parameters are32pinned RDKit facts. Conjugated charge split, attached vs separate
H and Python extrema NaN-order behavior are independently retained. No arbitrary
zero charges for unknown parameter modes; preserve the reference model's finite/
undefined behavior. See `../../descriptors/NEXT-CHARGE.md` and source notices.

## Genuine RED and independent goldens

Missing catalog API runtime RED; then2326atomic differences from the test-only
zero stub before model implementation. Pure sanitized oracle contexts expanded
to207molecules,3933scalar values,2372atomic q/qH/hybrid/radical controls and
conjugated-bond controls. All14PEOE bins have positive examples. Every boundary,
NaN and infinities checked. Mixed finite/unknown-parameter components test both
input orders. Defined ethanol charge→null mutation fails the independent subset
verifier, while the real native result passes. No missing value is called a
correct zero. Runtime RED and mutation receipts retained.

Rust workspace **750passed /0failed /1existing ignored**, no warnings. Swift
Core11passed. End-to-end workbench tests retain the original41reference inputs;
full207-source-SMILES C ABI comparison is separate:205full221+formula matches,
two known aromaticity gaps, all207both accepted. Strict gate fails.

## Hydride regression from the v15 full100k run

v15's completed202-field100k comparison:99,755both accepted,99,270all fields+
formula matches;485accepted differences+245unsupported.57old differing rows
resolved through FractionCSP3 (59old field differences→0),10new hydride rows.
The latter revealed C++ uint32 subtraction underflow in the pinned Chi delta
formula. v16 reproduces it explicitly with wrapping_sub, with runtime RED on
16scalar hydride values before fixing. All10actual source regressions now match
**all221fields+formula**, strict passed. Do not treat this targeted set as a
full197-standard100k pass.

## Actual App10k I/O/cache and independent new-field check

10,000distinct source inputs;221numeric/225CSV columns; **2,204,175numeric cells**
independently checked against actual Swift results. Original text/order,
unsupported25rows and CSV reimport retained. CSV20,115,599bytes. Disk cache
19,905,920bytes:all10k restored,0new native calls, exact values/order. Fits32MiB.
Mac cold43.555s, warm.359s, all10k reused; timings overlap correctness work and
are not a fair version/RDKit comparison or physical-device/UI speed claim.

Old202prefix/legacy unchanged on these10k except four Chi0n/v values on the two
hydride rows7359/7801; null→correct finite values. Prefix receipt retained.

Independent RDKit verification of **only the new19fields in actual Swift10k
rows**:9975both accepted,9964all selected values match,11accepted value rows
differ,25unsupported. All11are already differing structural/sanitization rows
in the v15full-corpus evidence; **no new differing row** in this subset. The
value and overall gates still fail. This is not an all221-field100k claim.
Examples include oxo-halogen/nitro cleanup and pre-existing aromaticity/surface
interpretation. The verifier retains every difference and computation error.

NCI4999inputs:4991both accepted,4986full221+formula matches;5value differences
+8acceptance differences remain. Strict fails. The old/new NCI difference-row sets are identical (new0,resolved0).

## Remaining verification / work

Current v16 iPad Air11M4/iOS26.5 Simulator build and four selected portrait UI
tests pass:partial charges,PEOE pages,positive radical count and FractionCSP3,
Graph pages and original explicit-H surface regression. Seven PNGs exported;
three post-fix images inspected at original resolution (charge values,PEOE
nonzero bins,radical count). An earlier functional run also passed;visual
inspection found Composition word wrapping,so a minimum picker width was added
and all four tests repeated. Header-before PNG retained. Search stays hittable,
no detail modal required. Seven languages119matching keys;translations draft.
App is installed and normally running;completed host/Apple build caches,
generated simulator.a and current xcresults removed after evidence export.
Only latest v16frozen Mac dylib and latest v16bulk raw output remain.

20standards remain:IPC2,BCUT8,Lipinski5,Morgan density3,
QED1,SPS1. Additional50must still have meaningful formulas/applicability/tests;
no aliases/unit conversions/profile padding. Numeric filters, proper chemical2D,
large-input UI,current iPhone,physical-device speed and prior rotation failure
remain open. IPC numerical-policy and additional50-use-case optional questions
are pending; other work continues.

## Selected current UI artifacts

- [Charge extrema](ipad-attachments/B05B5611-2390-4826-BE14-2B8368E2E901.png)
- [Nonzero PEOE surface bins](ipad-attachments/A02F578C-4447-419A-A1D5-8D1DD6FF866D.png)
- [Radical electron count](ipad-attachments/C52F8FFD-7BC8-48AA-94EE-B7D1E8FD0144.png)
