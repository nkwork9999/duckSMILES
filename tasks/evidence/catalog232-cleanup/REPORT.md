# Chemob v19 /chemical preparation checkpoint — 2026-10-18

Catalog unchanged:**208/217 standards +24 existing diagnostics =232**.
Additional50:**0/50**. The full goal remains unfinished. Revision:
`ducksmiles-workbench-v19-cleanup232-rdkit2026036`.
Frozen Mac dylib SHA:
`234f6246319cc7a5132f6a9d26fea06810fc68d2af3b86938ec6a8aca3778804`.

## Precision changes and genuine failing controls

Workbench-only N/P/oxo-halogen preparation,implicit-H refresh and shared-ring
aromaticity refresh. Raw native graph contract and original input text retained.
Strict rotor vertex-degree,total-H,charged-linkage and bond-type query rules
corrected. These implement computational definitions,not experimental property
accuracy. Details/source links:[NEXT-CLEANUP.md](../../descriptors/NEXT-CLEANUP.md).
No3D or new dependency. This remains a partial sanitizer.

Independent49 source inputs,**11,368 scalar values +49 formulas**. Runtime
RED1,303 differences→47→1→0. Separate73 strict-rotor inputs exposed9 differences;
all fixed. Intermediate failures retained. Full Rust workspace **760passed /
0failed /1existing ignored**,no new Rust warnings. Scoped new modules/tests and
rotor section formatted. One intermediate whole-suite run failed the still-open
rotor case;the subsequent full run passes. Raw graph/idempotence scope tested.

Full49-source C ABI independent strict comparison:**49 all232+formula matches**.
No native acceptance filtering,undefined→zero substitution or relaxed tolerance.

## Actual Swift10k and cache

10,000 distinct ChEBI source inputs,232 numeric/236 CSV columns. Python CSV
independently checked **2,312,540 numeric cells**,source names/SMILES/order,
25 unsupported rows with blank values,and successful CSV reimport.
CSV21,720,224bytes. Cold34.223s,warm.102s with all10k reused and zero extra native
calculations;export/snapshot11.527s. Single Mac host observation,no controlled
system load or equal-workload RDKit/device/UI speed claim.

Real Swift CacheIntegration now checks corrected nitromethane values as well
as duplicate/invalid-input reuse,row positions,graph reuse,exact disk reload,
stale-engine rejection and charge/isotope/stereo separation. Same test linked
to v18 fails the new literal,exit133;v19 passes,exit0. C ABI receipt records old
LogP0.4161 versus correct−0.10710000000000008. Corrected values survive disk
reload exactly. Full10k disk capacity was last tested inv18 at21,456,002bytes;
v19 tests10k in-memory reuse and focused disk behavior,not a new10k disk run.

## Changed-row independent audit — progress and a retained regression

Exactly15 of10k actual Swift result rows changed;9,985 are identical. Selection
was any result change **before** RDKit evaluation,not acceptance/agreement.
Original artifacts are untouched;copies reindex only for the existing validator
and retain original indices in projection.json.

Both versions'15 rows independently checked against all208 registered fields,
then24 profile fields+formula added. **Before:15 differing rows/363 field cells.
After:1 row/17 field cells.**347 field differences resolved,**1 new field
difference**. Fourteen rows now match all232 fields+formula.

Remaining:ChEBI-structure-23350,original index3573,a Ce nitrate coordination
complex. Missing dative/ring semantics affect17 fields. Its single-bond count
is newly different:the old wrong cleanup and missing dative behavior happened
to cancel. This is explicitly **not** a clean regression gate or full parity
claim;fix the underlying bond model rather than patch that one output count.
All details:changed-inputs/all232-plus-formula-delta.json.

## Broader independent checks

Actual Swift10k,**only recent30 fields** (charge19,BCUT8,structural3):9,975 common
accepted;9,958 all selected values match;17 value-differing rows+25 unsupported.
Same30 in actualv18 results:26 value-differing rows. Nine rows resolved,zero
new differing rows in that30-field scope. Both current value/overall gates fail.
The separate new profile-field regression above is outside this30-field scope.

NCI4,999:4,991 both accepted,**4,987 all232+formula matches**,4 value differences+
8 acceptance differences. Strict fails. Versusv17/catalog229:one source row871
resolved,no new differing source row,and common229 field error counts improve.
The baseline did not include the three structural fields;scope-delta says so.

These do not replace final217+50 whole-catalog100k evaluation. Completed full
100k remainsv15/catalog202/178 standards with its documented residuals.

## App verification

iPad Air11M4,arm64,iOS26.5 Simulator build and **two selected portrait UI tests
pass**:corrected nitromethane/explicit-H nitrate values with original source
text,and prior isotope-water surface/BCUT regression. Four PNGs exported;
two new chemistry images inspected at original resolution. Search remains
hittable,no detail modal or screen scrolling. Normal launch explicitly
terminated test mode first;Japanese UI restored the original six saved inputs.
Seven languages123 keys remain drafts. Same SDK AppIntents metadata warnings
as prior builds. Physical device/current iPhone/10k rendered UI/landscape and
the old rotation failure remain unverified. Connectivity preview is not a
chemical2D depiction;numeric filters remain unimplemented.

- [Corrected nitromethane](ipad-attachments/FFEB52F7-B620-45E5-AC21-FCA93483C161.png)
- [Explicit isotope H rotor](ipad-attachments/BE328560-1908-45FF-8289-181417D5A1C3.png)

Nine standards,meaningful additional50,coordination/aromaticity/ring/mass and
notation/size gaps plus remaining App gates are still required. Next precision
work is the dative bond model,starting from the retained Ce regression.
