# Chemob v15 topology checkpoint — superseded by v16

Catalog202 =178/217 standard descriptors +24 pre-existing profiles. Additional50
remains0/50. This is progress, **not objective completion or full RDKit parity**.
The frozen implementation is v15; v14 controls are superseded by a targeted CSP3
fix. v14's partial100k comparison was explicitly stopped and is not population
evidence. v15's full100k run is complete. Current development is v16/catalog221;
see `../catalog221-charges/REPORT.md`.

## Added19 / shared native implementation

Chi0/1; Chi0n..4n /Chi0v..4v; HallKierAlpha; Kappa1/2/3; Phi; BalabanJ; BertzCT.
Shared paths1..4, weighted distance matrix, atomic valence/conjugation context.
No new crate dependency or Python/RDKit runtime. Cache/C ABI revision bumped.
See `../../descriptors/NEXT-GRAPH.md` for formulas, sources and limitations.
Ipc/AvgIpc are **not in the catalog/count**: the numerical stability investigation
found severe errors even in the pinned reference on large linear graphs. The
experimental recurrence and independent analytic report remain reviewable.

## Tests completed so far

- v15 Rust workspace745 passed /0 failed /1 existing ignored; no warnings.
- Swift Core11 passed.
- Independent sanitized contexts:177 graphs ×19 =3363 scalar values;2280 atomic
  controls;708 path counts. The initial path implementation failed ring-with-tail
  cases; closure behavior was corrected from the pinned source.
- End-to-end registered Graph tests use the pre-existing41 reference inputs.
  This is distinct from algorithm testing on oracle-sanitized graphs.
- Full177-source-SMILES C ABI controls (not filtered to native successes):
  v14 had172/177 full matches,5 differing rows. v15 has175/177 matches and two
  aromaticity gaps. All177 are accepted by both. Strict gate still fails.
- Three radical-carbon CSP3 errors were fixed: the RDKit scalar counts four-
  coordinate carbon, even though its hybridization enum labels these radicals
  SP3. Runtime RED retained. The mixed radical control was corrected from
  `[CH3]C` (ordinary ethane) to `[CH2]C` after checking the oracle; this was a test
  authoring error, not a production fix hidden by changing a correct oracle.

## Actual App I/O / cache / selected portrait UI

v15 actual Swift Workbench10k distinct inputs:202numeric columns/206CSV fields,
2,014,946numeric cells independently compared to Swift results. Original text,
row order, invalid25rows and CSV reimport preserved. CSV17,550,384bytes.
Cold read/import/calculate37.304s, warm.202s; all10k reused,0additional native.
These Mac timings overlap full-corpus correctness work, not a fair RDKit/version
benchmark or physical-device performance. Disk cache17,474,201bytes, exact10k
restore/row order,0additional native calls. Fits32MiB.

Compared with previous v13 on these10k: every old183value and legacy field is
unchanged except two FractionCSP3 values; both corrected values independently
match RDKit. Prefix and I/O receipts retained.

v15 iPad Air11M4/iOS26.5 Simulator build and two selected UI tests pass.
Topology category, native values/two-page navigation, original explicit-H
surface regression, search hittability/no detail modal verified. Three original-
resolution PNGs inspected (including corrected zero format). The UI initially
showed exact negative-zero graph values as -0.0000; a UI runtime RED was recorded,
then display-only canonical-zero formatting fixed it. Native numeric values and
cache revision remain unchanged by that cosmetic fix. Test no longer taps Next
when the portrait second page is already the last page.

Seven languages117matching keys; scientific IDs remain canonical Latin names,
translations draft. App is installed and normally running; completed build
caches/generated simulator.a/xcresults removed after evidence export. The old
rotation-failure raw trace is still retained separately.

## Completed full100k v15 comparison

99,755both accepted;99,270all202+formula matches,99,274standard+formula matches.
485accepted value/formula differences +245unsupported =730differing rows,
6formula differences. Strict exit1: full compatibility is not achieved.
Against v13:57differing rows resolved,10new rows. The only changed legacy metric
is FractionCSP3:59mismatches→0, max error1→0. All other old183field mismatch
counts/max errors are unchanged. New10are unsigned-valence-delta hydrides.

v16 corrects the uint32 subtraction semantics and independently verifies all10
source SMILES against RDKit across all221fields+formula (strict passed). This
is a targeted correction, not a claim that v16 has completed a full221-field
100k run. Full197-standard/50additional compatibility remains required.

No physical-device speed claim. Numeric filters, proper chemical2D, current
iPhone, large-input UI and prior rotation failure remain open. Additional50
still0/50. Goal remains active.
