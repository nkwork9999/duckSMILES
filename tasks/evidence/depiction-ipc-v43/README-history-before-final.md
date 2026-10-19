# STOPPED — v43 follow-up, user requested safe stop

Stopped on 2026-10-18. All owned compiler/test processes have exited; no agent
or benchmark remains active. Changes are preserved, uncommitted. No push or release.

**Current source is WIP:** `crates/smiles/src/depiction_graph.rs` contains an
intentionally failing prepared-graph API stub. The App now calls that API.
Do not build/install this source as a finished App until the API is implemented.
The installed simulator App is the earlier verified phase1 v43 build.
Simulator boot restoration stalled; original iPad was safely shut down (exit 0).
Its six saved inputs match the original backup exactly.

Completed checkpoints:
- Ipc/AvgIpc diagnostics: stable exact-integer policy retained; legacy numerical
  differences remain 5,723 rows. No calculation arithmetic or tolerance changes.
- Swift Core: latest actual exit 0, 33 tests passed. Prepared-graph Rust RED:
  actual exit 101 (expected stub failures), three tests define pending behavior.
- Earlier 10,017-input drawing sweep: original10k had 9,365 drawings,
  610 quality rejections and 25 graph rejections. This predates the latest fan
  and prepared-graph changes and is not final current-source validation.
- One E/Z fan regression found and fixed; targeted reference round trip passed.
  Full-corpus post-fix verification is pending.
- Raw drawing-graph audit found 28 prepared-chemistry differences (10 charge
  differences), motivating the new isolated prepared-graph API.
- Physical iPhone execution remains unavailable; no physical speed claim.

Resume with the prepared-graph API implementation and its RED tests, then native
build/parity, regenerated prepared graph corpus, full stereo/layout checks,
cache proof and final App UI verification. Keep immutable v42 evidence unchanged.
Evidence: `tasks/evidence/depiction-ipc-v43/stop-checkpoint.json`.

--- Earlier checkpoints below are historical ---

# Chemob follow-up: 2D depiction and information diagnostics

Work in progress, 2026-10-18. Calculation revision remains v42; its immutable
100k audit and artifacts are in [the v42 report](../catalog289-sps/REPORT-v42.md).

## Ipc / AvgIpc

[Diagnostic](information-order-diagnostic.json) reuses all 5,723 recorded residual
rows and runs additional exact-integer FLINT controls on paths40/80 and original
source22, with the original string plus three atom-order variants. All source
identities are checked with fixed RDKit2026.03.6. Stable values match independent
coefficients and remain invariant. In source22 the pinned floating Ipc changes
by approximately3.88x across equivalent atom orders. The recorded residuals
include4,447 Ipc differences larger than100% relative to stable values.

These results do not establish a portable legacy numerical compatibility mode.
The fields remain explicitly named `Ipc_stable` / `AvgIpc_stable` in CSV.
No tolerance widening, new native arithmetic, or relabelling as strict numerical
RDKit interchangeability is performed.

## First 2D layout baseline

Swift core generates bounded ring templates, chain angles, bond-order strokes,
aromatic solid/dashed parallel strokes, and atom labels retaining hydrogen,
charge and isotope metadata. Crossings, overlapping vertices and bad bond
lengths reject the preview. No molecule or calculation input is modified.

The original first10,000 source positions plus10 curated cases were evaluated,
without selecting by drawing success. First baseline:8,569 drawings including
all10 curated,1,416 failed layout-quality gates,25 native graph API rejections.
The graph API remains limited to256 atoms/512 bonds, separately from the
descriptor engine's640/1280 limits. This baseline does **not** implement
specified stereochemical wedges or E/Z geometry; the App explicitly warns
when source stereo is omitted. Complex bridged/ring geometry needs further work.

[Baseline summary](drawing-probe/summary.json),
[unavailable sources](drawing-probe/unavailable.json),
[compressed native graph fixtures](graphs-10010.jsonl.gz).
The source/rights of these ChEBI release256 records are retained in
[v42 data notice](../catalog289-sps/v42-input-data/DATA_NOTICE.txt).

Actual Core tests:18 before changes;25 after the initial layout, all passed.
The six new positive cases were observed failing against the unavailable stub
before implementation. Actual native→Swift cache probe was observed failing
without memoization, then passed duplicate/disk/unsupported-result reuse,
hydrogen labels and source-stereo notice. Cache namespace separates drawing
input and coordinates from all descriptor values.

## Physical execution

Both paired phones remain unavailable, with no tunnel/DDI service. No physical
execution or speed result has been obtained. Pairing, trust and signing settings
were not changed. Simulator timing must not be presented as physical speed.

UI integration and further stereo/complex-layout work are pending in this file.
