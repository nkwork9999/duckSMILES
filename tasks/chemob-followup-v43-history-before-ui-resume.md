## ACTIVE — resumed by user 2026-10-18

停止指示は解除済み。保存済み Simulator Release 成果物で UI 検証を再開し、iPhone arm64 ビルドへ進む。以下の停止記録は履歴。

## STOPPED — 2026-10-18 user request

実行中の UI テストを中断（終了コード 73）し、今回作成した compact Simulator を停止（終了コード 0）。未完了の結果バンドルは `compact-runtime-attempt1.xcresult` に保存。UI テストの合格証拠として扱わない。

Prepared drawing graph API は実装済みで、スタブではない。Rust 811 件、Swift Core 33 件に合格。描画グラフ 10,003 件は RDKit 比較で差分 0（14 件は描画サイズ制限）。Simulator Release build-for-testing は成功。最終 App UI 検証、generic iPhone arm64 ビルド、実機動作・速度測定は未完了。元の iPad の保存入力 6 件はバックアップと一致。

ビルド成果物と作業差分を再開用に保持。コミット・push・リリースは行っていない。詳細は `tasks/evidence/depiction-ipc-v43/user-stop-current.json`。停止指示を解除するまで作業を再開しない。

---

# ACTIVE — resumed on 2026-10-18

Prepared drawing API implemented; verification in progress. Focused RED failures
reproduced on resume. Existing 71 Rust sources match frozen v42 except the new
module declaration; no descriptor arithmetic changes.

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

# Chemob follow-up — 2D, information compatibility, physical speed

2026-10-18. Authorized by user after v42 report. Preserve immutable v42 evidence.

- [x] Quantify Ipc/AvgIpc instability with input-order controls and independent coefficients; retain separate stable names and make compatibility limits explicit.
- [ ] Replace circular connectivity preview with cached chemical 2D geometry and correct bond/atom notation. Test rings/chains/fragments and reject misleading geometry; explicitly report stereo limitations.
- [ ] Exercise the actual iPad/iPhone UI and inspect original-resolution images.
- [ ] Prepare/run bounded physical-device measurements when a paired device is accessible. Never label simulator results physical.
- [ ] Preserve normal Japanese iPad inputs, document evidence, clean owned compiler products.

Initial physical-device read: both paired devices unavailable, tunnel unavailable,
DDI services false; developer mode enabled. No pairing/signing settings changed.
No numerical Rust change is planned unless diagnostics prove one is necessary.
Ryer catalog empty; no agents delegated; indexed discovery through Tsubame.
