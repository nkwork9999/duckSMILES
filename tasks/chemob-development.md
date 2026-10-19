# Current descriptor milestone — 2026-10-18

The v42 descriptor/App/data milestone has passed its declared gates:217RDKit-derived
2D definitions(215pinned numerical +2same-formula exact stable information),24existing
profiles,50additional defined nonregistered descriptors. Full100000 original inputs:
allaccepted/all291+formula matching declared references,0differences/referenceerrors.
PinnedIPCcompatibility differs5723rows and remains explicitly separate.

Actual App10000 input/output,295CSVcolumns,2806165finitecells/undefined/source/order,
diskreopen10000reuse0newcalculation,old-cache invalidation,iPad2/phone2UI,normalJA6
restoration,andunsigneddevicebuild pass. The observed fixed1000same-output291 workflow
is3.9754s native vs19.1658s RDKit+helperworkflow onMacM1,3alternatingpairs.

[Current evidence/scopes](evidence/catalog289-sps/REPORT-v42.md).
Remaining product work:physicaldevice runtime/speed when available; properchemical2D;
translation review; distribution/signing when authorized.3D is excluded. No commit/push/release.

--- Earlier work log ---

# Current checkpoint — 2026-10-18 / v37

291項目（RDKit由来217定義＋既存24＋追加50）をアプリに反映。10万入力のうち
受理99,752件は定義した参照・誤差内で全値一致。未対応248件が残り、入力受理の
完全互換基準は未達。IPC系2項目は安定計算を採用し、固定RDKitの数値結果とは
5,617件で差があるためCSV列名も区別する。実際の1万件入出力・保存復元・
iPad/iPhone Simulator画面検証は通過。実機向けビルドは成功、端末速度は未測定。

[最新の検証報告](evidence/catalog289-sps/REPORT.md) / [再開情報](descriptors/RESUME.md)

---

以下は初期方針と過去の実装計画。

# Chemob と Backbeat の開発

2026年10月17日の方針。ChemobをiPhone/iPadから使える化学ワークベンチとして作り、DuckSMILESの計算基盤とBackbeatの素材制作・書き出しを同じ実例で育てる。

## この実装単位

- DuckSMILES: 1分子1解析による基本12記述子の一括API、上限付き入力、行位置とエラー状態の保持、安全なC ABI。
- Chemob: SMILES/CSV/TSV取込、端末内計算、化合物一覧、計算値詳細、CSV出力、入力の復元。キャッチーなブランド画面と短いアニメーション。
- Backbeat: 編集可能なChemob原本、新しい透過PNGアトラス書き出し。編集→書き出し→Chemobへの同梱を再現できる形にする。

Chemobの初期の接続図は、原子と結合を確認するための簡易プレビュー。本格的な2Dレイアウトと立体表現の実装は次の段階。3Dは今回のアプリ機能には含めない。

## 次に測って改善するもの

1. 新APIと従来スカラー呼出しで、同じ12記述子の結果・入力数・時間を比較する。
2. 隣接リスト、環情報、原子型付けを共有し、依存計算の重複を減らす。
3. SDF取込、原本列の保持、記述子プランの選択、プロジェクト保存を増やす。
4. 本格2D描画の品質と可視範囲だけを生成するキャッシュを実装する。
5. Backbeatの素材書き出しを実例で改善し、アプリ内の操作からも利用しやすくする。
6. 部分構造・類似度検索とRDKitの機能別互換性を広げる。

## リポジトリ境界

- Chemob / DuckSMILES: `/Users/nk/Desktop/ducksmiles-workbench-20261017`、`feat/ios-workbench-20261017`。
- Backbeat: `/Users/nk/Desktop/renshu/backbeat`。既存の未コミット作業を維持し、追加機能は `tools/backbeat-cli` 内に限定。
- c-acheの実装は複製しない。今回はBackbeatの既存キャンバス描画経路を使う。

検証結果と未検証事項は `apps/Chemob/VERIFICATION.md` に記録する。

## 再開時の進捗

- サンプル→計算詳細→不正入力の自動UIテストが成功。
- CSV欄内の改行でエンジン行対応がずれる問題を修正し、Swiftの5テストが成功。
- 再現可能なMac上の一括計算ベンチマークを追加。初期基準値と測定条件は `apps/Chemob/VERIFICATION.md`。
- 本格2D描画は未実装。次は環・縮合環・分岐・塩の配置品質を独立したエンジンAPIとテストで育て、一覧の可視分子だけを生成する。

## 大量SMILES比較での追加進捗

- NCI 4,999行＋境界20入力をRDKit 2026.03.6と照合。互換性ゲートは未合格。修正優先は質量表・同位体・電荷補正、HBA/TPSA/LogPの差。
- LogP照合の隣接情報再構築を減らし、同条件1,000行で中央値3.874秒→0.592秒。RDKitは0.701秒、速度範囲は重なる。
- 実際のSwift Workbenchで10,000行の入出力を検証。結果CSVが1MiBを超え再読込できなかった問題を修正し、8MiB/10,000行へ。単一SMILESの上限は維持して他行への影響を防ぐ。
- Rust 727 tests、Swift 8 tests成功。詳細: `tasks/evidence/rdkit/REPORT.md`。

## キャッシュと質量整合の進捗

- 同一SMILESの記述子・接続図を版付きで保存。追加取込・再起動で再利用。異なる立体・同位体・電荷を統合しない。
- 質量表118元素/3,111同位体、電子質量補正、リンの暗黙水素、炭素なし分子式の順序を修正。共通受理5,002行で分子量・精密質量・分子式は一致。
- 値の不一致1,400→28行、受理状態差10行は残る。全互換とはしない。
- Rust 730 / Swift 10 / Simulator UI 1が成功。通常起動の保存・再起動も確認。
- 最新結果: `tasks/evidence/cache-mass/REPORT.md`。計算意味を変えるたび `workbench::REVISION` を更新する。
