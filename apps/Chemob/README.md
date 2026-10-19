# Chemob

iPhone／iPadでSMILESをまとめて計算する、端末内の化合物ワークベンチです。
DuckSMILESのRustエンジンを組み込み、アプリ実行時のPython・RDKit導入やネット接続は不要です。

## 現在の機能

- SMILES貼付け、CSV／TSV／SMIの読込。最大10,000行・32MiB。
- 291記述子: 固定RDKit由来217定義＋既存24プロファイル＋追加50。
- 名前／元SMILESの検索、状態カテゴリ、数値による絞込み、CSV共有。
- 入力の端末内保存、同じ元SMILESの結果・接続図の再利用。計算版を分けて保存します。
- iPadは検索・カテゴリ・ページ付き一覧を上半分、選択した化合物・値・状態を下半分に表示。
  通常文字サイズでは主画面にスクロール領域を置きません。iPhoneと大きい文字では適応した縦画面を使います。

v43では、化学計算と同じ前処理から2D構造式を生成し、電荷・同位体・付随水素・くさび／破線・二重結合の配置を表示します。
複雑な環や大きな構造など、描画できない場合も入力と計算結果を保持します。一部の立体配置を省く場合は表示で知らせます。
現行AppのiPad／iPhone Simulatorで2D表示と列配置・縦横切替を確認しています。3D、SDF取込、部分構造検索のUIは未実装です。
3Dは現在の開発対象から外しています。

## 検証済みの範囲

計算版 `ducksmiles-workbench-v42-cluster291-rdkit2026036`。
異なる元入力100,000件を全件受理し、全291項目・分子式が各検証基準内で一致しました。
元入力をネイティブの受理結果で選別していません。
実際のSwift／iPadで10,000件を取り込み、295列のCSV、元入力・順序・全数値を確認。
保存後の10,000件再利用は追加計算0で、全復元値が一致しました。
v42時点でiPad 2件＋iPhone 2件のUI試験、Simulator／実機向けarm64 Releaseビルドが合格しています。
v43はSimulator Release build-for-testingと署名なしiPhone arm64 Releaseビルド合格。
iPhone SE3で2D・E/Z・橋かけ構造・列配置、iPadで2D・立体省略表示・縦横の固定パネルの試験が合格しています。

215標準項目は固定RDKit 2026.03.6の数値と照合。`Ipc`／`AvgIpc` の2項目は同じ数式を
整数演算で安定化したもので、固定RDKitの数値とは5,723入力で差があります。
CSVでは `Ipc_stable`／`AvgIpc_stable` として出力します。RDKitで学習したモデルの入力ではこの違いを区別してください。
追加50は監査した標準記述子APIと重複しない既存の数式で、予測性能を保証するものではありません。

実機の実行・速度は、登録済み端末の接続トンネル／DDIサービスが利用できず未測定です。
[精度・データ・測定条件・全証拠](../../tasks/evidence/catalog289-sps/REPORT-v42.md)、
[検証の要約](VERIFICATION.md)、[v43描画・立体・キャッシュの検証](../../tasks/evidence/depiction-ipc-v43/REPORT-v43.md)。

## ビルド

Apple Siliconの開発環境、Xcode、XcodeGen、RustのAppleターゲットを使用します。

```sh
cd apps/Chemob
xcodegen generate
xcodebuild -project Chemob.xcodeproj -scheme Chemob -configuration Release \
  -destination 'generic/platform=iOS Simulator' \
  -derivedDataPath /private/tmp/chemob-descriptor-derived \
  CHEMOB_CARGO_TARGET_DIR=/private/tmp/chemob-descriptor-apple-target \
  CODE_SIGNING_ALLOWED=NO ARCHS=arm64 ONLY_ACTIVE_ARCH=YES build
```

pre-buildが同じリポジトリのRust staticlibを生成します。署名設定は変更しません。
同じエンジンを生成済みの場合だけ `CHEMOB_USE_PREBUILT_ENGINE=1` を指定して再利用できます。
プラットフォームを切り替える際は実行終了・証拠保存後に専用の生成物を片付け、各構成の出力を重複させないでください。

Coreのホスト検証は `swift test --package-path Core --scratch-path /private/tmp/chemob-descriptor-target/swift-host`。
UI試験は `Chemob` スキームにあります。1万件データと出典はUIテストのbundleだけに含め、本番Appへは入れません。

## キャッシュと表示

元SMILESをキーに、記述子10,000件・接続図256件を保存します。表記の異なる分子の同一性推定は行いません。
電荷・同位体・立体を保持し、エンジンの版が変わった結果は再計算します。
未定義値はダッシュ、CSVでは空欄です。不正・未対応の入力も元の文字列を保持します。

ブランド用の一度だけ再生するアニメーションは、Backbeat側の編集可能な
`tools/backbeat-cli/examples/chemob-intro.backbeat.json` とatlas書出しを原本にしています。
Appの `Resources/chemob-intro.png`／`.json` がその資産です。Reduce Motionでは静止表示します。

## International interface

English source strings and Japanese, Indonesian, Hindi, German, French and Italian dictionaries
share 224 keys and formatting arguments. Switzerland uses German/French/Italian with English fallback;
Romansh and other regional languages of India are not included. Wording is draft and needs native-speaker review.
The OS selects the app language. User names, SMILES, descriptor abbreviations and CSV machine column names are preserved.
No translation service is called at runtime.
