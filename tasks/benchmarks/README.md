# Chemob / RDKit 検証

対象はMac上のRustエンジンと、実際の `apps/Chemob/App/Chemistry.swift`。Simulator・実機の速度とは別の測定です。

## 環境とC ABI

```sh
uv venv /private/tmp/chemob-rdkit-env --python 3.13
uv pip install --python /private/tmp/chemob-rdkit-env/bin/python rdkit==2026.3.6
CARGO_TARGET_DIR=/private/tmp/chemob-rust-target cargo build --release -p ducksmiles_smiles
cc -dynamiclib -Wl,-force_load,/private/tmp/chemob-rust-target/release/libducksmiles_smiles.a \
  -framework Security -framework CoreFoundation -liconv \
  -o /private/tmp/chemob-rust-target/release/libchemob.dylib
/private/tmp/chemob-rdkit-env/bin/python tasks/benchmarks/compare_rdkit.py \
  --library /private/tmp/chemob-rust-target/release/libchemob.dylib \
  --output tasks/evidence/cache-mass --strict
```

`--strict` は数値・分子式・受理状態に不一致があれば終了コード1。現在の互換性ゲートは失敗しており、合格扱いにしません。`--trials 0` は速度測定を省略します。

RDKitのwheelに含まれる `Data/NCI/first_5K.smi`（4,999行・4,900種類）と20件の境界入力を使用。ファイルハッシュ、バージョン、誤差許容値はJSONに保存します。データの出典は [RDKit NCI data](https://github.com/rdkit/rdkit/tree/Release_2026_03_6/Data/NCI)。

RDKitは1行ごとに新しく分子を作り、sanitizeを有効、`removeHs=False` として明示H頂点を保持します。TPSAはN/O、回転可能結合はStrict。定義は [公式API](https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) を参照。塩除去・中和・重複キャッシュは行いません。

## 入出力

`apps/Chemob/Benchmarks/BulkIO.swift` を実際の `App/Chemistry.swift` と同じターゲットでコンパイルします。UIKitを使わないモデル部分をmacOSから実行します。Swift Coreを先に生成し、その後リンクします。

```sh
mkdir -p /private/tmp/chemob-swift-host
swiftc -O -parse-as-library -emit-library -emit-module -module-name ChemobCore \
  apps/Chemob/Core/Sources/ChemobCore/CompoundImport.swift \
  apps/Chemob/Core/Sources/ChemobCore/EngineInput.swift \
  apps/Chemob/Core/Sources/ChemobCore/MoleculeMemo.swift \
  -emit-module-path /private/tmp/chemob-swift-host/ChemobCore.swiftmodule \
  -o /private/tmp/chemob-swift-host/libChemobCore.dylib
swiftc -O -swift-version 6 -parse-as-library \
  -I /private/tmp/chemob-swift-host -L /private/tmp/chemob-swift-host -lChemobCore \
  -import-objc-header apps/Chemob/App/Engine.h \
  apps/Chemob/App/Chemistry.swift apps/Chemob/App/ChemistryStore.swift \
  apps/Chemob/Benchmarks/BulkIO.swift \
  /private/tmp/chemob-rust-target/release/libducksmiles_smiles.a \
  -framework Security -framework CoreFoundation -liconv \
  -Xlinker -rpath -Xlinker /private/tmp/chemob-swift-host \
  -o /private/tmp/chemob-swift-host/chemob-bulk-io
/private/tmp/chemob-swift-host/chemob-bulk-io \
  tasks/evidence/cache-mass/bulk-input.csv tasks/evidence/cache-mass --ui-testing
/private/tmp/chemob-rdkit-env/bin/python tasks/benchmarks/validate_cached_io.py \
  tasks/evidence/cache-mass /private/tmp/chemob-rust-target/release/libchemob.dylib
```

10,000行の取込・計算・出力・再読込を実行し、PythonのCSV実装でも全12列と名前・SMILES・順序を照合します。最新の検証ではキャッシュの結果を未キャッシュのC ABI結果と比較します。`baseline-native-results.json.gz` は前の高速化時の履歴です。`--ui-testing` はDocumentsへの保存を無効にします。出力CSVには未対応入力も残します。速度には画面描画は含みません。

## 同条件の速度比較

`paired_rdkit_timing.py <before.dylib> <after.dylib> <output.json>` は、両者が受理する先頭1,000行（997種類）を使い、修正前・修正後・RDKitの順番を3回入れ替えて測ります。他のビルド終了後に実行。読み込み済み分子を使い回さず、12記述子・分子式・JSONの直列化と復号を含みます。

修正前ライブラリは、現在の差分からPreparedMatcher導入だけを外したソースを別チェックアウトでビルドして再生成できます。測定元のソースハッシュと結果は `tasks/evidence/rdkit` を参照してください。ホストの他の作業は停止していないため、近い速度差を断定的に解釈しません。

## 質量参照データとキャッシュ

`generate_mass_data.py` はRDKit 2026.03.6の周期表データから、118元素と3,111同位体の質量をRust定数へ出力します。生成時だけRDKitを使い、配布アプリの実行時は不要です。数値データに付随する上流のライセンス通知は `crates/smiles/RDKIT-DATA-LICENSE.txt`。

`CacheIntegration.swift` は本物の `ChemistryStore` を使い、重複計算の回避、行位置の再割当、接続図の再利用、ディスク再読込、旧エンジン版のキャッシュ無効化、同位体・電荷・立体の区別を検証します。BulkIOと同じコンパイル手順で末尾のSwiftファイルを置き換えて実行できます。

計算仕様を更新するときは `workbench::REVISION` を更新してください。キャッシュのキーはこの版と完全一致SMILESです。異なるSMILES表記の同一性推定は行いません。現在保存する描画用データは接続図で、本格的な2D構造式ではありません。
