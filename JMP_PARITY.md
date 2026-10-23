# JMP / JMP Pro の計算機能差分と実装契約

初回調査日: 2026-10-21。既存APIとの対応を最後に照合した日: 2026-10-23。機能台帳の版基準は **JMP / JMP Pro 19.2** と取得した公開機能一覧。初期の24 SQL関数は限定的な解析基盤であり、JMPの計算機能全体を実装したものではない。現在版を含むJMP全体との完全同等性は主張しない。

全項目の台帳は [tasks/jmp-feature-matrix.json](tasks/jmp-feature-matrix.json)。14領域の公式ページに記載されたメソッド・機能の502出現を保持し、同じ名称を統合し、公式19.2リファレンスの具体的な下位機能も追加した。**552行**は粒度の異なる能力・手法・オプションを含む。件数をそのまま実装率に換算しない。公式リファレンスにはさらに細かなオプションがあるため、この台帳は追跡可能な下限であり、全オプションの完全な列挙ではない。

## 最新の照合結果（2026-10-23）

登録済みの解析APIは **267関数 / 拡張全体418 overloads**。台帳の現在の分類は、部分対応282行、未実装202行、DuckDB既存37行、従来の限定契約として検証済み6行、GUI除外17行、外部拡張4行、別製品4行、範囲未確定0行である。

既存APIを親項目にだけ関連付けていたため、下表の4項目が `absent` のままになっていた。ソース・Rustテスト・独立参照・過去のmatching native SQL検証を照合し、具体的なAPIと残差分を追記して `partial` に修正した。API側にも対応するfeature IDを追記した。アルゴリズムや公開関数の追加は0であり、JMPの方法群全体を完了にした変更ではない。

| 修正した項目 | 既存の限定実装 | 残る主な差分 |
| --- | --- | --- |
| GLMM binomial response | `chem_glmm_binomial` と `chem_glmm_predict`。Bernoulli-logit、Gaussian単一random interceptの積分MLと保存予測 | 集約binomial trials、他リンク、offset/重み、random slopes・複数random factors、JMP推定規約 |
| GLMM count response | `chem_glmm_poisson` と `chem_glmm_predict`。Poisson-log、Gaussian単一random interceptの積分MLと保存予測 | 他count family、offset/重み・overdispersion、複雑なrandom effectsと推論、JMP推定規約 |
| SEM fit indices | `chem_sem_fit_indices`。供給Gaussian covariance/meanモデルの適合指標 | 一般SEM推定、自動識別性判定、robust補正、modification indices、FIML・多群 |
| Text analysis | DTM・weighting・LSA・varimax・latent classと保存予測の既存8操作 | 文字列の言語処理、日本語辞書、sentiment、response-based term selection/text regression、巨大sparse corpus |

詳細とソースhashの照合結果は [対応付け監査記録](tasks/evidence/jmp-feature-mapping-audit-20261023.json)。GLMM・SEM・textの現行ソースは、それぞれ過去のnative SQL検証記録に保存されたhashと一致した。今回は台帳と比較文書の修正であり、Rust/C ABI/SQLテストの新しい実行結果は追加していない。[直近の実装検証記録](tasks/evidence/jmp-report-ordinal-validation.json)にあるRust 1,366件・SQL 603 assertions・418 catalog例の成功と、その当時の分類件数は履歴として保持する。

零過剰5系統、選択・正則化、標準誤差/信頼区間のJMP規約、図表の詳細オプションなどは引き続き未完了である。JMP本体との数値比較も未実施。台帳の残る `absent` 行も、APIが一切存在しないことを今回全件再検証した件数とは扱わない。

以下の各節の件数・「現在」「今回」は各調査・実装checkpoint当時の記録である。最新の集計は本節と台帳の `counting.status_counts` を参照する。

## 対象と扱い

- 以前の「GUI以外全部」という指示を踏まえ、純粋な画面操作は `excluded_gui` にした。
- 最新の「GPU以外」とGUIは別の条件。GPU専用でなければ、ニューラルネットを含む数値アルゴリズムはCPU実装の対象として残す。参照した公式一覧にGPU専用手法は見つからなかった。
- 静的図・帳票の出力は画面操作とは別である。再開前の `clarification_pending` 44項目は[非GUI受入契約](tasks/jmp-report-acceptance.json)へ具体化し、既存数値バックエンドの部分対応13行・未実装31行として対象に残した。図・帳票成果物の生成完了を新たに認定したものではなく、数値バックエンドも暗黙に除外しない。
- JMP Liveは公式に別製品。JMP ClinicalもJMP / JMP Proそのものと区別する。公式一覧にあるMarketplace拡張は `external_plugin` にした。
- DuckDBが既に提供する集計・相関・整形・結合・SQL入出力は `database_existing`。duckSMILESの新実装として数えない。外部形式・接続先は対象版で個別検証する。

版の区別は公式一覧の `icon-jmp-pro-tag` と凡例、または19.2リファレンスの明記に基づく。バッジが無い項目は、その一覧でPro制限が示されていないという意味であり、全てのオプションが標準版で使えることまで保証しない。ページ・プラットフォーム間で異なる場合は `option_or_platform_dependent`、根拠が不足する場合は `unknown` とした。

### 20.0 URLの確認と、版を混ぜないための記録

[20.0 helpのURL](https://www.jmp.com/support/help/en/20.0/jmp.html) は応答するが、調査時の取得本文は [19.2 help](https://www.jmp.com/support/help/en/19.2/jmp.html) と完全に同じだった。両方2,139,426 bytes、SHA-256 `0a83fa2ed5b842fba42fb6421f8b121b070a7d515591eff4a7d4094d279aa106`、固有 `.shtml` パス3,864件も一致した。この確認は目次の一致であり、全ページ内容や製品バイナリの一致を示すものではない。

一方、[公式ドキュメント入口](https://www.jmp.com/en/support/jmp-documentation) は19.2をCurrent Versionと表示し、[公式New Features / Release Notes一覧](https://community.jmp.com/t5/JMP-Knowledge-Base/JMP-New-Features-and-Release-Notes/ta-p/575202) の2026-09-23更新版にもJMP 19 / 19.2が掲載されている。20専用のWhat's Newは確認できなかった。したがって、URLが応答することだけから20.0の正式公開や新機能を認定せず、**20.0の追加差分は未確定**として追跡する。19.2基準の552行に推測の20.0項目を混ぜていない。

## 調査開始時の差分と受入条件

次表の「既存」は調査開始時の24関数の状態を指す。今回の追加実装は、[READMEの全SQL API](README.md#cpu-analytics-complete-sql-api-reference) と台帳の `current_sql_api` / 各項目の `current_implementation` に追記した。前回の登録は追加95関数、従来分と合わせて119関数で、その限定SQL契約の検証が完了した。8領域62操作、さらに6領域28操作を追加し、29操作の時点は254/405で、今回の打切り回帰・空間推定9操作を加え、現在の登録は **263 analytics関数 / 全414 overloads**。ソース・動作テスト・独立参照を関連付けた限定実装は `partial` のままとし、継続分のmatching native SQLゲートも合格した。公式の方法群全体の完了や、実行していない対象版/方法の同等性を推定しない。

| 領域 | 既存実装の限界と、残っている契約 | 公式根拠 |
| --- | --- | --- |
| 記述・分布 | 数値要約のみ。分布パラメータ推定、混合分布、適合度・正規性、重み/欠損を含む推定 | [基本解析](https://www.jmp.com/en/fragments/capabilities/basic-data-analysis-and-modeling) |
| 検定・区間 | 既存は両側Welch二群、一元配置ANOVA、BH/Holm。対応あり・一標本、同等性/非劣性、ノンパラメトリック、正確検定、多重比較、各種区間と効果量 | [基本解析](https://www.jmp.com/en/fragments/capabilities/basic-data-analysis-and-modeling) |
| 検出力・標本数 | 検定別の非心分布、割付比、最小整数標本数、信頼/許容区間・寿命試験・抜取検査の計画 | [基本解析](https://www.jmp.com/en/fragments/capabilities/basic-data-analysis-and-modeling) |
| 回帰・GLM | Ridge/PLS1と二次応答曲面は点推定が中心。カテゴリ/交互作用、重み、係数SE/CI、全体検定、診断、予測区間、分布・リンク・offset・Firth等 | [GLM 19.2](https://www.jmp.com/support/help/en/19.2/jmp/generalized-linear-models.shtml) |
| 一般化回帰・選択 | Lasso/Elastic Net、adaptive/二段階/Dantzig、選択経路、分布別尤度、正則化の選択と検証 | [推定法 19.2](https://www.jmp.com/support/help/en/19.2/jmp/estimation-method-options.shtml) |
| 混合・反復・空間 | 既存は均衡一元配置ランダム切片の分散成分だけ。ML/REML、非均衡・固定共変量・ランダム傾き、BLUP、反復/空間共分散、推論 | [混合モデル 19.2](https://www.jmp.com/support/help/en/19.2/jmp/mixed-models.shtml)、[共分散](https://www.jmp.com/support/help/en/19.2/jmp/statistical-details-for-repeated-covariance-structures.shtml) |
| GLMM | 非正規応答＋ランダム効果、推定法、条件付き/周辺予測と適切な不確実性 | [GLMM 19.2](https://www.jmp.com/support/help/en/19.2/jmp/generalized-linear-mixed-models.shtml) |
| 多変量・化学空間 | 既存PCA、k-means、PLS1。階層/SOM/混合、因子分析、MDS、判別、対応分析、潜在クラス、t-SNE/UMAP、多応答PLS/PLS-DA/VIP | [高度な解析](https://www.jmp.com/en/fragments/capabilities/advanced-statistical-modeling)、[PLS 19.2](https://www.jmp.com/support/help/en/19.2/jmp/partial-least-squares-models.shtml) |
| SEM | 測定誤差・媒介・潜在成長・多群、不変性、識別性と適合指標 | [SEM 19.2](https://www.jmp.com/support/help/en/19.2/jmp/structural-equation-models.shtml) |
| 非線形・曲線 | 一般非線形、溶出・シグモイド・指数・ピーク・薬物動態、スプライン、モデル推定と比較 | [高度な解析](https://www.jmp.com/en/fragments/capabilities/advanced-statistical-modeling) |
| 時系列 | AR/MA/ARIMA、季節差分・分解、状態空間/指数平滑、伝達関数、診断、時間順holdbackと予測区間 | [時系列 19.2](https://www.jmp.com/support/help/en/19.2/jmp/time-series-analysis.shtml)、[分解](https://www.jmp.com/support/help/en/19.2/jmp/differencing-and-decomposition.shtml) |
| 関数データ・スペクトル | 関数基底/FPCA、Fourier/wavelet、SG微分、SNV/MSC、背景補正、DTW、ピーク、functional DOE | [FDE 19.2](https://www.jmp.com/support/help/en/19.2/jmp/functional-data-explorer.shtml)、[前処理](https://www.jmp.com/support/help/en/19.2/jmp/data-processing-options.shtml) |
| DOE | 現状2水準完全要因、CCD、候補集合のD基準局所交換。部分実施、DSD、混合・カテゴリ、ブロック/split/strip、Taguchi、空間充填、covering、augmentation、alias・検出力 | [DOE](https://www.jmp.com/en/fragments/capabilities/design-of-experiments) |
| 最適化・因果 | GPは固定RBF・既知ハイパーパラメータで有限候補EIのみ。連続・バッチ・制約/多応答、シミュレーション・trade space、傾向スコア・マッチング・因果推論 | [高度な解析](https://www.jmp.com/en/fragments/capabilities/advanced-statistical-modeling) |
| 品質・MSA・抜取 | Cp/Cpkはwithin sigmaの外部指定。管理図のsigma/限界・ルール、変動/工程スクリーニング、測定者×部品のGage R&R、bias/linearity、agreement/kappa、OC | [品質](https://www.jmp.com/en/fragments/capabilities/quality-and-process-engineering) |
| 信頼性・生存 | KMは右打切りとGreenwood分散のみ。群比較、Cox、打切り尤度寿命分布、加速・累積損傷/疲労・劣化、成長/再発、RBD/RSS、試験計画 | [信頼性](https://www.jmp.com/en/fragments/capabilities/reliability-analysis) |
| 機械学習 | 木・forest・boosting、kNN/NB/SVM、CPU neural、分位点、モデル選択/平均/ensemble、確率・ROC/PR/lift、漏洩の無い検証・保存・score | [予測](https://www.jmp.com/en/fragments/capabilities/predictive-modeling-and-machine-learning)、[Neural 19.2](https://www.jmp.com/support/help/en/19.2/jmp/neural-networks.shtml) |
| テキスト・市場 | 文書特徴、LSA/潜在クラス、用語選択・感情・回帰、association、choice/conjoint/MaxDiff、uplift、調査分析 | [市場](https://www.jmp.com/en/fragments/capabilities/consumer-and-market-research) |
| ゲノム | 関連度・marker statistics/imputation、集団構造、GWAS/GBLUP、形式取込とcoding | [高度な解析](https://www.jmp.com/en/fragments/capabilities/advanced-statistical-modeling) |
| データ・自動化・連携 | 欠損パターン/補完・外れ値、モデルmetadata/score code、設定/前処理/seedの永続化、形式別入出力・外部システム連携 | [データ整形](https://www.jmp.com/en/fragments/capabilities/data-blending-and-cleanup)、[アクセス](https://www.jmp.com/en/fragments/capabilities/data-access)、[自動化](https://www.jmp.com/en/fragments/capabilities/automation-and-scripting) |

## TDDで進める単位

台帳の `acceptance_contracts` に、各領域で決める入力・仮定・出力・独立オラクルを記載した。全機能に一括で「テスト済み」を付けず、項目ごとに次を満たして昇格する。

1. 片側/両側、分布・リンク、欠損/重み/カテゴリ、打切り・同順位、標準化と学習範囲など、数学上の契約を先に固定する。
2. 閉形式の小例、公式例、独立した検証用実装に基づく期待値のテストを書く。未実装のために落ちるREDを実行し、コマンド・失敗理由を記録する。
3. 最小限の実装でGREENにし、数値境界・特異性・非有限値・入力形状・計算量上限を追加する。不適切な近似を方法名で隠さない。
4. 適切な不変条件も確認する。例: 線形モデルの平行移動、PCAの部分空間/再構成、確率和、設計の実行可能性、検出力の単調性、train/validation漏洩防止。
5. Rustだけで完了にしない。対応するビルド済み拡張をDuckDBにロードし、SQL型、NULL/選択ベクトル、JSON、保存/予測の往復、C ABIの境界まで検証する。
6. 証拠にrevision、OS/architecture、コマンド、数値許容差、独立参照の版、生成物hashを記録し、台帳のその限定契約だけを更新する。

例えば `ridge(lambda=0)` で係数が求まることと、OLSのSE/CI/F検定・カテゴリ効果・診断・予測区間が揃うことは別々に検証する。現行のKMと混合モデルについても、親領域の状態は `partial` のままにする。

この調査ではJMP本体を実行していない。JMPとの数値差分検証は未実施であり、同じ名称や独立参照との一致だけでJMPの全規約との一致を断定しない。

## 検証状態の更新

この台帳の既存実装評価は [前回の検証記録](tasks/evidence/sql-analytics-validation.json) と [独立参照](tasks/evidence/sql-analytics-reference.json)、実際のソース/テストに基づく。前回の865 Rustテスト・367 SQL assertions・18独立fixtureは、その記録のrevisionに対する証拠である。

今回の限定実装は、対象ソース・テストと利用可能な独立参照を関連付けたうえで `partial` とした。台帳552行のうち、関連限定実装/既存数値バックエンドがある項目は267行、未実装217行、ホスト既存37行、従来の限定契約として検証済み6行、GUI除外17行、範囲未確定0行、外部拡張4行、別製品4行。図・帳票44行は[非GUI受入契約](tasks/jmp-report-acceptance.json)から13行を既存数値バックエンドの部分対応、31行を未実装と具体化した。静的成果物の生成完了を認定した行は0で、実装済みSQL関数の件数は増えていない。これらは粒度が異なり、実装率やJMP同等性の割合ではない。

各追加SQL関数の入力・パラメータ・実行例・ソース/テスト・公式項目との関係を `current_sql_api` に保持した。直接Rust参照の証拠は [分布/ARIMA](tasks/evidence/jmp-distribution-arima-reference.json)、[回帰/多変量/混合](tasks/evidence/jmp-model-reference.json)、[検定/時系列/Cox](tasks/evidence/jmp-inference-timeseries-reference.json)。[追加DOE参照](tasks/evidence/jmp-advanced-design-reference.json) と [最終実SQL検証](tasks/evidence/jmp-sql-validation.json) も関連付けた。前回119関数のゲートはRust 1,006件（40 suites、失敗0、既存ignored 2）、native SQL 8 cases（skipなし）、catalog 270 overloadsと全270例の実行が合格。95追加操作のSQL oracleは199数値フィールド、保存モデル再利用・不正入力・query-label除外・選択/動的vectorも確認した。検証範囲はmacOS arm64 CPUの登録された限定契約であり、名称の一致をもとにJMPの方法群全体を完了にはしない。

版基準の確認元: [JMP documentation](https://www.jmp.com/en/support/jmp-documentation)。機能分類の入口: [JMP capabilities](https://www.jmp.com/en/software/capabilities)。ホストDBに既存の機能は [DuckDB aggregates](https://duckdb.org/docs/current/sql/functions/aggregates)、[PIVOT](https://duckdb.org/docs/current/sql/statements/pivot)、[UNPIVOT](https://duckdb.org/docs/current/sql/statements/unpivot)、[Friendly SQL](https://duckdb.org/docs/current/sql/dialect/friendly_sql) と区別した。

## 継続62操作の範囲と未完了の契約

次の追加実装と個別SQL例・引数・保存パラメータ・制限は、README内にまとめている。APIと証拠のリンクは再現用であり、READMEの関数一覧を置き換えない。直接Rust参照、公開C ABI、native DuckDB SQLは区別する。

| 領域 | 追加操作 | 現行の限定契約 | 残る主な契約 |
| --- | --- | --- | --- |
| 潜在変数/クラスタ | 8 | Classical MDS、単純対応分析、主因子法、CCA、full-covariance GMM、矩形batch SOMと保存予測 | optimized/nonmetric MDS、全covariance族/成分選択、SOM全topology/学習規約 |
| 評価/検証分割 | 12 | weighted regression/binary metrics、ROC/PR/lift/calibration、profit、固定ensemble、供給予測比較、holdout/rolling indices | 自動再学習/選択・nested tuning・前処理漏洩provenance・全推論 |
| 非線形/平滑化 | 7 | 10曲線族LM、clamped B-spline、係数2階差P-spline、tricube LOESSと保存予測 | 任意式/制約/重み/profile likelihood、全curve族、自動knots/lambda、JMP smoothing spline規約 |
| 関数/スペクトル | 10 | SNV/MSC、SG/差分、polynomial detrend、Haar/DFT、2系列DTW、共通grid FPCA、離散peak | irregular/sparse/平滑FPCA、自動基底選択/新曲線score、他wavelet/背景法、parametric peak、functional DOE |
| 多変量推論 | 4 | whole-model MANOVA、逐次CCA Bartlett近似、MCA調整慣性、ML factor+varimax | 任意effect/反復応答、regularized CCA推論、補助category/weights、斜交/全factor推論 |
| SVM/分位点/robust | 5 | 二値linear/RBF C-SVCのKKT付きSMO、check-loss ADMM、固定scale Huber IRLSと保存予測 | multiclass/SVR/calibration、罰則quantile/推論、scale同時推定・高breakdown S/MM/Cauchy |
| Diploid markers | 9 | 統計/HWE近似、imputation、VanRaden/dominance、距離、unlinked simulation/cross、固定ratio GBLUP | exact HWE、phased/recombination/polyploid、admixture、REML variance/GWAS/形式parser |
| 信頼性 | 7 | 右打切りAFT3族とsurvival/quantile、single-system time-truncated power NHPP、独立series/parallel/k-of-n blocks | 左/区間打切り・entry・frailty、全stress law/scale model、複数system/phase、任意graph依存/repair simulation |

181関数の過去checkpointでは台帳551行を維持し、関連した限定方法だけを `partial` に結び付けた。Canonical CCAの個別操作など、既存台帳に単独行がないものは `current_sql_api` の注記に残す。方法群の完了に見せるために、異なるオプションを暗黙に実装済みにはしない。

継続分のnative検証は [scripts/check_jmp_continuation.py](scripts/check_jmp_continuation.py) をmatching buildに対して実行して合格した。全62例の厳密JSON/SQL型境界、47数値・構造fixture操作の107結果field、11保存モデル往復、16不正入力NULL、6選択/dynamic行、7 interleaved NULL行、2 query-label除外を確認した。9queriesで134 JSONを復号した。181関数の過去checkpointの全体ゲートは **Rust 48 suites / 1,099 passed / failed 0 / 既存ignored 2**、native SQL **9 cases / 394 assertions / skip 0**、catalog **332 overloads / 全332例実行**。Strict Clippy・対象format・diff検査も合格した。[最終証拠](tasks/evidence/jmp-continuation-validation.json) とREADMEの独立参照表に実行層・数え方を分けて記録した。前回1006テスト/270overloadsは過去snapshotとして保持する。


## 2026-10-22: 209関数の検証checkpointと継続作業

28操作（混合/GLMM/因果/omics/古典DOE/market）を追加し、**209分析関数 / 全360 overloads**へ更新した。Rust **54 suites / 1,167 passed / failure0 / existing ignored2**、native SQL **11 cases / 420 assertions / skip0**、全360catalog例を実行済み。22-model checkerは22fixture操作57field、6保存model往復、15不正NULL。DOE checkerは6例、13構造case、14不正NULLに合格。検証層・制限・再現手順は[証拠](tasks/evidence/jmp-models-validation.json)と[README](README.md#cpu-analytics-complete-sql-api-reference)に記録。

公式MultivariateのItem ReliabilityをIRT Item Analysisと分けた追加下位項目として1行追加した。14capability断片の502出現は維持し、現在552行（partial 227, database_existing 37, absent 213, external_plugin 4, clarification_pending 44, implemented_verified 6, excluded_gui 17, separate_product 4）。粒度が異なるため率にしない。上の181関数/551行の記録は過去checkpointである。IRT、SEM、劣化/安定性、その他の未完契約は継続中で、**全体達成にはしていない**。


## 225関数checkpoint: IRT・SEM・劣化・ゲノム分散

16操作を追加し、225analysis関数/376overloadsをmatching nativeで検証した。Rust 1236tests/58suites、SQL 12cases/430assertions、全376catalog例が成功。保存モデル7往復とsubnormal予測の修正も実SQLで確認。552行の現在状態は partial 235, database_existing 37, absent 205, external_plugin 4, clarification_pending 44, implemented_verified 6, excluded_gui 17, separate_product 4。[証拠](tasks/evidence/jmp-specialized-validation.json)と[README](README.md#cpu-analytics-complete-sql-api-reference)に全契約・制限・別計算との照合を記録した。IRT、SEM、劣化、genomicモデルは限定部分対応で、全JMP同等は未達。時系列・残るfunctional/DOE等の必要作業を継続する。


## 254関数checkpoint（過去の停止時点）

今回の29操作を統合し、254analytics関数/全405overloadsをmatching nativeで検証。Rust1302tests/63suites、SQL13cases/439assertions、全405catalog例、29操作のproduction C ABI既知値、保存model/referenceの9SQL往復が成功した。数値レビューのGaussian covariance入力とspatial予測/半分散/物理共分散の不具合は実RED→GREEN回帰と独立計算で検証した。[最終証拠](tasks/evidence/jmp-streams-validation.json)、[全関数README](README.md#cpu-analytics-complete-sql-api-reference)。552行の状態は partial 254, database_existing 37, absent 186, external_plugin 4, clarification_pending 44, implemented_verified 6, excluded_gui 17, separate_product 4。19行をabsentからpartialへ変更したが、**全JMP同等は未達**。ユーザーの「それ完了したらpushでストップ」により今回の検証済み範囲をpush後、追加実装を停止する。choice/censored-regression等の次段階コードは追加していない。


## 2026-10-22: 全体目標の再開

ユーザーが全体目標を再開した。上のpush後停止は過去のcheckpointの記録であり、現在の作業状態ではない。main `442f54588d5882a3fc14f3f830c2d29d20701ace` への反映、Community CIの3対象、ソースCIの9対象の成功は、既存254関数の提出・互換性に対する証拠であり、JMP全体の達成を証明しない。再開時点の台帳は部分対応254行・未実装186行・出力範囲未確定44行だった。その後、44行の受入契約を具体化し、既存バックエンド部分対応13行・未実装31行へ分類したため、現在は部分対応267行・未実装217行・範囲未確定0行である。成果物の実装完了を追加認定したものではない。

継続中の契約は、Normal/Exponential/Gamma/Weibull/Lognormal/Betaの左・右・区間・正確観測を含む打切り回帰と、空間共分散パラメータのML/REML推定・推論・保存予測である。各契約のRED、独立数値参照、公開C ABI、実DuckDB SQLが揃うまでは検証済みの関数数を増やさない。限定的な推定法を追加しても、異方性・全推論オプションなどの残差分は保持する。

現目標が除外するのはGUI操作とGPU依存だけである。静的図・帳票のCPU数値バックエンドと出力は対象に残し、従来の44件の範囲未確定項目を一次資料から具体化する。未確定であることを除外や達成として扱わない。全体目標はactiveのまま継続する。

図・帳票の全44契約は [tasks/jmp-report-acceptance.json](tasks/jmp-report-acceptance.json) に、入力、計算/構造出力、静的artifact、非GUIオプション、一次資料、独立検証と残調査を記録した。範囲の具体化は未実装成果物の代わりにはならず、全体達成の根拠にしない。

## 再開後の検証済み9操作（全体は継続中）

6分布のexact/left/right/interval打切り回帰と保存予測7操作、および空間共分散ML/REML推定と保存予測2操作を追加した。契約・全引数・SQL例・制約はREADMEと `tasks/jmp-censored-api.json` / `tasks/jmp-spatial-estimation-api.json`。JMP本体との直接比較は行っておらず、該当方法群は `partial` を維持する。

- Rust: 65suite /1,340成功 /0失敗 /既存ignore2。厳格Clippy・scopedfmt・diff検査合格。
- 実Release拡張: matching DuckDB2 SDKで16 SQL cases /513assertions /skip0、全414カタログ実行例合格。
- 打切り回帰: 独立SciPy参照・実C ABI・実SQLの768検査。Gamma微小上側、有限CI、非対称保存共分散を含む。
- 空間推定: 実C ABI48calls /810fields、実SQL16独立fit /130fieldsと3保存再利用。

[集約証跡](tasks/evidence/jmp-parity-estimation-validation.json)。台帳は552行のまま、partial275 /absent209 /implemented_verified6 /database_existing37 /excluded_gui17 /external_plugin4 /separate_product4。44帳票項目の受入条件を具体化したが、このcheckpointでは新しい帳票成果物の完成は認定していない。全体達成率には換算しない。次の帳票契約は [histogram/box plotの調査](tasks/jmp-report-next-contracts.json) で具体化中。

化学DQの[別調査](tasks/chemistry-dq-research-20261022.json)は保存した。ユーザー指定により実装開始はJMP全体の完了後とし、DuckOrchの併用とDuckSMILESへの基盤組込みの両案を残す。
