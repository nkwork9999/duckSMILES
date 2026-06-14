# Phase 1: 3D Conformer Generation — 実装計画

branch: `feature/docking-vina`

## 目的

`smiles_to_3d(smiles VARCHAR) → STRUCT(x DOUBLE[], y DOUBLE[], z DOUBLE[])` を DuckDB 関数として公開し、後続のドッキング計算（Phase 2–4）の入力とする。

## アルゴリズム選択: Distance Geometry (DG)

### 選択理由

| 手法 | 精度 | 実装コスト | 外部依存 |
|---|---|---|---|
| ETKDG (RDKit default) | 高 | 高（実験的トーション知識DB必要） | なし（ただしDBが大きい） |
| Distance Geometry | 中 | 中 | なし |
| Fragment template assembly | 中 | 高（テンプレートDB必要） | なし |
| ランダム+FF最適化のみ | 低 | 低 | なし |

DG は外部 DB 不要・純 Rust 実装可能・ドッキング前処理として必要十分な精度。

### アルゴリズム概要

```
SMILES
  │
  ▼
[1] 距離制約行列構築 (bounds.rs)
    ├── 1-2 結合: 正確な結合長 (lower = upper = r_bond)
    ├── 1-3 ペア: 角度から計算 (law of cosines)
    ├── 1-4 ペア: トーション考慮 (cis/trans/ring)
    └── それ以外: VDW 下限 / 大きな上限
  │
  ▼
[2] 三角不等式平滑化 (smoothing.rs)
    Floyd-Warshall 変形:
    u[i][j] = min(u[i][j], u[i][k] + u[k][j])   ← 上限を縮める
    l[i][j] = max(l[i][j], l[i][k] - u[k][j])   ← 下限を広げる
  │
  ▼
[3] ランダム埋め込み (embed.rs)
    ├── l ≤ D ≤ u を満たすランダム距離行列をサンプル
    ├── 計量行列 G = -0.5 * H * D² * H  (H = I - 1/N * 11ᵀ)
    ├── G を固有値分解 → 上位3固有値/ベクトル
    └── 座標 = V₃ * sqrt(Λ₃)
  │
  ▼
[4] UFF-lite 力場最小化 (minimize.rs)
    ├── 結合伸縮: E = k_b(r - r₀)²
    ├── 角度変形: E = k_a(θ - θ₀)²
    ├── VDW (Lennard-Jones 12-6): E = ε[(r₀/r)¹² - 2(r₀/r)⁶]
    └── 最急降下法 100 step → 収束判定
  │
  ▼
Conformer { atoms: Vec<[f64; 3]> }
```

---

## ファイル構成

```
crates/smiles/src/
├── conformer/
│   ├── mod.rs        ← 公開 API: smiles_to_3d(), Conformer 構造体
│   ├── params.rs     ← 結合長・結合角・VDW 半径テーブル
│   ├── bounds.rs     ← [1] 距離制約行列構築
│   ├── smoothing.rs  ← [2] 三角不等式平滑化
│   ├── embed.rs      ← [3] ランダム埋め込み（固有値分解）
│   └── minimize.rs   ← [4] UFF-lite 力場最小化
└── lib.rs            ← conformer モジュールを追加
```

DuckDB FFI バインディングは `src/conformer/mod.rs` に `#[no_mangle] pub extern "C" fn duck_smiles_to_3d(...)` として追加。

---

## パラメータ仕様 (params.rs)

### 結合長 (Å) — BondOrder × 元素ペアの lookup

```
C-C  single: 1.540   C=C  double: 1.340   C≡C  triple: 1.200
C-N  single: 1.470   C=N  double: 1.270   C≡N  triple: 1.150
C-O  single: 1.430   C=O  double: 1.200
C-S  single: 1.820   C=S  double: 1.600
C-H         : 1.090   N-H: 1.010   O-H: 0.960
C-F         : 1.350   C-Cl: 1.770  C-Br: 1.940  C-I: 2.140
Aromatic C-C: 1.400   Aromatic C-N: 1.340
```

実装: `fn bond_length(a: &str, b: &str, order: BondOrder) -> f64`

### 結合角 (°) — ハイブリダイゼーション依存

```
sp3 中心: 109.5   (例: C-C-C, N-C-C, O-C-C)
sp2 中心: 120.0   (例: C=C-C, C-N=C)
sp  中心: 180.0   (例: C≡C-H)
水 (O-H): 104.5
アンモニア (N-H): 107.0
```

実装: `fn bond_angle(center: &str, hybridization: Hybridization) -> f64`

ハイブリダイゼーション判定:
- `sp3`: 4 重原子、または環外 N/O
- `sp2`: 二重結合を持つ、または芳香族環中
- `sp `: 三重結合を持つ

### VDW 半径 (Å)

```
H: 1.20   C: 1.70   N: 1.55   O: 1.52   S: 1.80
F: 1.47   Cl: 1.75  Br: 1.85  I: 1.98   P: 1.80
```

実装: `fn vdw_radius(symbol: &str) -> f64`

---

## 各モジュール詳細仕様

### bounds.rs

```rust
pub struct BoundsMatrix {
    pub n: usize,
    pub lower: Vec<f64>,  // n×n, lower[i*n+j]
    pub upper: Vec<f64>,  // n×n, upper[i*n+j]
}

pub fn build_bounds(mol: &Molecule) -> BoundsMatrix
```

**構築手順:**
1. 初期化: lower[i][j] = vdw_lower(i,j), upper[i][j] = 1000.0
2. 1-2 拘束: lower=upper= bond_length(i,j,order)
3. 1-3 拘束: 余弦定理で計算  
   `d_13 = sqrt(d_12² + d_23² - 2·d_12·d_23·cos(θ))`
   lower = upper = d_13
4. 1-4 拘束:
   - 二重結合周り (sp2): cis = d_cis, trans = d_trans (狭い範囲)
   - 単結合周り (sp3): lower = d_gauche_minus, upper = d_trans
   - 環内: 環サイズから距離を固定

### smoothing.rs

```rust
pub fn smooth_bounds(bounds: &mut BoundsMatrix)
```

Floyd-Warshall O(N³):
```
for k in 0..n:
  for i in 0..n:
    for j in 0..n:
      // 上限を縮める
      if upper[i][j] > upper[i][k] + upper[k][j]:
        upper[i][j] = upper[i][k] + upper[k][j]
      // 下限を広げる
      if lower[i][j] < lower[i][k] - upper[k][j]:
        lower[i][j] = lower[i][k] - upper[k][j]
```

N ≤ 100 の薬様化合物なら O(10⁶) 程度、問題なし。

### embed.rs

```rust
pub fn embed(bounds: &BoundsMatrix, seed: u64) -> Vec<[f64; 3]>
```

**手順:**
1. ランダム距離行列 D のサンプリング（LCG 乱数、lower ≤ D ≤ upper）
2. 計量行列変換:
   ```
   d²[i][j] = D[i][j]²
   G[i][j] = -0.5 * (d²[i][j] - mean_row[i] - mean_col[j] + grand_mean)
   ```
3. 固有値分解 (3×3 または Power Iteration):
   - Jacobi 法を実装（外部ライブラリ不使用）
   - 上位3固有値 λ₁,λ₂,λ₃ と固有ベクトル
   - 負固有値はゼロにクリップ
4. 座標: `pos[i] = [v1[i]*sqrt(λ₁), v2[i]*sqrt(λ₂), v3[i]*sqrt(λ₃)]`

固有値分解の実装: N×N 対称行列の Jacobi 法（収束基準 1e-10、最大 1000 反復）

### minimize.rs

```rust
pub fn minimize_uff(mol: &Molecule, coords: &mut Vec<[f64; 3]>, max_iter: usize)
```

**エネルギー関数:**
```
E_total = Σ E_bond + Σ E_angle + Σ E_vdw

E_bond  = 0.5 * k_b * (r - r₀)²          k_b = 700 kcal/mol/Å²
E_angle = 0.5 * k_a * (θ - θ₀)²          k_a = 100 kcal/mol/rad²
E_vdw   = ε * [(r_min/r)¹² - 2(r_min/r)⁶]  ε = 0.1 kcal/mol
```

**最適化:**
- 数値勾配 (中心差分、h=1e-5)
- 最急降下法、学習率 α=0.01
- 収束: |ΔE| < 1e-6 または max_iter 達成
- 水素原子は質量ゼロ扱いで座標固定オプション

---

## 公開 API (mod.rs)

```rust
pub struct Conformer {
    pub atom_symbols: Vec<String>,
    pub positions: Vec<[f64; 3]>,   // Å単位
}

/// SMILES から3D座標を生成する
/// seed: 乱数シード（再現性のため）
pub fn smiles_to_3d(smiles: &str, seed: u64) -> Option<Conformer>

/// DuckDB FFI
/// duck_smiles_to_3d(smiles, seed) → JSON文字列
/// "[{"sym":"C","x":1.2,"y":0.3,"z":-0.5}, ...]"
#[no_mangle]
pub extern "C" fn duck_smiles_to_3d(...)
```

DuckDB SQL での使用イメージ:
```sql
SELECT smiles, smiles_to_3d(smiles) AS coords
FROM compound_library;

-- 座標を展開
SELECT smiles, unnest(coords.x) AS x, unnest(coords.y) AS y
FROM (SELECT smiles, smiles_to_3d(smiles) AS coords FROM library);
```

---

## 実装順序とマイルストーン

### Step 1: params.rs (1〜2日)
- 結合長テーブル実装・単体テスト
- 結合角テーブル実装・単体テスト
- VDW 半径テーブル実装

完了条件: `bond_length("C", "C", BondOrder::Single) == 1.54`

### Step 2: ハイブリダイゼーション判定 (2〜3日)
- `Molecule` に `hybridization(idx: usize) -> Hybridization` を追加
- sp/sp2/sp3 の判定ロジック
- テスト: ベンゼン全炭素→sp2、エタン→sp3、アセチレン→sp

### Step 3: bounds.rs (3〜5日)
- BoundsMatrix 構造体
- 1-2 拘束の実装・テスト
- 1-3 拘束の実装・テスト
- 1-4 拘束の実装（単純化: trans のみ）
- VDW 下限の設定

完了条件: エタン (CC) の bounds が物理的に妥当

### Step 4: smoothing.rs (1日)
- Floyd-Warshall 実装
- lower ≤ upper の整合性チェック

### Step 5: embed.rs (4〜7日)
- LCG 乱数実装（外部不使用）
- 距離行列サンプリング
- 計量行列変換
- Jacobi 固有値分解（N×N 対称行列）← ここが最難関
- 座標生成

完了条件: 水分子 (O) 1原子で [0,0,0]、エタン (CC) で C-C ≈ 1.54Å

### Step 6: minimize.rs (3〜5日)
- 数値勾配計算
- エネルギー関数実装 (bond + angle + vdw)
- 最急降下法ループ
- 収束判定

完了条件: ベンゼンが平面六角形に収束する

### Step 7: 統合・DuckDB バインディング (2〜3日)
- mod.rs で全ステップを接続
- FFI 関数の実装
- JSON 出力フォーマット
- lib.rs への登録

### Step 8: テスト (2〜3日)
- 単原子分子 (C, N, O)
- 直線分子 (acetylene HCC#CCH → ほぼ直線)
- 平面分子 (benzene → all z ≈ 0)
- 3D 分子 (cyclohexane → chair 構造)
- 結合長・結合角の精度確認

---

## テスト戦略

### ユニットテスト (各モジュール内)

```rust
#[test]
fn test_bond_length_cc_single() {
    assert!((bond_length("C","C", BondOrder::Single) - 1.54).abs() < 0.01);
}

#[test]
fn test_bounds_ethane() {
    let mol = parse("CC").unwrap();
    let b = build_bounds(&mol);
    let d_cc = b.lower[0 * b.n + 1];
    assert!((d_cc - 1.54).abs() < 0.01);
}

#[test]
fn test_embed_bond_length_preserved() {
    let mol = parse("CC").unwrap();
    let conf = smiles_to_3d("CC", 42).unwrap();
    let d = dist(&conf.positions[0], &conf.positions[1]);
    assert!((d - 1.54).abs() < 0.05);
}

#[test]
fn test_benzene_planar() {
    let conf = smiles_to_3d("c1ccccc1", 42).unwrap();
    // 全炭素の z 座標が同一平面内 (±0.1Å)
    let zs: Vec<f64> = conf.positions.iter().map(|p| p[2]).collect();
    let z_range = zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                 - zs.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(z_range < 0.1);
}
```

### 精度ベンチマーク

RDKit で生成した座標との RMSD 比較（目標: RMSD < 1.0Å）。
ベンチマーク分子:
- aspirin, paracetamol, caffeine (既存テスト化合物)
- cyclohexane (立体化学)
- biphenyl (回転可能結合)

---

## 既知の制約・スコープ外

| 項目 | 対応 |
|---|---|
| 立体化学 (@, @@) | Phase 1 では無視。ランダム埋め込み後に反転チェックを追加予定 |
| 複数コンフォマー生成 | Phase 1 では最良の1コンフォマーのみ。ドッキングで必要なら Phase 3 で追加 |
| 巨大分子 (N > 200) | 固有値分解が O(N²)〜O(N³) になるため対象外 |
| 金属錯体 | VDW パラメータが不明。スキップ |
| 電荷の影響 | 静電項なし（UFF-lite の割り切り） |

---

## 依存追加なし方針

外部 crate を追加しない。必要な数学ルーティン:
- 行列演算: 手実装（N≤200 の薬様分子なら十分速い）
- 固有値分解: Jacobi 法（収束が遅くても許容）
- 乱数: LCG (Linear Congruential Generator)
- 三角関数: `f64::sin/cos/sqrt` のみ

これにより duckSMILES の「外部依存ゼロ」方針を維持する。
