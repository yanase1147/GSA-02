# PyPSA × LHS × PCE サロゲート グローバル感度分析 (GSA)

## 概要
単一ノードのPyPSA容量拡張モデル（太陽光・風力・定置型蓄電池・ディーゼル）の年間総費用（年換算CAPEX＋OPEX）を目的関数とし、
4つのコスト不確実性パラメータの影響を次のパイプラインで評価します。

```
LHS(320点) → PyPSA(LP) → PCEサロゲート学習 → サロゲート上GSA(10,240点) → Sobol S1/ST/Sij・Cij
```

- モデル: 3時間解像度 × 5ステップおき = **584スナップショット**（各15h重み、合計8760h）。負荷・太陽光/風力CFはスクリプト内で決定論的に生成する合成プロファイルです（実データではありません）。
- **CAPEXの年換算 (CRF)**: 初期建設費 [$/MW] を資本回収係数で年換算し `capital_cost` [$/MW/year] とします。

  `CRF = r(1+r)^n / ((1+r)^n − 1)`,  割引率 r = 5%

| 技術 | 耐用年数 n | CRF |
|---|---|---|
| 太陽光 solar | 25年 | 0.07095 |
| 風力 wind | 25年 | 0.07095 |
| 蓄電池 battery | 12年 | 0.11283 |
| ディーゼル diesel | 20年 | 0.08024 |

  ディーゼルの初期建設費は不確実性の対象外で 800,000 $/MW に固定しています。

- 不確実性パラメータ（変動範囲）:

| 名称 | 範囲 | 単位 |
|---|---|---|
| `solar_cost` | 500,000 – 900,000 | $/MW（初期建設費, CRFで年換算） |
| `wind_cost` | 800,000 – 1,450,000 | $/MW（同上） |
| `battery_cost` | 300,000 – 900,000 | $/MW（4時間蓄電池, MW基準, 同上） |
| `diesel_marginal_cost` | 200 – 300 | $/MWh |

  CAPEX範囲は、従来の年換算範囲（例: 太陽光 35,000–65,000 $/MW/year）を25年・5%のCRFで初期建設費に戻した値を丸めたものです。範囲を変えたい場合は `PROBLEM["bounds"]` を編集してください。

- サロゲート: `StandardScaler → PolynomialFeatures(2) → RidgeCV(alphas=logspace(-4,3,30), cv=10)`（基底15個）。
- 必須ライブラリが不足している場合は、不足名を表示して `sys.exit(1)` します（ダミーデータへのフォールバックなし）。PyPSAの最適化が失敗した場合も例外で停止します。

## ディレクトリ構成
```
GSA-02/
├── pypsa_pce_gsa.py            # 本体スクリプト
├── requirements.txt            # pip用
├── environment.yml             # conda用
├── README.md
└── (実行後に生成)
    ├── pypsa_lhs_320_results.csv / .nc
    ├── sobol_s2_matrix.csv
    ├── pce_interaction_matrix.csv
    └── pypsa_pce_gsa_results.png
```

## 環境構築
```bash
# pip
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# conda
conda env create -f environment.yml
conda activate pypsa-pce-gsa
```
求解器は HiGHS（`highspy`）を使用します。

## 実行方法
```bash
python pypsa_pce_gsa.py
```
所要時間は約5〜6分（LP 320回が大半。サロゲートによる10,240点の評価は数ミリ秒）。乱数シードは固定（42）です。

## 出力ファイル
| ファイル | 内容 |
|---|---|
| `pypsa_lhs_320_results.csv` / `.nc` | 320点の入力（CAPEX [$/MW], ディーゼル可変費 [$/MWh]）、最適総コスト `total_cost`[$/year]、最適容量 `solar_mw, wind_mw, diesel_mw, battery_mw` |
| `sobol_s2_matrix.csv` | Sobol 2次感度指標 Sij の4×4対称行列（対角=0） |
| `pce_interaction_matrix.csv` | PCE交差項係数 Cij の4×4対称行列（対角=2乗項係数 Cii） |
| `pypsa_pce_gsa_results.png` | 左: S1/ST棒グラフ、右: PyPSA vs PCE 1:1プロット（R², CV R²） |

## 結果の解釈ルール
- **S1**: 単独でのコスト分散への寄与割合。**ST**: 他パラメータとの交互作用を含む総寄与。
  `ST − S1` が大きいほど交互作用が強い。ΣS1≈1なら加法的なモデルです。
- **Sij (S2)**: 2パラメータ間の交互作用が出力分散に占める割合。信頼区間を考慮し、極端に小さい値はノイズの可能性があります。
- **Cij**: 標準化入力空間（平均0・標準偏差1）での2次交差項の回帰係数 [$/year]。符号は元のパラメータ空間の
  交差偏微分 ∂²C/∂θi∂θj と一致し、大きさは同一スケールで比較できます。

### Cij の符号解釈（総費用最小化 × 包絡線定理）
目的関数は総費用 C(θ) の最小化なので、包絡線定理により ∂C/∂θi = 最適量_i（CAPEXなら最適容量）。したがって

`Cij = ∂²C/∂θi∂θj = ∂(最適量_i)/∂θj`

- **Cij < 0 【補完関係 (容量連動・シナジー)】**: θj が上がると技術 i の最適容量も減る＝両技術の容量が連動して減少する関係。
  補完技術が同時に縮小するため、総費用の増加が抑えられます。
- **Cij > 0 【代替関係 (技術競合・置換)】**: θj が上がると技術 i の最適容量が増える＝コストが上がった技術から
  もう一方へ容量がシフト（置換）される競合関係。
- 注意: `diesel_marginal_cost` は容量ではなく発電量（エネルギー）に対するコストのため、ディーゼルを含むペアでは
  「最適量」はディーゼルの発電量を指します。
- Cij は符号付きの係数、Sij は分散寄与割合であり別の量です（Cijが大きくてもSijは小さいことがあります）。
- サロゲートの信頼性は 学習R²・CV R²・RMSE で確認してください（CV R²が低い場合はGSA結果を信頼しないでください）。

## 参考結果（seed=42 での実行例）
- 学習 R² = 0.9993、CV(10-fold) R² = 0.9992、CV RMSE ≈ 0.134 M$/year。
- S1 ≈ 風力CAPEX 0.73 / 太陽光 0.11 / ディーゼル可変費 0.09 / 蓄電池 0.05。ST は S1 とほぼ同じで、交互作用は小さい（Sij ≤ 0.004）。
- Cij: 太陽光×蓄電池のみ負（補完）、他の5ペアは正（代替）。
- 上記は本スクリプト内の合成プロファイルと固定パラメータに依存する結果で、一般的な結論ではありません。
- 代表日サンプリングは15h間隔の間引きのため、蓄電池のSOC連続性は近似されています。
