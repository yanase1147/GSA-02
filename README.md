# PyPSA × LHS × PCE サロゲート グローバル感度分析 (GSA)

## 概要
単一ノードのPyPSA容量拡張モデル（太陽光・風力・定置型蓄電池・ディーゼル）の年間総費用（CAPEX＋OPEX）を目的関数とし、
4つのコスト不確実性パラメータの影響を次のパイプラインで評価します。

```
LHS(320点) → PyPSA(LP) → PCEサロゲート学習 → サロゲート上GSA(10,240点) → Sobol S1/ST/Sij・Cij
```

- モデル: 3時間解像度 × 5ステップおき = **584スナップショット**（各15h重み、合計8760h）。負荷・太陽光/風力CFはスクリプト内で決定論的に生成する合成プロファイルです（実データではありません）。
- 不確実性パラメータ:

| 名称 | 範囲 | 単位 |
|---|---|---|
| `solar_cost` | 35,000 – 65,000 | $/MW/year |
| `wind_cost` | 56,000 – 104,000 | $/MW/year |
| `battery_cost` | 30,000 – 90,000 | $/MW/year（4時間蓄電池, MW基準） |
| `diesel_marginal_cost` | 60 – 140 | $/MWh |

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
所要時間は約6分（LP 320回が大半。サロゲートによる10,240点の評価は数ミリ秒）。乱数シードは固定（42）です。

## 出力ファイル
| ファイル | 内容 |
|---|---|
| `pypsa_lhs_320_results.csv` / `.nc` | 320点の入力コスト、最適総コスト `total_cost`[$/year]、最適容量 `solar_mw, wind_mw, diesel_mw, battery_mw` |
| `sobol_s2_matrix.csv` | Sobol 2次感度指標 Sij の4×4対称行列（対角=0） |
| `pce_interaction_matrix.csv` | PCE交差項係数 Cij の4×4対称行列（対角=2乗項係数 Cii） |
| `pypsa_pce_gsa_results.png` | 左: S1/ST棒グラフ、右: PyPSA vs PCE 1:1プロット（R², CV R²） |

## 結果の解釈ルール
- **S1**: 単独でのコスト分散への寄与割合。**ST**: 他パラメータとの交互作用を含む総寄与。
  `ST − S1` が大きいほど交互作用が強い。ΣS1≈1なら加法的なモデルです。
- **Sij (S2)**: 2パラメータ間の交互作用が出力分散に占める割合（0以上が目安）。値が小さくても、
  ばらつきが大きい場合は信頼区間内のノイズの可能性があります。
- **Cij**: 標準化入力空間（平均0・標準偏差1）での2次交差項の回帰係数 [$/year]。符号は方向を示し、大きさは同一スケールで比較できます。
  - **Cij > 0 【補完関係 (相乗効果 +)】**: 両パラメータが同時に高い（低い）とき、費用への影響が加算より大きい方向（一方の変化の効果がもう一方の高い水準で増幅）。
  - **Cij < 0 【代替関係 (競合相殺 -)】**: 一方が上がると他方の効果が相殺される方向（技術が互いの代替となる）。
  - Cij は符号付きの係数、Sij は分散寄与割合であり、両者は別の量です（Cijが大きくてもSijは小さいことがあります）。
- サロゲートの信頼性は 学習R²・CV R²・RMSE で確認してください（CV R²が低い場合はGSA結果を信頼しないでください）。

## 参考結果（seed=42 での実行例）
- 学習 R² = 0.9990、CV(10-fold) R² = 0.9988、CV RMSE ≈ 0.147 M$/year。
- S1 ≈ 風力CAPEX 0.61 / ディーゼル燃料費 0.30 / 太陽光 0.08 / 蓄電池 0.007。ST は S1 とほぼ同じで、交互作用は小さい（Sij ≤ 0.002）。
- 蓄電池は最適解で容量ゼロとなるケースが多く、感度が低い。これは本スクリプト内の合成プロファイルと固定パラメータ（ディーゼル設備費 90,000 $/MW/year等）に依存する結果で、一般的な結論ではありません。
- 代表日サンプリングは15h間隔の間引きのため、蓄電池のSOC連続性は近似されています。
