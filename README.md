# PyPSA × LHS × PCE サロゲート グローバル感度分析 (GSA)

## 概要
単一ノードのPyPSA容量拡張モデル（太陽光・風力・定置型蓄電池・ディーゼル、**8760時間 = 1年間フル時系列**）の
年間総費用（年換算CAPEX＋OPEX）を目的関数とし、4つのコスト不確実性パラメータの影響を次のパイプラインで評価します。

```
network_config.xlsx → LHS(320点) → PyPSA(LP, 8760h) → PCEサロゲート学習
                     → サロゲート上GSA(10,240点) → Sobol S1/ST/Sij・Cij
```

- ネットワーク構成（ノード・発電機・蓄電池・需要の基本定義）と8760時間の時系列データは、すべて
  `network_config.xlsx` 1ファイルから動的に読み込みます（コード埋め込みではありません）。
- **ファイルが存在しない場合はエラー終了せず**、リアリスティックな合成データを含むサンプル
  `network_config.xlsx` を自動生成して実行を継続します。
- **CAPEXの年換算 (CRF)**: 各技術シートの `capex` [$/MW] と `lifetime` [年] を用いて資本回収係数で年換算し
  `capital_cost` [$/MW/year] とします。

  `CRF = r(1+r)^n / ((1+r)^n − 1)`,  割引率 r = 5%（`DISCOUNT_RATE` で固定）

| 技術 | 耐用年数 n（既定値） | CRF |
|---|---|---|
| 太陽光 solar | 25年 | 0.07095 |
| 風力 wind | 25年 | 0.07095 |
| 蓄電池 battery | 12年 | 0.11283 |
| ディーゼル diesel | 20年 | 0.08024 |

  ディーゼルの初期建設費は不確実性の対象外で、Excelの `generators` シートの値（既定 800,000 $/MW）に固定です。

- 不確実性パラメータ（LHS/Sobolで振る4次元, `PROBLEM["bounds"]` で編集可能）:

| 名称 | 範囲 | 単位 |
|---|---|---|
| `solar_cost` | 500,000 – 900,000 | $/MW（初期建設費, CRFで年換算） |
| `wind_cost` | 800,000 – 1,450,000 | $/MW（同上） |
| `battery_cost` | 300,000 – 900,000 | $/MW（同上） |
| `diesel_marginal_cost` | 200 – 300 | $/MWh（燃料可変費） |

  ※ Excel側の `capex`/`marginal_cost` は「サンプル生成前のベース値」であり、実際の感度分析では上記範囲で
  独立にLHS/Sobolサンプリングした値で各技術の `capital_cost`（またはディーゼルの `marginal_cost`）を上書きします。

- サロゲート: `StandardScaler → PolynomialFeatures(2) → RidgeCV(alphas=logspace(-4,3,30), cv=10)`（基底15個）。
- 必須ライブラリが不足している場合は、不足名を表示して `sys.exit(1)` します（ダミーデータへのフォールバックなし）。
  PyPSAの最適化が失敗した場合も例外で停止します。

## ディレクトリ構成
```
GSA-02/
├── pypsa_pce_gsa.py             # 本体スクリプト
├── network_config.xlsx          # ネットワーク設定 + 8760h時系列（無ければ自動生成）
├── requirements.txt             # pip用
├── environment.yml              # conda用
├── README.md
└── results/
    └── YYYY-MM-DD_HH-MM/        # 実行ごとにタイムスタンプ付きフォルダを自動生成
        ├── pypsa_lhs_320_results.csv / .nc
        ├── sobol_s2_matrix.csv
        ├── pce_interaction_matrix.csv
        └── pypsa_pce_gsa_results.png
```

## `network_config.xlsx` のシート構成と編集方法
既存の値を変えたり技術構成を調整したりする場合は、このExcelファイルを直接編集してください。

| シート | 必須列 | 説明 |
|---|---|---|
| `buses` | `bus_name`, `v_nom` | ノード定義（既定は単一ノード `bus`） |
| `generators` | `name`, `bus`, `carrier`, `capex`, `marginal_cost`, `efficiency`, `lifetime`, `p_nom_extendable` | 発電機の基本定義。`carrier` が `solar`/`wind` の行は `timeseries` の対応する `*_p_max_pu` 列が自動で出力上限として適用されます |
| `storage_units` | `name`, `bus`, `carrier`, `capex`, `max_hours`, `efficiency_store`, `efficiency_dispatch`, `lifetime`, `p_nom_extendable` | 蓄電池の基本定義 |
| `loads` | `name`, `bus` | 需要の基本定義（実データは `timeseries.load_mw`） |
| `timeseries` | `timestamp`, `solar_p_max_pu`, `wind_p_max_pu`, `load_mw` | **8,760行**（1年・1時間刻み）。`solar_p_max_pu`/`wind_p_max_pu` は0.0〜1.0、`load_mw` は需要[MW] |

- 感度分析の対象4技術（`solar`, `wind`, `battery`, `diesel`）は、`generators`/`storage_units` シートの `name` 列で
  この名前と一致する行が使われます（各シートの `lifetime` がCRF年換算に使われます）。
- `timeseries` の行数が8760でない場合、または `*_p_max_pu` が0〜1の範囲外の場合はエラーで停止します。
- ファイルが存在しない場合、上記スキーマに従うサンプルデータ（決定論的な合成気象・負荷プロファイル、固定シード）
  を自動生成します。

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
オプション（動作確認・縮小実行用。既定値はいずれも要求仕様どおり）:
```bash
python pypsa_pce_gsa.py --config network_config.xlsx --n-lhs 320 --n-sobol 1024
```

**所要時間の目安**: 8760時間フル解像度のLPを320回解くため、PC性能に応じて数十分〜1時間程度かかることがあります
（サロゲート学習とSobol評価自体は数秒〜数十秒）。乱数シードは固定（42）です。動作確認だけ行いたい場合は
`--n-lhs 16 --n-sobol 8` のように小さい値を指定してください（`--n-lhs` は10-fold CVのため10以上を推奨）。

## 出力ファイル
| ファイル | 内容 |
|---|---|
| `pypsa_lhs_320_results.csv` / `.nc` | 各サンプルの入力（CAPEX [$/MW], ディーゼル可変費 [$/MWh]）、最適総コスト `total_cost`[$/year]、最適容量 `solar_mw, wind_mw, battery_mw, diesel_mw`[MW]、年間発電/充放電量 `solar_mwh, wind_mwh, battery_discharge_mwh, diesel_mwh`[MWh/year] |
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

## 注意事項
- `network_config.xlsx` を自動生成した場合の気象・負荷プロファイルはスクリプト内で決定論的に生成する合成データで、
  実測データではありません。実データに置き換える場合は `timeseries` シートを8760行のまま差し替えてください。
- `--n-lhs` を10未満にすると10-fold CVが機能しないため、動作確認目的でも10以上を指定してください。
