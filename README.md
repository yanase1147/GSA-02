# PyPSA × LHS × PCE サロゲート グローバル感度分析 (GSA)

## 概要
PyPSA容量拡張モデル（**8760時間 = 1年間フル時系列**）の年間総費用（年換算CAPEX＋OPEX）を目的関数とし、
Excelで自由に定義した任意個数・任意技術のコスト不確実性パラメータの影響を次のパイプラインで評価します。

```
network_config.xlsx → LHS(320点) → PyPSA(LP, 8760h) → PCEサロゲート学習
                     → サロゲート上GSA(10,240点) → Sobol S1/ST/Sij・Cij
```

- ネットワーク構成（ノード・発電機・蓄電池・リンク・需要の基本定義）、8760時間の時系列データ、および
  **不確実性パラメータの定義そのもの**が、すべて `network_config.xlsx` 1ファイルから動的に読み込まれます
  （コード埋め込みではありません）。
- **任意電源の追加はExcelのみで完結**: `generators`/`storage_units`（任意で`links`）シートに行を追加するだけで、
  コード変更なしに新しい技術（水素、水力、原子力等）がPyPSAネットワークへ組み込まれ最適化されます。
- **不確実性パラメータも任意個数**: `uncertainty_params` シートの行数がそのままLHS/Sobolの次元数、
  PCE多項式基底数（`C(n+2,2)`）、Sobol S1/ST/Sijマトリクスのサイズに自動反映されます（コード側の次元数は
  一切ハードコードされていません）。
- **ファイルが存在しない場合はエラー終了せず**、リアリスティックな合成データを含むサンプル
  `network_config.xlsx`（太陽光・風力・蓄電池・ディーゼルの4技術、4パラメータ）を自動生成して実行を継続します。
- **CAPEXの年換算 (CRF)**: 対象電源シートの `capex` [$/MW] と `lifetime` [年] を用いて資本回収係数で年換算し
  `capital_cost` [$/MW/year] とします（`target_attribute=capital_cost` のパラメータにのみ適用）。

  `CRF = r(1+r)^n / ((1+r)^n − 1)`,  割引率 r = 5%（`DISCOUNT_RATE` で固定、全電源共通）

  `target_attribute=marginal_cost`（可変費/OPEX）のパラメータはCRFを適用せず、サンプル値をそのまま使用します。

- 既定サンプル設定での不確実性パラメータ（4次元）:

| `param_name` | `component_type` | `component_name` | `target_attribute` | 範囲 | 単位 |
|---|---|---|---|---|---|
| `solar_capex` | generator | solar | capital_cost | 500,000 – 900,000 | $/MW（CRFで年換算） |
| `wind_capex` | generator | wind | capital_cost | 800,000 – 1,450,000 | $/MW（同上） |
| `battery_capex` | storage_unit | battery | capital_cost | 300,000 – 900,000 | $/MW（同上） |
| `diesel_opex` | generator | diesel | marginal_cost | 200 – 300 | $/MWh（燃料可変費） |

  ※ Excel側の各電源シートの `capex`/`marginal_cost` は「サンプリング前のベース値」で、`uncertainty_params` に
  登録された行についてはLHS/Sobolでサンプリングした値が対応する `capital_cost`/`marginal_cost` を上書きします。
  未登録の属性はExcelのベース値のまま固定されます。

- サロゲート: `StandardScaler → PolynomialFeatures(2) → RidgeCV(alphas=logspace(-4,3,30), cv=10)`
  （基底数は次元数nから `C(n+2,2)` で自動決定。既定4次元では15個）。
- 必須ライブラリが不足している場合は、不足名を表示して `sys.exit(1)` します（ダミーデータへのフォールバックなし）。
  PyPSAの最適化が失敗した場合も例外で停止します。

## ディレクトリ構成
```
GSA-02/
├── pypsa_pce_gsa.py             # 本体スクリプト (LHS → PyPSA → PCE → Sobol GSA)
├── plot_annual_dispatch.py      # 代表値での年間/月別ディスパッチ解析・作図スクリプト
├── network_config.xlsx          # ネットワーク設定 + 8760h時系列（無ければ自動生成）
├── requirements.txt             # pip用
├── environment.yml              # conda用
├── README.md
├── results/
│   └── YYYY-MM-DD_HH-MM/        # GSA実行ごとにタイムスタンプ付きフォルダを自動生成
│       ├── pypsa_lhs_320_results.csv / .nc
│       ├── sobol_s2_matrix.csv
│       ├── pce_interaction_matrix.csv
│       └── pypsa_pce_gsa_results.png
└── result_d/
    └── YYYY-MM-DD_HH-MM/        # plot_annual_dispatch.py 実行ごとに自動生成
        ├── annual_summary_metrics.csv
        ├── monthly_summary_metrics.csv
        ├── annual_dispatch_data.csv
        ├── annual_dispatch_balance.png
        ├── execution.log
        └── monthly_plots/01_Jan.png 〜 12_Dec.png
```

## `network_config.xlsx` のシート構成と編集方法
既存の値を変えたり技術構成・不確実性パラメータを調整したりする場合は、このExcelファイルを直接編集してください。
**コードの変更は不要です。**

| シート | 必須列 | 説明 |
|---|---|---|
| `buses` | `bus_name`, `v_nom` | ノード定義（既定は単一ノード `bus`） |
| `generators` | `name`, `bus`, `carrier`, `capex`, `marginal_cost`, `lifetime`, `p_nom_extendable` | 発電機の基本定義。行を追加するだけで新しい発電技術（水力・原子力・水素発電等）を組み込めます |
| `storage_units` | `name`, `bus`, `carrier`, `capex`, `max_hours`, `efficiency_store`, `efficiency_dispatch`, `lifetime`, `p_nom_extendable` | 蓄電池等の基本定義 |
| `loads` | `name`, `bus` | 需要の基本定義（実データは `timeseries` の対応列） |
| `timeseries` | `timestamp`, 各種`*_p_max_pu`, 各種`*_mw` | **8,760行**（1年・1時間刻み） |
| `uncertainty_params` | `param_name`, `component_type`, `component_name`, `target_attribute`, `lower_bound`, `upper_bound` | **不確実性パラメータの定義**。行数がそのままLHS/Sobol/PCEの次元数になります |
| `links`（任意） | `name`, `bus0`, `bus1`, `carrier`, `capex`, `marginal_cost`, `efficiency`, `lifetime`, `p_nom_extendable` | 電解槽・燃料電池等をリンクとして表現する場合にのみ追加（無ければ無視されます） |

### 命名規則（時系列列とコンポーネントの対応付け）
- 可変出力の発電機（太陽光・風力・水力等）は、`timeseries` に `<generators.name>_p_max_pu` 列（0.0〜1.0）を
  用意すると自動的に出力上限として適用されます（例: 名前 `hydro` なら列名 `hydro_p_max_pu`）。列が無い発電機は
  ディーゼル・原子力のように `p_nom` まで自由に出力できる電源として扱われます。
- 需要は `timeseries` に `<loads.name>_mw` 列（既定は `load_mw`）が必要です。

### `uncertainty_params` シート（不確実性パラメータの汎用登録）
| 列名 | 説明 |
|---|---|
| `param_name` | パラメータ識別名（一意。LHS/Sobol/PCEの変数名・結果列名として使用） |
| `component_type` | `generator` / `storage_unit` / `link` のいずれか |
| `component_name` | 対象電源の名前（対応シートの `name` 列と一致している必要あり） |
| `target_attribute` | `capital_cost`（固定費/CAPEX, CRFで自動年換算）または `marginal_cost`（可変費/OPEX, そのまま適用） |
| `lower_bound` / `upper_bound` | LHS/Sobolでサンプリングする範囲 |

- 行を追加・削除するだけで、LHS/PCE/Sobol/Cijのすべてが自動的に次元数を追従します（コード側は次元数を
  ハードコードしていません）。
- 同一電源に対して固定費と可変費の両方を不確実性パラメータとして同時登録することも可能です
  （例: `diesel_capex` と `diesel_opex` を両方登録）。
- `uncertainty_params` に登録されていない属性は、対応シートのベース値（`capex`/`marginal_cost`）に固定されます。
- `component_name` が対応シートに存在しない、`target_attribute` が不正、`lower_bound >= upper_bound` などの
  場合は起動時にエラーで停止します（サイレントに無視されることはありません）。
- `timeseries` の行数が8760でない場合、または `*_p_max_pu` 列が0〜1の範囲外の場合もエラーで停止します。
- ファイルが存在しない場合、上記スキーマに従うサンプルデータ（太陽光・風力・蓄電池・ディーゼルの4技術、
  決定論的な合成気象・負荷プロファイル、固定シード）を自動生成します。

### 拡張例: 新技術（原子力）の追加
1. `generators` シートに1行追加: `name=nuclear, bus=bus, carrier=nuclear, capex=6000000, marginal_cost=15, lifetime=40, p_nom_extendable=TRUE`
2. `uncertainty_params` シートに1行追加: `param_name=nuclear_capex, component_type=generator, component_name=nuclear, target_attribute=capital_cost, lower_bound=4000000, upper_bound=8000000`
3. `python pypsa_pce_gsa.py` を実行するだけで、5次元LHS/PCE/Sobol（基底数21）に自動的に拡張されます。
   出力CSVにも `nuclear_mw` / `nuclear_mwh` 列が自動的に追加されます。

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
| `pypsa_lhs_320_results.csv` / `.nc` | 各サンプルの入力（`uncertainty_params.param_name`列群）、最適総コスト `total_cost`[$/year]、**ネットワーク内の全発電機/蓄電池/リンクについて自動生成される** `<name>_mw`（最適容量）と `<name>_mwh`（年間発電量、蓄電池は`<name>_discharge_mwh`）列 |
| `sobol_s2_matrix.csv` | Sobol 2次感度指標 Sij のn×n対称行列（対角=0, n=不確実性パラメータ数） |
| `pce_interaction_matrix.csv` | PCE交差項係数 Cij のn×n対称行列（対角=2乗項係数 Cii） |
| `pypsa_pce_gsa_results.png` | 左: S1/ST棒グラフ、右: PyPSA vs PCE 1:1プロット（R², CV R²） |

## 代表値での年間ディスパッチ解析 (`plot_annual_dispatch.py`)
GSAとは別に、不確実性パラメータの**代表値（`lower_bound` と `upper_bound` の中央値）**で8,760時間の最適化を
1回だけ実行し、年間・月別の需給バランスを可視化します。ネットワーク構築・CRF年換算は `pypsa_pce_gsa.py` の
関数をそのまま再利用しています（電源名は `solar` / `wind` / `battery` / `diesel` を想定）。

```bash
python plot_annual_dispatch.py
python plot_annual_dispatch.py --config network_config.xlsx
```
所要時間の目安は数十秒です（LP 1回のみ）。

### 出力ファイル (`result_d/YYYY-MM-DD_HH-MM/`)
| ファイル | 内容 |
|---|---|
| `annual_summary_metrics.csv` | 電源別の適用コスト（代表値, CAPEXは`$/MW/year`・OPEXは`$/MWh`）、最適容量[MW]、容量比率[%]、年間発電量[MWh]、年間発電量比率[%]（蓄電池は年間放電量） |
| `monthly_summary_metrics.csv` | 1〜12月ごとの Solar/Wind/Diesel 発電量、蓄電池充放電量、需要[MWh]、再エネ比率[%] |
| `annual_dispatch_data.csv` | 8,760時間全時系列（各電源出力・蓄電池充放電・SOC・需要） |
| `annual_dispatch_balance.png` | 年間需給バランス（積層エリア＋需要線）と蓄電池SOC (16×9, 300dpi)。最適容量・構成比のインセット付き |
| `monthly_plots/01_Jan.png`〜`12_Dec.png` | 月別の需給バランスとSOC (14×8, 300dpi, 日単位の目盛り) |
| `execution.log` | コンソール出力の記録（PyPSAのログ・進捗バーを含む） |

- 配色は全図共通: Solar `#F1C40F` / Wind `#2980B9` / Battery Discharge `#2ECC71` / Battery Charge `#27AE60` /
  Diesel `#7F8C8D` / Load `#E74C3C`。グラフ内表記は英語で統一しています。
- 「年間発電量比率」「再エネ比率」は、年間（月別）総需要に対する割合です。蓄電池の充放電損失や出力の
  過剰分により、合計が100%を超えることがあります。
- `execution.log` には、HiGHSのネイティブなソルバーバナー（反復ログ）は保存されない場合があります。
  端末のOSレベル出力を付け替える方式は、Windowsの対話型コンソールでフリーズする恐れがあったため採用していません
  （ソルバーログは通常どおりコンソールに表示されます）。

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
- 注意: `target_attribute=marginal_cost` のパラメータ（可変費/OPEX）は容量ではなく発電量（エネルギー）に
  対するコストのため、そのパラメータを含むペアでは「最適量」は該当電源の発電量を指します。
- Cij は符号付きの係数、Sij は分散寄与割合であり別の量です（Cijが大きくてもSijは小さいことがあります）。
- サロゲートの信頼性は 学習R²・CV R²・RMSE で確認してください（CV R²が低い場合はGSA結果を信頼しないでください）。

## 注意事項
- `network_config.xlsx` を自動生成した場合の気象・負荷プロファイルはスクリプト内で決定論的に生成する合成データで、
  実測データではありません。実データに置き換える場合は `timeseries` シートを8760行のまま差し替えてください。
- `--n-lhs` を10未満にすると10-fold CVが機能しないため、動作確認目的でも10以上を指定してください。
