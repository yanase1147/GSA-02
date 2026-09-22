#!/usr/bin/env python
"""
Excel(network_config.xlsx) -> LHS(320点) -> PyPSA(8760h) -> PCEサロゲート
                            -> サロゲート上GSA(10,240点) -> Sobol / Cij 解析

不確実性パラメータ(固定費/可変費、任意の電源)は network_config.xlsx の
uncertainty_params シートで定義する。電源(generators/storage_units/links)も
Excelに行を追加するだけでコード変更なしに反映される。

実行:  python pypsa_pce_gsa.py
       python pypsa_pce_gsa.py --n-lhs 8 --n-sobol 16   (動作確認用の縮小実行)
"""
import importlib
import sys

# Windows既定コードページ(cp932等)では日本語出力が文字化けするため、UTF-8に固定する。
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

# ----------------------------------------------------------------------------
# 0. 必須ライブラリのチェック (不足時はフォールバックせず即停止)
# ----------------------------------------------------------------------------
REQUIRED = {  # import名: pip/conda名
    "numpy": "numpy",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "xarray": "xarray",
    "netCDF4": "netcdf4",
    "sklearn": "scikit-learn",
    "SALib": "salib",
    "pypsa": "pypsa",
    "highspy": "highspy",
    "openpyxl": "openpyxl",
}
missing = []
for mod, pkg in REQUIRED.items():
    try:
        importlib.import_module(mod)
    except ImportError:
        missing.append(pkg)
if missing:
    print("[ERROR] 必須ライブラリが不足しています: " + ", ".join(missing), file=sys.stderr)
    print("        pip install -r requirements.txt  または  conda env create -f environment.yml",
          file=sys.stderr)
    sys.exit(1)

import argparse
import logging
import os
import time
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from SALib.analyze import sobol as sobol_analyze
from SALib.sample import latin as latin_sample
from SALib.sample import sobol as sobol_sample
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

logging.getLogger("pypsa").setLevel(logging.WARNING)
logging.getLogger("linopy").setLevel(logging.ERROR)
logging.getLogger("highspy").setLevel(logging.ERROR)

BASE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = BASE_DIR / "network_config.xlsx"
SEED = 42
N_LHS = 320
N_SOBOL = 1024
N_HOURS = 8760
DISCOUNT_RATE = 0.05  # 割引率 (CAPEXの年換算 CRF 計算用, 全電源共通)

# component_type(uncertainty_params) -> PyPSA Network属性名
COMPONENT_ATTR = {"generator": "generators", "storage_unit": "storage_units", "link": "links"}
VALID_TARGET_ATTRS = {"capital_cost", "marginal_cost"}


def crf(rate, lifetime):
    """資本回収係数 CRF = r(1+r)^n / ((1+r)^n - 1)"""
    if rate == 0:
        return 1.0 / lifetime
    f = (1 + rate) ** lifetime
    return rate * f / (f - 1)


def annualized(capex, lifetime, rate=DISCOUNT_RATE):
    """初期建設費 [$/MW] -> 年換算固定費 capital_cost [$/MW/year]"""
    return capex * crf(rate, lifetime)


# ----------------------------------------------------------------------------
# 1. network_config.xlsx の読み込み / 自動生成
# ----------------------------------------------------------------------------
def _build_synthetic_timeseries():
    """決定論的な8760時間(1年・1時間刻み)の合成プロファイルを生成する。"""
    idx = pd.date_range("2025-01-01", periods=N_HOURS, freq="h")  # 非うるう年(365日x24h=8760h)
    h = np.asarray(idx.hour)
    doy = np.asarray(idx.dayofyear)
    rng = np.random.default_rng(2024)  # 固定シードの合成天候変動

    # 負荷 [MW]: 平均~100MW, 日変動 + 季節変動 + ノイズ
    load = (100 + 15 * np.sin(2 * np.pi * (h - 9) / 24)
            + 10 * np.cos(2 * np.pi * (doy - 20) / 365))
    load = load * (1 + 0.03 * rng.standard_normal(N_HOURS))
    load = np.clip(load, 5, None)

    # 太陽光出力比率: 日中の半正弦 x 季節振幅 x 曇天係数
    daylight = np.clip(np.sin(np.pi * (h - 6) / 12), 0, None)
    season = 0.75 + 0.25 * np.sin(2 * np.pi * (doy - 80) / 365)
    cloud_daily = np.clip(
        1 - 0.5 * np.abs(pd.Series(rng.standard_normal(366)).rolling(3, min_periods=1).mean().to_numpy()),
        0.2, 1,
    )
    cloud = cloud_daily[doy - 1]
    solar = np.clip(0.85 * daylight * season * cloud, 0, 1)

    # 風力出力比率: AR(1)的な風速変動 + 冬季強め + 夜間やや強め
    z = np.zeros(N_HOURS)
    e = rng.standard_normal(N_HOURS)
    for t in range(1, N_HOURS):
        z[t] = 0.97 * z[t - 1] + 0.25 * e[t]
    wind = np.clip(0.35 + 0.15 * z + 0.08 * np.cos(2 * np.pi * (doy - 10) / 365)
                   + 0.05 * np.cos(2 * np.pi * h / 24), 0, 1)

    return pd.DataFrame({
        "timestamp": idx,
        "solar_p_max_pu": solar,   # 列名は generators.name + "_p_max_pu" の命名規則
        "wind_p_max_pu": wind,
        "load_mw": load,           # 列名は loads.name + "_mw" の命名規則
    })


def generate_sample_config(path):
    """network_config.xlsx が存在しない場合に、リアリスティックなサンプルを自動生成する。"""
    timeseries = _build_synthetic_timeseries()
    buses = pd.DataFrame({"bus_name": ["bus"], "v_nom": [0.4]})
    generators = pd.DataFrame([
        {"name": "solar", "bus": "bus", "carrier": "solar", "capex": 700000.0,
         "marginal_cost": 0.0, "lifetime": 25, "p_nom_extendable": True},
        {"name": "wind", "bus": "bus", "carrier": "wind", "capex": 1100000.0,
         "marginal_cost": 0.0, "lifetime": 25, "p_nom_extendable": True},
        {"name": "diesel", "bus": "bus", "carrier": "diesel", "capex": 800000.0,
         "marginal_cost": 250.0, "lifetime": 20, "p_nom_extendable": True},
    ])
    storage_units = pd.DataFrame([
        {"name": "battery", "bus": "bus", "carrier": "battery", "capex": 600000.0,
         "max_hours": 4, "efficiency_store": 0.95, "efficiency_dispatch": 0.95,
         "lifetime": 12, "p_nom_extendable": True},
    ])
    loads = pd.DataFrame([{"name": "load", "bus": "bus"}])
    # 任意の電源の CAPEX(capital_cost)/OPEX(marginal_cost)を不確実性パラメータとして登録。
    # 行を追加/削除するだけでLHS/Sobol/PCEの次元数が自動追従する。
    uncertainty_params = pd.DataFrame([
        {"param_name": "solar_capex", "component_type": "generator", "component_name": "solar",
         "target_attribute": "capital_cost", "lower_bound": 500000, "upper_bound": 900000},
        {"param_name": "wind_capex", "component_type": "generator", "component_name": "wind",
         "target_attribute": "capital_cost", "lower_bound": 800000, "upper_bound": 1450000},
        {"param_name": "battery_capex", "component_type": "storage_unit", "component_name": "battery",
         "target_attribute": "capital_cost", "lower_bound": 300000, "upper_bound": 900000},
        {"param_name": "diesel_opex", "component_type": "generator", "component_name": "diesel",
         "target_attribute": "marginal_cost", "lower_bound": 200, "upper_bound": 300},
    ])

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        buses.to_excel(writer, sheet_name="buses", index=False)
        generators.to_excel(writer, sheet_name="generators", index=False)
        storage_units.to_excel(writer, sheet_name="storage_units", index=False)
        loads.to_excel(writer, sheet_name="loads", index=False)
        timeseries.to_excel(writer, sheet_name="timeseries", index=False)
        uncertainty_params.to_excel(writer, sheet_name="uncertainty_params", index=False)
    print(f"    [!] {path.name} が見つからないため、サンプル設定ファイルを自動生成しました: {path}")


_LINKS_COLUMNS = ["name", "bus0", "bus1", "carrier", "capex", "marginal_cost",
                   "efficiency", "lifetime", "p_nom_extendable"]


def _get_or_default(row, col, default):
    """row.get(col, default) と異なり、列は存在するがセルが空欄(NaN)の場合もdefaultを返す。"""
    val = row.get(col, default)
    return default if pd.isna(val) else val


def _to_bool(value):
    """Excel由来の値をboolへ変換する。"FALSE"/"0"等の文字列もFalseとして扱う。
    値がNaN(空欄セル)の場合はここでは判定せず、_validate_required_columns 側で
    事前にエラーとして検出する前提とする。"""
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


_REQUIRED_COLS = {
    "buses": ["bus_name"],
    "generators": ["name", "bus", "carrier", "capex", "marginal_cost", "lifetime", "p_nom_extendable"],
    "storage_units": ["name", "bus", "carrier", "capex", "max_hours",
                       "efficiency_store", "efficiency_dispatch", "lifetime", "p_nom_extendable"],
    "loads": ["name", "bus"],
    "links": ["name", "bus0", "bus1", "capex", "lifetime", "p_nom_extendable"],
    "uncertainty_params": ["param_name", "component_type", "component_name",
                            "target_attribute", "lower_bound", "upper_bound"],
}


def _validate_required_columns(df, sheet_name):
    """sheet_name の必須列について、列の欠落および空欄セル(NaN)を検出しエラーにする。
    値をfloat()/bool()変換する前にここで弾くことで、空欄がNaNとして
    サイレントにPyPSAへ渡ってしまう(または真偽値が意図せずTrueになる)事態を防ぐ。"""
    cols = _REQUIRED_COLS[sheet_name]
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"{sheet_name} シートに必須列がありません: {missing_cols}")
    for col in cols:
        na_rows = df.index[df[col].isna()]
        if len(na_rows):
            excel_rows = [int(i) + 2 for i in na_rows]  # 0-index -> Excel行番号(ヘッダー行+1)
            raise ValueError(
                f"{sheet_name}.{col} に空欄セルがあります (Excel行: {excel_rows})"
            )


def load_network_config(path):
    if not path.exists():
        generate_sample_config(path)

    xls = pd.ExcelFile(path)
    sheet_names = set(xls.sheet_names)
    required_sheets = {"buses", "generators", "storage_units", "loads", "timeseries", "uncertainty_params"}
    missing_sheets = required_sheets - sheet_names
    if missing_sheets:
        raise ValueError(f"{path.name} に必須シートがありません: {sorted(missing_sheets)}")

    cfg = {
        "buses": pd.read_excel(xls, sheet_name="buses"),
        "generators": pd.read_excel(xls, sheet_name="generators"),
        "storage_units": pd.read_excel(xls, sheet_name="storage_units"),
        "loads": pd.read_excel(xls, sheet_name="loads"),
        "timeseries": pd.read_excel(xls, sheet_name="timeseries"),
        "uncertainty_params": pd.read_excel(xls, sheet_name="uncertainty_params"),
    }
    # links は任意(オプション)シート: hydrogen electrolyzer/fuel cell等をリンクとして表現する場合に使用
    if "links" in sheet_names:
        cfg["links"] = pd.read_excel(xls, sheet_name="links")
    else:
        cfg["links"] = pd.DataFrame(columns=_LINKS_COLUMNS)

    for sheet_name in ("buses", "generators", "storage_units", "loads"):
        _validate_required_columns(cfg[sheet_name], sheet_name)
    if not cfg["links"].empty:
        _validate_required_columns(cfg["links"], "links")

    cfg["timeseries"]["timestamp"] = pd.to_datetime(cfg["timeseries"]["timestamp"])
    if len(cfg["timeseries"]) != N_HOURS:
        raise ValueError(
            f"timeseries シートの行数が{N_HOURS}(8760h)ではありません: {len(cfg['timeseries'])}行"
        )
    for col in cfg["timeseries"].columns:
        if col.endswith("_p_max_pu"):
            vals = cfg["timeseries"][col].to_numpy()
            if vals.min() < 0 or vals.max() > 1:
                raise ValueError(f"timeseries.{col} は0.0〜1.0の範囲である必要があります")

    _validate_uncertainty_params(cfg)
    return cfg


def _validate_uncertainty_params(cfg):
    udf = cfg["uncertainty_params"]
    _validate_required_columns(udf, "uncertainty_params")
    if udf.empty:
        raise ValueError("uncertainty_params シートに行がありません（最低1パラメータが必要です）")
    if udf["param_name"].duplicated().any():
        dup = udf.loc[udf["param_name"].duplicated(), "param_name"].tolist()
        raise ValueError(f"uncertainty_params.param_name に重複があります: {dup}")

    name_index = {
        "generator": set(cfg["generators"]["name"].astype(str)),
        "storage_unit": set(cfg["storage_units"]["name"].astype(str)),
        "link": set(cfg["links"]["name"].astype(str)) if not cfg["links"].empty else set(),
    }
    for _, row in udf.iterrows():
        pname = row["param_name"]
        ctype = str(row["component_type"]).strip().lower()
        cname = str(row["component_name"]).strip()
        attr = str(row["target_attribute"]).strip()
        if ctype not in COMPONENT_ATTR:
            raise ValueError(
                f"uncertainty_params: 不正な component_type='{ctype}' (param='{pname}'). "
                f"有効値: {sorted(COMPONENT_ATTR)}"
            )
        if attr not in VALID_TARGET_ATTRS:
            raise ValueError(
                f"uncertainty_params: 不正な target_attribute='{attr}' (param='{pname}'). "
                f"有効値: {sorted(VALID_TARGET_ATTRS)}"
            )
        if cname not in name_index[ctype]:
            raise ValueError(
                f"uncertainty_params: component_name='{cname}' (param='{pname}') が "
                f"{COMPONENT_ATTR[ctype]} シートに見つかりません"
            )
        if float(row["lower_bound"]) >= float(row["upper_bound"]):
            raise ValueError(f"uncertainty_params: lower_bound >= upper_bound (param='{pname}')")


def get_lifetime(cfg, component_type, name):
    sheet = COMPONENT_ATTR[component_type]
    df = cfg[sheet]
    row = df.loc[df["name"].astype(str) == name]
    if row.empty:
        raise KeyError(f"'{name}' が {sheet} シートに見つかりません")
    return float(row.iloc[0]["lifetime"])


# ----------------------------------------------------------------------------
# 2. 不確実性パラメータ (uncertainty_params シートから動的に構成)
# ----------------------------------------------------------------------------
def build_problem(cfg):
    udf = cfg["uncertainty_params"]
    names = udf["param_name"].astype(str).tolist()
    bounds = udf[["lower_bound", "upper_bound"]].astype(float).values.tolist()
    problem = {"num_vars": len(names), "names": names, "bounds": bounds}
    return problem, udf


def build_labels(udf):
    labels = []
    for _, row in udf.iterrows():
        tech = str(row["component_name"]).replace("_", " ").title()
        kind = "CAPEX" if str(row["target_attribute"]).strip() == "capital_cost" else "OPEX"
        labels.append(f"{tech} {kind}")
    return labels


def build_lifetimes(cfg, udf):
    """capital_cost を対象とするパラメータについてのみ (component_type, name) -> lifetime を引く"""
    lifetimes = {}
    for _, row in udf.iterrows():
        if str(row["target_attribute"]).strip() != "capital_cost":
            continue
        ctype = str(row["component_type"]).strip().lower()
        cname = str(row["component_name"]).strip()
        lifetimes[(ctype, cname)] = get_lifetime(cfg, ctype, cname)
    return lifetimes


# ----------------------------------------------------------------------------
# 3. PyPSA モデル構築 (Excel設定を動的に反映, 8760スナップショット, 任意電源対応)
# ----------------------------------------------------------------------------
def build_network(cfg):
    n = pypsa.Network()
    ts = cfg["timeseries"].set_index("timestamp")
    n.set_snapshots(ts.index)  # 1時間刻み x 8760 -> snapshot_weightings は既定で1

    carriers = set(cfg["generators"]["carrier"]) | set(cfg["storage_units"]["carrier"])
    if not cfg["links"].empty:
        carriers |= set(cfg["links"]["carrier"].dropna().astype(str))
    carriers |= {"AC"}
    n.add("Carrier", sorted(carriers))

    for _, row in cfg["buses"].iterrows():
        n.add("Bus", str(row["bus_name"]), v_nom=float(_get_or_default(row, "v_nom", 1.0)))

    for _, row in cfg["loads"].iterrows():
        name = str(row["name"])
        col = f"{name}_mw"
        if col not in ts.columns:
            raise KeyError(f"timeseries に負荷列 '{col}' がありません (loads.name='{name}')")
        n.add("Load", name, bus=str(row["bus"]), p_set=ts[col])

    for _, row in cfg["generators"].iterrows():
        name = str(row["name"])
        # 注: Generatorのefficiencyは(Linkと異なり)PyPSAのLOPFでは使われず
        # 燃料費/CO2排出換算のロジックもこのスクリプトには無いため、意味を持たない。
        # 誤解を避けるためGeneratorには渡さない(Link.efficiencyはp1=-p0*efficiencyの
        # フロー計算に実際に使われるため引き続き渡す)。
        kwargs = dict(
            bus=str(row["bus"]),
            carrier=str(row["carrier"]),
            p_nom_extendable=_to_bool(row["p_nom_extendable"]),
            capital_cost=annualized(float(row["capex"]), float(row["lifetime"])),
            marginal_cost=float(row["marginal_cost"]),
        )
        pu_col = f"{name}_p_max_pu"  # 命名規則: generators.name + "_p_max_pu"
        if pu_col in ts.columns:
            kwargs["p_max_pu"] = ts[pu_col]
        n.add("Generator", name, **kwargs)

    for _, row in cfg["storage_units"].iterrows():
        n.add(
            "StorageUnit", str(row["name"]),
            bus=str(row["bus"]), carrier=str(row["carrier"]),
            p_nom_extendable=_to_bool(row["p_nom_extendable"]),
            max_hours=float(row["max_hours"]),
            capital_cost=annualized(float(row["capex"]), float(row["lifetime"])),
            efficiency_store=float(row["efficiency_store"]),
            efficiency_dispatch=float(row["efficiency_dispatch"]),
            cyclic_state_of_charge=True,
        )

    for _, row in cfg["links"].iterrows():
        n.add(
            "Link", str(row["name"]),
            bus0=str(row["bus0"]), bus1=str(row["bus1"]),
            carrier=str(_get_or_default(row, "carrier", "")),
            p_nom_extendable=_to_bool(row["p_nom_extendable"]),
            capital_cost=annualized(float(row["capex"]), float(row["lifetime"])),
            marginal_cost=float(_get_or_default(row, "marginal_cost", 0.0)),
            efficiency=float(_get_or_default(row, "efficiency", 1.0)),
        )
    return n


def apply_uncertainty(n, udf, params, lifetimes):
    """LHS/Sobolでサンプリングされた params を、uncertainty_params定義に従って
    対応コンポーネントの capital_cost(CRF年換算) または marginal_cost に適用する。"""
    for value, (_, row) in zip(params, udf.iterrows()):
        ctype = str(row["component_type"]).strip().lower()
        cname = str(row["component_name"]).strip()
        attr = str(row["target_attribute"]).strip()
        comp_df = getattr(n, COMPONENT_ATTR[ctype])
        if attr == "capital_cost":
            comp_df.loc[cname, attr] = annualized(float(value), lifetimes[(ctype, cname)])
        else:
            comp_df.loc[cname, attr] = float(value)


def run_pypsa(cfg, params, udf, lifetimes):
    """1サンプル分のLP(8760h)を解き、最適費用・全電源の容量・発電/充放電量を返す。"""
    n = build_network(cfg)
    apply_uncertainty(n, udf, params, lifetimes)

    status, cond = n.optimize(solver_name="highs", log_to_console=False,
                              include_objective_constant=False)
    if status != "ok":
        raise RuntimeError(f"PyPSA最適化に失敗: status={status}, condition={cond}, params={params}")

    w = n.snapshot_weightings.generators
    result = {"total_cost": float(n.objective)}

    for name in n.generators.index:
        result[f"{name}_mw"] = float(n.generators.at[name, "p_nom_opt"])
        result[f"{name}_mwh"] = float((n.generators_t.p[name] * w).sum())

    for name in n.storage_units.index:
        result[f"{name}_mw"] = float(n.storage_units.at[name, "p_nom_opt"])
        disp = n.storage_units_t.p[name].clip(lower=0)
        result[f"{name}_discharge_mwh"] = float((disp * w).sum())

    for name in n.links.index:
        result[f"{name}_mw"] = float(n.links.at[name, "p_nom_opt"])
        flow = n.links_t.p0[name].clip(lower=0)
        result[f"{name}_mwh"] = float((flow * w).sum())

    return result


# ----------------------------------------------------------------------------
# 4. メイン
# ----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="LHS+PCEサロゲートによるPyPSA GSA (任意電源・任意パラメータ対応)")
    p.add_argument("--config", default=str(CONFIG_PATH), help="network_config.xlsx のパス")
    p.add_argument("--n-lhs", type=int, default=N_LHS, help="LHSサンプル数 (既定320)")
    p.add_argument("--n-sobol", type=int, default=N_SOBOL, help="Sobolベースサンプル数 (既定1024)")
    return p.parse_args()


def main():
    args = parse_args()
    n_lhs = args.n_lhs
    n_sobol = args.n_sobol
    config_path = Path(args.config)

    t0 = time.time()
    print("=" * 78)
    print(f" network_config.xlsx -> LHS({n_lhs}) -> PyPSA(8760h) -> PCEサロゲート"
          f" -> Sobol GSA({n_sobol})")
    print("=" * 78)

    # --- 出力先ディレクトリ (タイムスタンプ付き) -----------------------------------
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    output_dir = os.path.join(BASE_DIR, "results", timestamp)
    os.makedirs(output_dir, exist_ok=True)
    print(f"[0] 出力先フォルダ: {output_dir}")

    # --- ネットワーク設定読み込み (無ければ自動生成) --------------------------------
    print(f"[1] ネットワーク設定を読み込み中: {config_path}")
    cfg = load_network_config(config_path)
    problem, udf = build_problem(cfg)
    names = problem["names"]
    labels = build_labels(udf)
    lifetimes = build_lifetimes(cfg, udf)
    num_vars = problem["num_vars"]
    print(f"    電源: generators={cfg['generators']['name'].tolist()}, "
          f"storage_units={cfg['storage_units']['name'].tolist()}, "
          f"links={cfg['links']['name'].tolist() if not cfg['links'].empty else []}")
    print(f"    不確実性パラメータ({num_vars}次元): {names}")
    print(f"    スナップショット数: {len(cfg['timeseries'])} (8760h, 各1時間重み)")

    # --- LHS -----------------------------------------------------------------
    X = latin_sample.sample(problem, n_lhs, seed=SEED)
    print(f"[2] LHSサンプル生成: {X.shape}")

    # --- PyPSA ----------------------------------------------------------------
    print("[3] PyPSA最適化(8760h)を実行中 ...")
    rows = []
    report_every = max(1, n_lhs // 8)
    try:
        for i, x in enumerate(X):
            rows.append(run_pypsa(cfg, x, udf, lifetimes))
            if (i + 1) % report_every == 0 or (i + 1) == n_lhs:
                print(f"    {i + 1}/{n_lhs} done  ({time.time() - t0:.0f}s)")
    except Exception:
        # 失敗時、それまでに完了した分だけでも退避する(8760h LPは1点あたり高コストなため)。
        # フォールバックはせず、失敗自体はここで揉み消さずに再送出する。
        if rows:
            partial_res = pd.DataFrame(rows)
            partial_df = pd.concat(
                [pd.DataFrame(X[:len(rows)], columns=names), partial_res], axis=1
            )
            partial_df.index.name = "sample"
            partial_path = os.path.join(
                output_dir, f"pypsa_lhs_{n_lhs}_results_partial_{len(rows)}.csv"
            )
            partial_df.to_csv(partial_path)
            print(f"    [!] {len(rows)}/{n_lhs} 件の完了分を退避しました: {partial_path}",
                  file=sys.stderr)
        raise
    res = pd.DataFrame(rows)
    df = pd.concat([pd.DataFrame(X, columns=names), res], axis=1)
    df.index.name = "sample"
    lhs_csv_name = f"pypsa_lhs_{n_lhs}_results.csv"
    lhs_nc_name = f"pypsa_lhs_{n_lhs}_results.nc"
    df.to_csv(os.path.join(output_dir, lhs_csv_name))

    ds = xr.Dataset(
        {c: ("sample", df[c].to_numpy()) for c in df.columns},
        coords={"sample": np.arange(n_lhs)},
        attrs={
            "title": "PyPSA LHS results",
            "uncertainty_params": ", ".join(names),
            "total_cost_unit": "USD/year", "capacity_unit": "MW", "energy_unit": "MWh/year",
        },
    )
    ds.to_netcdf(os.path.join(output_dir, lhs_nc_name), engine="netcdf4")
    print(f"    保存: {lhs_csv_name} / {lhs_nc_name}")

    y = df["total_cost"].to_numpy()

    # --- PCEサロゲート -----------------------------------------------------------
    print("[4] PCEサロゲート学習 (StandardScaler -> Poly(2) -> RidgeCV)")
    pce = Pipeline([
        ("scaler", StandardScaler()),
        # RidgeCVがfit_intercept=True(既定)で自前の切片を推定するため、
        # PolynomialFeatures側のバイアス列(定数項)は不要(二重の切片を避ける)。
        ("poly", PolynomialFeatures(degree=2, include_bias=False)),
        ("ridge", RidgeCV(alphas=np.logspace(-4, 3, 30), cv=10)),
    ])
    pce.fit(X, y)
    y_fit = pce.predict(X)
    kf = KFold(n_splits=10, shuffle=True, random_state=SEED)
    y_cv = cross_val_predict(pce, X, y, cv=kf)
    r2_train = r2_score(y, y_fit)
    r2_cv = r2_score(y, y_cv)
    rmse_train = float(np.sqrt(mean_squared_error(y, y_fit)))
    rmse_cv = float(np.sqrt(mean_squared_error(y, y_cv)))
    n_basis = pce.named_steps["poly"].n_output_features_
    print(f"    次元数 = {num_vars},  基底数 = {n_basis},  最適alpha = {pce.named_steps['ridge'].alpha_:.4g}")
    print(f"    学習 R^2 = {r2_train:.5f}   CV(10-fold) R^2 = {r2_cv:.5f}")
    print(f"    学習 RMSE = {rmse_train:,.1f}   CV RMSE = {rmse_cv:,.1f}  [USD/year]")

    # --- Sobol on surrogate ------------------------------------------------------
    print(f"[5] サロゲート上でSobol解析 (N={n_sobol} -> {n_sobol * (2 * num_vars + 2)}点)")
    Xs = sobol_sample.sample(problem, n_sobol, calc_second_order=True, seed=SEED)
    ts_sobol = time.time()
    Ys = pce.predict(Xs)
    print(f"    サロゲート評価 {len(Xs)}点: {(time.time() - ts_sobol) * 1000:.1f} ms")
    Si = sobol_analyze.analyze(problem, Ys, calc_second_order=True, seed=SEED,
                               print_to_console=False)
    sens = pd.DataFrame({"S1": Si["S1"], "S1_conf": Si["S1_conf"],
                         "ST": Si["ST"], "ST_conf": Si["ST_conf"]}, index=names)
    print("\n--- 第1次 (S1) / 総 (ST) 感度指標 ---")
    print(sens.round(4).to_string())

    s2 = np.nan_to_num(np.asarray(Si["S2"], dtype=float), nan=0.0)
    s2 = s2 + s2.T  # SALibは上三角のみ
    s2_df = pd.DataFrame(s2, index=names, columns=names)
    print("\n--- Sobol 2次感度指標 (Sij) マトリクス (対角=0) ---")
    print(s2_df.round(4).to_string())
    s2_df.to_csv(os.path.join(output_dir, "sobol_s2_matrix.csv"))

    # --- Cij -------------------------------------------------------------------
    poly = pce.named_steps["poly"]
    coef = pce.named_steps["ridge"].coef_
    fnames = poly.get_feature_names_out([f"x{i}" for i in range(num_vars)])
    cij = np.zeros((num_vars, num_vars))
    for fname, c in zip(fnames, coef):
        if "^2" in fname:  # 対角: 2乗項 Cii
            i = int(fname.split("^")[0][1:])
            cij[i, i] = c
        elif " " in fname:  # 交差項 Cij
            i, j = (int(t[1:]) for t in fname.split())
            cij[i, j] = cij[j, i] = c
    cij_df = pd.DataFrame(cij, index=names, columns=names)
    print("\n--- PCE交差項係数 (Cij) マトリクス [標準化入力空間, 単位: USD/year] "
          "(対角=2乗項Cii) ---")
    print(cij_df.round(1).to_string())
    cij_df.to_csv(os.path.join(output_dir, "pce_interaction_matrix.csv"))

    print("\n--- パラメータペアの関係判定 (総費用最小化: 包絡線定理 dC/dθi = 最適量_i) ---")
    print("    Cij = d2C/dθi dθj = d(最適量_i)/dθj")
    print("    Cij < 0: 補完 (θj上昇で両技術の最適容量が同時に減少 / 総費用増を抑制)")
    print("    Cij > 0: 代替 (θj上昇で技術iへ容量がシフト・置換)")
    for i in range(num_vars):
        for j in range(i + 1, num_vars):
            c = cij[i, j]
            tag = ("【補完関係 (容量連動・シナジー)】" if c < 0
                   else "【代替関係 (技術競合・置換)】" if c > 0 else "【無相互作用】")
            print(f"  {names[i]:>21s} x {names[j]:<21s}  Cij = {c:>12,.1f}  {tag}  "
                  f"(Sij = {s2_df.iloc[i, j]:.4f})")

    # --- プロット ----------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(max(9, 2.6 * num_vars + 4), 5.5))
    ax = axes[0]
    pos = np.arange(num_vars)
    w = 0.38
    ax.bar(pos - w / 2, sens["S1"], w, yerr=sens["S1_conf"], label="S1", color="#4c72b0", capsize=3)
    ax.bar(pos + w / 2, sens["ST"], w, yerr=sens["ST_conf"], label="ST", color="#dd8452", capsize=3)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("Sobol index")
    ax.set_title("Sobol sensitivity (on PCE surrogate)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1]
    ax.scatter(y / 1e6, y_fit / 1e6, s=14, alpha=0.7, color="#4c72b0")
    lo, hi = min(y.min(), y_fit.min()) / 1e6, max(y.max(), y_fit.max()) / 1e6
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
    ax.set_xlabel("PyPSA total cost [M$/year]")
    ax.set_ylabel("PCE predicted cost [M$/year]")
    ax.set_title("PyPSA vs PCE surrogate")
    ax.text(0.05, 0.95, f"$R^2$ = {r2_train:.4f}\nCV $R^2$ = {r2_cv:.4f}",
            transform=ax.transAxes, va="top", bbox=dict(boxstyle="round", fc="white", ec="gray"))
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "pypsa_pce_gsa_results.png"), dpi=200)
    plt.close(fig)

    print("\n" + "=" * 78)
    print(f"完了 ({time.time() - t0:.0f}s)")
    print(f"成果物の保存先: {output_dir}")
    print(f"  - {lhs_csv_name}")
    print(f"  - {lhs_nc_name}")
    print("  - sobol_s2_matrix.csv")
    print("  - pce_interaction_matrix.csv")
    print("  - pypsa_pce_gsa_results.png")
    print("=" * 78)


if __name__ == "__main__":
    main()
