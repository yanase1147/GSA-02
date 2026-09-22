#!/usr/bin/env python
"""
Excel(network_config.xlsx) -> LHS(320点) -> PyPSA(8760h) -> PCEサロゲート
                            -> サロゲート上GSA(10,240点) -> Sobol / Cij 解析

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

# ----------------------------------------------------------------------------
# 1. 不確実性パラメータ (LHS/Sobolで振る4つのコストパラメータ)
# ----------------------------------------------------------------------------
PROBLEM = {
    "num_vars": 4,
    "names": ["solar_cost", "wind_cost", "battery_cost", "diesel_marginal_cost"],
    # 初期建設費 CAPEX [$/MW] (太陽光/風力/蓄電池) と ディーゼル可変費 [$/MWh]
    "bounds": [[500000, 900000], [800000, 1450000], [300000, 900000], [200, 300]],
}
NAMES = PROBLEM["names"]
LABELS = ["Solar CAPEX", "Wind CAPEX", "Battery CAPEX", "Diesel fuel cost"]

DISCOUNT_RATE = 0.05  # 割引率 (CAPEXの年換算 CRF 計算用)


def crf(rate, lifetime):
    """資本回収係数 CRF = r(1+r)^n / ((1+r)^n - 1)"""
    if rate == 0:
        return 1.0 / lifetime
    f = (1 + rate) ** lifetime
    return rate * f / (f - 1)


def annualized(capex, lifetime, rate=DISCOUNT_RATE):
    """初期建設費 [$/MW] -> 年換算コスト capital_cost [$/MW/year]"""
    return capex * crf(rate, lifetime)


# ----------------------------------------------------------------------------
# 2. network_config.xlsx の読み込み / 自動生成
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
        "solar_p_max_pu": solar,
        "wind_p_max_pu": wind,
        "load_mw": load,
    })


def generate_sample_config(path):
    """network_config.xlsx が存在しない場合に、リアリスティックなサンプルを自動生成する。"""
    timeseries = _build_synthetic_timeseries()
    buses = pd.DataFrame({"bus_name": ["bus"], "v_nom": [0.4]})
    generators = pd.DataFrame([
        {"name": "solar", "bus": "bus", "carrier": "solar", "capex": 700000.0,
         "marginal_cost": 0.0, "efficiency": 1.0, "lifetime": 25, "p_nom_extendable": True},
        {"name": "wind", "bus": "bus", "carrier": "wind", "capex": 1100000.0,
         "marginal_cost": 0.0, "efficiency": 1.0, "lifetime": 25, "p_nom_extendable": True},
        {"name": "diesel", "bus": "bus", "carrier": "diesel", "capex": 800000.0,
         "marginal_cost": 250.0, "efficiency": 0.40, "lifetime": 20, "p_nom_extendable": True},
    ])
    storage_units = pd.DataFrame([
        {"name": "battery", "bus": "bus", "carrier": "battery", "capex": 600000.0,
         "max_hours": 4, "efficiency_store": 0.95, "efficiency_dispatch": 0.95,
         "lifetime": 12, "p_nom_extendable": True},
    ])
    loads = pd.DataFrame([{"name": "load", "bus": "bus"}])

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        buses.to_excel(writer, sheet_name="buses", index=False)
        generators.to_excel(writer, sheet_name="generators", index=False)
        storage_units.to_excel(writer, sheet_name="storage_units", index=False)
        loads.to_excel(writer, sheet_name="loads", index=False)
        timeseries.to_excel(writer, sheet_name="timeseries", index=False)
    print(f"    [!] {path.name} が見つからないため、サンプル設定ファイルを自動生成しました: {path}")


def load_network_config(path):
    if not path.exists():
        generate_sample_config(path)
    cfg = {
        "buses": pd.read_excel(path, sheet_name="buses"),
        "generators": pd.read_excel(path, sheet_name="generators"),
        "storage_units": pd.read_excel(path, sheet_name="storage_units"),
        "loads": pd.read_excel(path, sheet_name="loads"),
        "timeseries": pd.read_excel(path, sheet_name="timeseries"),
    }
    cfg["timeseries"]["timestamp"] = pd.to_datetime(cfg["timeseries"]["timestamp"])
    if len(cfg["timeseries"]) != N_HOURS:
        raise ValueError(
            f"timeseries シートの行数が{N_HOURS}(8760h)ではありません: {len(cfg['timeseries'])}行"
        )
    for col in ("solar_p_max_pu", "wind_p_max_pu"):
        vals = cfg["timeseries"][col].to_numpy()
        if vals.min() < 0 or vals.max() > 1:
            raise ValueError(f"timeseries.{col} は0.0〜1.0の範囲である必要があります")
    return cfg


def get_lifetime(cfg, name):
    for sheet in ("generators", "storage_units"):
        df = cfg[sheet]
        row = df.loc[df["name"] == name]
        if not row.empty:
            return float(row.iloc[0]["lifetime"])
    raise KeyError(f"'{name}' が generators / storage_units シートに見つかりません")


# ----------------------------------------------------------------------------
# 3. PyPSA モデル構築 (Excel設定を動的に反映, 8760スナップショット)
# ----------------------------------------------------------------------------
def build_network(cfg):
    n = pypsa.Network()
    ts = cfg["timeseries"].set_index("timestamp")
    n.set_snapshots(ts.index)  # 1時間刻み x 8760 -> snapshot_weightings は既定で1

    carriers = sorted(set(cfg["generators"]["carrier"]) | set(cfg["storage_units"]["carrier"]) | {"AC"})
    n.add("Carrier", carriers)

    for _, row in cfg["buses"].iterrows():
        n.add("Bus", str(row["bus_name"]), v_nom=float(row.get("v_nom", 1.0)))

    for _, row in cfg["loads"].iterrows():
        n.add("Load", str(row["name"]), bus=str(row["bus"]), p_set=ts["load_mw"])

    for _, row in cfg["generators"].iterrows():
        carrier = str(row["carrier"])
        kwargs = dict(
            bus=str(row["bus"]),
            carrier=carrier,
            p_nom_extendable=bool(row["p_nom_extendable"]),
            capital_cost=annualized(float(row["capex"]), float(row["lifetime"])),
            marginal_cost=float(row["marginal_cost"]),
            efficiency=float(row.get("efficiency", 1.0)),
        )
        if carrier == "solar":
            kwargs["p_max_pu"] = ts["solar_p_max_pu"]
        elif carrier == "wind":
            kwargs["p_max_pu"] = ts["wind_p_max_pu"]
        n.add("Generator", str(row["name"]), **kwargs)

    for _, row in cfg["storage_units"].iterrows():
        n.add(
            "StorageUnit", str(row["name"]),
            bus=str(row["bus"]), carrier=str(row["carrier"]),
            p_nom_extendable=bool(row["p_nom_extendable"]),
            max_hours=float(row["max_hours"]),
            capital_cost=annualized(float(row["capex"]), float(row["lifetime"])),
            efficiency_store=float(row["efficiency_store"]),
            efficiency_dispatch=float(row["efficiency_dispatch"]),
            cyclic_state_of_charge=True,
        )
    return n


def run_pypsa(cfg, params, lifetimes):
    """1サンプル分のLP(8760h)を解き、最適費用・容量・発電/充放電量を返す。"""
    n = build_network(cfg)
    n.generators.loc["solar", "capital_cost"] = annualized(params[0], lifetimes["solar"])
    n.generators.loc["wind", "capital_cost"] = annualized(params[1], lifetimes["wind"])
    n.storage_units.loc["battery", "capital_cost"] = annualized(params[2], lifetimes["battery"])
    n.generators.loc["diesel", "marginal_cost"] = params[3]

    status, cond = n.optimize(solver_name="highs", log_to_console=False,
                              include_objective_constant=False)
    if status != "ok":
        raise RuntimeError(f"PyPSA最適化に失敗: status={status}, condition={cond}, params={params}")

    w = n.snapshot_weightings.generators
    solar_mwh = float((n.generators_t.p["solar"] * w).sum())
    wind_mwh = float((n.generators_t.p["wind"] * w).sum())
    diesel_mwh = float((n.generators_t.p["diesel"] * w).sum())
    batt_p = n.storage_units_t.p["battery"]
    battery_discharge_mwh = float((batt_p.clip(lower=0) * w).sum())

    return {
        "total_cost": float(n.objective),
        "solar_mw": float(n.generators.at["solar", "p_nom_opt"]),
        "wind_mw": float(n.generators.at["wind", "p_nom_opt"]),
        "battery_mw": float(n.storage_units.at["battery", "p_nom_opt"]),
        "diesel_mw": float(n.generators.at["diesel", "p_nom_opt"]),
        "solar_mwh": solar_mwh,
        "wind_mwh": wind_mwh,
        "battery_discharge_mwh": battery_discharge_mwh,
        "diesel_mwh": diesel_mwh,
    }


# ----------------------------------------------------------------------------
# 4. メイン
# ----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="LHS+PCEサロゲートによるPyPSA GSA")
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
    lifetimes = {name: get_lifetime(cfg, name) for name in ("solar", "wind", "battery", "diesel")}
    print(f"    耐用年数: {lifetimes}  (割引率={DISCOUNT_RATE})")
    print(f"    スナップショット数: {len(cfg['timeseries'])} (8760h, 各1時間重み)")

    # --- LHS -----------------------------------------------------------------
    X = latin_sample.sample(PROBLEM, n_lhs, seed=SEED)
    print(f"[2] LHSサンプル生成: {X.shape}")

    # --- PyPSA ----------------------------------------------------------------
    print("[3] PyPSA最適化(8760h)を実行中 ...")
    rows = []
    report_every = max(1, n_lhs // 8)
    for i, x in enumerate(X):
        rows.append(run_pypsa(cfg, x, lifetimes))
        if (i + 1) % report_every == 0 or (i + 1) == n_lhs:
            print(f"    {i + 1}/{n_lhs} done  ({time.time() - t0:.0f}s)")
    res = pd.DataFrame(rows)
    df = pd.concat([pd.DataFrame(X, columns=NAMES), res], axis=1)
    df.index.name = "sample"
    df.to_csv(os.path.join(output_dir, "pypsa_lhs_320_results.csv"))

    ds = xr.Dataset(
        {c: ("sample", df[c].to_numpy()) for c in df.columns},
        coords={"sample": np.arange(n_lhs)},
        attrs={
            "title": "PyPSA LHS results",
            "units_costs": "CAPEX USD/MW (overnight), diesel USD/MWh",
            "total_cost_unit": "USD/year", "capacity_unit": "MW", "energy_unit": "MWh/year",
        },
    )
    ds.to_netcdf(os.path.join(output_dir, "pypsa_lhs_320_results.nc"), engine="netcdf4")
    print("    保存: pypsa_lhs_320_results.csv / pypsa_lhs_320_results.nc")

    y = df["total_cost"].to_numpy()

    # --- PCEサロゲート -----------------------------------------------------------
    print("[4] PCEサロゲート学習 (StandardScaler -> Poly(2) -> RidgeCV)")
    pce = Pipeline([
        ("scaler", StandardScaler()),
        ("poly", PolynomialFeatures(degree=2, include_bias=True)),
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
    print(f"    基底数 = {n_basis},  最適alpha = {pce.named_steps['ridge'].alpha_:.4g}")
    print(f"    学習 R^2 = {r2_train:.5f}   CV(10-fold) R^2 = {r2_cv:.5f}")
    print(f"    学習 RMSE = {rmse_train:,.1f}   CV RMSE = {rmse_cv:,.1f}  [USD/year]")

    # --- Sobol on surrogate ------------------------------------------------------
    print(f"[5] サロゲート上でSobol解析 (N={n_sobol} -> {n_sobol * (2 * 4 + 2)}点)")
    Xs = sobol_sample.sample(PROBLEM, n_sobol, calc_second_order=True, seed=SEED)
    ts_sobol = time.time()
    Ys = pce.predict(Xs)
    print(f"    サロゲート評価 {len(Xs)}点: {(time.time() - ts_sobol) * 1000:.1f} ms")
    Si = sobol_analyze.analyze(PROBLEM, Ys, calc_second_order=True, seed=SEED,
                               print_to_console=False)
    sens = pd.DataFrame({"S1": Si["S1"], "S1_conf": Si["S1_conf"],
                         "ST": Si["ST"], "ST_conf": Si["ST_conf"]}, index=NAMES)
    print("\n--- 第1次 (S1) / 総 (ST) 感度指標 ---")
    print(sens.round(4).to_string())

    s2 = np.nan_to_num(np.asarray(Si["S2"], dtype=float), nan=0.0)
    s2 = s2 + s2.T  # SALibは上三角のみ
    s2_df = pd.DataFrame(s2, index=NAMES, columns=NAMES)
    print("\n--- Sobol 2次感度指標 (Sij) マトリクス (対角=0) ---")
    print(s2_df.round(4).to_string())
    s2_df.to_csv(os.path.join(output_dir, "sobol_s2_matrix.csv"))

    # --- Cij -------------------------------------------------------------------
    poly = pce.named_steps["poly"]
    coef = pce.named_steps["ridge"].coef_
    fnames = poly.get_feature_names_out([f"x{i}" for i in range(4)])
    cij = np.zeros((4, 4))
    for name, c in zip(fnames, coef):
        if "^2" in name:  # 対角: 2乗項 Cii
            i = int(name.split("^")[0][1:])
            cij[i, i] = c
        elif " " in name:  # 交差項 Cij
            i, j = (int(t[1:]) for t in name.split())
            cij[i, j] = cij[j, i] = c
    cij_df = pd.DataFrame(cij, index=NAMES, columns=NAMES)
    print("\n--- PCE交差項係数 (Cij) マトリクス [標準化入力空間, 単位: USD/year] "
          "(対角=2乗項Cii) ---")
    print(cij_df.round(1).to_string())
    cij_df.to_csv(os.path.join(output_dir, "pce_interaction_matrix.csv"))

    print("\n--- 技術ペアの関係判定 (総費用最小化: 包絡線定理 dC/dθi = 最適量_i) ---")
    print("    Cij = d2C/dθi dθj = d(最適量_i)/dθj")
    print("    Cij < 0: 補完 (θj上昇で両技術の最適容量が同時に減少 / 総費用増を抑制)")
    print("    Cij > 0: 代替 (θj上昇で技術iへ容量がシフト・置換)")
    for i in range(4):
        for j in range(i + 1, 4):
            c = cij[i, j]
            tag = ("【補完関係 (容量連動・シナジー)】" if c < 0
                   else "【代替関係 (技術競合・置換)】" if c > 0 else "【無相互作用】")
            print(f"  {NAMES[i]:>21s} x {NAMES[j]:<21s}  Cij = {c:>12,.1f}  {tag}  "
                  f"(Sij = {s2_df.iloc[i, j]:.4f})")

    # --- プロット ----------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    pos = np.arange(4)
    w = 0.38
    ax.bar(pos - w / 2, sens["S1"], w, yerr=sens["S1_conf"], label="S1", color="#4c72b0", capsize=3)
    ax.bar(pos + w / 2, sens["ST"], w, yerr=sens["ST_conf"], label="ST", color="#dd8452", capsize=3)
    ax.set_xticks(pos)
    ax.set_xticklabels(LABELS, rotation=15)
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
    print("  - pypsa_lhs_320_results.csv")
    print("  - pypsa_lhs_320_results.nc")
    print("  - sobol_s2_matrix.csv")
    print("  - pce_interaction_matrix.csv")
    print("  - pypsa_pce_gsa_results.png")
    print("=" * 78)


if __name__ == "__main__":
    main()
