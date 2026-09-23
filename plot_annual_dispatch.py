#!/usr/bin/env python
"""
network_config.xlsx の代表値(不確実性パラメータの中央値/ベースライン値)を用いて
PyPSA の 8,760時間最適化を1回実行し、以下を出力する。

  - 年間集計指標(適用コスト・最適容量・容量比率・年間発電量・発電量比率)
  - 年間8,760時間の需給バランス/SOC推移グラフ
  - 月別(1〜12月)の需給バランス/SOC推移グラフ (12枚)
  - 年間/月別サマリーCSV, 全時系列CSV, 実行ログ

出力先: result_d/YYYY-MM-DD_HH-MM/ (実行の都度、タイムスタンプ付きで自動作成)

実行:  python plot_annual_dispatch.py
       python plot_annual_dispatch.py --config network_config.xlsx
"""
import sys

# Windows既定コードページ(cp932等)では日本語出力が文字化けするため、UTF-8に固定する。
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

import argparse
import logging
import os
import time
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pypsa_pce_gsa import (
    BASE_DIR,
    CONFIG_PATH,
    apply_uncertainty,
    build_lifetimes,
    build_network,
    build_problem,
    load_network_config,
)

# pypsa_pce_gsa 側でWARNING抑制されているログレベルを、本スクリプトでは
# ソルバーログもexecution.logへ残すためINFOへ引き上げる。
logging.getLogger("pypsa").setLevel(logging.INFO)
logging.getLogger("linopy").setLevel(logging.WARNING)
logging.getLogger("highspy").setLevel(logging.WARNING)

MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
RESOURCES = ["Solar", "Wind", "Battery", "Diesel"]

# ----------------------------------------------------------------------------
# 統一カラーパレット(全プロット共通)
# ----------------------------------------------------------------------------
COLORS = {
    "solar": "#F1C40F",
    "wind": "#2980B9",
    "battery_discharge": "#2ECC71",
    "battery_charge": "#27AE60",
    "diesel": "#7F8C8D",
    "load": "#E74C3C",
}


class Tee:
    """標準出力/エラー出力を複数ストリーム(コンソール + ログファイル)へ同時書き込みする。"""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
        return len(data)

    def flush(self):
        for s in self.streams:
            s.flush()


# ----------------------------------------------------------------------------
# 出力先ディレクトリ
# ----------------------------------------------------------------------------
def make_output_dirs():
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    output_dir = os.path.join(BASE_DIR, "result_d", timestamp)
    monthly_dir = os.path.join(output_dir, "monthly_plots")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(monthly_dir, exist_ok=True)
    return output_dir, monthly_dir


# ----------------------------------------------------------------------------
# 代表値(中央値)ネットワークの構築
# ----------------------------------------------------------------------------
def representative_params(udf):
    """不確実性パラメータの代表値 = (lower_bound + upper_bound) / 2 (中央値)。"""
    lo = udf["lower_bound"].astype(float)
    hi = udf["upper_bound"].astype(float)
    return ((lo + hi) / 2.0).to_numpy()


def build_representative_network(cfg):
    problem, udf = build_problem(cfg)
    lifetimes = build_lifetimes(cfg, udf)
    params = representative_params(udf)
    n = build_network(cfg)
    apply_uncertainty(n, udf, params, lifetimes)
    return n, udf, params


def print_applied_costs(n):
    """代表値として適用された各電源のコスト(CAPEX年換算 or OPEX)を表示する。"""
    rows = []
    if "solar" in n.generators.index:
        rows.append(("Solar", "CAPEX", float(n.generators.at["solar", "capital_cost"]), "$/MW/year"))
    if "wind" in n.generators.index:
        rows.append(("Wind", "CAPEX", float(n.generators.at["wind", "capital_cost"]), "$/MW/year"))
    if "battery" in n.storage_units.index:
        rows.append(("Battery", "CAPEX", float(n.storage_units.at["battery", "capital_cost"]), "$/MW/year"))
    if "diesel" in n.generators.index:
        rows.append(("Diesel", "OPEX", float(n.generators.at["diesel", "marginal_cost"]), "$/MWh"))

    print("\n[代表値コスト] 適用コスト (中央値/ベースライン)")
    for name, kind, val, unit in rows:
        print(f"    {name:8s} {kind:5s} = {val:14,.2f} {unit}")
    return rows


def run_optimization(n):
    """n.optimize()を実行する。ソルバーログはコンソールへライブ表示しつつ、Tee経由で
    execution.logにも書き出す(pypsaロガー/進捗バー経由分)。
    注: HiGHSネイティブのバナー行はOSファイルディスクリプタへ直接出力されるため、
    Pythonのsys.stdout差し替え(Tee)だけでは一部captureされない場合があるが、
    実端末のfd 1/2を書き換える方式は環境によってハングするリスクがあるため採用しない。"""
    status, cond = n.optimize(solver_name="highs", log_to_console=True,
                               include_objective_constant=False)
    if status != "ok":
        raise RuntimeError(f"PyPSA最適化に失敗しました: status={status}, condition={cond}")
    return status, cond


# ----------------------------------------------------------------------------
# 結果集計
# ----------------------------------------------------------------------------
def compute_dispatch_frame(n):
    idx = n.snapshots
    w = n.snapshot_weightings.generators

    df = pd.DataFrame(index=idx)
    df["solar_mw"] = n.generators_t.p["solar"] if "solar" in n.generators.index else 0.0
    df["wind_mw"] = n.generators_t.p["wind"] if "wind" in n.generators.index else 0.0
    df["diesel_mw"] = n.generators_t.p["diesel"] if "diesel" in n.generators.index else 0.0

    if "battery" in n.storage_units.index:
        batt_p = n.storage_units_t.p["battery"]
        df["battery_discharge_mw"] = batt_p.clip(lower=0)
        df["battery_charge_mw"] = -batt_p.clip(upper=0)
        df["battery_soc_mwh"] = n.storage_units_t.state_of_charge["battery"]
    else:
        df["battery_discharge_mw"] = 0.0
        df["battery_charge_mw"] = 0.0
        df["battery_soc_mwh"] = 0.0

    if len(n.loads_t.p_set.columns) > 0:
        df["load_mw"] = n.loads_t.p_set.sum(axis=1)
    else:
        df["load_mw"] = n.loads_t.p.sum(axis=1)

    df["weight_h"] = w.to_numpy()
    return df


def compute_capacity_table(n):
    caps = {}
    if "solar" in n.generators.index:
        caps["Solar"] = float(n.generators.at["solar", "p_nom_opt"])
    if "wind" in n.generators.index:
        caps["Wind"] = float(n.generators.at["wind", "p_nom_opt"])
    if "battery" in n.storage_units.index:
        caps["Battery"] = float(n.storage_units.at["battery", "p_nom_opt"])
    if "diesel" in n.generators.index:
        caps["Diesel"] = float(n.generators.at["diesel", "p_nom_opt"])
    total = sum(caps.values())
    ratio = {k: (v / total * 100.0 if total > 0 else 0.0) for k, v in caps.items()}
    return caps, ratio, total


def compute_annual_generation(df):
    w = df["weight_h"]
    gen = {
        "Solar": float((df["solar_mw"] * w).sum()),
        "Wind": float((df["wind_mw"] * w).sum()),
        "Diesel": float((df["diesel_mw"] * w).sum()),
        "Battery": float((df["battery_discharge_mw"] * w).sum()),
    }
    battery_charge_mwh = float((df["battery_charge_mw"] * w).sum())
    total_demand = float((df["load_mw"] * w).sum())
    ratio = {k: (v / total_demand * 100.0 if total_demand > 0 else 0.0) for k, v in gen.items()}
    return gen, ratio, total_demand, battery_charge_mwh


def build_annual_summary(applied_costs, caps, cap_ratio, gen, gen_ratio):
    cost_map = {name: (val, unit) for name, _kind, val, unit in applied_costs}
    rows = []
    for name in RESOURCES:
        cost_val, cost_unit = cost_map.get(name, (np.nan, ""))
        rows.append({
            "電源/リソース": name,
            "適用コスト(代表値)": cost_val,
            "コスト単位": cost_unit,
            "最適容量[MW]": caps.get(name, 0.0),
            "容量比率[%]": cap_ratio.get(name, 0.0),
            "年間発電量[MWh]": gen.get(name, 0.0),
            "年間発電量比率[%]": gen_ratio.get(name, 0.0),
        })
    return pd.DataFrame(rows)


def build_monthly_summary(df):
    tmp = df.copy()
    tmp["month"] = tmp.index.month
    rows = []
    for m in range(1, 13):
        sub = tmp[tmp["month"] == m]
        w = sub["weight_h"]
        solar_mwh = float((sub["solar_mw"] * w).sum())
        wind_mwh = float((sub["wind_mw"] * w).sum())
        diesel_mwh = float((sub["diesel_mw"] * w).sum())
        batt_dis_mwh = float((sub["battery_discharge_mw"] * w).sum())
        batt_chg_mwh = float((sub["battery_charge_mw"] * w).sum())
        demand_mwh = float((sub["load_mw"] * w).sum())
        renewable_ratio = ((solar_mwh + wind_mwh) / demand_mwh * 100.0) if demand_mwh > 0 else 0.0
        rows.append({
            "月": MONTH_NAMES[m - 1],
            "Solar[MWh]": solar_mwh,
            "Wind[MWh]": wind_mwh,
            "Diesel[MWh]": diesel_mwh,
            "Battery放電[MWh]": batt_dis_mwh,
            "Battery充電[MWh]": batt_chg_mwh,
            "需要[MWh]": demand_mwh,
            "再エネ比率[%]": renewable_ratio,
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# インセットテキスト
# ----------------------------------------------------------------------------
def _capacity_mix_text(caps, cap_ratio, gen_ratio=None, header="Optimal capacity:"):
    lines = [header]
    for name in RESOURCES:
        if name in caps:
            lines.append(f"  {name}: {caps[name]:,.2f} MW ({cap_ratio[name]:.1f}%)")
    if gen_ratio is not None:
        lines.append("Annual generation mix:")
        for name in RESOURCES:
            if name in gen_ratio:
                lines.append(f"  {name}: {gen_ratio[name]:.1f}%")
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# 年間8,760時間 需給バランス/SOCグラフ
# ----------------------------------------------------------------------------
def plot_annual_dispatch(df, caps, cap_ratio, gen_ratio, out_path):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(16, 9), dpi=300,
        gridspec_kw={"height_ratios": [2.2, 1]}, sharex=True,
    )

    idx = df.index
    solar = df["solar_mw"].to_numpy()
    wind = df["wind_mw"].to_numpy()
    diesel = df["diesel_mw"].to_numpy()
    batt_dis = df["battery_discharge_mw"].to_numpy()
    batt_chg = -df["battery_charge_mw"].to_numpy()

    ax1.stackplot(
        idx, solar, wind, diesel, batt_dis,
        colors=[COLORS["solar"], COLORS["wind"], COLORS["diesel"], COLORS["battery_discharge"]],
        labels=["Solar", "Wind", "Diesel", "Battery Discharge"],
    )
    ax1.stackplot(idx, batt_chg, colors=[COLORS["battery_charge"]], labels=["Battery Charge"])
    ax1.plot(idx, df["load_mw"], color=COLORS["load"], lw=1.0, label="Load")
    ax1.axhline(0, color="black", lw=0.6)
    ax1.set_ylabel("Power [MW]")
    ax1.set_title("Annual Supply-Demand Balance (8,760 hours, Representative Case)")
    ax1.legend(loc="upper right", ncol=3, fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))

    ax1.text(
        0.01, 0.98, _capacity_mix_text(caps, cap_ratio, gen_ratio),
        transform=ax1.transAxes, va="top", ha="left", fontsize=8.5,
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.85),
    )

    soc = df["battery_soc_mwh"]
    ax2.fill_between(idx, 0, soc, color=COLORS["battery_discharge"], alpha=0.35)
    ax2.plot(idx, soc, color=COLORS["battery_charge"], lw=0.8)
    ax2.set_ylabel("Battery SOC [MWh]")
    ax2.set_xlabel("Month")
    ax2.grid(alpha=0.3)
    ax2.xaxis.set_major_locator(mdates.MonthLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b"))

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")


# ----------------------------------------------------------------------------
# 月別 需給バランス/SOCグラフ
# ----------------------------------------------------------------------------
def plot_monthly_dispatch(df, month, caps, cap_ratio, out_path):
    sub = df[df.index.month == month]
    idx = sub.index

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(14, 8), dpi=300,
        gridspec_kw={"height_ratios": [2.2, 1]}, sharex=True,
    )

    solar = sub["solar_mw"].to_numpy()
    wind = sub["wind_mw"].to_numpy()
    diesel = sub["diesel_mw"].to_numpy()
    batt_dis = sub["battery_discharge_mw"].to_numpy()
    batt_chg = -sub["battery_charge_mw"].to_numpy()

    ax1.stackplot(
        idx, solar, wind, diesel, batt_dis,
        colors=[COLORS["solar"], COLORS["wind"], COLORS["diesel"], COLORS["battery_discharge"]],
        labels=["Solar", "Wind", "Diesel", "Battery Discharge"],
    )
    ax1.stackplot(idx, batt_chg, colors=[COLORS["battery_charge"]], labels=["Battery Charge"])
    ax1.plot(idx, sub["load_mw"], color=COLORS["load"], lw=1.0, label="Load")
    ax1.axhline(0, color="black", lw=0.6)
    ax1.set_ylabel("Power [MW]")
    ax1.set_title(f"Monthly Supply-Demand Balance - {MONTH_NAMES[month - 1]}")
    ax1.legend(loc="upper right", ncol=3, fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d"))

    w = sub["weight_h"]
    demand_mwh = float((sub["load_mw"] * w).sum())
    monthly_gen = {
        "Solar": float((sub["solar_mw"] * w).sum()),
        "Wind": float((sub["wind_mw"] * w).sum()),
        "Battery": float((sub["battery_discharge_mw"] * w).sum()),
        "Diesel": float((sub["diesel_mw"] * w).sum()),
    }
    lines = [_capacity_mix_text(caps, cap_ratio, header="Optimal capacity (annual):"),
             "Monthly generation:"]
    for name in RESOURCES:
        val = monthly_gen.get(name, 0.0)
        pct = (val / demand_mwh * 100.0) if demand_mwh > 0 else 0.0
        lines.append(f"  {name}: {val:,.1f} MWh ({pct:.1f}%)")
    ax1.text(
        0.01, 0.98, "\n".join(lines), transform=ax1.transAxes, va="top", ha="left",
        fontsize=8, bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.85),
    )

    soc = sub["battery_soc_mwh"]
    ax2.fill_between(idx, 0, soc, color=COLORS["battery_discharge"], alpha=0.35)
    ax2.plot(idx, soc, color=COLORS["battery_charge"], lw=0.8)
    ax2.set_ylabel("Battery SOC [MWh]")
    ax2.set_xlabel("Day")
    ax2.grid(alpha=0.3)
    ax2.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d"))

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")


# ----------------------------------------------------------------------------
# メイン
# ----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="代表値(中央値)によるPyPSA 8760h最適化 -> 年間/月別ディスパッチグラフ出力"
    )
    p.add_argument("--config", default=str(CONFIG_PATH), help="network_config.xlsx のパス")
    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()
    output_dir, monthly_dir = make_output_dirs()

    log_path = os.path.join(output_dir, "execution.log")
    log_file = open(log_path, "w", encoding="utf-8")
    orig_stdout, orig_stderr = sys.stdout, sys.stderr
    sys.stdout = Tee(orig_stdout, log_file)
    sys.stderr = Tee(orig_stderr, log_file)

    try:
        print("=" * 78)
        print(" network_config.xlsx (代表値) -> PyPSA(8760h) -> 年間/月別ディスパッチ解析")
        print("=" * 78)
        print(f"[0] 出力先フォルダ: {output_dir}")
        print(f"    月別グラフ: {monthly_dir}")

        print(f"\n[1] ネットワーク設定を読み込み中: {args.config}")
        cfg = load_network_config(Path(args.config))

        print("[2] 代表値(中央値)ネットワークを構築中 ...")
        n, udf, params = build_representative_network(cfg)
        applied_costs = print_applied_costs(n)

        print("\n[3] PyPSA最適化(8760h)を実行中 ...")
        status, cond = run_optimization(n)
        print(f"\n    最適化ステータス: status={status}, condition={cond}")

        print("\n[4] 結果を集計中 ...")
        df = compute_dispatch_frame(n)
        caps, cap_ratio, total_cap = compute_capacity_table(n)
        gen, gen_ratio, total_demand, battery_charge_mwh = compute_annual_generation(df)

        print("\n--- 最適容量 [MW] / 容量比率 [%] ---")
        for name in RESOURCES:
            if name in caps:
                print(f"    {name:8s}: {caps[name]:10.2f} MW  ({cap_ratio[name]:5.1f}%)")
        print(f"    合計容量: {total_cap:10.2f} MW")

        print("\n--- 年間発電量 [MWh] / 年間発電量比率 [%] (対 年間総需要) ---")
        for name in RESOURCES:
            if name in gen:
                print(f"    {name:8s}: {gen[name]:12,.1f} MWh  ({gen_ratio[name]:5.1f}%)")
        print(f"    年間総需要:       {total_demand:12,.1f} MWh")
        print(f"    (参考)Battery充電量: {battery_charge_mwh:12,.1f} MWh")

        annual_summary = build_annual_summary(applied_costs, caps, cap_ratio, gen, gen_ratio)
        annual_summary_path = os.path.join(output_dir, "annual_summary_metrics.csv")
        annual_summary.to_csv(annual_summary_path, index=False, encoding="utf-8-sig")
        print(f"\n[5] 年間サマリー:\n{annual_summary.to_string(index=False)}")
        print(f"    保存: {annual_summary_path}")

        monthly_summary = build_monthly_summary(df)
        monthly_summary_path = os.path.join(output_dir, "monthly_summary_metrics.csv")
        monthly_summary.to_csv(monthly_summary_path, index=False, encoding="utf-8-sig")
        print(f"\n    月別サマリー:\n{monthly_summary.to_string(index=False)}")
        print(f"    保存: {monthly_summary_path}")

        dispatch_csv_path = os.path.join(output_dir, "annual_dispatch_data.csv")
        df.drop(columns=["weight_h"]).to_csv(dispatch_csv_path, index_label="timestamp",
                                              encoding="utf-8-sig")
        print(f"\n    保存: {dispatch_csv_path}")

        print("\n[6] 年間ディスパッチグラフを作成中 ...")
        annual_png_path = os.path.join(output_dir, "annual_dispatch_balance.png")
        plot_annual_dispatch(df, caps, cap_ratio, gen_ratio, annual_png_path)
        print(f"    保存: {annual_png_path}")

        print("\n[7] 月別ディスパッチグラフ(12枚)を作成中 ...")
        for m in range(1, 13):
            fname = f"{m:02d}_{MONTH_NAMES[m - 1]}.png"
            fpath = os.path.join(monthly_dir, fname)
            plot_monthly_dispatch(df, m, caps, cap_ratio, fpath)
            print(f"    保存: {fpath}")

        print("\n" + "=" * 78)
        print(f"完了 ({time.time() - t0:.0f}s)")
        print(f"成果物の保存先: {output_dir}")
        print("=" * 78)
    finally:
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        log_file.close()


if __name__ == "__main__":
    main()
