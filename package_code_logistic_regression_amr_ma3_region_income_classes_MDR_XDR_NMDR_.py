#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AMR MA(3) figures (vector PDF) with pale group lines AND trend stars based on classic logistic regression.

Stars are placed at the LAST VALID POINT of the central MA(3) line:
  - Blue (#2f89f5): Significant global temporal increase (p < 0.05)
  - Red  (#ff0000): Significant global temporal decrease (p < 0.05)

Trend tests:
  - Binary outcomes (12 classes, MDR, XDR): classic individual-level logistic regression
      GLM(Binomial, logit): logit(P(Y=1)) = beta0 + beta1*year
  - N_MDR (continuous): WLS on yearly means (weights = yearly isolate count)

This script loads each Excel file ONCE and precomputes group-year proportions for speed.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, to_hex
from matplotlib.lines import Line2D
import statsmodels.api as sm

# ================= USER SETTINGS =================
DATA_12 = "Final_Data_for_graphs_and_analysis.xlsx"
DATA_MDRXDR = "3 MDR and XDR.xlsx"
DATA_NMDR = "N_mdr data.xlsx"

CLASSES_12 = [
    'Aminoglycosides','Carbapenems','Third-generation cephalosporins','Polymyxins',
    'Fosfomycins','Fluoroquinolones','Macrolides','Phenicols','Rifampicin',
    'Trimethoprim + Sulfamethoxazoles','Tetracyclines','Tigecycline'
]
CLASSES_MDRXDR = ["MDR","XDR"]

CENTER_LINE_HEX = "#2f89f5"
DECREASE_STAR_HEX = "#ff0000"
WINDOW = 3
Z = 1.96
ALPHA = 0.05
PALE_MIX_WITH_WHITE = 0.65

REGION_ORDER = ["Europe","Africa","Asia","Oceania","Latin America and the Caribbean","US / Canada"]
INCOME_ORDER = ["High income","Upper middle income","Lower middle income","Low income"]

FIGSIZE_4x3 = (18, 10.5)
FIGSIZE_1x2 = (9, 3.8)
FIGSIZE_NMDR = (7, 5.5)
TITLE_FS = 10
AXIS_FS = 8

OUT_4x3_REGION = "global_MA3_panel_4x3_VECTOR_region_pale_package_LOGITCLASSIC_TRENDSTARS.pdf"
OUT_4x3_INCOME = "global_MA3_panel_4x3_VECTOR_income_pale_package_LOGITCLASSIC_TRENDSTARS.pdf"
OUT_1x2_REGION = "MDR_XDR_panel_1x2_VECTOR_region_pale_package_LOGITCLASSIC_TRENDSTARS.pdf"
OUT_1x2_INCOME = "MDR_XDR_panel_1x2_VECTOR_income_pale_package_LOGITCLASSIC_TRENDSTARS.pdf"
OUT_NMDR_REGION = "N_MDR_global_MA3_pale_region_package_LOGITCLASSIC_TRENDSTARS.pdf"
OUT_NMDR_INCOME = "N_MDR_global_MA3_pale_income_package_LOGITCLASSIC_TRENDSTARS.pdf"

# ================= HELPERS =================
cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', list(plt.cm.tab10.colors))

def pale_color(hex_color: str, mix_with_white: float = PALE_MIX_WITH_WHITE) -> str:
    r, g, b = to_rgb(hex_color)
    return to_hex((r*(1-mix_with_white)+mix_with_white,
                   g*(1-mix_with_white)+mix_with_white,
                   b*(1-mix_with_white)+mix_with_white))

def build_pale_map(order_list):
    base_hex = [to_hex(to_rgb(cycle[i % len(cycle)])) for i in range(len(order_list))]
    return {k: pale_color(h) for k, h in zip(order_list, base_hex)}

def ensure_fallback(groups, cmap):
    for i, g in enumerate(groups):
        if g not in cmap:
            cmap[g] = pale_color(to_hex(to_rgb(cycle[i % len(cycle)])))
    return cmap

def centered_ma_ci(series: pd.Series):
    s = pd.to_numeric(series, errors="coerce")
    ma = s.rolling(WINDOW, center=True, min_periods=2).mean()
    sd = s.rolling(WINDOW, center=True, min_periods=2).std(ddof=1)
    n = s.rolling(WINDOW, center=True, min_periods=2).count()
    se = sd / np.sqrt(n)
    lo = ma - Z*se
    hi = ma + Z*se
    return ma.astype(float), lo.astype(float), hi.astype(float)

def last_valid_xy(x_vals, y_vals):
    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if not m.any():
        return None
    j = np.where(m)[0][-1]
    return float(x[j]), float(y[j])

def trend_logit_individual(df: pd.DataFrame, ycol: str):
    d = df[["year", ycol]].dropna().copy()
    if d.shape[0] < 10:
        return None, None
    y = pd.to_numeric(d[ycol], errors="coerce")
    d = d.loc[y.notna()].copy()
    y = y.loc[y.notna()].astype(int)
    if y.nunique() < 2:
        return None, None
    X = sm.add_constant(d["year"].astype(float), has_constant="add")
    try:
        m = sm.GLM(y.astype(float), X, family=sm.families.Binomial()).fit(disp=0)
        return float(m.params["year"]), float(m.pvalues["year"])
    except Exception:
        return None, None

def trend_wls_yearmeans(df: pd.DataFrame, ycol: str):
    d = df[["year", ycol]].dropna().copy()
    if d.shape[0] < 10:
        return None, None
    yearly = d.groupby("year")[ycol].agg(["mean","count"]).reset_index()
    yearly = yearly[yearly["count"] > 0]
    if yearly.shape[0] < 3:
        return None, None
    X = sm.add_constant(yearly["year"].astype(float), has_constant="add")
    try:
        m = sm.WLS(yearly["mean"].astype(float), X, weights=yearly["count"].astype(float)).fit()
        return float(m.params["year"]), float(m.pvalues["year"])
    except Exception:
        return None, None

def add_trend_star(ax, ma_series, beta, pval):
    if beta is None or pval is None:
        return False
    try:
        b = float(beta); p = float(pval)
    except Exception:
        return False
    if p >= ALPHA:
        return False
    xy = last_valid_xy(ma_series.index.values, ma_series.values)
    if xy is None:
        return False
    x_last, y_last = xy
    col = CENTER_LINE_HEX if b > 0 else DECREASE_STAR_HEX
    ax.scatter([x_last], [y_last], marker="*", s=180, color=col, edgecolors="none", zorder=10)
    return True

def add_star_legend(fig, anchor_y: float, fontsize: int):
    inc = Line2D([0], [0], marker="*", linestyle="none", markersize=11,
                 markerfacecolor=CENTER_LINE_HEX, markeredgecolor="none",
                 label="Significant global temporal increase (p < 0.05)")
    dec = Line2D([0], [0], marker="*", linestyle="none", markersize=11,
                 markerfacecolor=DECREASE_STAR_HEX, markeredgecolor="none",
                 label="Significant global temporal decrease (p < 0.05)")
    fig.legend(handles=[inc, dec], loc="lower right", bbox_to_anchor=(0.98, anchor_y),
               frameon=False, fontsize=fontsize)

def make_group_legend(fig, groups, cmap, title, ncol, fontsize, title_fontsize):
    handles = [Line2D([0],[0], color=cmap[g], lw=3, label=g) for g in groups]
    fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False,
               fontsize=fontsize, title=title, title_fontsize=title_fontsize)

def precompute_props(df, group_col, cls_list):
    years = np.sort(df["year"].unique())
    out = {}
    for cls in cls_list:
        agg = df.groupby([group_col,"year"])[cls].agg(["sum","count"]).reset_index()
        agg = agg[agg["count"]>0]
        agg["prop"] = (agg["sum"]/agg["count"])*100.0
        out[cls] = agg.pivot(index="year", columns=group_col, values="prop").reindex(years)
    return out

# ================= LOAD DATA (ONCE) =================
def load_12(path: str) -> pd.DataFrame:
    usecols12 = ["year","region","income_group"] + CLASSES_12
    df = pd.read_excel(path, usecols=lambda c: c in set(usecols12))
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int)
    df["region"] = df["region"].astype(str).str.strip()
    df["income_group"] = df["income_group"].astype(str).str.strip()
    for c in CLASSES_12:
        df[c] = df[c].astype(bool).astype(int)
    return df

def load_mdrxdr(path: str) -> pd.DataFrame:
    usecolsm = ["year","region","income_group"] + CLASSES_MDRXDR
    df = pd.read_excel(path, usecols=lambda c: c in set(usecolsm))
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int)
    df["region"] = df["region"].astype(str).str.strip()
    df["income_group"] = df["income_group"].astype(str).str.strip()
    for c in CLASSES_MDRXDR:
        df[c] = df[c].astype(bool).astype(int)
    return df

def load_nmdr(path: str):
    df = pd.read_excel(path)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int)
    df["region"] = df["region"].astype(str).str.strip()
    if "income_group" in df.columns:
        inc_col = "income_group"
    elif "income" in df.columns:
        inc_col = "income"
    else:
        raise ValueError("No income column found in N_MDR file.")
    df[inc_col] = df[inc_col].astype(str).str.strip()
    nmdr_col = None
    for c in df.columns:
        if c.lower() in ["n_mdr","nmdr"]:
            nmdr_col = c
            break
    if nmdr_col is None:
        raise ValueError("N_MDR column not found.")
    df[nmdr_col] = pd.to_numeric(df[nmdr_col], errors="coerce")
    df = df.dropna(subset=[nmdr_col]).copy()
    return df, inc_col, nmdr_col

# ================= PLOTS =================
def plot_4x3(df, group_col, order, out_pdf, legend_title, legend_ncol, star_anchor_y=0.13):
    present = set([g for g in df[group_col].dropna().unique().tolist() if str(g).strip()!=""])
    groups = [g for g in order if g in present] + sorted(list(present - set(order)))
    cmap = ensure_fallback(groups, build_pale_map(order))
    props = precompute_props(df, group_col, CLASSES_12)
    yearly_all = {cls: df.groupby("year")[cls].agg(["sum","count"]).assign(prop=lambda x:(x["sum"]/x["count"])*100.0)
                  for cls in CLASSES_12}
    trend = {cls: trend_logit_individual(df, cls) for cls in CLASSES_12}

    fig, axes = plt.subplots(3,4, figsize=FIGSIZE_4x3)
    axes = axes.flatten()
    has_star = False

    for i, cls in enumerate(CLASSES_12):
        ax = axes[i]
        piv = props[cls]
        for g in groups:
            if g not in piv.columns:
                continue
            y = piv[g].dropna()
            if y.empty:
                continue
            ax.plot(y.index.astype(float), y.values.astype(float), color=cmap[g], linewidth=1.2, alpha=0.95, zorder=1)

        ya = yearly_all[cls]
        ax.plot(ya.index.astype(float), ya["prop"].astype(float), marker="o", linestyle="none", alpha=0.6, zorder=3)

        ma, lo, hi = centered_ma_ci(ya["prop"])
        ax.fill_between(ma.index.astype(float), lo.values, hi.values, alpha=0.2, zorder=2)
        ax.plot(ma.index.astype(float), ma.values, linewidth=2.5, color=CENTER_LINE_HEX, zorder=4)

        beta, pval = trend[cls]
        has_star = add_trend_star(ax, ma, beta, pval) or has_star

        ax.set_title(cls, fontsize=TITLE_FS)
        ax.set_xlabel("Year", fontsize=AXIS_FS)
        ax.set_ylabel("Resistant proportion (%)", fontsize=AXIS_FS)
        ax.set_ylim(0,100)

    make_group_legend(fig, groups, cmap, legend_title, legend_ncol, fontsize=9, title_fontsize=9)
    if has_star:
        add_star_legend(fig, anchor_y=star_anchor_y, fontsize=9)
    plt.tight_layout(rect=[0,0.08,1,1])
    fig.savefig(out_pdf)
    plt.close(fig)

def plot_1x2(df, group_col, order, out_pdf, legend_title, legend_ncol, star_anchor_y=0.28):
    present = set([g for g in df[group_col].dropna().unique().tolist() if str(g).strip()!=""])
    groups = [g for g in order if g in present] + sorted(list(present - set(order)))
    cmap = ensure_fallback(groups, build_pale_map(order))
    props = precompute_props(df, group_col, CLASSES_MDRXDR)
    yearly_all = {cls: df.groupby("year")[cls].agg(["sum","count"]).assign(prop=lambda x:(x["sum"]/x["count"])*100.0)
                  for cls in CLASSES_MDRXDR}
    trend = {cls: trend_logit_individual(df, cls) for cls in CLASSES_MDRXDR}

    fig, axes = plt.subplots(1,2, figsize=FIGSIZE_1x2)
    has_star = False

    for ax, cls in zip(axes, CLASSES_MDRXDR):
        piv = props[cls]
        for g in groups:
            if g not in piv.columns:
                continue
            y = piv[g].dropna()
            if y.empty:
                continue
            ax.plot(y.index.astype(float), y.values.astype(float), color=cmap[g], linewidth=1.2, alpha=0.95, zorder=1)

        ya = yearly_all[cls]
        ax.plot(ya.index.astype(float), ya["prop"].astype(float), marker="o", linestyle="none", alpha=0.6, zorder=3)

        ma, lo, hi = centered_ma_ci(ya["prop"])
        ax.fill_between(ma.index.astype(float), lo.values, hi.values, alpha=0.2, zorder=2)
        ax.plot(ma.index.astype(float), ma.values, linewidth=2.5, color=CENTER_LINE_HEX, zorder=4)

        beta, pval = trend[cls]
        has_star = add_trend_star(ax, ma, beta, pval) or has_star

        ax.set_title(cls, fontsize=TITLE_FS)
        ax.set_xlabel("Year", fontsize=AXIS_FS)
        ax.set_ylabel("Resistant proportion (%)", fontsize=AXIS_FS)
        ax.set_ylim(0,100)

    make_group_legend(fig, groups, cmap, legend_title, legend_ncol, fontsize=8, title_fontsize=8)
    if has_star:
        add_star_legend(fig, anchor_y=star_anchor_y, fontsize=8)
    plt.tight_layout(rect=[0,0.18,1,1])
    fig.savefig(out_pdf)
    plt.close(fig)

def plot_nmdr(df, group_col, order, nmdr_col, out_pdf, legend_title, legend_ncol, star_anchor_y=0.35):
    present = set([g for g in df[group_col].dropna().unique().tolist() if str(g).strip()!=""])
    groups = [g for g in order if g in present] + sorted(list(present - set(order)))
    cmap = ensure_fallback(groups, build_pale_map(order))

    fig, ax = plt.subplots(figsize=FIGSIZE_NMDR)
    has_star = False

    for g in groups:
        sub = df[df[group_col]==g][["year", nmdr_col]].dropna()
        if sub.empty:
            continue
        yr = sub.groupby("year")[nmdr_col].mean()
        ax.plot(yr.index.astype(float), yr.values.astype(float), color=cmap[g], linewidth=1.2, alpha=0.95, zorder=1)

    yearly = df.groupby("year")[nmdr_col].mean()
    ax.plot(yearly.index.astype(float), yearly.values.astype(float), marker="o", linestyle="none", alpha=0.6, zorder=3)

    ma, lo, hi = centered_ma_ci(yearly)
    ax.fill_between(ma.index.astype(float), lo.values, hi.values, alpha=0.2, zorder=2)
    ax.plot(ma.index.astype(float), ma.values, linewidth=2.5, color=CENTER_LINE_HEX, zorder=4)

    beta, pval = trend_wls_yearmeans(df, nmdr_col)
    has_star = add_trend_star(ax, ma, beta, pval) or has_star

    ax.set_title("N_MDR (Global mean)", fontsize=TITLE_FS)
    ax.set_xlabel("Year", fontsize=AXIS_FS)
    ax.set_ylabel("Mean N_MDR", fontsize=AXIS_FS)
    ax.set_ylim(0,5)
    ax.set_yticks(np.arange(0,6,1))

    make_group_legend(fig, groups, cmap, legend_title, legend_ncol, fontsize=8, title_fontsize=8)
    if has_star:
        add_star_legend(fig, anchor_y=star_anchor_y, fontsize=8)
    plt.tight_layout(rect=[0,0.25,1,1])
    fig.savefig(out_pdf)
    plt.close(fig)

def main():
    df12 = load_12(DATA_12)
    dfm = load_mdrxdr(DATA_MDRXDR)
    dfn, inc_col_n, nmdr_col = load_nmdr(DATA_NMDR)

    plot_4x3(df12, "region", REGION_ORDER, OUT_4x3_REGION, "Regions", 3)
    plot_4x3(df12, "income_group", INCOME_ORDER, OUT_4x3_INCOME, "Income groups", 2)

    plot_1x2(dfm, "region", REGION_ORDER, OUT_1x2_REGION, "Regions", 3)
    plot_1x2(dfm, "income_group", INCOME_ORDER, OUT_1x2_INCOME, "Income groups", 2)

    plot_nmdr(dfn, "region", REGION_ORDER, nmdr_col, OUT_NMDR_REGION, "Regions", 3)
    plot_nmdr(dfn, inc_col_n, INCOME_ORDER, nmdr_col, OUT_NMDR_INCOME, "Income groups", 2)

    print("Saved PDFs:")
    for p in [OUT_4x3_REGION, OUT_4x3_INCOME, OUT_1x2_REGION, OUT_1x2_INCOME, OUT_NMDR_REGION, OUT_NMDR_INCOME]:
        print(" -", p)

if __name__ == "__main__":
    main()
