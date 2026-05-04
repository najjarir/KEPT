import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, to_hex
from matplotlib.lines import Line2D
import statsmodels.api as sm

# ================= SETTINGS =================
DATA_PATH = r"/mnt/data/Virulence Data_correct.xlsx"

VIR_PARAMS = ['Potential for Virulence', 'High potential for virulence', 'Potential for hypermucoidity', 'Potential for hypervirulence', 'Convergence of MDR and hypermucoidity', 'Convergence of MDR and hypervirulence']

CENTER_LINE_HEX = "#2f89f5"
DECREASE_STAR_HEX = "#ff0000"
WINDOW = 3
Z = 1.96
ALPHA = 0.05

FIGSIZE = (9.2, 10.5)
TITLE_FS = 10
AXIS_FS = 8

REGION_ORDER = ['Europe', 'Africa', 'Asia', 'Oceania', 'Latin America and the Caribbean', 'US / Canada']
INCOME_ORDER = ['High income', 'Upper middle income', 'Lower middle income', 'Low income']

OUT_PDF_REGION = r"/mnt/data/virulence_MA3_panel_3x2_VECTOR_region_pale_LOGITCLASSIC_TRENDSTARS.pdf"
OUT_PDF_INCOME = r"/mnt/data/virulence_MA3_panel_3x2_VECTOR_income_pale_LOGITCLASSIC_TRENDSTARS.pdf"

PALE_MIX_WITH_WHITE = 0.65

# ================= HELPERS =================
cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", list(plt.cm.tab10.colors))

def pale_color(hex_color: str, mix_with_white: float = PALE_MIX_WITH_WHITE) -> str:
    r, g, b = to_rgb(hex_color)
    return to_hex((
        r*(1-mix_with_white) + mix_with_white,
        g*(1-mix_with_white) + mix_with_white,
        b*(1-mix_with_white) + mix_with_white
    ))

def build_pale_map(order_list):
    orig_hex = [to_hex(to_rgb(cycle[i % len(cycle)])) for i in range(len(order_list))]
    pale_hex = [pale_color(h, PALE_MIX_WITH_WHITE) for h in orig_hex]
    return {k: v for k, v in zip(order_list, pale_hex)}

def ensure_fallback(groups, cmap):
    for i, g in enumerate(groups):
        if g not in cmap:
            base = to_hex(to_rgb(cycle[i % len(cycle)]))
            cmap[g] = pale_color(base, PALE_MIX_WITH_WHITE)
    return cmap

def centered_ma_ci(series: pd.Series):
    series = pd.to_numeric(series, errors="coerce")
    ma = series.rolling(WINDOW, center=True, min_periods=2).mean()
    sd = series.rolling(WINDOW, center=True, min_periods=2).std(ddof=1)
    n = series.rolling(WINDOW, center=True, min_periods=2).count()
    se = sd / np.sqrt(n)
    lo = ma - Z * se
    hi = ma + Z * se
    return ma.astype(float), lo.astype(float), hi.astype(float)

def last_valid_xy(x_vals, y_vals):
    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if not m.any():
        return None
    idx = np.where(m)[0][-1]
    return float(x[idx]), float(y[idx])

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

def add_trend_star(ax, ma_series: pd.Series, beta, pval):
    if beta is None or pval is None:
        return False
    try:
        p = float(pval); b = float(beta)
    except Exception:
        return False
    if p >= ALPHA:
        return False
    xy = last_valid_xy(ma_series.index.values, ma_series.values)
    if xy is None:
        return False
    x_last, y_last = xy
    color = CENTER_LINE_HEX if b > 0 else DECREASE_STAR_HEX
    ax.scatter([x_last], [y_last], marker="*", s=180, color=color, edgecolors="none", zorder=10)
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

def load_df(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int)

    if "region" not in df.columns:
        raise ValueError("Expected a 'region' column.")
    df["region"] = df["region"].astype(str).str.strip()

    if "income_group" in df.columns:
        df["income_group"] = df["income_group"].astype(str).str.strip()
    elif "income" in df.columns:
        df["income"] = df["income"].astype(str).str.strip()
    else:
        raise ValueError("Income column not found (expected 'income_group' or 'income').")

    for c in VIR_PARAMS:
        if c not in df.columns:
            raise ValueError(f"Missing column: {c}")
        s = pd.to_numeric(df[c], errors="coerce")
        uniq = set(s.dropna().unique().tolist())
        df[c] = s.astype(float) if uniq.issubset({0, 1}) else (s > 0).astype(int)

    return df

def plot_3x2(data: pd.DataFrame, group_col: str, group_order: list[str], out_pdf: str, legend_title: str, legend_ncol: int, star_anchor_y: float):
    present = set([g for g in data[group_col].dropna().unique().tolist() if str(g).strip() != ""])
    groups = [g for g in group_order if g in present] + sorted(list(present - set(group_order)))
    cmap = ensure_fallback(groups, build_pale_map(group_order))

    fig, axes = plt.subplots(3, 2, figsize=FIGSIZE)
    axes = axes.flatten()
    has_star = False

    trends = {var: trend_logit_individual(data, var) for var in VIR_PARAMS}

    for i, var in enumerate(VIR_PARAMS):
        ax = axes[i]

        for g in groups:
            sub = data[data[group_col] == g][["year", var]].dropna()
            if sub.empty:
                continue
            agg = sub.groupby("year")[var].agg(["sum", "count"])
            agg = agg[agg["count"] > 0]
            if agg.empty:
                continue
            prop_g = (agg["sum"] / agg["count"]) * 100.0
            ax.plot(prop_g.index.astype(float), prop_g.values.astype(float),
                    color=cmap[g], linewidth=1.2, alpha=0.95, zorder=1)

        yearly = data.groupby("year")[var].agg(["sum", "count"])
        yearly = yearly[yearly["count"] > 0]
        yearly["prop"] = (yearly["sum"] / yearly["count"]) * 100.0

        ax.plot(yearly.index.astype(float), yearly["prop"].astype(float),
                marker="o", linestyle="none", alpha=0.6, zorder=3)

        ma, lo, hi = centered_ma_ci(yearly["prop"])
        ax.fill_between(ma.index.astype(float), lo.values, hi.values, alpha=0.2, zorder=2)
        ax.plot(ma.index.astype(float), ma.values, linewidth=2.5, color=CENTER_LINE_HEX, zorder=4)

        beta, pval = trends[var]
        has_star = add_trend_star(ax, ma, beta, pval) or has_star

        ax.set_title(var, fontsize=TITLE_FS)
        ax.set_xlabel("Year", fontsize=AXIS_FS)
        ax.set_ylabel("Proportion (%)", fontsize=AXIS_FS)
        ax.set_ylim(0, 100)

    handles = [Line2D([0], [0], color=cmap[g], lw=3, label=g) for g in groups]
    fig.legend(handles=handles, loc="lower center", ncol=legend_ncol, frameon=False,
               fontsize=9, title=legend_title, title_fontsize=9)

    if has_star:
        add_star_legend(fig, anchor_y=star_anchor_y, fontsize=9)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(out_pdf)
    plt.close(fig)

def main():
    df = load_df(DATA_PATH)
    income_col = "income_group" if "income_group" in df.columns else "income"
    plot_3x2(df, "region", REGION_ORDER, OUT_PDF_REGION, "Regions", 3, star_anchor_y=0.13)
    plot_3x2(df, income_col, INCOME_ORDER, OUT_PDF_INCOME, "Income groups", 2, star_anchor_y=0.13)
    print("Saved:", OUT_PDF_REGION)
    print("Saved:", OUT_PDF_INCOME)

if __name__ == "__main__":
    main()
