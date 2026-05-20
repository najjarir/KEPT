import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "axes.linewidth":    0.8,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size":  3,
    "ytick.major.size":  3,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
})

OUTPUT_DIR = "Graficos_Distribuicao"
os.makedirs(OUTPUT_DIR, exist_ok=True)

CMAP = plt.get_cmap("tab20")
BASE_COLORS = CMAP.colors

# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------
df = pd.read_csv("C:\Users\mathe\OneDrive\Área de Trabalho\TESTE\tabela_genes_limpa.csv")

# ---------------------------------------------------------------------------
# Figure 1 — sample distribution by geographic region (pie chart)
# ---------------------------------------------------------------------------
region_pct = df["region"].dropna().value_counts(normalize=True) * 100

fig, ax = plt.subplots(figsize=(8, 8))
wedges, _ = ax.pie(
    region_pct,
    colors=BASE_COLORS[: len(region_pct)],
    startangle=140,
)
ax.legend(
    wedges,
    [f"{idx} ({val:.1f}%)" for idx, val in region_pct.items()],
    title="Region",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    fontsize=12,
    title_fontsize=13,
)
ax.set_title("")
plt.savefig(f"{OUTPUT_DIR}/Distribuicao_Regiao.pdf", dpi=600, bbox_inches="tight")
plt.show()
plt.close()

# ---------------------------------------------------------------------------
# Figure 2 — sample distribution by income group (pie chart)
# ---------------------------------------------------------------------------
income_pct = df["income"].dropna().value_counts(normalize=True) * 100

fig, ax = plt.subplots(figsize=(8, 8))
wedges, _ = ax.pie(
    income_pct,
    colors=BASE_COLORS[: len(income_pct)],
    startangle=140,
)
ax.legend(
    wedges,
    [f"{idx} ({val:.1f}%)" for idx, val in income_pct.items()],
    title="Income Group",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    fontsize=12,
    title_fontsize=13,
)
ax.set_title("")
plt.savefig(f"{OUTPUT_DIR}/Distribuicao_Income.pdf", dpi=600, bbox_inches="tight")
plt.show()
plt.close()

# ---------------------------------------------------------------------------
# Figure 3 — sample distribution by year, simplified layout
# ---------------------------------------------------------------------------
year_pct = df["year"].dropna().value_counts(normalize=True).sort_index() * 100

fig, ax = plt.subplots(figsize=(12, 6))
bar_objects = ax.bar(
    year_pct.index.astype(int).astype(str),
    year_pct.values,
    color=BASE_COLORS[: len(year_pct)],
    edgecolor="black",
    linewidth=0.5,
)
for bar in bar_objects:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        h + 0.5,
        f"{h:.1f}%",
        ha="center",
        va="bottom",
        fontsize=7,
    )
ax.set_xlabel("Year", fontsize=14, labelpad=10)
ax.set_ylabel("Percentage of Samples (%)", fontsize=14, labelpad=10)
ax.set_title("")
ax.set_ylim(0, year_pct.max() * 1.15)
ax.tick_params(axis="x", rotation=45, labelsize=12)
ax.tick_params(axis="y", labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.savefig(f"{OUTPUT_DIR}/Distribuicao_Periodo.pdf", dpi=600, bbox_inches="tight")
plt.show()
plt.close()

# ---------------------------------------------------------------------------
# Figure 4 — sample distribution by year, publication layout
#             (absolute counts above bars, single fill colour)
# ---------------------------------------------------------------------------
BAR_COLOR = "#2166AC"

year_pct2  = df["year"].dropna().value_counts(normalize=True).sort_index() * 100
year_abs   = df["year"].dropna().value_counts().sort_index()
years      = year_pct2.index.astype(int).astype(str)
pct_values = year_pct2.values
abs_values = year_abs.values

fig, ax = plt.subplots(figsize=(10, 5), dpi=300, facecolor="white")
fig.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.15)

bars = ax.bar(
    years,
    pct_values,
    color=BAR_COLOR,
    edgecolor="white",
    linewidth=0.4,
    width=0.68,
    zorder=3,
)
for bar, n in zip(bars, abs_values):
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        h + pct_values.max() * 0.013,
        f"{n:,}",
        ha="center",
        va="bottom",
        fontsize=5.2,
        color="#222222",
        zorder=4,
    )

ax.set_xlabel("Year", fontsize=8, labelpad=6)
ax.set_ylabel("Percentage of Samples (%)", fontsize=8, labelpad=6)
ax.set_ylim(0, pct_values.max() * 1.20)
ax.tick_params(axis="x", rotation=45, labelsize=6.5)
ax.tick_params(axis="y", labelsize=6.5)
ax.yaxis.grid(True, linewidth=0.4, color="#cccccc", linestyle="--", zorder=0)
ax.set_axisbelow(True)

plt.savefig(
    f"{OUTPUT_DIR}/barras_por_ano.pdf",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
)
plt.show()
plt.close()
