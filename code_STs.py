import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import math
from matplotlib.gridspec import GridSpec

# =========================
# PATHS
# =========================

input_file = r"D:\Python\AMR_Iris\Input\ST_data.xlsx"
output_folder = r"D:\Python\AMR_Iris\Output"

os.makedirs(output_folder, exist_ok=True)

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.size"] = 11

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# =========================
# LOAD DATA
# =========================

df = pd.read_excel(input_file)
df["year"] = df["year"].astype(int)

regions = sorted(df["region"].dropna().unique())

# =========================
# FUNCTION TO IDENTIFY HIGHLIGHT ST
# =========================

def identify_highlight_ST(data):

    yearly_counts = data.groupby(["year", "ST"]).size().reset_index(name="n")
    yearly_total = data.groupby("year").size().reset_index(name="total")

    merged = yearly_counts.merge(yearly_total, on="year")
    merged["rel_prev"] = merged["n"] / merged["total"]

    pivot = merged.pivot(index="year", columns="ST", values="rel_prev").fillna(0)

    top5 = data["ST"].value_counts().head(5).index.tolist()

    recent = pivot.loc[(pivot.index >= 2013) & (pivot.index <= 2022)]

    rule_ST = []

    for st in recent.columns:
        if (recent[st] > 0.05).sum() >= 2:
            rule_ST.append(st)

    return sorted(list(set(top5 + rule_ST)))

# =========================
# IDENTIFY ALL ST TO COLOR
# =========================

all_highlight_ST = set()
all_highlight_ST.update(identify_highlight_ST(df))

for region in regions:
    subset = df[df["region"] == region]
    all_highlight_ST.update(identify_highlight_ST(subset))

all_highlight_ST = sorted(list(all_highlight_ST))

# =========================
# DISTINCT COLOR PALETTE
# =========================

distinct_colors = [
    "#ff7f0e", "#1f77b4", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22",
    "#7570b3", "#d95f02", "#1b9e77", "#e7298a", "#66a61e",
    "#e6ab02", "#a6761d", "#a6cee3", "#fb9a99",
    "#b2df8a", "#fdbf6f", "#cab2d6", "#ffff99", "#b15928",
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
    "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
    "#ccebc5", "#ffed6f", "#3b528b", "#5ec962", "#440154"
]

color_map = {}
for i, st in enumerate(all_highlight_ST):
    color_map[st] = distinct_colors[i % len(distinct_colors)]

if "ST323" in all_highlight_ST:
    color_map["ST323"] = "#e31a1c"

# =========================
# FUNCTION TO PREPARE PANEL DATA
# =========================

def prepare_panel_data(data):

    yearly_counts = data.groupby(["year", "ST"]).size().reset_index(name="n")
    yearly_total = data.groupby("year").size().reset_index(name="total")

    merged = yearly_counts.merge(yearly_total, on="year")
    merged["rel_prev"] = merged["n"] / merged["total"]

    pivot = merged.pivot(index="year", columns="ST", values="rel_prev").fillna(0)

    highlight_ST = identify_highlight_ST(data)

    data_year_min = int(pivot.index.min())
    data_year_max = int(pivot.index.max())

    left_limit = min(2000, data_year_min)
    right_limit = max(2022, data_year_max)

    return {
        "pivot": pivot,
        "highlight_ST": highlight_ST,
        "left_limit": left_limit,
        "right_limit": right_limit
    }

# =========================
# FUNCTION TO DRAW EACH PANEL
# =========================

def draw_panel(ax, panel_data, title=None, title_size=12):

    pivot = panel_data["pivot"]
    highlight_ST = panel_data["highlight_ST"]
    left_limit = panel_data["left_limit"]
    right_limit = panel_data["right_limit"]

    # =========================
    # 1. PLOT OTHERS (FUNDO)
    # =========================
    for st in pivot.columns:
        if st not in highlight_ST:
            ax.plot(
                pivot.index,
                pivot[st] * 100,
                color="lightgrey",
                linewidth=1,
                alpha=0.6,
                zorder=1   # <- fundo
            )

    # =========================
    # 2. PLOT HIGHLIGHT (FRENTE)
    # =========================
    for st in pivot.columns:
        if st in highlight_ST:
            ax.plot(
                pivot.index,
                pivot[st] * 100,
                linewidth=2.2,
                color=color_map[st],
                zorder=3   # <- frente
            )

    # =========================
    # AXES
    # =========================
    ax.set_xticks([2000, 2005, 2010, 2015, 2020])
    ax.set_xticklabels([2000, 2005, 2010, 2015, 2020])
    ax.set_xlim(left_limit - 1.0, right_limit + 1.0)

    ax.set_xlabel("Year")
    ax.set_ylabel("Relative prevalence (%)")

    if title is not None:
        ax.set_title(title, fontsize=title_size, pad=6)

# =========================
# PREPARE ALL PANELS
# =========================

panel_info = {}
panel_info["Global"] = prepare_panel_data(df)

for region in regions:
    subset = df[df["region"] == region]
    panel_info[region] = prepare_panel_data(subset)

# =========================
# FIGURE LAYOUT
# =========================

fig = plt.figure(figsize=(17, 28))
gs = GridSpec(
    nrows=5,
    ncols=2,
    height_ratios=[1.75, 1.15, 1.15, 1.15, 1.15],
    hspace=0.28,
    wspace=0.18,
    figure=fig
)

# =========================
# GLOBAL PANEL
# =========================

ax_global = fig.add_subplot(gs[0, :])
draw_panel(
    ax_global,
    panel_info["Global"],
    title="Global",
    title_size=14
)

# =========================
# REGIONAL PANELS
# =========================

region_positions = [
    (1, 0), (1, 1),
    (2, 0), (2, 1),
    (3, 0), (3, 1)
]

for i, region in enumerate(regions[:6]):
    r, c = region_positions[i]
    ax = fig.add_subplot(gs[r, c])
    draw_panel(
        ax,
        panel_info[region],
        title=region,
        title_size=12
    )

for j in range(len(regions[:6]), 6):
    r, c = region_positions[j]
    ax_empty = fig.add_subplot(gs[r, c])
    ax_empty.axis("off")

# =========================
# BOTTOM LEGEND AREA
# =========================

legend_ax = fig.add_subplot(gs[4, :])
legend_ax.set_xlim(0, 1)
legend_ax.set_ylim(0, 1)
legend_ax.axis("off")

legend_entries = [("Global", panel_info["Global"]["highlight_ST"])]
for region in regions[:6]:
    legend_entries.append((region, panel_info[region]["highlight_ST"]))

display_names = {
    "Global": "Global",
    "Africa": "Africa",
    "Asia": "Asia",
    "Europe": "Europe",
    "Latin America and the Caribbean": "Latin America and\nthe Caribbean",
    "North America": "North America",
    "Oceania": "Oceania"
}

# one single row of legend blocks
n_blocks = len(legend_entries)
left_margin = 0.015
right_margin = 0.015
usable_width = 1 - left_margin - right_margin
block_width = usable_width / n_blocks

for i, (panel_name, st_list) in enumerate(legend_entries):
    x0 = left_margin + i * block_width
    y0 = 1.10

    title_text = display_names.get(panel_name, panel_name)

    legend_ax.text(
        x0, y0,
        title_text,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        transform=legend_ax.transAxes
    )

    # more spacing between title and ST list
    n_title_lines = title_text.count("\n") + 1
    y = y0 - (0.09 * n_title_lines)

    line_x1 = x0
    line_x2 = x0 + 0.018
    text_x = x0 + 0.025

    for st in st_list:
        legend_ax.plot(
            [line_x1, line_x2],
            [y, y],
            transform=legend_ax.transAxes,
            color=color_map[st],
            linewidth=2.0,
            solid_capstyle="round",
            clip_on=False
        )
        legend_ax.text(
            text_x, y,
            st,
            ha="left",
            va="center",
            fontsize=10,
            transform=legend_ax.transAxes
        )
        y -= 0.055

    legend_ax.plot(
        [line_x1, line_x2],
        [y, y],
        transform=legend_ax.transAxes,
        color="lightgrey",
        linewidth=2.0,
        solid_capstyle="round",
        clip_on=False
    )
    legend_ax.text(
        text_x, y,
        "Others",
        ha="left",
        va="center",
        fontsize=10,
        transform=legend_ax.transAxes
    )

# =========================
# SAVE
# =========================

outfile = os.path.join(output_folder, "ST_trends_combined_panel.pdf")
plt.savefig(outfile, bbox_inches="tight")
plt.close()

print(f"Saved: {outfile}")