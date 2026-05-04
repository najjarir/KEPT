import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import circlify
from matplotlib.colors import to_rgb

# =====================
# PATHS
# =====================

input_file = r"D:\Python\AMR_Iris\Input\ST_and_Enzymes_KL.xlsx"
output_folder = r"D:\Python\AMR_Iris\Output"

os.makedirs(output_folder, exist_ok=True)

# =====================
# MATPLOTLIB SETTINGS
# =====================

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.size"] = 10

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# =====================
# LOAD DATA
# =====================

df_full = pd.read_excel(input_file)

df_full["ST"] = df_full["ST"].astype(str).str.strip()
df_full["K_locus"] = df_full["K_locus"].astype(str).str.strip()

df_full = df_full[
    df_full["ST"].notna() &
    df_full["K_locus"].notna() &
    (df_full["ST"] != "") &
    (df_full["K_locus"] != "")
].copy()

# =====================
# FIRST FIGURE
# KL within top ST
# =====================

target_st = ["ST258", "ST11", "ST307", "ST147", "ST15"]

df = df_full[df_full["ST"].isin(target_st)].copy()

counts = (
    df.groupby(["ST", "K_locus"])
    .size()
    .reset_index(name="count")
)

filtered = []
pct_lookup = {}

for st in target_st:

    sub = counts[counts["ST"] == st].copy()

    total = sub["count"].sum()

    sub["pct"] = sub["count"]/total*100

    sub = sub[sub["pct"] > 1]

    for _, r in sub.iterrows():
        pct_lookup[(st, r["K_locus"])] = r["pct"]

    filtered.append(sub)

counts_filt = pd.concat(filtered)

st_totals = (
    counts.groupby("ST")["count"]
    .sum()
)

st_colors = {
    "ST258": "#f4a3a3",
    "ST11": "#a8c6ff",
    "ST307": "#a8e6a1",
    "ST147": "#d2b6ff",
    "ST15": "#ffd2a6"
}


def darker(color):

    r, g, b = to_rgb(color)

    return (r*0.8, g*0.8, b*0.8)

# =====================
# BUILD HIERARCHY
# =====================


hier = []

for st in target_st:

    sub = counts_filt[counts_filt["ST"] == st]

    children = []

    for _, row in sub.iterrows():

        children.append({
            "id": row["K_locus"],
            "datum": row["count"]
        })

    hier.append({
        "id": st,
        "datum": st_totals[st],
        "children": children
    })

# =====================
# CIRCLES
# =====================

circles = circlify.circlify(
    hier,
    show_enclosure=False,
    target_enclosure=circlify.Circle(0, 0, 1)
)

level1 = [c for c in circles if c.level == 1]
level2 = [c for c in circles if c.level == 2]

# =====================
# FIGURE
# =====================

fig, ax = plt.subplots(figsize=(14, 10))

ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)

ax.set_aspect("equal")

ax.axis("off")

# =====================
# ST circles
# =====================

for c in level1:

    st = c.ex["id"]

    x, y, r = c.x, c.y, c.r

    ax.add_patch(
        plt.Circle(
            (x, y),
            r,
            facecolor=darker(st_colors[st]),
            edgecolor="white",
            linewidth=0,
            alpha=0.65
        )
    )

    n_total = int(st_totals[st])

    label_x = x
    label_y = y + r + 0.01
    ha = "center"
    va = "bottom"

    if st == "ST258":
        label_x = x + r + 0.1
        label_y = y
        ha = "center"
        va = "center"

    elif st == "ST147":
        label_x = x
        label_y = y - r - 0.01
        ha = "center"
        va = "top"

    elif st == "ST11":
        label_x = x - r * 0.25
        label_y = y + r + 0.01
        ha = "center"
        va = "bottom"

    ax.text(
        label_x,
        label_y,
        f"{st}\n(n = {n_total})",
        ha=ha,
        va=va,
        fontsize=13,
        fontweight="bold"
    )

# =====================
# KL circles
# =====================

for c in level2:

    x, y, r = c.x, c.y, c.r

    kl = c.ex["id"]

    parent = None

    for stc in level1:

        dx = (x-stc.x)
        dy = (y-stc.y)

        if (dx*dx+dy*dy)**0.5 + r <= stc.r:

            parent = stc.ex["id"]

            break

    if parent is None:
        continue

    if parent == "ST307" and kl == "KL102":
        edge = "none"
        lw = 0
    else:
        edge = "black"
        lw = 0.6

    ax.add_patch(
        plt.Circle(
            (x, y),
            r,
            facecolor=st_colors[parent],
            edgecolor=edge,
            linewidth=lw,
            alpha=0.75
        )
    )

    pct = pct_lookup[(parent, kl)]

    ax.text(
        x,
        y,
        f"{kl}\n{pct:.1f}%",
        ha="center",
        va="center",
        fontsize=max(3, 12*r*10)
    )

plt.tight_layout()

outfile = os.path.join(
    output_folder,
    "PackedCircles_KL_ST.pdf"
)

plt.savefig(outfile, bbox_inches="tight")

plt.close()

print("Saved:", outfile)

# =========================================================
# SECOND FIGURE
# ST within selected KL
# =========================================================

target_kl = ["KL107", "KL64", "KL102", "KL2", "KL106", "KL24", "KL51", "KL1"]

df_kl = df_full[df_full["K_locus"].isin(target_kl)].copy()

counts_kl = (
    df_kl.groupby(["K_locus", "ST"])
    .size()
    .reset_index(name="count")
)

filtered_kl = []
pct_lookup_kl = {}

for kl in target_kl:

    sub = counts_kl[counts_kl["K_locus"] == kl].copy()

    if sub.empty:
        continue

    total = sub["count"].sum()
    sub["pct"] = sub["count"] / total * 100
    sub = sub[sub["pct"] > 1].copy()

    for _, r in sub.iterrows():
        pct_lookup_kl[(kl, r["ST"])] = r["pct"]

    filtered_kl.append(sub)

counts_filt_kl = pd.concat(filtered_kl, ignore_index=True)

kl_totals = (
    counts_kl.groupby("K_locus")["count"]
    .sum()
)

kl_colors = {
    "KL107": "#ffd6e0",
    "KL64": "#d9edff",
    "KL102": "#e2f7d4",
    "KL2": "#efe0ff",
    "KL106": "#ffe7cc",
    "KL24": "#dff6f0",
    "KL51": "#fff3bf",
    "KL1": "#66ccda"
}

hier_kl = []

for kl in target_kl:

    sub = counts_filt_kl[counts_filt_kl["K_locus"] == kl]

    children = []

    for _, row in sub.iterrows():
        children.append({
            "id": row["ST"],
            "datum": row["count"]
        })

    hier_kl.append({
        "id": kl,
        "datum": int(kl_totals.get(kl, 0)),
        "children": children
    })

circles_kl = circlify.circlify(
    hier_kl,
    show_enclosure=False,
    target_enclosure=circlify.Circle(0, 0, 1)
)

level1_kl = [c for c in circles_kl if c.level == 1]
level2_kl = [c for c in circles_kl if c.level == 2]

fig, ax = plt.subplots(figsize=(14, 10))
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.18, 1.15)
ax.set_aspect("equal")
ax.axis("off")

for c in level1_kl:

    kl = c.ex["id"]
    x, y, r = c.x, c.y, c.r

    ax.add_patch(
        plt.Circle(
            (x, y),
            r,
            facecolor=darker(kl_colors[kl]),
            edgecolor="none",
            linewidth=0,
            alpha=0.65
        )
    )

    n_total = int(kl_totals.get(kl, 0))

    # posição padrão
    label_x = x
    label_y = y + r + 0.01
    ha = "center"
    va = "bottom"

    # ajustes finos
    if kl == "KL102":
        label_x = x + r * 0.15
        label_y = y + r + 0.01

    elif kl == "KL106":
        label_x = x - r - 0.11
        label_y = y
        ha = "center"
        va = "center"

    elif kl == "KL24":
        label_x = x + r + 0.11
        label_y = y
        ha = "center"
        va = "center"

    elif kl == "KL107":
        label_x = x + r + 0.11
        label_y = y
        ha = "center"
        va = "center"

    elif kl == "KL64":
        label_x = x - r - 0.11
        label_y = y
        ha = "center"
        va = "center"

    elif kl == "KL2":
        label_x = x
        label_y = y - r - 0.02
        ha = "center"
        va = "top"

    ax.text(
        label_x,
        label_y,
        f"{kl}\n(n = {n_total})",
        ha=ha,
        va=va,
        fontsize=13,
        fontweight="bold"
    )


for c in level2_kl:

    x, y, r = c.x, c.y, c.r
    st = c.ex["id"]

    parent = None

    for klc in level1_kl:
        dx = x - klc.x
        dy = y - klc.y

        if (dx * dx + dy * dy) ** 0.5 + r <= klc.r + 1e-9:
            parent = klc.ex["id"]
            break

    if parent is None:
        continue

    ax.add_patch(
        plt.Circle(
            (x, y),
            r,
            facecolor=kl_colors[parent],
            edgecolor="black",
            linewidth=0.6,
            alpha=0.75
        )
    )

    pct = pct_lookup_kl[(parent, st)]

    ax.text(
        tx,
        ty,
        f"{st}\n{pct:.1f}%",
        ha=center,
        va=center,
        fontsize=max(3, 12 * r * 10)
    )

plt.tight_layout()

outfile_kl = os.path.join(
    output_folder,
    "PackedCircles_ST_KL.pdf"
)

plt.savefig(outfile_kl, bbox_inches="tight")
plt.close()

print("Saved:", outfile_kl)

# =========================================================
# THIRD AND FOURTH FIGURES
# KL mostly encountered within selected enzymes
# =========================================================

# =====================
# CLEAN ENZYME COLUMNS
# =====================


def clean_enzymes(cell):
    if pd.isna(cell):
        return []
    enzymes = [x.strip() for x in str(cell).split(";") if str(x).strip() != ""]
    return list(dict.fromkeys(enzymes))  # remove duplicates preserving order


df_full["Carb_list"] = df_full["Bla_Carb_acquired"].apply(clean_enzymes)
df_full["ESBL_list"] = df_full["Bla_ESBL_acquired"].apply(clean_enzymes)

# =====================
# COLORS
# =====================

carb_colors = {
    "KPC-2": "#d73027",
    "KPC-3": "#fc8d59",
    "NDM-1": "#4575b4",
    "NDM-5": "#74add1",
    "OXA-48": "#66bd63"
}

esbl_colors = {
    "CTX-M-15": "#7b3294",
    "CTX-M-14": "#c2a5cf",
    "CTX-M-65": "#008837",
    "CTX-M-2": "#fdb863",
    "SHV-12": "#5ab4ac"
}

dark_fill = {
    "KPC-2", "NDM-1",
    "CTX-M-15", "CTX-M-65"
}

# =====================
# GENERIC FUNCTION
# =====================


def plot_enzyme_to_kl_packed(
    data,
    enzyme_list_column,
    focus_enzymes,
    color_map,
    outfile_name
):

    # explode enzyme column
    tmp = data.copy().explode(enzyme_list_column)
    tmp = tmp[tmp[enzyme_list_column].isin(focus_enzymes)].copy()

    # keep valid KL
    tmp = tmp[
        tmp["K_locus"].notna() &
        (tmp["K_locus"] != "") &
        (tmp["K_locus"].astype(str).str.lower() != "nan")
    ].copy()

    # count KL within each enzyme
    counts = (
        tmp.groupby([enzyme_list_column, "K_locus"])
        .size()
        .reset_index(name="count")
    )

    # filter for plotting
    # Untypeable is kept in denominator, but never plotted
    filtered = []
    pct_lookup = {}

    for enz in focus_enzymes:
        sub = counts[counts[enzyme_list_column] == enz].copy()

        if sub.empty:
            continue

        total = sub["count"].sum()
        sub["pct"] = sub["count"] / total * 100

        sub_plot = sub[
            (sub["pct"] > 1) &
            (sub["K_locus"].astype(str).str.strip().str.lower() != "untypeable")
        ].copy()

        for _, r in sub_plot.iterrows():
            pct_lookup[(enz, r["K_locus"])] = r["pct"]

        if not sub_plot.empty:
            filtered.append(sub_plot)

    counts_filt = pd.concat(filtered, ignore_index=True)

    # totals per enzyme (includes Untypeable)
    enz_totals = counts.groupby(enzyme_list_column)["count"].sum()

    # build hierarchy
    hier = []

    for enz in focus_enzymes:
        sub = counts_filt[counts_filt[enzyme_list_column] == enz].copy()

        children = []

        for _, row in sub.iterrows():
            children.append({
                "id": row["K_locus"],
                "datum": row["count"]
            })

        hier.append({
            "id": enz,
            "datum": int(enz_totals.get(enz, 0)),
            "children": children
        })

    # run circlify
    circles = circlify.circlify(
        hier,
        show_enclosure=False,
        target_enclosure=circlify.Circle(0, 0, 1)
    )

    level1 = [c for c in circles if c.level == 1]  # outer enzyme circles
    level2 = [c for c in circles if c.level == 2]  # inner KL circles

    # figure
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.18, 1.15)
    ax.set_aspect("equal")
    ax.axis("off")

    # draw outer enzyme circles
    for c in level1:
        enz = c.ex["id"]
        x, y, r = c.x, c.y, c.r

        ax.add_patch(
            plt.Circle(
                (x, y),
                r,
                facecolor=darker(color_map[enz]),
                edgecolor="none",
                linewidth=0,
                alpha=0.65
            )
        )

        n_total = int(enz_totals.get(enz, 0))

        # default label position
        label_x = x
        label_y = y + r + 0.01
        ha = "center"
        va = "bottom"
        fs = 13

        # fine adjustments
        if enz == "KPC-2":
            label_x = x + r + 0.11
            label_y = y
            ha = "center"
            va = "center"

        elif enz == "KPC-3":
            label_x = x - r - 0.11
            label_y = y
            ha = "center"
            va = "center"

        elif enz == "OXA-48":
            label_x = x
            label_y = y - r - 0.02
            ha = "center"
            va = "top"

        elif enz == "CTX-M-14":
            label_x = x - r - 0.11
            label_y = y
            ha = "center"
            va = "center"

        elif enz == "SHV-12":
            label_x = x - r - 0.11
            label_y = y
            ha = "center"
            va = "center"

        elif enz == "CTX-M-65":
            label_x = x - r - 0.11
            label_y = y
            ha = "center"
            va = "center"

        ax.text(
            label_x,
            label_y,
            f"{enz}\n(n = {n_total})",
            ha=ha,
            va=va,
            fontsize=fs,
            fontweight="bold"
        )

    # draw inner KL circles
    for c in level2:
        x, y, r = c.x, c.y, c.r
        kl = c.ex["id"]

        parent = None

        for ec in level1:
            dx = x - ec.x
            dy = y - ec.y

            if (dx * dx + dy * dy) ** 0.5 + r <= ec.r + 1e-9:
                parent = ec.ex["id"]
                break

        if parent is None:
            continue

        ax.add_patch(
            plt.Circle(
                (x, y),
                r,
                facecolor=color_map[parent],
                edgecolor="black",
                linewidth=0.6,
                alpha=0.75
            )
        )

        pct = pct_lookup[(parent, kl)]
        text_color = "white" if parent in dark_fill else "black"

        # condição de tamanho do círculo
        if r < 0.015:
            label = f"{kl}"          # apenas KL
        else:
            label = f"{kl}\n{pct:.1f}%"   # KL + %

        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=max(3, 12 * r * 10),
            color=text_color
        )

    plt.tight_layout()

    outfile = os.path.join(output_folder, outfile_name)
    plt.savefig(outfile, bbox_inches="tight")
    plt.close()

    print("Saved:", outfile)

# =====================
# FIGURE 3
# CARBAPENEMASE -> KL
# =====================


carb_focus = ["KPC-2", "KPC-3", "NDM-1", "NDM-5", "OXA-48"]

plot_enzyme_to_kl_packed(
    data=df_full,
    enzyme_list_column="Carb_list",
    focus_enzymes=carb_focus,
    color_map=carb_colors,
    outfile_name="PackedCircles_KL_Carbapenemases.pdf"
)

# =====================
# FIGURE 4
# ESBL -> KL
# =====================

esbl_focus = ["CTX-M-15", "CTX-M-14", "CTX-M-65", "CTX-M-2", "SHV-12"]

plot_enzyme_to_kl_packed(
    data=df_full,
    enzyme_list_column="ESBL_list",
    focus_enzymes=esbl_focus,
    color_map=esbl_colors,
    outfile_name="PackedCircles_KL_ESBL.pdf"
)
