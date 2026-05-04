import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import squarify

# =====================
# PATHS
# =====================

input_file = r"D:\Python\AMR_Iris\Input\ST_and_Enzymes.xlsx"
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

df = pd.read_excel(input_file)

# =====================
# CLEAN ENZYMES
# =====================

def clean_enzymes(cell):

    if pd.isna(cell):
        return []

    enzymes = [x.strip() for x in str(cell).split(";") if str(x).strip() != ""]

    # remove duplicates while preserving order
    return list(dict.fromkeys(enzymes))

df["Carb_list"] = df["Bla_Carb_acquired"].apply(clean_enzymes)
df["ESBL_list"] = df["Bla_ESBL_acquired"].apply(clean_enzymes)

# =====================
# ENZYMES OF INTEREST
# =====================

carb_focus = [
"KPC-2",
"KPC-3",
"NDM-1",
"NDM-5",
"OXA-48"
]

esbl_focus = [
"CTX-M-15",
"CTX-M-14",
"CTX-M-65",
"CTX-M-2",
"SHV-12"
]

# =====================
# COLORS
# =====================

carb_colors = {
"KPC-2":"#d73027",
"KPC-3":"#fc8d59",
"NDM-1":"#4575b4",
"NDM-5":"#74add1",
"OXA-48":"#66bd63"
}

esbl_colors = {
"CTX-M-15":"#7b3294",
"CTX-M-14":"#c2a5cf",
"CTX-M-65":"#008837",
"CTX-M-2":"#fdb863",
"SHV-12":"#5ab4ac"
}

dark_fill = {
"KPC-2","NDM-1",
"CTX-M-15","CTX-M-65"
}

# =====================
# FONT SIZE FUNCTION
# =====================

def get_fontsize(w, h, st_text, pct_text, enzyme=None):

    area = w * h
    min_side = min(w, h)

    # casos extremamente pequenos
    if min_side < 2.2:
        return 4

    if area > 2200:
        return 34

    elif area > 1800:
        return 30

    elif area > 1500:
        return 27

    elif area > 1200:
        return 24

    elif area > 950:
        return 22

    elif area > 750:
        return 20

    elif area > 600:
        return 18

    elif area > 470:
        return 16

    elif area > 360:
        return 15

    elif area > 280:
        return 14

    elif area > 220:
        return 13

    elif area > 170:
        return 12

    elif area > 130:
        return 11

    elif area > 100:
        return 10

    elif area > 75:
        return 9

    elif area > 55:
        return 8

    elif area > 40:
        return 7

    elif area > 30:
        return 6

    elif area > 22:
        return 5.5

    elif area > 15:
        return 5

    else:
        return 4.5

# =====================
# TREEMAP FUNCTION
# =====================

def make_enzyme_treemap(
    data,
    enzyme_column,
    focus_list,
    color_map,
    legend_title,
    output_filename,
    filter_threshold=0.01
):

    tmp = data.explode(enzyme_column).copy()
    tmp = tmp[tmp[enzyme_column].isin(focus_list)].copy()

    counts = (
        tmp.groupby([enzyme_column,"ST"])
        .size()
        .reset_index(name="count")
    )

    filtered_parts = []

    for enzyme in focus_list:

        sub = counts[counts[enzyme_column]==enzyme].copy()

        if sub.empty:
            continue

        total_enzyme = sub["count"].sum()

        sub["freq_within_enzyme"] = sub["count"]/total_enzyme

        sub = sub[sub["freq_within_enzyme"]>filter_threshold].copy()

        if not sub.empty:
            filtered_parts.append(sub)

    counts_filt = pd.concat(filtered_parts)

    enzyme_totals = (
        counts.groupby(enzyme_column)["count"]
        .sum()
        .reindex(focus_list)
        .dropna()
    )

    dataset_total = enzyme_totals.sum()

    enzyme_pct = enzyme_totals/dataset_total*100

    fig, ax = plt.subplots(figsize=(16,9))

    ax.set_xlim(0,100)
    ax.set_ylim(0,100)
    ax.axis("off")

    top_sizes = enzyme_totals.values.tolist()
    top_labels = enzyme_totals.index.tolist()

    top_rects = squarify.squarify(
        squarify.normalize_sizes(top_sizes,100,100),
        0,0,100,100
    )

    for rect,enzyme in zip(top_rects,top_labels):

        x = rect["x"]
        y = rect["y"]
        dx = rect["dx"]
        dy = rect["dy"]

        ax.add_patch(
            plt.Rectangle(
                (x,y),
                dx,
                dy,
                facecolor="none",
                edgecolor="white",
                linewidth=2
            )
        )

        sub = counts_filt[counts_filt[enzyme_column]==enzyme].copy()

        sub = sub.sort_values("count",ascending=False)

        st_sizes = sub["count"].tolist()

        st_rects = squarify.squarify(
            squarify.normalize_sizes(st_sizes,dx,dy),
            x,y,dx,dy
        )

        for st_rect,(_,row) in zip(st_rects,sub.iterrows()):

            sx = st_rect["x"]
            sy = st_rect["y"]
            sdx = st_rect["dx"]
            sdy = st_rect["dy"]

            ax.add_patch(
                plt.Rectangle(
                    (sx,sy),
                    sdx,
                    sdy,
                    facecolor=color_map[enzyme],
                    edgecolor="white",
                    linewidth=1.2
                )
            )

            pct_label = f"{row['freq_within_enzyme']*100:.1f}%"

            label = f"{row['ST']}\n{pct_label}"

            fs = get_fontsize(
                sdx,
                sdy,
                row["ST"],
                pct_label,
                enzyme=enzyme
            )

            text_color = "white" if enzyme in dark_fill else "black"

            ax.text(
                sx+sdx/2,
                sy+sdy/2,
                label,
                ha="center",
                va="center",
                fontsize=fs,
                color=text_color,
                linespacing=0.92,
                clip_on=True
            )

    handles = []

    for enzyme in top_labels:

        label = f"{enzyme} ({enzyme_pct[enzyme]:.1f}%)"

        handles.append(
            plt.Rectangle(
                (0,0),
                1,
                1,
                facecolor=color_map[enzyme],
                edgecolor="none",
                label=label
            )
        )

    ax.legend(
        handles=handles,
        title=legend_title,
        loc="upper center",
        bbox_to_anchor=(0.5,-0.03),
        ncol=len(handles),
        frameon=False
    )

    plt.tight_layout()

    outfile = os.path.join(output_folder,output_filename)

    plt.savefig(outfile,bbox_inches="tight")

    plt.close()

    print(f"Saved: {outfile}")

# =====================
# CARBAPENEMASE TREEMAP
# =====================

make_enzyme_treemap(
    data=df,
    enzyme_column="Carb_list",
    focus_list=carb_focus,
    color_map=carb_colors,
    legend_title="Carbapenemase",
    output_filename="Treemap_Carbapenemases_Hierarchical_final.pdf"
)

# =====================
# ESBL TREEMAP
# =====================

make_enzyme_treemap(
    data=df,
    enzyme_column="ESBL_list",
    focus_list=esbl_focus,
    color_map=esbl_colors,
    legend_title="ESBL",
    output_filename="Treemap_ESBL_Hierarchical_final.pdf"
)