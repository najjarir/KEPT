import os
import warnings

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Font configuration
# ---------------------------------------------------------------------------
_arial_path = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"
if os.path.exists(_arial_path):
    fm.fontManager.addfont(_arial_path)
    plt.rcParams["font.family"] = "Arial"
else:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Liberation Sans", "Arial"]

# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------
df = pd.read_csv("")

df = df.reset_index(drop=True)
df["genome_id"] = df.index
df["ST_formatted"] = (
    df["ST"].astype(str).apply(lambda x: x if x.startswith("ST") else f"ST{x}")
)

# ---------------------------------------------------------------------------
# Gene panel (ordered by resistance class)
# ---------------------------------------------------------------------------
GENES = [
    # Aminoglycoside
    "strA", "strB", "rmtB", "armA", "aph3-Ia", "aph(3')-VI", "aadA2", "aadA",
    "aac(6')-Ib-cr", "aac(6')-Ib'", "aac(3)-IIa", "aac(3)-IId",
    # Carbapenem
    "KPC-2", "KPC-3", "NDM-1", "NDM-5", "OXA-48", "VIM-1", "IMP-4",
    # Colistin
    "mcr-9.1", "mcr-1.1",
    # Third-generation cephalosporins
    "CTX-M-15", "CTX-M-14", "CTX-M-65", "CTX-M-2", "SHV-12",
    # Fosfomycin
    "fosA3", "fosA7",
    # Fluoroquinolone
    "qnrA1", "qnrS1", "qnrB1",
    # MLS
    "ermB", "mphA", "msrE", "mphE",
    # Phenicol
    "catA1", "catII.2", "catB3", "floR", "cmlA1",
    # Sulfonamide
    "sul1", "sul2", "sul3",
    # Trimethoprim
    "dfrA1", "dfrA12", "dfrA14", "dfrA5",
    # Tetracycline
    "tet(A)", "tet(D)", "tet(B)",
    # Tigecycline
    "tmexCD1-toprJ1", "Tet(X4)",
]

AMR_COLUMNS = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired", "Gly_acquired",
    "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired",
    "Tgc_acquired", "Tmt_acquired", "Bla_acquired", "Bla_inhR_acquired",
    "Bla_ESBL_acquired", "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired", "Bla_chr",
]

TARGET_STS = ["ST258", "ST11", "ST307", "ST147", "ST15"]

# Resistance class boundaries (number of genes per class, same order as GENES)
RESISTANCE_CLASSES = {
    "Aminoglycoside": 12,
    "Carbapenem": 7,
    "Colistin": 2,
    "3rd-gen. cephalosporins": 5,
    "Fosfomycin": 2,
    "Fluoroquinolone": 3,
    "MLS": 4,
    "Phenicol": 5,
    "Sulfonamide": 3,
    "Trimethoprim": 4,
    "Tetracycline": 3,
    "Tigecycline": 2,
}

# ---------------------------------------------------------------------------
# Reshape: one row per genome × gene
# ---------------------------------------------------------------------------
keep_cols = ["genome_id", "year", "region", "ST_formatted"] + AMR_COLUMNS
df_long = (
    df[keep_cols]
    .melt(
        id_vars=["genome_id", "year", "region", "ST_formatted"],
        value_vars=AMR_COLUMNS,
        value_name="raw_genes",
    )
    .dropna(subset=["raw_genes"])
)
df_long["gene"] = df_long["raw_genes"].str.split(";")
df_long = df_long.explode("gene")
df_long["gene"] = df_long["gene"].str.strip()
df_long = df_long[df_long["gene"].isin(GENES)].drop_duplicates(
    subset=["genome_id", "gene"]
)

# ---------------------------------------------------------------------------
# Block 1 — prevalence by year
# ---------------------------------------------------------------------------
genomes_per_year = df.groupby("year")["genome_id"].nunique()
all_years = sorted(df["year"].dropna().unique())

gene_year = (
    df_long.groupby(["year", "gene"])["genome_id"]
    .nunique()
    .reset_index(name="count")
)
gene_year["total"] = gene_year["year"].map(genomes_per_year)
gene_year["prevalence"] = gene_year["count"] / gene_year["total"] * 100

pivot_year = (
    gene_year.pivot(index="gene", columns="year", values="prevalence")
    .reindex(columns=all_years, fill_value=0)
    .reindex(GENES)
    .fillna(0)
)
pivot_year.columns = pivot_year.columns.astype(str)

# ---------------------------------------------------------------------------
# Block 2 — prevalence by region (with global summary)
# ---------------------------------------------------------------------------
total_global = df["genome_id"].nunique()
genomes_per_region = df.groupby("region")["genome_id"].nunique()

gene_global = (
    df_long.groupby("gene")["genome_id"]
    .nunique()
    .reset_index(name="count")
)
gene_global["prevalence"] = gene_global["count"] / total_global * 100
gene_global["region"] = "Global"

gene_region = (
    df_long.groupby(["region", "gene"])["genome_id"]
    .nunique()
    .reset_index(name="count")
)
gene_region["total"] = gene_region["region"].map(genomes_per_region)
gene_region["prevalence"] = gene_region["count"] / gene_region["total"] * 100

pivot_region = (
    pd.concat(
        [gene_global[["region", "gene", "prevalence"]],
         gene_region[["region", "gene", "prevalence"]]]
    )
    .pivot(index="gene", columns="region", values="prevalence")
)
region_order = ["Global"] + [c for c in pivot_region.columns if c != "Global"]
pivot_region = pivot_region[region_order].reindex(GENES).fillna(0)

# ---------------------------------------------------------------------------
# Block 3 — prevalence by sequence type
# ---------------------------------------------------------------------------
df_st = df[df["ST_formatted"].isin(TARGET_STS)]
genomes_per_st = df_st.groupby("ST_formatted")["genome_id"].nunique()

gene_st = (
    df_long[df_long["ST_formatted"].isin(TARGET_STS)]
    .groupby(["ST_formatted", "gene"])["genome_id"]
    .nunique()
    .reset_index(name="count")
)
gene_st["total"] = gene_st["ST_formatted"].map(genomes_per_st)
gene_st["prevalence"] = gene_st["count"] / gene_st["total"] * 100

pivot_st = gene_st.pivot(index="gene", columns="ST_formatted", values="prevalence")
st_order = [st for st in TARGET_STS if st in pivot_st.columns]
pivot_st = pivot_st.reindex(columns=st_order, fill_value=0).reindex(GENES).fillna(0)

# ---------------------------------------------------------------------------
# Composite pivot (year | region | ST)
# ---------------------------------------------------------------------------
pivot_master = pd.concat([pivot_year, pivot_region, pivot_st], axis=1)

n_genes = len(GENES)
div1 = len(pivot_year.columns)
div2 = div1 + len(pivot_region.columns)

# Horizontal separator positions (between resistance classes)
h_lines = []
pos = 0
for count in RESISTANCE_CLASSES.values():
    pos += count
    if pos < n_genes:
        h_lines.append(pos)

# ---------------------------------------------------------------------------
# Figure 1 — canvas version (PNG, no square cells)
# ---------------------------------------------------------------------------
fig_w = max(15, len(pivot_master.columns) * 0.7)
fig_h = n_genes * 0.35

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

sns.heatmap(
    pivot_master,
    ax=ax,
    cmap="flare",
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    vmin=0,
    vmax=100,
    annot_kws={"size": 9},
    cbar_kws={"shrink": 0.8, "aspect": 30, "pad": 0.02},
)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=10)
cbar.set_label("Prevalence (%)", fontsize=12)

ax.vlines(x=[div1, div2], ymin=0, ymax=n_genes, colors="white", lw=5)

ax.set_title("")
ax.set_xlabel("Year  |  Region  |  Sequence Type", fontsize=14, labelpad=10)
ax.set_ylabel("Gene", fontsize=14, labelpad=10)
ax.tick_params(axis="x", rotation=45, labelsize=12)
ax.tick_params(axis="y", rotation=0, labelsize=11)

plt.savefig("Heatmap_Master_Canvas.png", dpi=600, bbox_inches="tight")
plt.show()
plt.close()

# ---------------------------------------------------------------------------
# Figure 2 — publication version (PDF, square cells, class separators)
# ---------------------------------------------------------------------------
fig_w2 = max(12, len(pivot_master.columns) * 0.6)
fig_h2 = n_genes * 0.4

fig, ax = plt.subplots(figsize=(fig_w2, fig_h2))

sns.heatmap(
    pivot_master,
    ax=ax,
    cmap="flare",
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    vmin=0,
    vmax=100,
    square=True,
    annot_kws={"size": 9},
    cbar_kws={"shrink": 0.6, "aspect": 40, "pad": 0.02},
)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=11)
cbar.set_label("Prevalence (%)", fontsize=12)

ax.vlines(x=[div1, div2], ymin=0, ymax=n_genes, colors="white", lw=4)
ax.hlines(y=h_lines, xmin=0, xmax=len(pivot_master.columns), colors="white", lw=6)

ax.set_title("")
ax.set_xlabel("Year  |  Region  |  Sequence Type", fontsize=14, labelpad=15)
ax.set_ylabel("Gene", fontsize=14, labelpad=10)
ax.tick_params(axis="x", rotation=45, labelsize=11)
ax.tick_params(axis="y", rotation=0, labelsize=11)

plt.savefig("Heatmap_Master.pdf", dpi=600, bbox_inches="tight", format="pdf")
plt.show()
plt.close()
