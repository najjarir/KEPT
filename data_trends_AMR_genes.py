import pandas as pd
import numpy as np
import statsmodels.api as sm
import os
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# =========================
# INPUT / OUTPUT
# =========================

input_file = r"D:\Python\AMR_Iris\Input\Data_GENES_AMR.xlsx"
output_file = r"D:\Python\AMR_Iris\Output\ARG_trends_region_income.xlsx"

# =========================
# SETTINGS
# =========================

genes_selecionados = [
    # Aminoglycoside
    "strA", "strB", "rmtB", "armA", "aph3-Ia", "aph(3')-VI", "aadA2", "aadA",
    "aac(6')-Ib-cr", "aac(6')-Ib'", "aac(3)-IIa", "aac(3)-IId",
    # Carbapenem
    "KPC-2", "KPC-3", "NDM-1", "NDM-5", "OXA-48", "VIM-1", "IMP-4",
    # Colistin
    "mcr-9.1", "mcr-1.1",
    # Third generation cephalosporins
    "CTX-M-15", "CTX-M-14", "CTX-M-65", "CTX-M-2", "SHV-12",
    # Fosfomycin
    "fosA3", "fosA7",
    # Fluoroquinolone
    "qnrA1", "qnrS1", "qnrB1",
    # MLS
    "ermB", "mphA", "msrE", "mphE",
    # Phenicol
    "catA1", "catII.2", "catB3", "floR", "cmlA1",
    # Sulfametaxazole
    "sul1", "sul2", "sul3",
    # Trimethoprim
    "dfrA1", "dfrA12", "dfrA14", "dfrA5",
    # Tetracycline
    "tet(A)", "tet(D)", "tet(B)",
    # Tigecycline
    "tmexCD1-toprJ1", "Tet(X4)"
]

colunas_genes = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired", "Gly_acquired",
    "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired",
    "Tgc_acquired", "Tmt_acquired", "Bla_acquired", "Bla_inhR_acquired",
    "Bla_ESBL_acquired", "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired", "Bla_chr"
]

MIN_POSITIVE = 5

# =========================
# LOAD DATA
# =========================

df = pd.read_excel(input_file)
df.columns = df.columns.str.strip()

# detectar coluna de income
if "income_group" in df.columns:
    income_col = "income_group"
elif "income" in df.columns:
    income_col = "income"
else:
    raise ValueError("No income column found. Expected 'income_group' or 'income'.")

# checagens mínimas
required_cols = ["year", "region", income_col] + [c for c in colunas_genes if c in df.columns]
missing_required = [c for c in ["year", "region", income_col] if c not in df.columns]
if missing_required:
    raise ValueError(f"Missing required columns: {missing_required}")

# limpar colunas principais
df["year"] = pd.to_numeric(df["year"], errors="coerce")
df["region"] = df["region"].astype(str).str.strip()
df[income_col] = df[income_col].astype(str).str.strip()

df = df[
    df["year"].notna() &
    df["region"].notna() &
    df[income_col].notna() &
    (df["region"] != "") &
    (df[income_col] != "") &
    (df["region"].str.lower() != "nan") &
    (df[income_col].str.lower() != "nan")
].copy()

df["year"] = df["year"].astype(int)

# preencher colunas de genes ausentes com vazio, se alguma não existir no arquivo
for col in colunas_genes:
    if col not in df.columns:
        df[col] = ""

# =========================
# EXTRACT UNIQUE GENES PER ISOLATE
# =========================

def extract_genes(row):
    genes = []
    for col in colunas_genes:
        val = row[col]
        if pd.isna(val):
            continue
        val = str(val).strip()
        if val == "" or val.lower() == "nan":
            continue
        parts = [g.strip() for g in val.split(";") if g.strip() != ""]
        genes.extend(parts)

    # remove duplicatas preservando ordem
    seen = set()
    unique_genes = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)

    return unique_genes

df["all_genes"] = df.apply(extract_genes, axis=1)

# =========================
# CREATE BINARY PRESENCE COLUMNS
# =========================

for gene in genes_selecionados:
    df[gene] = df["all_genes"].apply(lambda x: 1 if gene in x else 0)

# =========================
# LOGISTIC TRENDS FUNCTION
# =========================

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

def run_trends(data, group_var):
    results = []
    yearly_rows = []

    groups = sorted(data[group_var].dropna().unique())

    for group in groups:
        df_group = data[data[group_var] == group].copy()

        if df_group.empty:
            continue

        n_total_group = len(df_group)
        year_min = int(df_group["year"].min())
        year_max = int(df_group["year"].max())

        for gene in genes_selecionados:
            y = df_group[gene].astype(int)

            n_positive = int(y.sum())
            n_negative = int((1 - y).sum())

            # guardar tabela anual
            yearly_summary = (
                df_group.groupby("year")[gene]
                .agg(N_positive="sum", N_total="count")
                .reset_index()
            )
            yearly_summary["Prevalence_%"] = (yearly_summary["N_positive"] / yearly_summary["N_total"] * 100).round(4)
            yearly_summary[group_var] = group
            yearly_summary["Gene"] = gene

            yearly_rows.append(
                yearly_summary[[group_var, "Gene", "year", "N_positive", "N_total", "Prevalence_%"]]
                .rename(columns={"year": "Year"})
            )

            # critérios mínimos para regressão
            if n_positive < MIN_POSITIVE:
                results.append({
                    group_var: group,
                    "Gene": gene,
                    "N_total": n_total_group,
                    "N_positive": n_positive,
                    "N_negative": n_negative,
                    "Prevalence_%": round(n_positive / n_total_group * 100, 4),
                    "Year_min": year_min,
                    "Year_max": year_max,
                    "beta_year": np.nan,
                    "OR_per_year": np.nan,
                    "CI95_low_year": np.nan,
                    "CI95_high_year": np.nan,
                    "OR_per_decade": np.nan,
                    "CI95_low_decade": np.nan,
                    "CI95_high_decade": np.nan,
                    "p_value": np.nan,
                    "Status": f"< {MIN_POSITIVE} positives"
                })
                continue

            if y.nunique() < 2:
                results.append({
                    group_var: group,
                    "Gene": gene,
                    "N_total": n_total_group,
                    "N_positive": n_positive,
                    "N_negative": n_negative,
                    "Prevalence_%": round(n_positive / n_total_group * 100, 4),
                    "Year_min": year_min,
                    "Year_max": year_max,
                    "beta_year": np.nan,
                    "OR_per_year": np.nan,
                    "CI95_low_year": np.nan,
                    "CI95_high_year": np.nan,
                    "OR_per_decade": np.nan,
                    "CI95_low_decade": np.nan,
                    "CI95_high_decade": np.nan,
                    "p_value": np.nan,
                    "Status": "No variation"
                })
                continue

            try:
                # centralizar ano para melhorar estabilidade numérica
                year_centered = df_group["year"] - df_group["year"].mean()
                X = sm.add_constant(year_centered)

                model = sm.Logit(y, X).fit(disp=0)

                coef = model.params.iloc[1]
                pval = model.pvalues.iloc[1]

                ci = model.conf_int()
                ci_low = ci.iloc[1, 0]
                ci_high = ci.iloc[1, 1]

                or_year = np.exp(coef)
                or_year_low = np.exp(ci_low)
                or_year_high = np.exp(ci_high)

                or_decade = np.exp(coef * 10)
                or_decade_low = np.exp(ci_low * 10)
                or_decade_high = np.exp(ci_high * 10)

                # limpar infinitos
                if not np.isfinite(or_year):
                    or_year = np.nan
                if not np.isfinite(or_year_low):
                    or_year_low = np.nan
                if not np.isfinite(or_year_high):
                    or_year_high = np.nan
                if not np.isfinite(or_decade):
                    or_decade = np.nan
                if not np.isfinite(or_decade_low):
                    or_decade_low = np.nan
                if not np.isfinite(or_decade_high):
                    or_decade_high = np.nan

                results.append({
                    group_var: group,
                    "Gene": gene,
                    "N_total": n_total_group,
                    "N_positive": n_positive,
                    "N_negative": n_negative,
                    "Prevalence_%": round(n_positive / n_total_group * 100, 4),
                    "Year_min": year_min,
                    "Year_max": year_max,
                    "beta_year": coef,
                    "OR_per_year": or_year,
                    "CI95_low_year": or_year_low,
                    "CI95_high_year": or_year_high,
                    "OR_per_decade": or_decade,
                    "CI95_low_decade": or_decade_low,
                    "CI95_high_decade": or_decade_high,
                    "p_value": pval,
                    "Status": "OK"
                })

            except Exception:
                results.append({
                    group_var: group,
                    "Gene": gene,
                    "N_total": n_total_group,
                    "N_positive": n_positive,
                    "N_negative": n_negative,
                    "Prevalence_%": round(n_positive / n_total_group * 100, 4),
                    "Year_min": year_min,
                    "Year_max": year_max,
                    "beta_year": np.nan,
                    "OR_per_year": np.nan,
                    "CI95_low_year": np.nan,
                    "CI95_high_year": np.nan,
                    "OR_per_decade": np.nan,
                    "CI95_low_decade": np.nan,
                    "CI95_high_decade": np.nan,
                    "p_value": np.nan,
                    "Status": "Model failed"
                })

    results_df = pd.DataFrame(results)
    yearly_df = pd.concat(yearly_rows, ignore_index=True) if yearly_rows else pd.DataFrame()

    return results_df, yearly_df

# =========================
# RUN TRENDS
# =========================

region_results, region_yearly = run_trends(df, "region")
income_results, income_yearly = run_trends(df, income_col)

# =========================
# SAVE EXCEL
# =========================

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    region_results.to_excel(writer, sheet_name="region_trends", index=False)
    income_results.to_excel(writer, sheet_name="income_trends", index=False)
    region_yearly.to_excel(writer, sheet_name="region_yearly_counts", index=False)
    income_yearly.to_excel(writer, sheet_name="income_yearly_counts", index=False)

print(f"Saved: {output_file}")