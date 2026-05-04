import pandas as pd
import numpy as np
import statsmodels.api as sm
import os
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# =========================
# INPUT / OUTPUT
# =========================

input_file = r"D:\Python\AMR_Iris\Input\ST_and_Enzymes_KL.xlsx"
output_file = r"D:\Python\AMR_Iris\Output\ST_KL_association_full.xlsx"

# =========================
# LOAD DATA
# =========================

df = pd.read_excel(input_file)

df.columns = df.columns.str.strip()

df["ST"] = df["ST"].astype(str).str.strip()
df["region"] = df["region"].astype(str).str.strip()
df["K_locus"] = df["K_locus"].astype(str).str.strip()

df = df[
    df["region"].notna() &
    df["ST"].notna() &
    df["K_locus"].notna() &
    (df["region"] != "") &
    (df["ST"] != "") &
    (df["K_locus"] != "") &
    (df["K_locus"].str.lower() != "untypeable")
].copy()

# =========================
# SETTINGS
# =========================

target_ST = ["ST258", "ST11", "ST307", "ST147", "ST15"]

MIN_ST_KL = 5

regions = df["region"].dropna().unique()

results = []

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# =========================
# LOOP
# =========================

for region in regions:

    df_region = df[df["region"] == region].copy()

    if df_region.empty:
        continue

    total_n = len(df_region)

    valid_kl = df_region["K_locus"].dropna().unique().tolist()

    for st in target_ST:

        df_region["ST_bin"] = (df_region["ST"] == st).astype(int)

        n_st = int(df_region["ST_bin"].sum())
        if n_st == 0:
            continue

        if df_region["ST_bin"].nunique() < 2:
            continue

        for kl in valid_kl:

            df_region["KL_bin"] = (df_region["K_locus"] == kl).astype(int)

            n_kl = int(df_region["KL_bin"].sum())
            n_st_and_kl = int(((df_region["ST_bin"] == 1) & (df_region["KL_bin"] == 1)).sum())

            if n_kl == 0:
                continue
            if n_kl == len(df_region):
                continue
            if n_st_and_kl < MIN_ST_KL:
                continue
            if df_region["KL_bin"].nunique() < 2:
                continue

            # =========================
            # PERCENTUAIS
            # =========================

            pct_st_with_kl = (n_st_and_kl / n_st) * 100
            pct_kl_with_st = (n_st_and_kl / n_kl) * 100

            # =========================
            # REGRESSÃO
            # =========================

            try:
                X = sm.add_constant(df_region["KL_bin"])
                y = df_region["ST_bin"]

                model = sm.Logit(y, X).fit(disp=0)

                coef = model.params["KL_bin"]
                pval = model.pvalues["KL_bin"]

                OR = np.exp(coef)
                CI_low, CI_high = np.exp(model.conf_int().loc["KL_bin"])

                if not np.isfinite(OR):
                    OR = np.nan
                if not np.isfinite(CI_low):
                    CI_low = np.nan
                if not np.isfinite(CI_high):
                    CI_high = np.nan

            except Exception:
                OR, CI_low, CI_high, pval = np.nan, np.nan, np.nan, np.nan

            results.append({
                "Region": region,
                "ST": st,
                "KL": kl,
                "N_region": total_n,
                "N_ST": n_st,
                "N_KL": n_kl,
                "N_ST_and_KL": n_st_and_kl,
                "%_ST_with_KL": round(pct_st_with_kl, 2),
                "%_KL_with_ST": round(pct_kl_with_st, 2),
                "OR": OR,
                "CI_low": CI_low,
                "CI_high": CI_high,
                "p_value": pval
            })

# =========================
# OUTPUT
# =========================

results_df = pd.DataFrame(results)

if results_df.empty:
    print("No valid associations were generated.")
else:
    results_df = results_df.sort_values(
        ["Region", "ST", "%_ST_with_KL"],
        ascending=[True, True, False]
    )

    results_df.to_excel(output_file, index=False)

    print("Saved:", output_file)