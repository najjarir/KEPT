import os
import pandas as pd

# =====================
# PATHS
# =====================

input_file = r"D:\Python\AMR_Iris\Input\Data_Classes_MDR_XDR_virulence.xlsx"
output_file = r"D:\Python\AMR_Iris\Output\Region_Income_vs_Classes_MDR_XDR_virulence.xlsx"

os.makedirs(os.path.dirname(output_file), exist_ok=True)

# =====================
# LOAD DATA
# =====================

df = pd.read_excel(input_file)

# =====================
# COLUMNS OF INTEREST
# =====================

class_cols = [
    "Aminoglycosides",
    "Carbapenems",
    "Third-generation cephalosporins",
    "Polymyxins",
    "Fosfomycins",
    "Fluoroquinolones",
    "Macrolides",
    "Phenicols",
    "Rifampicin",
    "Trimethoprim + Sulfamethoxazoles",
    "Tetracyclines",
    "Tigecycline",
    "MDR",
    "XDR",
    "Potential for Virulence",
    "High potential for virulence",
    "Potential for hypermucoidity",
    "Potential for hypervirulence",
    "Convergence of MDR and hypermucoidity",
    "Convergence of MDR and hypervirulence"
]

# =====================
# ENSURE BOOLEAN
# =====================

for col in class_cols:
    df[col] = df[col].astype(bool)

# =====================
# COUNT POSITIVE CASES BY REGION
# =====================

region_counts = df.groupby("region")[class_cols].sum().reset_index()

region_total = df.groupby("region").size().reset_index(name="total_isolates")

region_counts = region_total.merge(region_counts, on="region", how="left")

# =====================
# COUNT POSITIVE CASES BY INCOME
# =====================

income_counts = df.groupby("income_group")[class_cols].sum().reset_index()

income_total = df.groupby("income_group").size().reset_index(name="total_isolates")

income_counts = income_total.merge(income_counts, on="income_group", how="left")

# =====================
# OPTIONAL: LONG FORMAT
# =====================

region_long = region_counts.melt(
    id_vars=["region", "total_isolates"],
    value_vars=class_cols,
    var_name="class",
    value_name="positive_cases"
)

income_long = income_counts.melt(
    id_vars=["income_group", "total_isolates"],
    value_vars=class_cols,
    var_name="class",
    value_name="positive_cases"
)

# =====================
# SAVE TO EXCEL
# =====================

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    region_counts.to_excel(writer, sheet_name="region_counts_wide", index=False)
    income_counts.to_excel(writer, sheet_name="income_counts_wide", index=False)
    region_long.to_excel(writer, sheet_name="region_counts_long", index=False)
    income_long.to_excel(writer, sheet_name="income_counts_long", index=False)

print(f"Saved: {output_file}")