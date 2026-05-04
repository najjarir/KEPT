# Script KEPT dataset analysis
# Author : Iris Najjar
# Date : 22.10.2025

library(tidyverse)
library(rio)
library(here)


#import--------------------------------------------------------------------------------------------
df <- import(here("data" , "data_clean_Cat_2_CRISPR.csv"))

n_total <- nrow(df)

#cleaning---------------------------------------------------------------------------------------------
# Reassign "" to Unknown in sample_type
df <- df |>
mutate(sample_type = case_when(sample_type == "" ~ "Unknown", TRUE ~ sample_type)) 

df |>
  count(sample_type)

# Reassign Ethiopia and venezuela to correct income groups

df <- df |>
mutate(income = case_when(country == "Venezuela" ~ "Upper middle income", country == "Ethiopia" ~ "Low income", TRUE ~ income))

df |>
  filter(country == "Ethiopia" | country == "Venezuela") |>
  select(country, income)

df |>
  count(income) |>
  mutate(percent = n / n_total * 100)

# Create new AMR columns ESBL_all_acquired, Tmt_sul_acquired, Col_mut_acquired, Flq_mut_acquired

df <- df |>
  mutate(
    ESBL_all_acquired = ifelse(Bla_ESBL_acquired != "" | Bla_ESBL_inhR_acquired != "", TRUE, "")) 

df|>
  count(ESBL_all_acquired)

df <- df |>
  mutate(
    Tmt_Sul_acquired = ifelse(Tmt_acquired != "" | Sul_acquired != "", TRUE, "")) 

df|>
  count(Tmt_Sul_acquired)

df <- df |>
  mutate(
    Col_mut_acquired = ifelse(Col_acquired != "" | Col_mutations != "", TRUE, "")) 

df|>
  count(Col_mut_acquired)

df <- df |>
  mutate(
    Flq_mut_acquired = ifelse(Flq_acquired != "" | Flq_mutations != "", TRUE, "")) 

df|>
  count(Flq_mut_acquired)

#analyze----------------------------------------------------------------------------------------------

# Dataset description-----------------------------------------------------------------------------------

# number of KL types

df <- df |> 
mutate(K_new = case_when(K_locus_confidence == "Untypeable" ~ "untypeable",
                          .default = K_locus))

KL_percentage <- df |>
  count(K_new) |>
  arrange(desc(n)) |>
  mutate(prop = n / n_total * 100)

KL_percentage$K_new[15]

df |>
  filter(K_new == "KL107") |>
  
  count(ST) |>
    mutate(KL107_ST = n / 100)

df |>
  pull(K_new) |>
  unique()

# Distribution of ST
df |> 
  count(ST) |>
  mutate(percent = n / n_total * 100) |>
  arrange(desc(n)) |>
  head(20)

# Distribution of ST by region

df |>
  count(region, ST) |>
  group_by(region) |>
  mutate(percent = n / sum(n) * 100) |>
  arrange(desc(n), .by_group = TRUE) |>
  slice_head(n = 20) |>
  print(n=120) 

df |>
  count(region, ST) |>
  group_by(region) |>
    summarise(total_percent = sum(n / sum(n) * 100 ))


# Distribution of ST by sample type
df |> 
  filter(sample_type == "Animal") |>
  count(ST) |>
  mutate(percent = n / n_total * 100) |>
  arrange(desc(n)) |>
  head(20)

df |> 
  filter(sample_type == "Human") |>
  count(ST) |>
  mutate(percent = n / n_total * 100) |>
  arrange(desc(n)) |>
  head(20)

df |> 
  filter(sample_type == "Environment") |>
  count(ST) |>
  mutate(percent = n / n_total * 100) |>
  arrange(desc(n)) |>
  head(20)

# unique ST
df |>
  pull(ST) |>
  unique() |>
  length()

# Mean and Median od ST
df |>
  count(ST) |>
  summarise(mean_count = mean(n), median_count = median(n), sd_count =sd (n))

# Mean and Median od ST by region
df |>
  count(ST,region) |>
  group_by(region) |>
  summarise(mean_count = mean(n), median_count = median(n), sd_count =sd (n))

# unique ST by sample type
df |>
  filter(sample_type == "Animal") |>
  pull(ST) |>
  unique() |>
  length()

df |>
  filter(sample_type == "Human") |>
  pull(ST) |>
  unique() |>
  length()

df |>
  filter(sample_type == "Environment") |>
  pull(ST) |>
  unique() |>
  length()


# Distribution by region
df |> 
  count(region) |>
    mutate(percent = n / n_total * 100)

# Count the number of countries and percentage
df |>
  count(country) |>
  mutate(percent = n / n_total * 100)


# Distribution by sample type 
df |> 
  count(sample_type) |>
    mutate(percent = n / n_total * 100)

# Distribution by year
df |> 
  count(year) |>
    mutate(percent = n / n_total * 100)

# AMR description----------------------------------------------------------------------------------------
 # General dataset
df |>
  summarise(across(ends_with("_acquired"), ~ sum( . != ""), .names = "n_{.col}")) |>
   pivot_longer(everything(), names_to = "Resistance", values_to = "n resistant") |>
    print(n = Inf)

df |>
  summarise(across(ends_with("_acquired"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(everything(), names_to = "Resistance", values_to = "Percent resistant") |>
    print(n = Inf)

df |>
  summarise(across(ends_with("_mutations"), ~ sum( . != ""), .names = "n_{.col}")) |>
   pivot_longer(everything(), names_to = "Resistance_mutations", values_to = "n resistant") |>
    print(n = Inf)

df |>
  summarise(across(ends_with("_mutations"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(everything(), names_to = "Resistance_mutations", values_to = "Percent resistant") |>
    print(n = Inf)

 # AMR by region
df |>
  group_by(region) |>
  summarise(across(ends_with("_acquired"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-region, names_to = "Resistance", values_to = "n resistant") |>
  print(n = Inf)

df |>
  group_by(region) |>
  summarise(across(ends_with("_acquired"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-region, names_to = "Resistance", values_to = "Percent resistant") |>
   print(n = Inf)
 
df |>
  group_by(region) |>
  summarise(across(ends_with("_mutations"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-region, names_to = "Resistance_mutations", values_to = "n resistant") |>
  print(n = Inf)

df |>
  group_by(region) |>
  summarise(across(ends_with("_mutations"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-region, names_to = "Resistance_mutations", values_to = "Percent resistant") |>
   print(n = Inf)

 # AMR by Income group
df |>
  group_by(income) |>
  summarise(across(ends_with("_acquired"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-income, names_to = "Resistance", values_to = "n resistant") |>
  print(n = Inf)

df |>
  group_by(income) |>
  summarise(across(ends_with("_acquired"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-income, names_to = "Resistance", values_to = "Percent resistant") |>
   print(n = Inf)
 
df |>
  group_by(income) |>
  summarise(across(ends_with("_mutations"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-income, names_to = "Resistance_mutations", values_to = "n resistant") |>
  print(n = Inf)

df |>
  group_by(income) |>
  summarise(across(ends_with("_mutations"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-income, names_to = "Resistance_mutations", values_to = "Percent resistant") |>
   print(n = Inf)
 
 # AMR by sample type
df |>
  group_by(sample_type) |>
  summarise(across(ends_with("_acquired"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-sample_type, names_to = "Resistance", values_to = "n resistant") |>
  print(n = Inf)

df |>
  group_by(sample_type) |>
  summarise(across(ends_with("_acquired"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-sample_type, names_to = "Resistance", values_to = "Percent resistant") |>
   print(n = Inf)
 
df |>
  group_by(sample_type) |>
  summarise(across(ends_with("_mutations"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-sample_type, names_to = "Resistance_mutations", values_to = "n resistant") |>
  print(n = Inf)

df |>
  group_by(sample_type) |>
  summarise(across(ends_with("_mutations"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-sample_type, names_to = "Resistance_mutations", values_to = "Percent resistant") |>
   print(n = Inf)

 # AMR by sequence type
df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
  summarise(across(ends_with("_acquired"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-ST, names_to = "Resistance", values_to = "n resistant") |>
  print(n = Inf)

df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
  summarise(across(ends_with("_acquired"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-ST, names_to = "Resistance", values_to = "Percent resistant") |>
   print(n = Inf)
 
df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
  summarise(across(ends_with("_mutations"), ~ sum(. != ""), .names = "n_{.col}")) |>
   pivot_longer(-ST, names_to = "Resistance_mutations", values_to = "n resistant") |>
  print(n = Inf)

df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
  summarise(across(ends_with("_mutations"), ~ mean( . != "") * 100, .names = "percent_resistant_{.col}")) |>
   pivot_longer(-ST, names_to = "Resistance_mutations", values_to = "Percent resistant") |>
   print(n = Inf)

 # MDR and XDR ------------------------------------------------------------------------------------------
df |>
   summarise(n_total = n(), n_MDR = sum(MDR, na.rm = TRUE), percent_MDR = 100 * n_MDR / n_total)

df |>
   summarise(n_total = n(), n_XDR = sum(XDR, na.rm = TRUE), percent_XDR = 100 * n_XDR / n_total)

# MDR and XDR by income group
df |>
  group_by(income) |>
 summarise(n_total = n(), n_MDR = sum(MDR, na.rm = TRUE), percent_MDR = 100 * n_MDR / n_total)

df |>
  group_by(income) |>
 summarise(n_total = n(), n_XDR = sum(XDR, na.rm = TRUE), percent_XDR = 100 * n_XDR / n_total)

# MDR and XDR by region
df |>
  group_by(region) |>
 summarise(n_total = n(), n_MDR = sum(MDR, na.rm = TRUE), percent_MDR = 100 * n_MDR / n_total)

df |>
  group_by(region) |>
 summarise(n_total = n(), n_XDR = sum(XDR, na.rm = TRUE), percent_XDR = 100 * n_XDR / n_total)

# MDR and XDR by sample type
df |>
  group_by(sample_type) |>
 summarise(n_total = n(), n_MDR = sum(MDR, na.rm = TRUE), percent_MDR = 100 * n_MDR / n_total)

df |>
  group_by(sample_type) |>
 summarise(n_total = n(), n_XDR = sum(XDR, na.rm = TRUE), percent_XDR = 100 * n_XDR / n_total)

# MDR and XDR by ST
df |>
   filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
 summarise(n_total = n(), n_MDR = sum(MDR, na.rm = TRUE), percent_MDR = 100 * n_MDR / n_total)

df |>
    filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
 summarise(n_total = n(), n_XDR = sum(XDR, na.rm = TRUE), percent_XDR = 100 * n_XDR / n_total)
 
# MDR by ST
df |>
  group_by(ST) |>
  summarise(n_total = n(), n_MDR = sum(MDR, na.rm = TRUE), percent_MDR = 100 * n_MDR / n_total) |>
    arrange(desc(percent_MDR)) |>
  head(20)

# Mean number of ARG-------------------------------------------------------------------------------------

df |>
  count(n_mdr) |>
  summarise(
    mean_n_mdr = sum(n_mdr * n) / sum(n))

# mean by region
df |>
  group_by(region) |>
  count(n_mdr) |>
  summarise(mean_n_mdr = sum(n_mdr * n) / sum(n))

# mean by income group
df |>
  group_by(income) |>
  count(n_mdr) |>
  summarise(mean_n_mdr = sum(n_mdr * n) / sum(n))

# mean by sample type
df |>
  group_by(sample_type) |>
  count(n_mdr) |>
  summarise(mean_n_mdr = sum(n_mdr * n) / sum(n))

# mean by ST
df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  count(ST, n_mdr) |>
  group_by(ST) |>
  summarise(mean_n_mdr = sum(n_mdr * n) / sum(n))

# Virulence --------------------------------------------------------------------------------------------
df |>
  mutate(virulence = ifelse(virulence_score > 0, "Yes", "No")) |>
  count(virulence) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

df |>
  mutate(virulence = ifelse(virulence_score > 2, "Yes", "No")) |>
  count(virulence) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

# Virulence by region
df |>
  group_by(region) |>
  mutate(virulence = ifelse(virulence_score > 0, "Yes", "No")) |>
  count(virulence) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

# Virulence by income
df |>
  group_by(income) |>
  mutate(virulence = ifelse(virulence_score > 0, "Yes", "No")) |>
  count(virulence) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

# Virulence by sample_type
df |>
  group_by(sample_type) |>
  mutate(virulence = ifelse(virulence_score > 0, "Yes", "No")) |>
  count(virulence) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

# Virulence by ST
df |>
    filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
  group_by(ST) |>
  mutate(virulence = ifelse(virulence_score > 0, "Yes", "No")) |>
  count(virulence) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity----------------------------------------------------------------------------------------

df |>
  mutate(hypermucoidity = if_else(!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated")), "Yes", "No")) |>
    count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

df |>
  mutate(hypermucoidity = if_else(!(rmpA %in% c("-", "truncated")), "Yes", "No")) |>
  count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

df |>
  mutate(hypermucoidity = if_else(!(rmpA2 %in% c("-", "truncated")), "Yes", "No")) |>
  count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity by region
df |>
  mutate(hypermucoidity = if_else(!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated")), "Yes", "No")) |>
  group_by(region) |>
    count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity by income group
df |>
  mutate(hypermucoidity = if_else(!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated")), "Yes", "No")) |>
  group_by(income) |>
    count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity by sample type
df |>
  mutate(hypermucoidity = if_else(!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated")), "Yes", "No")) |>
  group_by(sample_type) |>
    count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity by ST
df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
      mutate(hypermucoidity = if_else(!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated")), "Yes", "No")) |>
    group_by(ST) |>
    count(hypermucoidity) |>
    mutate(n_total = sum(n), percent = n / n_total * 100) |>
  head(20)

# Hypermucoidity AND virulence --------------------------------------------------------------------------

df |>
  mutate(virulent_mucoid = if_else((!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated"))) & virulence_score > 0, "Yes", "No")) |>
  count(virulent_mucoid) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

df |> # among hypermucoid samples, how many are virulent
  filter(!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated"))) |>
  mutate(virulent = ifelse(virulence_score > 0, "Yes", "No" )) |>
  count(virulent) |>
   mutate(n_total = sum(n), percent = n / n_total * 100)


# Hypermucoidity AND virulence by region

  df |>
    group_by(region) |>
  mutate(virulent_mucoid = if_else((!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated"))) & virulence_score > 0, "Yes", "No")) |>
  count(virulent_mucoid) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity AND virulence by income

  df |>
    group_by(income) |>
  mutate(virulent_mucoid = if_else((!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated"))) & virulence_score > 0, "Yes", "No")) |>
  count(virulent_mucoid) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity AND virulence by sample type

  df |>
    group_by(sample_type) |>
  mutate(virulent_mucoid = if_else((!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated"))) & virulence_score > 0, "Yes", "No")) |>
  count(virulent_mucoid) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Hypermucoidity AND virulence by ST

  df |>
  filter(ST %in% c("ST258", "ST11", "ST147", "ST307", "ST15")) |>
   mutate(virulent_mucoid = if_else((!(rmpA %in% c("-", "truncated")) 
   | !(rmpA2 %in% c("-", "truncated"))) & virulence_score > 0, "Yes", "No")) |>
    group_by(ST) |>
  count(virulent_mucoid) |>
    mutate(n_total = sum(n), percent = n / n_total * 100)

# Virulence AND MDR--------------------------------------------------------------------------

df |>
  summarise(n_MDR_virulent = sum(MDR == "true" & virulence_score > 2, na.rm = TRUE)) |>
    count(n_MDR_virulent) |>
  mutate(n_total = sum(n), percent = n / n_total * 100)

# Carbapenemase and ESBL --------------------------------------------------------------------------
df |>
  separate_rows(Bla_Carb_acquired, sep = ";") |>
  distinct(strain, Bla_Carb_acquired, .keep_all = TRUE) |>
  filter(Bla_Carb_acquired == "KPC-3") |>
  group_by(ST) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(n_total = sum(n), percent = n / n_total * 100) |>
  arrange(desc(n)) |>
print(n=120)

df |>
  separate_rows(Bla_Carb_acquired, sep = ";") |>
  distinct(strain, Bla_Carb_acquired, .keep_all = TRUE) |>
  filter(Bla_Carb_acquired == "KPC-2") |>
  group_by(ST) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(n_total = sum(n), percent = n / n_total * 100) |>
  arrange(desc(n)) |>
print(n=120)

df |>
  separate_rows(Bla_Carb_acquired, sep = ";") |>
  distinct(strain, Bla_Carb_acquired, .keep_all = TRUE) |>
  filter(Bla_Carb_acquired == "NDM-1") |>
  group_by(ST) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(n_total = sum(n), percent = n / n_total * 100) |>
  arrange(desc(n)) |>
print(n=120)

df |>
  separate_rows(Bla_Carb_acquired, sep = ";") |>
  distinct(strain, Bla_Carb_acquired, .keep_all = TRUE) |>
  filter(Bla_Carb_acquired == "NDM-5") |>
  group_by(ST) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(n_total = sum(n), percent = n / n_total * 100) |>
  arrange(desc(n)) |>
print(n=120)

df |>
  separate_rows(Bla_Carb_acquired, sep = ";") |>
  distinct(strain, Bla_Carb_acquired, .keep_all = TRUE) |>
  filter(Bla_Carb_acquired == "OXA-48") |>
  group_by(ST) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(n_total = sum(n), percent = n / n_total * 100) |>
  arrange(desc(n)) |>
print(n=120)


df |>
  separate_rows(Bla_ESBL_acquired, sep = ";") |>
  distinct(strain, Bla_ESBL_acquired, .keep_all = TRUE) |>
  filter(Bla_ESBL_acquired == "CTX-M-15") |>
  group_by(ST) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(n_total = sum(n), percent = n / n_total * 100) |>
  arrange(desc(n)) |>
print(n=120)

df |>
  filter(str_detect(Bla_ESBL_acquired, "CTX-M-15")) |>
  summarise(missing_ST = sum(is.na(ST)))

df |>
  count(Bla_ESBL_inhR_acquired) |>
  mutate(n/n_total*100)

# ARG ------------------------------------------------------------------------------------------------
# Bla_Chr
df |>
  separate_rows(Bla_chr, sep = ";") |>
  distinct(strain, Bla_chr) |>
  group_by(Bla_chr) |>
  summarise(n_isolates = n(), .groups = "drop") |>
  mutate(n_total = n_distinct(df$strain), percent = n_isolates / n_total * 100) |>
  group_by(region)
  arrange(desc(n_isolates)) |>
print(n=120)

df |>
  separate_rows(Bla_chr, sep = ";") |>
  distinct(strain, region, Bla_chr) |>
  group_by(region, Bla_chr) |>
  summarise(n_isolates = n(), .groups = "drop") |>
  group_by(region) |>
  mutate(
    n_total = sum(n_isolates),
    percent = n_isolates / n_total * 100) |>
  arrange(region, desc(n_isolates)) |>
  print(n = 500)

# Bla_acquired
df |>
  separate_rows(Bla_acquired, sep = ";") |>
  distinct(strain, Bla_acquired) |>
  group_by(Bla_acquired) |>
  summarise(n_isolates = n(), .groups = "drop") |>
  mutate(n_total = n_distinct(df$strain), percent = n_isolates / n_total * 100) |>
  arrange(desc(n_isolates)) |>
print(n=120)

df |>
  separate_rows(Bla_acquired, sep = ";") |>
  distinct(strain, region, Bla_acquired) |>
  group_by(region, Bla_acquired) |>
  summarise(n_isolates = n(), .groups = "drop") |>
  group_by(region) |>
  mutate(
    n_total = sum(n_isolates),
    percent = n_isolates / n_total * 100) |>
  arrange(region, desc(n_isolates)) |>
  print(n = 500)

df |>
  separate_rows(Bla_inhR_acquired, sep = ";") |>
  distinct(strain, Bla_inhR_acquired) |>
  group_by(Bla_inhR_acquired) |>
  summarise(n_isolates = n(), .groups = "drop") |>
  mutate(n_total = n_distinct(df$strain), percent = n_isolates / n_total * 100) |>
  arrange(desc(n_isolates)) |>
print(n=120)

df |>
  separate_rows(Bla_inhR_acquired, sep = ";") |>
  distinct(strain, region, Bla_inhR_acquired) |>
  group_by(region, Bla_inhR_acquired) |>
  summarise(n_isolates = n(), .groups = "drop") |>
  group_by(region) |>
  mutate(
    n_total = sum(n_isolates),
    percent = n_isolates / n_total * 100) |>
  arrange(region, desc(n_isolates)) |>
  print(n = 500)

# Colistin total resistance

df |>
  summarise(
    n_total = n(),
    n_colR = sum(
      (!is.na(Col_mutations) & Col_mutations != "") |
      (!is.na(Col_acquired) & Col_acquired != "")
    ), percent_colR = n_colR / n_total * 100)

# Fluoroquinolone total resistance

df |>
  summarise(
    n_total = n(),
    n_flqR = sum(
      (!is.na(Flq_mutations) & Flq_mutations != "") |
      (!is.na(Flq_acquired) & Flq_acquired != "")
    ), percent_flqR = n_flqR / n_total * 100)

# Capsule anaylsis ----------------------------------------------------------------------------------

df |> 
  count(K_locus) |>
    mutate(percent = n / n_total * 100) |>
   arrange(desc(n)) |>
  print(all())

# O type anaylsis ----------------------------------------------------------------------------------

df |> 
  count(O_type) |>
    mutate(percent = n / n_total * 100) |>
   arrange(desc(n)) |>
  print(all())

# CRISPR CAS Analysis -------------------------------------------------------------------------------

df |>
  dplyr::distinct(Cas_Types)


#visualize----------------------------
# Visualize ST distribution by region

df_top20_region  <- df |>
  count(region, ST) |>
  group_by(region) |>
  arrange(desc(n), .by_group = TRUE) |>
  slice_head(n = 20) |>
  mutate(percent = n / sum(n) * 100, ST_reordered = reorder_within(ST, percent, region)) |>
  ungroup() 

df_top20_region <- df_top20_region |>
  group_by(region) |>
   mutate(ST = fct_reorder(ST, percent, .desc = TRUE)) |>
  ungroup() |>
  mutate(ST = fct_rev(ST))  

df_top20_region <- df_top20_region |>
  mutate(region = if_else(region == "Latin America and the Caribbean",
                          "Latin America\nand the Caribbean",
                          region))

ggplot(df_top20_region, aes(x = ST, y = percent, fill = ST)) +
  geom_col() + geom_text(aes(label = round(percent, 1)), hjust = -0.1, size = 1) +
  facet_wrap(~region, scales = "free_y") + coord_flip() +
  labs(
    x = "none",
    y = "Proportion of samples (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 2), axis.text.y = element_text(size = 2), strip.text = element_text(size = 6), legend.position = "none")

