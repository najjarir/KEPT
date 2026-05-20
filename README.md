Repository with scripts used for temporal, descriptive, and visualization analyses of Klebsiella pneumoniae genomic data.

Scripts included:
- Clean_main.R : cleaning of data as described in the methods section of the article (quality filtering, removal of incomplete or truncated genes, exclusion of tetR and Mrx gene, tmexCD1-toprJ1 gene clustering);  creation of new columns (monobactam resistance, fluoroquinolone resistance, colistin resistance, trimethoprim-sulfamethoxazole reisstance, MDR, XDR); cleaning of K-locus and O-locus columns. (R)
- Clean_metadat.R: grouping into regions, income groups, sample type. (R)
- extract_unique.R: extraction of unique values for downstream amalysis. (R)
- Table1_data.R: generates data for Table_1. (R)
- package_code_logistic_regression_amr_ma3_region_income_for_classes_MDR_XDR_NMDR.py: generates temporal trend analyses and summary figures for 12 antimicrobial resistance classes, MDR, XDR, and N_MDR by region and income group. (Python)
- package_code_virulence_logistic_regression_ma3_region_income.py: generates temporal trend analyses for virulence-related outcomes by region and income group.(Python)
- data_trends_AMR_genes.py: builds isolate-level ARG presence/absence matrices and performs temporal trend analyses for selected ARGs by region and income group.(Python)
- code_association_ST_KL.py: evaluates regional associations between major sequence types and K-loci using logistic regression and descriptive overlap metrics.(Python)
- code_STs.py: generates temporal prevalence plots for sequence types, highlighting dominant and recently expanding STs.(Python)
- code_Sts_and_KL.py: generates packed-circle visualizations for relationships among STs, K-loci, carbapenemases, and ESBLs.(Python)
- code_treemap.py: generates treemap visualizations of ST distributions within selected carbapenemases and ESBLs.(Python)
- code_counts.py: produces descriptive count tables of resistance, MDR/XDR, and virulence outcomes by region and income group.(Python)
- analizyer.py: Calculation and graphic generation for sample distribution figures (region, income group, year). (Python)
- heatmaps: Multi-panel AMR gene prevalence heatmaps. (Python)
