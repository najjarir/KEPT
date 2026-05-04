# --------------------------------------------------------------------------------------------------
# side script to extract unique values from df
#
# author : cat eisenhauer
# date : january 2026
# --------------------------------------------------------------------------------------------------

library(tidyverse)

extract_unique_genes <- function(column) {
  column <- column[!is.na(column)]
  if (length(column) == 0) return(NA)
  
  # using strsplit with fixed=TRUE for speed
  split_genes <- unlist(strsplit(column, ";", fixed = TRUE))
  unique_genes <- sort(unique(split_genes))
  
  return(unique_genes)
}


# main ---------------------------------------------------------------------------------------------

df <- rio::import(here::here('data', 'clean', 'data_clean.rds'))

tmp <- df |>
  select(O_locus, O_type, K_locus, K_type, Yersiniabactin,  YbST, CbST, Colibactin, AbST, 
         Aerobactin, SmST, Salmochelin, RmST, RmpADC, rmpA, rmpA2, AGly_acquired, Col_acquired, 
         Fcyn_acquired, Flq_acquired, Gly_acquired, MLS_acquired, Phe_acquired, Rif_acquired, 
         Sul_acquired, Tet_acquired, Tmt_acquired, Bla_acquired, Bla_inhR_acquired, 
         Bla_ESBL_acquired, Bla_ESBL_inhR_acquired, Bla_Carb_acquired, Bla_chr, Col_mutations, 
         Flq_mutations)

# apply to all columns efficiently
gene_lists <- vector("list", ncol(tmp))
names(gene_lists) <- names(tmp)

for (i in seq_along(tmp)) {
  gene_lists[[i]] <- extract_unique_genes(tmp[[i]])
}

max_length <- max(lengths(gene_lists))

# predefine matrix for speed
result_matrix <- matrix(NA_character_, nrow = max_length, ncol = ncol(tmp))
colnames(result_matrix) <- names(tmp)

for (i in seq_along(gene_lists)) {
  len <- length(gene_lists[[i]])
  if (len > 0) {
    result_matrix[1:len, i] <- gene_lists[[i]]
  }
}

result <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
result |>
  rio::export(here::here('data', 'clean', 'unique_values.csv'))

