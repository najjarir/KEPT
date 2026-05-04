# --------------------------------------------------------------------------------------------------
# cleaning of raw data
#
# author : cat eisenhauer
# date : january 2026
# --------------------------------------------------------------------------------------------------

library(tidyverse)


# import -------------------------------------------------------------------------------------------
meta <- rio::import(here::here('data', 'ref', 'meta_clean.rds'))

df_raw <- data.frame()
for (i in 1:2) {
  tmp <- rio::import(here::here('data', 'raw', paste0('48k_parte', i, '.csv')), 
                    format = 'tsv',
                    header = FALSE,
                    blank.lines.skip = TRUE,
                    comment.char = '',
                    fill = TRUE) |>
    filter(!grepl('^strain', V1))

  colnames(tmp) <- scan(here::here('data', 'raw', paste0('48k_parte', i, '.csv')),
                       what = character(), 
                       nlines = 1, 
                       sep = "\t", 
                       quiet = TRUE)

  df_raw <- bind_rows(df_raw, tmp)
}

df_raw <- left_join(df_raw, meta)



# wrangle ------------------------------------------------------------------------------------------

df <- df_raw |>
  mutate(
    total_size = as.numeric(total_size),
    n_ambiguous = case_when(
      ambiguous_bases == 'no' ~ 0,
      .default = as.numeric(str_extract(ambiguous_bases, "\\d+"))
     )
  ) |>

  # filter by exclusion criteria
  filter(
    n_ambiguous < 40000,            # 2 removed
    !is.na(year),                   # 5728 removed
    !is.na(country),                # 524 removed
  ) |>
  filter(
    year >= 2000 & year != 2023     # 518 removed
  ) |>

  # cleanup genes
  mutate(

    # remove problem genes (%, ?, $) and remove [^, *] characters and harmonize NA style
    across(ends_with('_acquired'), remove_problem_genes),
    across(ends_with('_mutations'), ~na_if(., "-")),
    Bla_chr = remove_problem_genes(Bla_chr),

    # clean tetR and Mrx
    Tet_acquired = exclude_genes(Tet_acquired, 'tetR'),
    MLS_acquired = exclude_genes(MLS_acquired, 'Mrx'),

    # create gene clusters
    Tgc_acquired = create_gene_cluster(Tgc_acquired,
                                       cluster_genes = c('tmexC1', 'tmexD1', 'toprJ1'),
                                       cluster_name = 'tmexCD1-toprJ1'),

    # add resistance columns for antibiotics dependant on multiple enzymees
    monobactam = !is.na(Bla_ESBL_acquired) | !is.na(Bla_ESBL_inhR_acquired) |
                 grepl('^KPC|^GES', Bla_Carb_acquired),
    n_mdr = as.integer(!is.na(AGly_acquired)) + 
            as.integer(!is.na(Flq_acquired) | !is.na(Flq_mutations)) +
            as.integer(!is.na(Sul_acquired) | !is.na(Tmt_acquired)) + 
            as.integer(!is.na(Tgc_acquired)) +
            as.integer(!is.na(Phe_acquired)) + 
            as.integer(!is.na(Tet_acquired)) + 
            as.integer(!is.na(Fcyn_acquired)) +
            as.integer(!is.na(Col_acquired) | !is.na(Col_mutations)),
    MDR_betalactamase = !is.na(Bla_ESBL_acquired) | !is.na(Bla_ESBL_inhR_acquired) |
                        !is.na(Bla_Carb_acquired),
    MDR_wo_betalactamase = !MDR_betalactamase &
                           n_mdr >= 3,
    MDR = MDR_betalactamase | MDR_wo_betalactamase,
    XDR = ((!is.na(Bla_ESBL_acquired) | !is.na(Bla_ESBL_inhR_acquired)) & 
           is.na(Bla_Carb_acquired) & 
           n_mdr >= 7) |
          (!is.na(Bla_Carb_acquired) & (as.integer(monobactam) + n_mdr) >= 7),

    # fix K_locus and O_locus
    K_locus = case_when(K_locus_confidence == 'Untypeable' ~ 'Untypeable',
                        .default = K_locus),
    K_type = case_when(K_locus_confidence == 'Untypeable' ~ 'Untypeable',
                        .default = K_type),
    O_locus = case_when(O_locus_confidence == 'Untypeable' ~ 'Untypeable',
                        .default = O_locus),
    O_type = case_when(O_locus_confidence == 'Untypeable' ~ 'Untypeable',
                        .default = O_type),

    # remove truncation on rmpA(2)
    across(c(rmpA, rmpA2), mark_truncated_genes),
         
    # remove versions (.v1 / .v2) 
    across(c(AGly_acquired, Flq_acquired, MLS_acquired, Phe_acquired, Tet_acquired, 
            Tmt_acquired, Bla_acquired, Bla_Carb_acquired, Bla_ESBL_acquired),
          ~ gsub("\\.v[12]", "", .)),

    # remove -[1-9]+LV tail from ST
    ST = gsub('-\\d+LV', '', ST)
  ) |>

  # remove unnecessary columns
  select(-species_match, -starts_with('clb'), -starts_with('iro'), -starts_with('iuc'),
         -starts_with('ybt'), -ends_with('hits'), -contig_count, -largest_contig, -gapA, 
         -infB, -mdh, -pgi, -phoE, -rpoB, -tonB, -fyuA, -irp2, -irp1, -iutA, n_mdr, n_ambiguous)



# export -------------------------------------------------------------------------------------------
df |>
  rio::export(here::here('data', 'clean', 'data_clean.rds'))

df |>
  rio::export(here::here('data', 'clean', 'data_clean.csv'))

