# --------------------------------------------------------------------------------------------------
# cleaning of meta data
#
# author : cat eisenhauer
# date : january 2026
# --------------------------------------------------------------------------------------------------

library(tidyverse)



# import joinables ---------------------------------------------------------------------------------
regions <- rio::import(here::here('data', 'ref', 'world.rds')) |> 
  as.data.frame() |>
  select(country = 'name', region = 'continent') |>
  mutate(
    region = case_when(
      country %in% c('United States', 'Canada') ~ 'US / Canada',
      region %in% c('North America', 'South America') ~ 'Latin America and the Caribbean',
      .default = region
    )
  )

incomes <- rio::import(here::here('data', 'ref', 'WBD_clean.xlsx')) |>
  select(country = 'Economy', income = 'Income group')

iso_cats <- rio::import(here::here('data', 'ref', 'isolation-source_categorized.csv')) |>
  select(iso2, iso2_type) |>
  mutate(
    iso2_type = case_when(
      iso2_type == 'environmental' ~ 'Environment',
      iso2_type == 'human' ~ 'Human',
      iso2_type == 'animal' ~ 'Animal',
      iso2_type == 'na' ~ NA
    )
  )

  


# import main metadata and clean -------------------------------------------------------------------
meta <- rio::import(here::here('data', 'ref', 'metadados.csv')) |>
  select(
    # id for joining
    strain = 'displayname',
    # for sample collection year
    year, 
    # for sample type
    base_type = 'sample type',
    host1 = 'host',
    host2 = 'Host',
    iso1 = 'isolation source',
    iso2 = 'Isolation source',
    env1 = 'environmental sample',
    env2 = 'Environmental sample',
    # for country
    country1 = 'country',
    country2 = 'Country',
    state = 'usa state (if applicable)'
  ) |>
  mutate(
    # sample type ----------------------------------------------------------------------------------
    sample_type = case_when(
      # human 
      host1 == 'Human' | 
      host2 %in% c('Human', 'human', 'patient', 'Human a', 'Human b', 'human being',
                   'homosapiens', 'humans') |
      iso1 != ''
      ~ 'Human',

      # animal
      host2 %in% c('Phocarctos hookeri', 'Porcine', 'chicken', 'Bos taurus', 'dog',
                   'cat', 'Primate', 'Mus musculus', 'Avian', 'Equine', 'Canis lupus familiaris',
                   'Chickens', 'horse', 'bovine', 'Meleagris gallopavo', 'gallus gallus',
                   'cattle', 'Canine', 'Stercorarius antarcticus', 'Megadyptes antipodes',
                   'Chicken', 'Swine', 'Reptiles', 'Mytilus edulis', 'Sus scrofa',
                   'Pteropus poliocephalus', 'Anopheline darlingi', 'Anastrepha ludens',
                   'Phoca vitulina', 'Sus domesticus', 'bird', 'Sus scrofa domesticus',
                   'feral swine', 'Felis catus', 'Tegillarca granosa', 'Misgurnus anguillicaudatus',
                   'Chroicocephalus ridibundus', 'Dog', 'Cat', 'Rodentia', 'duck', 'Gallus gallus',
                   'swine', 'Equus caballus', 'livestock', 'animal', 'Tilapia',
                   'Gallus gallus domesticus')
      ~ 'Animal',

      # environment
      host2 %in% c('Water', 'Environmental', 'Environment', 'wastewater', 'environmental',
                   'environment', 'hosptial surface', 'well water', 'Food', 'Milk powder') |
      base_type == 'Environmental' 
      ~ 'Environment',

      # other
      host2 %in% c('nan', 'missing', 'Lab', 'Laboratory', 'Lab strain') ~ 'Unknown',

      # everything else is unknown
      .default = NA
    ),

    # country --------------------------------------------------------------------------------------
    country = case_when(
      country1 != '' ~ country1,
      country2 != '' ~ country2,
      .default = NA
    ),
    country = case_when(
      country == 'Other (International Space Station)' ~ NA,
      country == 'Viet Nam' ~ 'Vietnam',
      country == 'USA' ~ 'United States',
      country %in% c('UK', 'United Kingdom', 'United Kingdom (England/Wales/N. Ireland)',
                     'United Kingdom (Scotland)') ~ 'United Kingdom',
      country == 'UAE' ~ 'United Arab Emirates',
      country == 'Gambia' ~ 'The Gambia',
      country == 'missing' ~ NA,
      .default = country
    )
  ) |>
  left_join(regions) |>
  left_join(incomes) |>
  left_join(iso_cats) |>
  mutate(
    sample_type = case_when(
      is.na(sample_type) & !is.na(iso2_type) ~ iso2_type,
      .default = sample_type
    ),
    region = case_when(
      country %in% c('South Korea', 'Russia', 'KSA', 'Laos', 'West Bank') ~ 'Asia',
      country == 'Macedonia' ~ 'Europe',
      country %in% c('Curacao', 'Martinique', 'Guadeloupe') ~ 'Latin America and the Caribbean',
      country == 'Gambia' ~ 'Africa',
      .default = region
    ),
    income = case_when(
      country %in% c('Curacao', 'Czech Republic', 'KSA', 'Guadeloupe', 'Hong Kong',
                     'Martinique', 'Russia', 'Slovakia', 'South Korea', 'Taiwan') ~ 'High income',
      country %in% c('Iran', 'Macedonia', 'Turkey') ~ 'Upper middle income',
      country %in% c('Egypt', 'Laos', 'West Bank') ~ 'Lower middle income',
      country %in% c('The Gambia') ~ 'Low income',
      country %in% c('Ethiopia', 'Venezuela') ~ 'Unclassified',
      .default = income
    )
  ) |>

  # remove full duplicates
  distinct() |>

  # select only final vars
  select(strain, year, sample_type, country, region, income)


# some strains have duplicated meta-data -- remove them entirely
dupes <- meta |>
  count(strain) |>
  filter(n > 1) |>
  pull(strain)

meta <- meta |>
  filter(!(strain %in% dupes))

# export -------------------------------------------------------------------------------------------
meta |>
  rio::export(here::here('data', 'ref', 'meta_clean.rds'))
