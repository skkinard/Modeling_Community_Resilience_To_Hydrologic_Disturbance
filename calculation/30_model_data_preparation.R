# model_data_preparation
# Sean Kinard
# 2023-07-04
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')

d_density <- read_csv("data/06_fill_data/fish_density_fill.csv")
d_biomass <- read_csv("data/06_fill_data/fish_biomass.csv")
d_diversity <- read_csv('exploration/output/fish_diversity_density.csv')
d_ordination <- read_csv('exploration/output/rda_all_extract_site.csv')
d_traits <- read_csv(
  "data/02_clean_data/fish_biological_scales_of_comparison.csv")
d_lte <- read_csv('exploration/output/environment_longterm_with_ln.csv')
d_ste <- read_csv('exploration/output/environment_shortterm_with_ln.csv')

# -----------------------------------------------------------------------------
# Data prep: Biological outcomes
# -----------------------------------------------------------------------------
# biological outcome variables
slim_density <- d_density %>% 
  group_by(site_code, collection_period) %>%
  dplyr:: summarize(density = sum(density, na.rm=T))

slim_biomass <- d_biomass %>% 
  group_by(site_code, collection_period) %>%
  dplyr:: summarize(biomass = sum(AFDMg_m2, na.rm=T))

slim_ordination <- d_ordination %>% 
  select(site_code, collection_period, centroid_dist)

# merge checks
anti_join(slim_density, slim_biomass) # missing biomass spring 2017
anti_join(slim_density, slim_ordination) # no missing
anti_join(slim_density, d_diversity) # no missing

# join biological outcomes
d_outcome <- slim_density %>% 
  left_join(slim_biomass) %>%
  left_join(d_diversity) %>%
  left_join(slim_ordination) %>%
  ungroup()

# impute missing biomass spring 2017 with baseline average
d_outcome <- d_outcome %>%
  mutate(qtr=quarter(collection_period)) %>%
  mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm = TRUE))), 
         .by = c(site_code, qtr))

# -----------------------------------------------------------------------------
# Data prep: Environmental Predictors
# -----------------------------------------------------------------------------
# fill missing environmental predictors with linear interpolation
d_ste_fill <- d_ste %>%
  rowwise() %>%
  mutate(total_algae = sum(sqrt_green_algae^2, 
                           sqrt_diatoms^2, sqrt_bluegreen_cyano^2)) %>%
  unique() %>%
  as_tsibble(key=site_code,
             index=collection_period) %>%  
  na_ma(weighting="linear") %>%
  ungroup() %>%
  as_tibble() %>%
  pivot_longer(cols=-any_of(c('site_code', 'collection_period')),
               names_to='xname',
               values_to = 'x') %>%
  group_by(site_code, collection_period, xname) %>%
  dplyr::summarize(x=mean(x ,na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from = x) %>%
  ungroup()

# merge biological outcomes with environmental predictors
d_predictor <- left_join(slim_density, d_ste_fill) %>%
  select(-density) %>%
  create_date_vars() %>%
  add_rain() %>%
  ungroup()

# check for merge errors
d_predictor %>% filter(is.na(sqrt_conductivity)) # missing GC

# impute site average to fill missing e-predictors at GC fall 2019
d_predictor <- d_predictor %>%
  mutate(qtr=quarter(collection_period)) %>%
  mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm = TRUE))), 
         .by = c(site_code, qtr))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
d_combined <- left_join(d_outcome, d_predictor) %>% 
  add_rain() %>%
  select(site_code, collection_period, year, qtr, month, day, nday, year_week, year_month, everything())

write_csv(d_outcome, 
          'exploration/output/model_outcome_vars.csv')
write_csv(d_predictor, 
          'exploration/output/model_predictor_vars.csv')
write_csv(d_combined, 
          'exploration/output/model_combined_vars.csv')

# -----------------------------------------------------------------------------
# End model_data_preparation