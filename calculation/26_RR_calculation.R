# RR_calculation
# Sean Kinard
# 2023-07-01
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')

d_density <- read_csv("data/06_fill_data/fish_density_fill.csv")
d_biomass <- read_csv("data/06_fill_data/fish_biomass.csv")
d_diversity <- read_csv('exploration/output/fish_diversity_density.csv')
d_ordination <- read_csv('exploration/output/rda_all_extract_site.csv')
d_traits <- read_csv("data/02_clean_data/fish_biological_scales_of_comparison.csv")
d_lte <- read_csv('exploration/output/environment_longterm_with_ln.csv')
d_ste <- read_csv('exploration/output/environment_shortterm_with_ln.csv')

# -----------------------------------------------------------------------------
# Response Ratio Function
# -----------------------------------------------------------------------------
create_is_baseline <- function(x) {
  x %>%
    mutate(is_baseline=
             ifelse(collection_period > ymd('2017-06-01') &
                      collection_period < ymd("2018-09-01"), 'test',
                    'baseline')) }

calculate_RR <- function(my_data, my_var, my_level) { 
  # data prep 
  df <- my_data %>%
    create_date_vars() %>%
    mutate(com='community') %>%
    create_is_baseline() %>%
    dplyr::rename(XX={{my_var}}) %>% # XX=variable
    dplyr::rename(YY={{my_level}}) %>% # YY=level of comparison
    select(any_of(c('site_code', 'collection_period', 'qtr', 'is_baseline',
                    'XX', 'YY'))) %>%
    group_by(site_code, qtr, collection_period, YY, is_baseline) %>%
    dplyr::summarize(XX=sum(XX,na.rm=T)) %>%
    unique() %>%
    ungroup()
  
  # 2020 quarterly mean
  b_qtr <- df %>%
    filter(is_baseline == 'baseline') %>%
    group_by(site_code, qtr, YY) %>%
    dplyr::summarize(XX_b_qtr=mean(XX,na.rm=T)) %>%
    unique() %>%
    ungroup()
  
  # 2020 annual mean
  b_all <- df %>%
    filter(is_baseline == 'baseline') %>%
    group_by(site_code, YY) %>%
    dplyr::summarize(XX_b_all=mean(XX, na.rm=T),
                     XX_sd20=sd(XX, na.rm=T)) %>%
    ungroup() 
  
  # Response Ratios=(x-x_base) / (x_base)
  RR <- df %>%
    filter(is_baseline == 'test') %>%
    left_join(b_all) %>%
    mutate(XX_b_all=ifelse(is.na(XX_b_all), 0, XX_b_all)) %>%
    left_join(b_qtr) %>%
    mutate( # gapfill baselines with 2020 mean
      XX_b_qtr_fill=ifelse(is.na(XX_b_qtr), XX_b_all, XX_b_qtr)) %>%
    mutate( # calculate response ratio (RR)
      XX_RR_b_qtr=(XX-XX_b_qtr_fill) / (XX_b_qtr_fill + .01),
      XX_RR_b_all=(XX-XX_b_all) / (XX_b_all + .01)) %>%
    mutate( # calculate log response ratio (LRR)
      XX_LRR_b_qtr=log(abs(XX_RR_b_qtr)),
      XX_LRR_b_all=log(abs(XX_RR_b_all))) %>%
    mutate( # replace infinity with -100000 LRR
      XX_LRR_b_qtr=ifelse(is.infinite(XX_LRR_b_qtr), -100000, XX_LRR_b_qtr),
      XX_LRR_b_all=ifelse(is.infinite(XX_LRR_b_all), -100000, XX_LRR_b_all)) %>%
    ungroup() %>%
    fix_site_order()
  
  # Revert Column names
  colnames(RR) <- str_replace_all(colnames(RR), 'XX', substr(my_var,1,4))
  colnames(RR) <- str_replace_all(colnames(RR), 'YY', my_level)
  
  return(RR) }

# -----------------------------------------------------------------------------
# Calculate Response Ratios
# -----------------------------------------------------------------------------
# Abundance: biomass
RR_biomass_spe <- d_biomass %>%
  calculate_RR(my_var='AFDMg_m2', my_level='lowest_taxon')
RR_biomass_fam <- d_biomass %>% 
  left_join(d_traits) %>%
  calculate_RR(my_var='AFDMg_m2', my_level='family')
RR_biomass_com <- d_biomass %>%
  calculate_RR(my_var='AFDMg_m2', my_level='com')

# Abundance: density
RR_density_spe <- d_density %>%
  calculate_RR(my_var='density', my_level='lowest_taxon')
RR_density_fam <- d_density %>% 
  left_join(d_traits) %>%
  calculate_RR(my_var='density', my_level='family')
RR_density_com <- d_density %>% 
  calculate_RR(my_var='density', my_level='com')

# Diversity: richness, SW-index, Simpson-index
RR_shan_com <- d_diversity %>% 
  calculate_RR(my_var='shannon', my_level='com')
RR_simp_com <- d_diversity %>% 
  calculate_RR(my_var='simpson', my_level='com')
RR_rich_com <- d_diversity %>% 
  calculate_RR(my_var='richness', my_level='com')

# Composition: centroid-distance
RR_cdis_com <- d_ordination %>%
  create_is_baseline() %>%
  filter(is_baseline == 'baseline') %>%
  group_by(site_code) %>%
  dplyr::summarize(
    axis1_mu20 = mean(axis1, na.rm=T),
    axis2_mu20 = mean(axis2, na.rm=T) ) %>%
  ungroup() %>%
  right_join(d_ordination) %>%
  mutate(
    cdis = sqrt((axis1-axis1_mu20)^2 + 
                  (axis2-axis2_mu20)^2) ) %>%
  calculate_RR(my_var = 'cdis', 
               my_level = 'com')

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
my_objects <- ls()
my_csv_names <- my_objects[str_detect(my_objects, 'RR_')]

my_tables<-list(RR_biomass_com, RR_biomass_fam, RR_biomass_spe, RR_cdis_com, 
                RR_density_com, RR_density_fam, RR_density_spe, RR_rich_com, 
                RR_shan_com, RR_simp_com)

names(my_tables) <- my_csv_names

for (i in 1:length(my_tables)) {
  my_place <- paste('exploration/output/', names(my_tables[i]), ".csv", sep='')
  my_object <- my_tables[[i]]
  write_csv(my_object, my_place) }

# -----------------------------------------------------------------------------
# End RR_calculation