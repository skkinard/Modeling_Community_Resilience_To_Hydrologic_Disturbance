# 06_site_table
# Sean Kinard
# 2023-06-23

# summarize environmental variables
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

dl <- read_csv("exploration/output/environment_longterm.csv")
dw <- read_csv("data/02_clean_data/watershed.csv")
ds <- read_csv("exploration/output/environment_shortterm.csv")


pretty_titles <- function(df) {
  x <- df
  colnames(x) <- str_replace(colnames(x), '_', ' ')
  colnames(x) <- str_replace(colnames(x), 'land', ' land')
  colnames(x) <- str_replace(colnames(x), 'annual', 'annual ')
  colnames(x) <- str_replace(colnames(x), 'q q', 'q ')
  colnames(x) <- str_to_title(colnames(x))
  colnames(x) <- str_replace(colnames(x), 'Staid', 'USGS STAID')
  colnames(x) <- str_replace(colnames(x), 'Doc', 'DOC')
  colnames(x) <- str_replace(colnames(x), 'Ph', 'pH')
  colnames(x) <- str_replace(colnames(x), 'pHo', 'Pho')
  return(x) }

#------------------------------------------------------------------------------
# summarize: long term environment
#------------------------------------------------------------------------------
lte_watershed <- dl %>%
  select(site_code, any_of(location), any_of(climate), any_of(landuse)) %>%
  right_join(dw %>% select(site_code, staid) ) %>% 
  pretty_titles()
  
lte_flow <- dl %>% select(site_code, contains('q_')) %>%
  pretty_titles()
colnames(lte_flow) <- c('Site Code', 
                        colnames(select(lte_flow, contains('Q'))) %>%
  str_to_upper())

lte_water_quality <- dl %>% select(site_code, any_of(water_quality)) %>%
  pretty_titles()

lte_geomorph <- dl %>% select(site_code, any_of(geomorph)) %>% pretty_titles()

lte_algae <- dl %>% select(site_code, any_of(algae)) %>% pretty_titles()

lte_table_long <- left_join(lte_watershed, lte_flow) %>% 
  left_join(lte_water_quality) %>%
  left_join(lte_geomorph) %>%
  left_join(lte_algae) %>%
  select("USGS STAID", everything()) %>%
  column_to_rownames(var = "Site Code") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var='Variable') %>%
  as_tibble()

caption_lte_summary <- 'Table of long-term environmental variables for sampling locations. Values represent 20 year averages. Flow is average annual discharge, HFPP is the proportion of the annual discharge that is 3x higher than the Flow, LFPP is the proportion of discharge below the 25th percentile, Flashiness is the cumulative changes in day to day discharge divided by cumulative annual discharge, Season approximates the degree to which the flow varies during the course of a single year.'

# # Export lte Tables
# lte_summary_tables <- list(lte_table_long,
#                            lte_watershed,
#                            lte_flow,
#                            lte_water_quality,
#                            lte_geomorph,
#                            lte_algae)
# 
# names(lte_summary_tables) <- c('lte_table_long',
#                                'lte_watershed',
#                                'lte_flow',
#                                'lte_water_quality',
#                                'lte_geomorph',
#                                'lte_algae')
# 
# for (i in 1:length(lte_summary_tables)) {
#   my_place <- paste('exploration/output/', 
#                     names(lte_summary_tables[i]), 
#                     "_summary_table.csv", sep='')
#   my_object <- lte_summary_tables[[i]]
#   write_csv(my_object, my_place) }

#------------------------------------------------------------------------------
# Short Term Quarterly
#------------------------------------------------------------------------------
find_colmeans <- function(df) {
  df %>% dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) }

ste_qtr_long <- ds %>% 
  mutate(year_qtr = yearquarter(collection_period)) %>%
  select(-collection_period) %>%
  group_by(site_code, year_qtr) %>%
  nest() %>%
  mutate(qtr_mean = map(data, find_colmeans)) %>%
  select(-data) %>%
  unnest(qtr_mean) %>%
  arrange(year_qtr, site_code)

test <- ds %>% 
  mutate(year_qtr = yearquarter(collection_period)) %>%
  select(-collection_period) %>%
  group_by(site_code, year_qtr) %>%
  nest()

test$data[[10]] %>% 
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))



ste_qtr_flow <- ste_qtr_long %>% select(site_code, contains('q2wk')) %>%
  pretty_titles()
colnames(lte_flow) <- c('Site Code', 
                        colnames(select(lte_flow, contains('Q'))) %>%
                          str_to_upper())

ste_qtr_water_quality <- ste_qtr_long %>% select(site_code, any_of(water_quality)) %>%
  pretty_titles()

ste_qtr_geomorph <- ste_qtr_long %>% select(site_code, any_of(geomorph)) %>% pretty_titles()

ste_qtr_algae <- ste_qtr_long %>% select(site_code, any_of(algae)) %>% pretty_titles()

# # Export ste Tables
# ste_qtr_summary_tables <- list(ste_qtr_long,
#                            ste_qtr_flow,
#                            ste_qtr_water_quality,
#                            ste_qtr_geomorph,
#                            ste_qtr_algae)
# 
# names(ste_qtr_summary_tables) <- c('ste_qtr_long',
#                                'ste_qtr_flow',
#                                'ste_qtr_water_quality',
#                                'ste_qtr_geomorph',
#                                'ste_qtr_algae')
# 
# for (i in 1:length(ste_qtr_summary_tables)) {
#   my_place <- paste('exploration/output/', 
#                     names(ste_qtr_summary_tables[i]), 
#                     "_summary_table.csv", sep='')
#   my_object <- ste_qtr_summary_tables[[i]]
#   write_csv(my_object, my_place) }

#------------------------------------------------------------------------------
# End 06_site_table