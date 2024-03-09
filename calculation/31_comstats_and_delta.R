# com_stats_and_delta
# Sean Kinard
# 7-25-2022
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
library(tidyverse) # dplyr and ggplot
library(lubridate) # handling dates
library(tsibble) # time series
library(imputeTS) # time series imputation

my_sites <- c('TR', 'SF', 'AR', 'MR', 'PD', 'PL', 'GC', 'WM', 'EM')

good_vars <- c('site_code', 'collection_period', 'lowest_taxon', 'biomass', 
               'density')
d_den <- read_csv('data/06_fill_data/fish_density_fill.csv') %>% 
  select(any_of(good_vars))

d_bio <- read_csv('data/06_fill_data/fish_biomass.csv') %>% 
  rename(biomass = AFDMg_m2) %>% 
  select(any_of(good_vars))

#------------------------------------------------------------------------------
# function: comstats creates (mean,median,n) x (alltime, year, quarter)
#------------------------------------------------------------------------------
comstats <- function(x, xvar) {
  
  df <- x %>%
    mutate(collection_year = year(collection_period)) %>%
    mutate(collection_qtr = quarter(collection_period))
  
  colnames(df) <- str_replace_all(colnames(df), xvar, 'XXX')
  
  d_sum <- df %>%
    group_by(site_code, collection_period, collection_year, collection_qtr) %>%
    dplyr::summarize(XXX_com=sum(XXX,na.rm=T)) 
  
  d_stats_alltime <- d_sum %>%
    group_by(site_code) %>%
    dplyr::summarize(XXX_com_mu_at = mean(XXX_com, na.rm=T),
                     XXX_com_me_at = median(XXX_com, na.rm=T),
                     XXX_com_n_at = length(XXX_com))
  
  d_stats_year <- d_sum %>%
    group_by(site_code, collection_year) %>%
    dplyr::summarize(XXX_com_mu_yr = mean(XXX_com, na.rm=T),
                     XXX_com_me_yr = median(XXX_com, na.rm=T),
                     XXX_com_n_at_yr = length(XXX_com))
  
  d_stats_quarter <- d_sum %>%
    group_by(site_code, collection_qtr) %>%
    dplyr::summarize(XXX_com_mu_qt = mean(XXX_com, na.rm=T),
                     XXX_com_me_qt = median(XXX_com, na.rm=T),
                     XXX_com_n_qt = length(XXX_com))
  
  my_output <- df %>%
    left_join(d_sum) %>%
    left_join(d_stats_alltime) %>%
    left_join(d_stats_year) %>%
    left_join(d_stats_quarter) %>%
    arrange(site_code, collection_period)
  
  colnames(my_output) <- str_replace_all(colnames(my_output), 'XXX', xvar)
  
  return(my_output)}
#------------------------------------------------------------------------------
# Iterate comstats
#------------------------------------------------------------------------------
density_comstats <- d_den %>% comstats('density')

biomass_comstats <- d_bio %>% comstats('biomass')

#------------------------------------------------------------------------------
# Function: create_delta (delta_k1, delta_mu, delta_me) * (at, yr, qt)
#------------------------------------------------------------------------------
create_delta <- function(x, xvar) {
  df <- x
  colnames(df) <- str_replace_all(colnames(df), xvar, 'YYY')
  
  df <- df %>%
    select(-lowest_taxon, -YYY) %>%
    mutate(ymonth=yearmonth(collection_period)) %>%
    ungroup() %>%
    unique()
  
  # create tibble with all sites by all dates
  date_sequence <- seq.Date(as_date("2017-03-01"), 
                            length.out = 46L, 
                            by = "month") %>% yearmonth()
  site_x_date_cross <- tibble(site_code = list(my_sites),
                              ymonth = list(date_sequence)) %>%
    unnest(site_code) %>%
    unnest(ymonth)
  
  # join data to site_x_date_cross
  d_gaps <- left_join(site_x_date_cross, df)
  
  # use tsibble to calculate lag-difference (k=1)
  d_tsibble <- d_gaps %>%
    as_tsibble(key=site_code,
               index=ymonth)%>%
    mutate(YYY_delta_k1 = difference(YYY_com, lag = 1)) %>%
    as_tibble() %>%
    ungroup()
  
  # calculate all_time difference
  d_delta <- d_gaps %>%
    mutate(YYY_com_delta_at_mu = YYY_com-YYY_com_mu_at,
           YYY_com_delta_at_me = YYY_com-YYY_com_me_at,
           YYY_com_delta_yr_mu = YYY_com-YYY_com_mu_yr,
           YYY_com_delta_yr_me = YYY_com-YYY_com_me_yr,
           YYY_com_delta_qt_mu = YYY_com-YYY_com_mu_qt,
           YYY_com_delta_qt_me = YYY_com-YYY_com_me_qt) %>%
    select(site_code, ymonth, YYY_com, contains('delta'))
  
  colnames(d_delta) <- str_replace_all(colnames(d_delta), 'YYY', xvar)
  
  return(d_delta) }

#------------------------------------------------------------------------------
# Iterate create_delta
#------------------------------------------------------------------------------
delta_density <- d_den %>%
  comstats('density') %>%
  create_delta('density')

delta_biomass <- d_den %>%
  comstats('density') %>%
  create_delta('density')

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(density_comstats, 'exploration/output/fish_density_comstats.csv')
write_csv(biomass_comstats, 'exploration/output/fish_biomass_comstats.csv')

write_csv(delta_density, 'exploration/output/fish_density_delta.csv')
write_csv(delta_biomass, 'exploration/output/fish_biomass_delta.csv')

#------------------------------------------------------------------------------
# com_stats_and_delta 