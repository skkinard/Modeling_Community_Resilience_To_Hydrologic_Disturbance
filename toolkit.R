# resilience_tools
# Sean Kinard
# 2023-08-06

#------------------------------------------------------------------------------
# load packages
#------------------------------------------------------------------------------
library(tidyverse) # dplyr and ggplot
library(lubridate) # handling dates
library(ggdark) # plot themes
library(ggrepel) # plot labels
library(patchwork) # combo plots
library(corrr) # autocorrelation
library(tsibble) # time series
library(imputeTS) # time series imputation
library(broom.mixed) # tidy mixed models
library(performance) # SEM peformance statistics
library(ggpmisc)
library(ggConvexHull)
library(ggpubr)

#------------------------------------------------------------------------------
# vectors
#------------------------------------------------------------------------------
my_sites <- c('TR', 'SF', 'AR', 'MR', 'PD', 'PL', 'GC', 'WM', 'EM')

my_colors <- c('#F5191CFF', '#A54E21FF', '#EC7014FF',
               '#FFC72CFF', '#FFF7BCFF', '#66C2A4FF',
               '#006D2CFF', '#6BAED6FF', '#08519CFF')

# group predictor colnames for easier indexing
location <- c('lat', 'lon')
climate <- c('annualtemp', 'annualrain')
landuse <- c('developedland', 'forestland', 'cropland', 'otherland')
lt_flow <- colnames(read_csv('exploration/output/environment_longterm.csv') %>%
                   dplyr::select(contains('q_')))
st_flow <- read_csv("data/06_fill_data/site_flow_4week_flood_duration.csv") %>%
  select(-site_period, -collection_period, -lf_min_return,
         -contains('is_'), -contains('date') ) %>%
  colnames()
flow <- c(lt_flow, st_flow)
water_quality <- c('ammonia', 'conductivity', 'd_oxygen',
                   'doc', 'nitrate', 'ph', 'phosphate',
                   'temperature', 'turbidity')
geomorph <- c('gravel', 'sand', 'silt',
              'canopy', 'depth_mx', 'width')
algae <- c('bluegreen_cyano', 'green_algae', 'diatoms')

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
fix_site_order <- function(x) { 
  mutate(x, site_code = fct_relevel(site_code, my_sites)) }

add_rain <- function(x) { 
  read_csv('exploration/output/environment_longterm.csv') %>%
    dplyr::select(site_code, annualrain) %>%
    right_join(x) }

add_hurricane <- function(x) { 
  read_csv('data/04_match_data/site_hurricane.csv') %>%
    right_join(x) }

# determine baseline
is_baseline <- function(x, end_test) {
  x %>%
    mutate(is_baseline=ifelse(
      collection_period > ymd("2017-04-20") &
        collection_period <= ymd(end_test), "test", "baseline")) }

restore_category <- function(x) { 
  case_when(
    x %in% ln_vars ~ paste('ln_', x, sep=''),
    x %in% sqrt_vars ~ paste('sqrt_', x, sep=''),
    TRUE ~ x) }

r_friendly_colnames <- function(d) {
  output <- d
  colnames(output) <- str_to_lower(colnames(output))
  colnames(output) <- str_replace_all(colnames(output), ' ', '_')
  colnames(output) <- str_replace_all(colnames(output), '\\.', '_') 
  return(output) }

pretty_x <- function(x) {
  x <- str_replace_all(x, '_', ' ')
  x <- str_to_title(x)
  return(x) }

grid_date <- function(x) {
  x %>% group_by(site_code, collection_period) %>%
    dplyr:: summarize(n=length(site_code)) %>%
    arrange(collection_period) %>%
    pivot_wider(names_from = site_code,
                values_from = n) }

remove_na <- function(df) {
  df %>% 
    rowid_to_column("ID") %>%
    pivot_longer(
      cols=-c('site_code', 'ID'), 
      names_to='predictor', 
      values_to = 'x') %>%
    filter(!is.na(x)) %>%
    pivot_wider(names_from=predictor, values_from=x) %>%
    dplyr::select(-ID) %>%
    na.omit() }

create_site_group <- function(df) {
  df %>% mutate(site_group = case_when(
    site_code %in% c('TR', 'SF', 'AR') ~ 'Semi-Arid',
    site_code %in% c('MR', 'PD', 'PL') ~ 'Transition',
    site_code %in% c('GC', 'WM', 'EM') ~ 'Sub-Humid'))%>%
    mutate(site_group = fct_relevel(site_group, c('Semi-Arid',
                                                  'Transition',
                                                  'Sub-Humid')))}
add_season <- function(df) {
  Winter <- c(11, 12, 1, 2, 3) %>% month()
  Spring <- c(4,5,6) %>% month()
  Summer <- c(7,8) %>% month()
  Fall <- c(9,10) %>% month()
  
  df %>% mutate(
    month=month(collection_period)%>%as.numeric(),
    season = case_when(
      month %in% Winter ~ 'Winter',
      month %in% Spring ~ 'Spring',
      month %in% Summer ~ "Summer",
      month %in% Fall ~ "Fall"),
    season = fct_relevel(season, 'Winter', 'Spring', 'Summer', 'Fall')) %>%
    select(-month) }

remove_na_without_site <- function(df) {
  df %>% 
    rowid_to_column("ID") %>%
    pivot_longer(
      cols=-c('ID'), 
      names_to='predictor', 
      values_to = 'x') %>%
    filter(!is.na(x)) %>%
    pivot_wider(names_from=predictor, values_from=x) %>%
    dplyr::select(-ID) %>%
    na.omit() }

scale_by_site <- function(df, xvar) {
  my_data <- df
  colnames(my_data) <- str_replace_all(colnames(my_data), xvar, 'XXX')
  
  output <- my_data %>%
    pivot_wider(values_from=XXX, names_from = site_code) %>%
    mutate(AR = scale(AR),
           EM = scale(EM),
           GC = scale(GC),
           MR = scale(MR),
           PD = scale(PD),
           PL = scale(PL),
           SF = scale(SF),
           TR = scale(TR),
           WM = scale(WM)) %>%
    pivot_longer(cols=AR:WM, names_to='site_code', values_to='XXX_scaled') %>%
    na.omit() 
  
  colnames(output) <- str_replace_all(colnames(output), 'XXX', xvar)
  return(output) }

create_date_vars <- function(my_data) {
  my_data %>%
    mutate(year=year(collection_period),
           month=month(collection_period),
           day=day(collection_period),
           nday=yday(collection_period),
           qtr = quarter(collection_period),
           qtr = case_when( # consolidate 2019 quarters
             collection_period > ymd("2019-01-01") & 
               collection_period < ymd("201-06-01") ~ 1,
             collection_period > ymd("2019-06-01") & 
               collection_period < ymd("201-12-01") ~ 3,
             TRUE ~ qtr),
           year_week = yearweek(collection_period),
           year_month = yearmonth(collection_period)) }
# site_period = lumping sampling events + or - 2 weeks to beginning of month
create_site_period <- function(d) {
  my_d <- d %>%
    mutate(my_date = Collection_Date) %>%
    separate(my_date, c("Year", "Month", "Day"), sep = '-') %>%
    mutate(Day = as.numeric(Day),
           Month = as.numeric(Month)) %>%
    mutate(Collection_Month = case_when(
      Site_Code == 'AR' & Collection_Date == ymd("2019-03-21") ~ as.character(Month),
      Day <= 20 ~ as.character(Month),
      TRUE ~ as.character(Month + 1))) %>%
    mutate(Collection_Period = ym(paste(Year, as.character(Collection_Month), 
                                        sep = '-'))) %>%
    mutate(site_period = paste(Site_Code, Collection_Period, sep = '_')) %>%
    dplyr::select( - Year, -Month, -Day, -Collection_Month) 
  
  return(my_d)  }

create_terrg_period <- function(d) {
  my_d <- d %>%
    mutate(my_date = Collection_Date) %>%
    separate(my_date, c("Year", "Month", "Day"), sep = '-') %>%
    mutate(Day = as.numeric(Day),
           Month = as.numeric(Month)) %>%
    mutate(Collection_Month = ifelse(Day < 20,
                                     Month, as.character(Month + 1))) %>%
    mutate(Collection_Period = ym(paste(Year, as.character(Collection_Month), 
                                        sep = '-'))) %>%
    mutate(site_period = paste(Site_Code, Collection_Period, sep = '_')) %>%
    dplyr::select( - Year, -Month, -Day, -Collection_Month) 
  
  return(my_d) }

combo_create_period <- function(df) {
  pre_2020 <- filter(df, collection_date < ymd('2020-01-01')) %>%
    rename(Site_Code = site_code,
           Collection_Date = collection_date) %>%
    create_site_period() %>%
    r_friendly_colnames()
  
  post_2020 <- filter(df, collection_date >= ymd('2020-01-01')) %>%
    rename(Site_Code = site_code,
           Collection_Date = collection_date) %>%
    create_terrg_period() %>%
    r_friendly_colnames()
  
  output <- full_join(pre_2020, post_2020) 
  
  return(output) }