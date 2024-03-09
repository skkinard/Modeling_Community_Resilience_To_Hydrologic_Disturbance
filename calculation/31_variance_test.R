# 31_variance_test
# Sean Kinard
# 2023-07-09
#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')
library(tidymodels)
library(tsibble)

trash_vars <- c('year', 'yrqtr', 'month', 'day', 'nday', 'year_week', 
                'year_month', 'is_baseline', 'annualrain')
id_vars <- c('site_code', 'collection_period')
hurricane_date <- ymd('2017-08-27')
end_RAPID_date <- ymd('2018-12-31')

make_rapid <- function(df) {
  df %>% mutate(
    is_rapid = ifelse(collection_period > hurricane_date &
                        collection_period <= end_RAPID_date,
                      'yes', 'no')) %>%
    select(is_rapid, everything()) }

d <- read_csv('exploration/output/model_combined_vars.csv') %>% 
  select(-any_of(trash_vars)) %>%
  fix_site_order() %>%
  arrange(collection_period, site_code) %>%
  mutate(ln_density = log(density)) %>%
  make_rapid()

d_scaled <- read_csv('exploration/output/model_combined_vars_SCALED.csv') %>% 
  select(-any_of(trash_vars)) %>%
  fix_site_order() %>%
  arrange(collection_period, site_code) %>%
  make_rapid()

#------------------------------------------------------------------------------
# monthly variance
#------------------------------------------------------------------------------
response_vars <- c('ln_density', 'ln_biomass', 'richness', 
                   'shannon', 'simpson', 'ln_centroid_dist')

# visualize regional variance over time
d_scaled %>%
  mutate(yrqtr= yearquarter(collection_period + 35)) %>%
  select(site_code, collection_period, yrqtr, any_of(response_vars)) %>%
  pivot_longer(cols=any_of(response_vars), names_to='xvar', values_to='xval') %>%
  group_by(yrqtr, xvar) %>%
  dplyr::summarise(x_variance = var(xval, na.rm=T),
                   x_events = length(site_code)) %>%
  ungroup() %>%
  mutate(qtr=quarter(yrqtr)%>%as.factor()) %>%
  ggplot(aes(x=yrqtr, y=x_variance, color = qtr)) +
  facet_wrap(~xvar, scales='free_y') +
  dark_theme_grey(base_size=12) +
  geom_path(lty=2,lwd=.3) +
  geom_point(size=3) +
  scale_color_manual(values=c('green', 'orange', 'red', 'skyblue'))

# Regional biomass variation misses a seasonal peak in Q4 of 2017
# Regional density variation remains high in Q4 of 2017
# Implication: size of fish is homogenized in q4 of 2017
