# 31_baseline_t_test
# Sean Kinard
# 2023-07-09
#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')
library(tidymodels)

trash_vars <- c('year', 'qtr', 'month', 'day', 'nday', 'year_week', 
                'year_month', 'is_baseline', 'annualrain')
id_vars <- c('site_code', 'collection_period')
hurricane_date <- ymd('2017-08-27')
end_RAPID_date <- ymd('2018-03-31')

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
# t_test density (ln-transformed for normality & scaled by site)
#------------------------------------------------------------------------------
visualize_x_vs_time <- function(df, my_variable) {
  colnames(df) <- str_replace_all(colnames(df), my_variable, 'XXX')
  df %>%
  ggplot(aes(collection_period, XXX, fill=is_rapid, color=is_rapid)) +
    facet_wrap(~site_code) +
    geom_hline(yintercept=0, lty=2, lwd=.3) +
    geom_point(size = 3, shape = 21, alpha=.4) +
    geom_point(size = 3, shape = 21, fill=NA) +
    dark_theme_grey(base_size = 12)+
    ylab(my_variable)}

plot_density_vs_time <- visualize_x_vs_time(d_scaled, 'ln_density')

make_tt_density <- function(df) {
  t.test(x = df %>%pull(ln_density),
         y = d_scaled %>% filter(is_rapid == 'no') %>% pull(ln_density),
         alternative = "two.sided", paired = FALSE, var.equal = FALSE) %>%
    broom::tidy()}

tt_density <- d_scaled %>% filter(is_rapid == 'yes') %>% select(-is_rapid) %>%
  mutate(qtr = yearquarter(collection_period)) %>%
  filter(qtr < yearquarter("2018 Q4")) %>%
  group_by(collection_period) %>% nest() %>%
  mutate(tt_density = map(data, make_tt_density)) %>%
  unnest(tt_density)

# site_scaled_ln_density t-test is negative for months 1-2 after hurricane.
# However, positive t-tests for November through February indicate delayed

#------------------------------------------------------------------------------
# t_test biomass (ln-transformed for normality & scaled by site)
#------------------------------------------------------------------------------
plot_ln_biomass_vs_time <- visualize_x_vs_time(d_scaled, 'ln_biomass')

make_tt_biomass <- function(df) {
  t.test(x = df %>%pull(ln_biomass),
         y = d_scaled %>% filter(is_rapid == 'no') %>% pull(ln_biomass),
         alternative = "two.sided", paired = FALSE, var.equal = FALSE) %>%
    broom::tidy()}

tt_biomass <- d_scaled %>% filter(is_rapid == 'yes') %>% select(-is_rapid) %>%
  mutate(qtr = yearquarter(collection_period)) %>%
  filter(qtr < yearquarter("2018 Q4")) %>%
  group_by(collection_period) %>% nest() %>%
  mutate(tt_density = map(data, make_tt_biomass)) %>%
  unnest(tt_density)

# biomass t-tests are negative, except in October after hurricane
#------------------------------------------------------------------------------
# t_test richness (scaled by site)
#------------------------------------------------------------------------------
plot_richness_vs_time <- visualize_x_vs_time(d_scaled, 'richness')

make_tt_richness <- function(df) {
  t.test(x = df %>%pull(richness),
         y = d_scaled %>% filter(is_rapid == 'no') %>% pull(richness),
         alternative = "two.sided", paired = FALSE, var.equal = FALSE) %>%
    broom::tidy()}

tt_richness <- d_scaled %>% filter(is_rapid == 'yes') %>% select(-is_rapid) %>%
  mutate(qtr = yearquarter(collection_period)) %>%
  filter(qtr < yearquarter("2018 Q4")) %>%
  group_by(collection_period) %>% nest() %>%
  mutate(tt_density = map(data, make_tt_richness)) %>%
  unnest(tt_density)

# richness t-test is only positive for September after Hurricane
#------------------------------------------------------------------------------
# t test shannon
#------------------------------------------------------------------------------
plot_shannon_vs_time <- visualize_x_vs_time(d_scaled, 'shannon')

make_tt_shannon <- function(df) {
  t.test(x = df %>%pull(shannon),
         y = d_scaled %>% filter(is_rapid == 'no') %>% pull(shannon),
         alternative = "two.sided", paired = FALSE, var.equal = FALSE) %>%
    broom::tidy()}

tt_shannon <- d_scaled %>% filter(is_rapid == 'yes') %>% select(-is_rapid) %>%
  mutate(qtr = yearquarter(collection_period)) %>%
  filter(qtr < yearquarter("2018 Q4")) %>%
  group_by(collection_period) %>% nest() %>%
  mutate(tt_density = map(data, make_tt_shannon)) %>%
  unnest(tt_density)

# shannon t-test oscillates between positive and negative until February
# Extreme difference is detected in Februrary

#------------------------------------------------------------------------------
# t_test simpson
#------------------------------------------------------------------------------
plot_simpson_vs_time <- visualize_x_vs_time(d_scaled, 'simpson')

make_tt_simpson <- function(df) {
  t.test(x = df %>%pull(simpson),
         y = d_scaled %>% filter(is_rapid == 'no') %>% pull(simpson),
         alternative = "two.sided", paired = FALSE, var.equal = FALSE) %>%
    broom::tidy()}

tt_simpson <- d_scaled %>% filter(is_rapid == 'yes') %>% select(-is_rapid) %>%
  mutate(qtr = yearquarter(collection_period)) %>%
  filter(qtr < yearquarter("2018 Q4")) %>%
  group_by(collection_period) %>% nest() %>%
  mutate(tt_density = map(data, make_tt_simpson)) %>%
  unnest(tt_density)

# simpson t-test oscillates between positive and negative until February
# Extreme difference is detected in Februrary
#------------------------------------------------------------------------------
# Combine t_test results
#------------------------------------------------------------------------------
dc <- select(tt_density, collection_period, p.value) %>% rename(density = p.value) %>%
  left_join(select(tt_biomass, collection_period, p.value) %>% rename(biomass = p.value)) %>%
  left_join(select(tt_richness, collection_period, p.value) %>% rename(richness = p.value)) %>%
  left_join(select(tt_shannon, collection_period, p.value) %>% rename(shannon = p.value)) %>%
  left_join(select(tt_simpson, collection_period, p.value) %>% rename(simpson = p.value)) %>%
  pivot_longer(cols=-collection_period, names_to='metric', values_to='p_value')

critical_alpha <- 1-(1-0.05)^(1/6)

plot_t_test_pvalues <- dc %>%
  ggplot(aes(collection_period, p_value)) +
  scale_x_date(breaks='1 month', date_labels='%m') +
  facet_wrap(~pretty_x(metric)) +
  geom_hline(yintercept= critical_alpha, lty=2, lwd=.3, color='red') +
  geom_point(color='black', fill='grey50', size=3, shape=21) +
  geom_point(data = filter(dc, p_value<critical_alpha),
             aes(collection_period, p_value),
             shape=21, size=7, color='red', fill=NA) +
  scale_y_log10() +
  theme_bw(base_size=12) +
  labs(x='Month', y='p value')





  

#------------------------------------------------------------------------------
# Export table
#------------------------------------------------------------------------------
write_csv(dc, 'exploration/output/t_test_baseline_results.csv')

#------------------------------------------------------------------------------
# Export figures
#------------------------------------------------------------------------------
my_objects<-ls()
my_figure_names <- my_objects[str_detect(my_objects,'plot_')]

my_figures <- list(
plot_density_vs_time, plot_ln_biomass_vs_time, plot_richness_vs_time, 
plot_shannon_vs_time, plot_simpson_vs_time, plot_t_test_pvalues)

names(my_figures) <- my_figure_names

for (i in 1:length(my_figures)) {
  my_place <- paste('exploration/visualization/31_', names(my_figures[i]), ".png", sep='')
  my_object <- my_figures[[i]]
  ggsave(my_place,
         plot = my_object,
         width = 10,
         height = 10,
         units = c("in")) }

#------------------------------------------------------------------------------
# End 31_baseline_t_test