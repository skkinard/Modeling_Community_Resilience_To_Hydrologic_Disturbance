# 33_site_models
# Sean Kinard
# 2023-08-06
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')
library(tidymodels)
library(patchwork)
library(tsibble)
library(ggpmisc)
tidymodels_prefer()
theme_set(dark_theme_grey(base_size = 20))

d <- read_csv('exploration/output/31_combined_vars_SCALED.csv') %>%
  add_rain() %>%
  fix_site_order()

colnames(d) <- str_replace_all(colnames(d), 'z_', '') # all variables are site scaled

# -----------------------------------------------------------------------------
# min q: Visual Exploration
# -----------------------------------------------------------------------------
plot_den_minq <- d %>%
  ggplot(aes(x=ln_wk4_min_to_annual_med, y = ln_density)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # no obvious pattern

plot_bio_minq <- d %>%
  ggplot(aes(x=ln_wk4_min_to_annual_med, y = ln_biomass)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # consistent weak (-)

plot_ric_minq <- d %>%
  ggplot(aes(x=ln_wk4_min_to_annual_med, y = richness)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # consistent weak (-)

plot_sha_minq <- d %>%
  ggplot(aes(x=ln_wk4_min_to_annual_med, y = shannon)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen4') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # somewhat weak (-)

# -----------------------------------------------------------------------------
# min_q: regression functions
# -----------------------------------------------------------------------------
# fits linear regression: response ~ min_q
lm_den_minq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(ln_density ~ ln_wk4_min_to_annual_med, data = x)
  return(lm_x_fit) }
lm_bio_minq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(ln_biomass ~ ln_wk4_min_to_annual_med, data = x)
  return(lm_x_fit) }
lm_ric_minq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(richness ~ ln_wk4_min_to_annual_med, data = x)
  return(lm_x_fit) }
lm_sha_minq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(shannon ~ ln_wk4_min_to_annual_med, data = x)
  return(lm_x_fit) }

# -----------------------------------------------------------------------------
# min_q regression x site: coefficients
# -----------------------------------------------------------------------------
d_minq <- d %>%
  fix_site_order() %>%
  arrange(site_code, collection_period) %>%
  group_by(site_code) %>%
  nest() %>%
  mutate(my_lm_den = map(data, lm_den_minq)) %>%
  mutate(my_lm_bio = map(data, lm_bio_minq)) %>%
  mutate(my_lm_ric = map(data, lm_ric_minq)) %>%
  mutate(my_lm_sha = map(data, lm_sha_minq))

visualize_coef <- function(df) {
  
  d_prep <- df %>%
    unnest(my_tidy) %>%
    filter(term == 'ln_wk4_min_to_annual_med') %>%
    select(-term) %>%
    add_rain()
  
  output <- d_prep %>%
    ggplot(aes(x=annualrain, y = estimate)) +
    geom_point(size=3, shape=21, fill = 'palegreen4') +
    geom_smooth(method='lm', se=F, lty=2, lwd=.8) +
    labs(x=element_blank(), y=element_blank()) +
    stat_poly_eq(color='white', size = 3)
  
  return(output) }

plot_coef_den_minq <- d_minq %>%
  mutate(my_tidy = map(my_lm_den, broom::tidy)) %>%
  visualize_coef() +
  ggtitle('Density')

plot_coef_bio_minq <- d_minq %>%
  mutate(my_tidy = map(my_lm_bio, broom::tidy)) %>%
  visualize_coef()+
  ggtitle('Biomass')

plot_coef_ric_minq <- d_minq %>%
  mutate(my_tidy = map(my_lm_ric, broom::tidy)) %>%
  visualize_coef()+
  ggtitle('Richness')

plot_coef_sha_minq <- d_minq %>%
  mutate(my_tidy = map(my_lm_ric, broom::tidy)) %>%
  visualize_coef()+
  ggtitle('Shannon')

plot_minq_combo <- (plot_coef_den_minq + ylab('min_q slope') + 
    plot_coef_bio_minq) /
  (plot_coef_ric_minq + ylab('min_q slope') + xlab('Annual Rain') +
     plot_coef_sha_minq + xlab('Annual Rain'))

# -----------------------------------------------------------------------------
# max_q: Visual Exploration
# -----------------------------------------------------------------------------
plot_den_maxq <- d %>%
  ggplot(aes(x=ln_wk4_max_to_annual_med, y = ln_density)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # consistent (-)

plot_bio_maxq <- d %>%
  ggplot(aes(x=ln_wk4_max_to_annual_med, y = ln_biomass)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # consistent (-)

plot_ric_maxq <- d %>%
  ggplot(aes(x=ln_wk4_max_to_annual_med, y = richness)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) # 2 Arid (+) vs  7x (-)

plot_sha_maxq <- d %>%
  ggplot(aes(x=ln_wk4_max_to_annual_med, y = shannon)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(size=3, shape=21, fill='palegreen', alpha=.3) +
  geom_point(size=3, shape=21, color='palegreen') +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5) + # TR,SF,PD (+) vs 6x (-)
  labs( x = 'Antecedent Max Discharge',
        y = 'Diversity')

d %>%
  filter(site_code %in% c('TR', 'PL', 'GC')) %>%
  mutate(site_code = case_when(
    site_code == 'TR' ~ 'Tranquitas Creek',
    site_code == 'PL' ~ 'Placedo Creek', 
    site_code == 'GC' ~ 'Garcitas Creek')) %>%
  mutate(site_code = fct_relevel(site_code, 
                                 c('Semi-Humid', 'Mesic', 'Semi-Arid'))) %>%
  ggplot(aes(x=ln_wk4_max_to_annual_med, y = richness)) +
  facet_wrap(~site_code, scales = 'free') +
  geom_point(aes(fill=site_code, color = site_code),
             size=3, shape=21, alpha=.3) +
  geom_point(aes(color = site_code), fill=NA,
             size=3, shape=21) +
  scale_color_manual(values=c('skyblue2', 'green3', 'gold')) +
  scale_fill_manual(values=c('skyblue2', 'green3', 'gold')) +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5, color = 'white') +
  labs(y='Species Richness', x = 'Q-Max') +
  theme(legend.position = 'none')

d %>%
  create_site_group() %>%
  mutate(site_group = str_replace_all(site_group, 'Sub', 'Semi')) %>%
  mutate(site_group = str_replace_all(site_group, 'Transition', 'Mesic')) %>%
  mutate(site_code = fct_relevel(site_group, 
                                 c('Semi-Humid', 'Mesic', 'Semi-Arid'))) %>%
  ggplot(aes(x=ln_wk4_max_to_annual_med, y = richness)) +
  facet_wrap(~site_group) +
  geom_point(aes(fill=site_group, color = site_group),
             size=3, shape=21, alpha=.3) +
  geom_point(aes(color = site_group), fill=NA,
             size=3, shape=21) +
  scale_color_manual(values=c('skyblue2', 'green3', 'gold')) +
  scale_fill_manual(values=c('skyblue2', 'green3', 'gold')) +
  geom_smooth(method='lm', se=F, lty=2, lwd=.5, color = 'white') +
  labs(y='Species Richness', x = 'Q-Max') +
  theme(legend.position = 'none')


# -----------------------------------------------------------------------------
# max_q: regression functions
# -----------------------------------------------------------------------------
# fits linear regression: response ~ max_q
lm_den_maxq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(ln_density ~ ln_wk4_max_to_annual_med, data = x)
  return(lm_x_fit) }
lm_bio_maxq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(ln_biomass ~ ln_wk4_max_to_annual_med, data = x)
  return(lm_x_fit) }
lm_ric_maxq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(richness ~ ln_wk4_max_to_annual_med, data = x)
  return(lm_x_fit) }
lm_sha_maxq <- function(x) {
  lm_x <- linear_reg() %>% set_engine("lm")
  lm_x_fit <- lm_x %>% fit(shannon ~ ln_wk4_max_to_annual_med, data = x)
  return(lm_x_fit) }

# -----------------------------------------------------------------------------
# max_q regression x site: coefficients
# -----------------------------------------------------------------------------
d_maxq <- d %>%
  fix_site_order() %>%
  arrange(site_code, collection_period) %>%
  group_by(site_code) %>%
  nest() %>%
  mutate(my_lm_den = map(data, lm_den_maxq)) %>%
  mutate(my_lm_bio = map(data, lm_bio_maxq)) %>%
  mutate(my_lm_ric = map(data, lm_ric_maxq)) %>%
  mutate(my_lm_sha = map(data, lm_sha_maxq))

visualize_coef <- function(df) {
  
  d_prep <- df %>%
    unnest(my_tidy) %>%
    filter(term == 'ln_wk4_max_to_annual_med') %>%
    select(-term) %>%
    add_rain()
  
  output <- d_prep %>%
    ggplot(aes(x=annualrain, y = estimate)) +
    geom_point(size=3, shape=21, fill = 'palegreen4') +
    geom_smooth(method='lm', se=F, lty=2, lwd=.5) +
    labs(x=element_blank(), y=element_blank()) +
    stat_poly_eq(use_label(c("R2", "P")), color='white', size = 6, 
                 label.x = "right")
  
  return(output) }

plot_coef_den_maxq <- d_maxq %>%
  mutate(my_tidy = map(my_lm_den, broom::tidy)) %>%
  visualize_coef() +
  ggtitle('Density')

plot_coef_bio_maxq <- d_maxq %>%
  mutate(my_tidy = map(my_lm_bio, broom::tidy)) %>%
  visualize_coef()+
  ggtitle('Biomass')

plot_coef_ric_maxq <- d_maxq %>%
  mutate(my_tidy = map(my_lm_ric, broom::tidy)) %>%
  visualize_coef()+
  ggtitle('Richness')

plot_coef_sha_maxq <- d_maxq %>%
  mutate(my_tidy = map(my_lm_ric, broom::tidy)) %>%
  visualize_coef()+
  ggtitle('Shannon')

plot_maxq_combo <- (plot_coef_den_maxq + ylab('max_q slope') + 
                      plot_coef_bio_maxq) /
  (plot_coef_ric_maxq + ylab('max_q slope') + xlab('Annual Rain') +
     plot_coef_sha_maxq + xlab('Annual Rain'))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
my_objects<-ls()
my_figure_names <- my_objects[str_detect(my_objects,'plot_')]

my_figures<-list(
plot_bio_maxq, plot_bio_minq, plot_coef_bio_maxq, plot_coef_bio_minq, 
plot_coef_den_maxq, plot_coef_den_minq, plot_coef_ric_maxq, plot_coef_ric_minq, 
plot_coef_sha_maxq, plot_coef_sha_minq, plot_den_maxq, plot_den_minq, 
plot_maxq_combo, plot_minq_combo, plot_ric_maxq, plot_ric_minq, plot_sha_maxq, 
plot_sha_minq )

names(my_figures) <- str_replace_all(my_figure_names, 'plot_', '33_site_mdoels_')

for (i in 1:length(my_figures)) {
  my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
  my_object <- my_figures[[i]]
  ggsave(my_place,
         plot = my_object,
         width = 10,
         height = 10,
         units = c("in")) }
# -----------------------------------------------------------------------------
# End 33_site_models

