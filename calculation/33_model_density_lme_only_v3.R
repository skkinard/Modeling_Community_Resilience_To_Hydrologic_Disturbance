# 33_model_density_lme_only_v3
# Sean Kinard
# 2023-08-06
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')
library(MASS)
library(tidymodels)
library(multilevelmod)
library(patchwork)
library(tsibble)
library(vip)
library(ggpubr)
library(nlme)
tidymodels_prefer()
theme_set(theme_bw(base_size = 14))

d <- read_csv('exploration/output/31_combined_vars_SCALED.csv') %>%
  add_rain()

# -----------------------------------------------------------------------------
# Data modifications
# -----------------------------------------------------------------------------
colnames(d) <- str_replace_all(colnames(d), 'z_', '')

d <- d %>% rename(q_min = ln_wk4_min_to_annual_med,
                  q_max = ln_wk4_max_to_annual_med) 

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------
my_qq <- function(x_lme_extract) {
  broom::augment(x = x_lme_extract,
                 data = d) %>%
    ggqqplot(x='.resid') }

my_shapiro <- function(x_lme_extract) {
  broom::augment(x = x_lme_extract,
                 data = d) %>%
    pull(.resid) %>%
    shapiro.test()%>%
    tidy() %>%
    pull(p.value) }

# -----------------------------------------------------------------------------
# random intercept
# -----------------------------------------------------------------------------
make_ri_fit <- function(df) {
  ri_spec <- 
    linear_reg() %>% 
    set_engine("lme", 
               random = ~ 1 | site_code, 
               corr = corAR1(form = ~1|site_code))
  
  x_fit <- ri_spec %>% fit(as.formula(df$x_formula), data = d)
  return(x_fit) }

outcome <- c('ln_density', 'ln_biomass', 'richness', 'shannon', 
              'sqrt_centroid_dist')
predictor <- c('q_min', 'q_max')

# create crossed data frame
d_ri <- crossing(outcome, predictor) %>%
  mutate(outcome_abr = outcome %>%
           str_replace_all('ln_', '') %>%
           str_replace_all('sqrt_', '') %>%
           substr(1,3),
         model_type = 'ri', 
         autoregression = 'yes',
         model_name = paste(model_type, outcome_abr, predictor, sep='_'),
         x_formula = 
           paste(outcome, predictor, sep=' ~ '))

# fit models and diagnostics
d_ri_fitted <- d_ri %>%
  group_by(model_name) %>%
  nest() %>%
  mutate(lme_fit = map(data, make_ri_fit),
         lme_extract = map(lme_fit, extract_fit_engine),
         lme_glance = map(lme_extract, broom::glance),
         res_qq = map(lme_extract, my_qq),
         shapiro_p = map_dbl(lme_extract, my_shapiro),
         normal_res = ifelse(shapiro_p<0.05, 'no', 'yes')) %>%
  unnest(lme_glance) %>%
  select(-nobs, - sigma, -logLik, -BIC) %>%
  arrange(AIC)

# -----------------------------------------------------------------------------
# Random Intercept with interaction (rainfall)
# -----------------------------------------------------------------------------
make_rx_fit <- function(df) {
  rx_spec <- 
    linear_reg() %>% 
    set_engine("lme", 
               random = ~ 1 | site_code, 
               corr = corAR1(form = ~1|site_code))
  
  x_fit <- ri_spec %>% fit(as.formula(df$x_formula), data = d)
  return(x_fit) }

# create crossed data frame
d_rx <- crossing(outcome, predictor) %>%
  mutate(outcome_abr = outcome %>%
           str_replace_all('ln_', '') %>%
           str_replace_all('sqrt_', '') %>%
           substr(1,3),
         model_type = 'rx', 
         autoregression = 'yes',
         model_name = paste(model_type, outcome_abr, predictor, sep='_'),
         x_formula = 
           paste(outcome, paste(predictor, "*annualrain", sep=''), sep=' ~ '))

# fit models and diagnostics
d_rx_fitted <- d_rx %>%
  group_by(model_name) %>%
  nest() %>%
  mutate(lme_fit = map(data, make_rx_fit),
         lme_extract = map(lme_fit, extract_fit_engine),
         lme_glance = map(lme_extract, broom::glance),
         res_qq = map(lme_extract, my_qq),
         shapiro_p = map_dbl(lme_extract, my_shapiro),
         normal_res = ifelse(shapiro_p<0.05, 'no', 'yes')) %>%
  unnest(lme_glance) %>%
  select(-nobs, - sigma, -logLik, -BIC) %>%
  arrange(AIC)

# -----------------------------------------------------------------------------
# Random Slope
# -----------------------------------------------------------------------------
rs_min_spec <- 
  linear_reg() %>% 
  set_engine("lme",
             random = ~ 1 + q_min | site_code, 
             corr = corAR1(form = ~ 1 | site_code))

rs_max_spec <- 
  linear_reg() %>% 
  set_engine("lme",
             random = ~ 1 + q_min | site_code, 
             corr = corAR1(form = ~ 1 | site_code))

# Density
rs_den_qmin <- rs_min_spec %>% fit(ln_density ~ q_min, data = d)
rs_den_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_den_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_den_qmin %>% extract_fit_engine()) # normal-fail

rs_den_qmax <- rs_max_spec %>% fit(ln_density ~ q_max, data = d)
rs_den_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, rs_den_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_den_qmax %>% extract_fit_engine()) # normal-fail

# Biomass
rs_bio_qmin <- rs_min_spec %>% fit(ln_biomass ~ q_min, data = d)
rs_bio_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_bio_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_bio_qmin %>% extract_fit_engine()) # normal-fail

rs_bio_qmax <- rs_max_spec %>% fit(ln_biomass ~ q_max, data = d)
rs_bio_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, rs_bio_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_bio_qmax %>% extract_fit_engine()) # normal-fail

# Richness
rs_ric_qmin <- rs_min_spec %>% fit(richness ~ q_min, data = d)
rs_ric_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_ric_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_ric_qmin %>% extract_fit_engine()) # normal-fail

rs_ric_qmax <- rs_max_spec %>% fit(richness ~ q_max, data = d)
rs_ric_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, rs_ric_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_ric_qmax %>% extract_fit_engine()) # normal-fail

# shannon
rs_sha_qmin <- rs_min_spec %>% fit(shannon ~ q_min, data = d)
rs_sha_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_sha_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_sha_qmin %>% extract_fit_engine()) # normal-fail

rs_sha_qmax <- rs_max_spec %>% fit(shannon ~ q_max, data = d)
rs_sha_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_sha_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_sha_qmax %>% extract_fit_engine()) # normal-fail

# sqrt_centroid_dist
rs_cdi_qmin <- rs_min_spec %>% fit(sqrt_centroid_dist ~  q_min, data = d)
rs_cdi_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_cdi_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_cdi_qmin %>% extract_fit_engine()) # normal-fail

rs_cdi_qmax <- rs_max_spec %>% fit(sqrt_centroid_dist ~ q_max, data = d)
rs_cdi_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_cdi_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_cdi_qmax %>% extract_fit_engine()) # normal-fail

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------