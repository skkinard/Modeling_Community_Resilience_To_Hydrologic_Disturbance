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
library(rsample)
library(ggrepel)
library(tidyposterior)
library(rstanarm)
library(patchwork)
library(rules)
library(baguette)
library(finetune)
library(xgboost)
library(kknn)
library(tsibble)
library(vip)
tidymodels_prefer()
theme_set(dark_theme_grey(base_size = 14))

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
my_qq <- function(my_data, my_function) {
  broom::augment(x = my_function,
                 data = my_data) %>%
    ggqqplot(x='.fitted') }

my_shapiro <- function(my_data, my_function) {
  broom::augment(x = my_function,
                 data = my_data) %>%
    pull(.fitted) %>%
    shapiro.test() }

# -----------------------------------------------------------------------------
# random intercept
# -----------------------------------------------------------------------------
ri_spec <- 
  linear_reg() %>% 
  set_engine("lme", 
             random = ~ 1 | site_code, 
             corr = corAR1(form = ~1|site_code))

# Density
ri_den_qmin <- ri_spec %>% fit(ln_density ~ q_min, data = d)
ri_den_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_den_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri_den_qmin %>% extract_fit_engine()) # normal-fail

ri_den_qmax <- ri_spec %>% fit(ln_density ~ q_min, data = d)
ri_den_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_den_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri_den_qmax %>% extract_fit_engine()) # normal-fail

# Biomass
ri_bio_qmin <- ri_spec %>% fit(ln_biomass ~ q_min, data = d)
ri_bio_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_bio_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri_bio_qmin %>% extract_fit_engine()) # normal-fail

ri_bio_qmax <- ri_spec %>% fit(ln_biomass ~ q_max, data = d)
ri_bio_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, ri_bio_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri_bio_qmax %>% extract_fit_engine()) # normal-fail

# Richness
ri_ric_qmin <- ri_spec %>% fit(richness ~ q_min, data = d)
ri_ric_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_ric_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri_ric_qmin %>% extract_fit_engine()) # normal-fail

ri_ric_qmax <- ri_spec %>% fit(richness ~ q_max, data = d)
ri_ric_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, ri_ric_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri_ric_qmax %>% extract_fit_engine()) # normal-fail

# shannon
ri_sha_qmin <- ri_spec %>% fit(shannon ~ q_min, data = d)
ri_sha_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_sha_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri_sha_qmin %>% extract_fit_engine()) # normal-fail

ri_sha_qmax <- ri_spec %>% fit(shannon ~ q_max, data = d)
ri_sha_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_sha_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri_sha_qmax %>% extract_fit_engine()) # normal-fail

# sqrt_centroid_dist
ri_cdi_qmin <- ri_spec %>% fit(sqrt_centroid_dist ~ q_min, data = d)
ri_cdi_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_cdi_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri_cdi_qmin %>% extract_fit_engine()) # normal-fail

ri_cdi_qmax <- ri_spec %>% fit(sqrt_centroid_dist ~ q_max, data = d)
ri_cdi_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri_cdi_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri_cdi_qmax %>% extract_fit_engine()) # normal-fail

# -----------------------------------------------------------------------------
# Random Intercept with interaction (rainfall)
# -----------------------------------------------------------------------------
ri2_spec <- 
  linear_reg() %>% 
  set_engine("lme", 
             random = ~ annualrain*1 | site_code, 
             corr = corAR1(form = ~1|site_code))

# Density
ri2_den_qmin <- ri2_spec %>% fit(ln_density ~ annualrain*q_min, data = d)
ri2_den_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_den_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri2_den_qmin %>% extract_fit_engine()) # normal-pass

ri2_den_qmax <- ri2_spec %>% fit(ln_density ~ annualrain*q_max, data = d)
ri2_den_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_den_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri2_den_qmax %>% extract_fit_engine()) # normal-pass

# Biomass
ri2_bio_qmin <- ri2_spec %>% fit(ln_biomass ~ annualrain*q_min, data = d)
ri2_bio_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_bio_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri2_bio_qmin %>% extract_fit_engine()) # normal-fail

ri2_bio_qmax <- ri2_spec %>% fit(ln_biomass ~ annualrain*q_max, data = d)
ri2_bio_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT, not interaction
my_qq(d, ri2_bio_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri2_bio_qmax %>% extract_fit_engine()) # normal-fail

# Richness
ri2_ric_qmin <- ri2_spec %>% fit(richness ~ annualrain*q_min, data = d)
ri2_ric_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_ric_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri2_ric_qmin %>% extract_fit_engine()) # normal-fail

ri2_ric_qmax <- ri2_spec %>% fit(richness ~ annualrain*q_max, data = d)
ri2_ric_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT, not interaction
my_qq(d, ri2_ric_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri2_ric_qmax %>% extract_fit_engine()) # normal-fail

# shannon
ri2_sha_qmin <- ri2_spec %>% fit(shannon ~ annualrain*q_min, data = d)
ri2_sha_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_sha_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri2_sha_qmin %>% extract_fit_engine()) # normal-fail

ri2_sha_qmax <- ri2_spec %>% fit(shannon ~ annualrain*q_max, data = d)
ri2_sha_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_sha_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri2_sha_qmax %>% extract_fit_engine()) # normal-fail

# sqrt_centroid_dist
ri2_cdi_qmin <- ri2_spec %>% fit(sqrt_centroid_dist ~ annualrain*q_min, data = d)
ri2_cdi_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_cdi_qmin %>% extract_fit_engine()) 
my_shapiro(d, ri2_cdi_qmin %>% extract_fit_engine()) # normal-fail

ri2_cdi_qmax <- ri2_spec %>% fit(sqrt_centroid_dist ~ annualrain*q_max, data = d)
ri2_cdi_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, ri2_cdi_qmax %>% extract_fit_engine()) 
my_shapiro(d, ri2_cdi_qmax %>% extract_fit_engine()) # normal-fail

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
             random = ~ 1 + q_max | site_code, 
             corr = corAR1(form = ~ 1 | site_code))

make_slope_plot <- function(fit_model_engine) {
  fit_model_engine$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column('site_code') %>%
    add_rain() %>%
    ggplot(aes(x=annualrain, y=site_code.q_min)) +
    geom_point(size=4, color = 'white') +
    geom_smooth(method='lm', lty=2, lwd=.5)}

make_slope_plot2 <- function(fit_model_engine) {
  fit_model_engine$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column('site_code') %>%
    add_rain() %>%
    ggplot(aes(x=annualrain, y=site_code.q_max)) +
    geom_point(size=4, color = 'white') +
    geom_smooth(method='lm', lty=2, lwd=.5)}

# Density
rs_den_qmin <- rs_min_spec %>% fit(ln_density ~ q_min, data = d)
rs_den_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_den_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_den_qmin %>% extract_fit_engine()) # normal-fail
make_slope_plot(rs_den_qmin %>% extract_fit_engine())

rs_den_qmax <- rs_max_spec %>% fit(ln_density ~ q_max, data = d)
rs_den_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, rs_den_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_den_qmax %>% extract_fit_engine()) # normal-fail
make_slope_plot2(rs_den_qmax %>% extract_fit_engine())

# Biomass
rs_bio_qmin <- rs_min_spec %>% fit(ln_biomass ~ q_min, data = d)
rs_bio_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_bio_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_bio_qmin %>% extract_fit_engine()) # normal-fail
make_slope_plot(rs_bio_qmin %>% extract_fit_engine())


rs_bio_qmax <- rs_max_spec %>% fit(ln_biomass ~ q_max, data = d)
rs_bio_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, rs_bio_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_bio_qmax %>% extract_fit_engine()) # normal-fail
make_slope_plot2(rs_bio_qmax %>% extract_fit_engine())

# Richness
rs_ric_qmin <- rs_min_spec %>% fit(richness ~ q_min, data = d)
rs_ric_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_ric_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_ric_qmin %>% extract_fit_engine()) # normal-fail
make_slope_plot(rs_ric_qmin %>% extract_fit_engine())

rs_ric_qmax <- rs_max_spec %>% fit(richness ~ q_max, data = d)
rs_ric_qmax %>% extract_fit_engine() %>% anova() # SIGNIFICANT
my_qq(d, rs_ric_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_ric_qmax %>% extract_fit_engine()) # normal-fail
make_slope_plot2(rs_ric_qmax %>% extract_fit_engine())

# shannon
rs_sha_qmin <- rs_min_spec %>% fit(shannon ~ q_min, data = d)
rs_sha_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_sha_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_sha_qmin %>% extract_fit_engine()) # normal-fail
make_slope_plot(rs_sha_qmin %>% extract_fit_engine())

rs_sha_qmax <- rs_max_spec %>% fit(shannon ~ q_max, data = d)
rs_sha_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_sha_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_sha_qmax %>% extract_fit_engine()) # normal-fail
make_slope_plot2(rs_sha_qmax %>% extract_fit_engine())

# sqrt_centroid_dist
rs_cdi_qmin <- rs_min_spec %>% fit(sqrt_centroid_dist ~  q_min, data = d)
rs_cdi_qmin %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_cdi_qmin %>% extract_fit_engine()) 
my_shapiro(d, rs_cdi_qmin %>% extract_fit_engine()) # normal-fail
make_slope_plot(rs_cdi_qmin %>% extract_fit_engine())

rs_cdi_qmax <- rs_max_spec %>% fit(sqrt_centroid_dist ~ q_max, data = d)
rs_cdi_qmax %>% extract_fit_engine() %>% anova() # NS
my_qq(d, rs_cdi_qmax %>% extract_fit_engine()) 
my_shapiro(d, rs_cdi_qmax %>% extract_fit_engine()) # normal-fail
make_slope_plot2(rs_cdi_qmax %>% extract_fit_engine())

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------