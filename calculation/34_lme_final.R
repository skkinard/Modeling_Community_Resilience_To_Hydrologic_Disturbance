# 34_lme_final
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
library(ggpmisc)
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
# random intercept: y = x + (1 | site_code)
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
# Random Intercept with interaction: y = (annualrain * x) + (1 | site_code)
# -----------------------------------------------------------------------------
make_rx_fit <- function(df) {
  rx_spec <- 
    linear_reg() %>% 
    set_engine("lme", 
               random = ~ 1 | site_code, 
               corr = corAR1(form = ~1|site_code))
  
  x_fit <- rx_spec %>% fit(as.formula(df$x_formula), data = d)
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
# Random Slope: : y= x + (1+ qmin | site_code)
# -----------------------------------------------------------------------------
# create crossed data frame
d_rs <- crossing(outcome, predictor) %>%
  mutate(outcome_abr = outcome %>%
           str_replace_all('ln_', '') %>%
           str_replace_all('sqrt_', '') %>%
           substr(1,3),
         model_type = 'rs', 
         autoregression = 'yes',
         model_name = paste(model_type, outcome_abr, predictor, sep='_'),
         x_formula = 
           paste(outcome, predictor, sep=' ~ '))

make_rx_fit_qmin <- function(df) {
  rs_min_spec <- 
    linear_reg() %>% 
    set_engine("lme",
               random = ~ 1 + q_min | site_code, 
               corr = corAR1(form = ~ 1 | site_code))
  
  x_fit <- rs_min_spec %>% fit(as.formula(df$x_formula), data = d)
  return(x_fit) }

make_slope_plot_qmin <- function(fit_model_engine) {
  my_label <- paste(fit_model_engine$terms[[2]], 'vs', 'Q-Min', sep=' ') %>%
    str_replace_all('ln_', '') %>%
    str_replace_all('sqrt_', '') %>%
    str_to_title() %>%
    str_replace_all('_', ' ') %>%
    str_replace_all('Vs', 'vs')
  
  my_output <- fit_model_engine$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column('site_code') %>%
    add_rain() %>%
    ggplot(aes(x=annualrain, y=site_code.q_min)) +
    geom_point(size=4, color = 'black') +
    stat_poly_eq(label.x = "right",
                 label.y = "bottom",
                 use_label(c("R2"))) +
    labs(title=my_label,
         y='Slope',
         x='Annual Rainfall (cm)')
  
  # use R^2 to add or leave out line
  my_r2 <- make_slope_lm_qmin(fit_model_engine) %>% pull(lms_r.squared)
  
  if (my_r2 > 0.1) {
    my_output <- my_output + geom_smooth(method='lm', se=F, lty=2, lwd=.5)
    return(my_output) } 
  else 
  {return(my_output)}   }

extract_anova_p <- function(fit_model_engine) {
  fit_model_engine %>%
    anova() %>% 
    rownames_to_column('coef') %>%
    filter(coef!='(Intercept)') %>% 
    pull(`p-value`) }

make_slope_lm_qmin <- function(fit_model_engine) {
  d_slope <- fit_model_engine$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column('site_code') %>%
    add_rain()
  
  my_lm <- linear_reg() %>% 
    set_engine("lm") %>%
    fit(site_code.q_min ~ annualrain, data = d_slope )%>% 
    extract_fit_engine()
  
  my_est <- my_lm %>% tidy() %>% filter(term=='annualrain') %>% pull(estimate)
  
  my_output <- my_lm %>% glance() %>% select(r.squared, AIC, p.value) %>%
    mutate(estimate=my_est)
  colnames(my_output) <- paste('lms', colnames(my_output), sep='_')
  return(my_output)   }

# fit models and diagnostics
d_rs_fitted_qmin <- d_rs %>%
  filter(predictor == 'q_min') %>%
  group_by(model_name) %>%
  nest() %>%
  mutate(lme_fit = map(data, make_rx_fit_qmin),
         lme_extract = map(lme_fit, extract_fit_engine),
         lme_anova_p = map_dbl(lme_extract, extract_anova_p),
         lme_anova_sig = ifelse(lme_anova_p < 0.05, 'sig', 'ns'),
         lme_glance = map(lme_extract, broom::glance),
         res_qq = map(lme_extract, my_qq),
         shapiro_p = map_dbl(lme_extract, my_shapiro),
         normal_res = ifelse(shapiro_p<0.05, 'no', 'yes'),
         plot_slope = map(lme_extract, make_slope_plot_qmin),
         lm_slope = map(lme_extract, make_slope_lm_qmin)) %>%
  unnest(lme_glance) %>%
  select(-nobs, - sigma, -logLik, -BIC) %>%
  arrange(AIC) %>%
  unnest(lm_slope)

# -----------------------------------------------------------------------------
# Random Slope: : y= x + (1+ qmax | site_code)
# -----------------------------------------------------------------------------
# create crossed data frame
d_rs <- crossing(outcome, predictor) %>%
  mutate(outcome_abr = outcome %>%
           str_replace_all('ln_', '') %>%
           str_replace_all('sqrt_', '') %>%
           substr(1,3),
         model_type = 'rs', 
         autoregression = 'yes',
         model_name = paste(model_type, outcome_abr, predictor, sep='_'),
         x_formula = 
           paste(outcome, predictor, sep=' ~ '))

make_rx_fit_qmax <- function(df) {
  rs_min_spec <- 
    linear_reg() %>% 
    set_engine("lme",
               random = ~ 1 + q_max | site_code, 
               corr = corAR1(form = ~ 1 | site_code))
  
  x_fit <- rs_min_spec %>% fit(as.formula(df$x_formula), data = d)
  return(x_fit) }

make_slope_plot_qmax <- function(fit_model_engine) {
  my_label <- paste(fit_model_engine$terms[[2]], 'vs', 'Q-Max', sep=' ') %>%
    str_replace_all('ln_', '') %>%
    str_replace_all('sqrt_', '') %>%
    str_to_title() %>%
    str_replace_all('_', ' ') %>%
    str_replace_all('Vs', 'vs')
  
  my_output <- fit_model_engine$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column('site_code') %>%
    add_rain() %>%
    ggplot(aes(x=annualrain, y=site_code.q_max)) +
    geom_hline(color='grey40', yintercept=0, lwd=.3) +
    geom_point(size=4, color = 'black') +
    stat_poly_eq(label.x = "right",
                 label.y = "bottom",
                 use_label(c("R2", "P"))) +
    labs(y=my_label,
         x='Annual Rainfall (cm)')
  
  # use R^2 to add or leave out line
  my_r2 <- make_slope_lm_qmax(fit_model_engine) %>% pull(lms_r.squared)
  
  if (my_r2 > 0.9) {
    my_output <- my_output + geom_smooth(method='lm', se=F, lty=2, lwd=.5)
    return(my_output) } 
  else 
    {return(my_output)}   }

extract_anova_p <- function(fit_model_engine) {
  fit_model_engine %>%
    anova() %>% 
    rownames_to_column('coef') %>%
    filter(coef!='(Intercept)') %>% 
    pull(`p-value`) }

make_slope_lm_qmax <- function(fit_model_engine) {
  d_slope <- fit_model_engine$coefficients$random %>%
    as.data.frame() %>%
    rownames_to_column('site_code') %>%
    add_rain()
  
  my_lm <- linear_reg() %>% 
    set_engine("lm") %>%
    fit(site_code.q_max ~ annualrain, data = d_slope )%>% 
    extract_fit_engine()
  
  my_est <- my_lm %>% tidy() %>% filter(term=='annualrain') %>% pull(estimate)
  
  my_output <- my_lm %>% glance() %>% select(r.squared, AIC, p.value) %>%
    mutate(estimate=my_est)
  colnames(my_output) <- paste('lms', colnames(my_output), sep='_')
  return(my_output)   }

# fit models and diagnostics
d_rs_fitted_qmax <- d_rs %>%
  filter(predictor == 'q_max') %>%
  group_by(model_name) %>%
  nest() %>%
  mutate(lme_fit = map(data, make_rx_fit_qmax),
         lme_extract = map(lme_fit, extract_fit_engine),
         lme_anova_p = map_dbl(lme_extract, extract_anova_p),
         lme_anova_sig = ifelse(lme_anova_p < 0.05, 'sig', 'ns'),
         lme_glance = map(lme_extract, broom::glance),
         res_qq = map(lme_extract, my_qq),
         shapiro_p = map_dbl(lme_extract, my_shapiro),
         normal_res = ifelse(shapiro_p<0.05, 'no', 'yes'),
         plot_slope = map(lme_extract, make_slope_plot_qmax),
         lm_slope = map(lme_extract, make_slope_lm_qmax)) %>%
  unnest(lme_glance) %>%
  select(-nobs, - sigma, -logLik, -BIC) %>%
  arrange(AIC) %>%
  unnest(lm_slope)

# merge random slope tables
d_rs_fitted <- full_join(d_rs_fitted_qmin, d_rs_fitted_qmax)

# -----------------------------------------------------------------------------
# final plot
# -----------------------------------------------------------------------------
d_rs_fitted %>%
  filter(lme_anova_sig == 'sig') %>%
  pull(plot_slope)

sig_plots <- d_rs_fitted %>%
  filter(lme_anova_sig == 'sig')

plot_lme_slope <- sig_plots$plot_slope[[1]] +
  sig_plots$plot_slope[[2]]  +
  sig_plots$plot_slope[[3]]

table_lme <- 
  d_rs_fitted %>%
    unnest(data) %>%
    ungroup() %>%
    mutate(corr=ifelse(autoregression=='yes', 'autoregression', 'none')) %>%
    select(-model_name, -outcome_abr, -model_type, -lme_fit, -lme_extract, 
           -res_qq, - plot_slope, - x_formula, -autoregression) %>%
    rename(lms_r2 = lms_r.squared,
           lms_p = lms_p.value) %>%
    select(outcome, predictor, corr, AIC, contains('anova'), shapiro_p, normal_res, lms_estimate, lms_r2, lms_AIC, lms_p) 
    
# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
write_csv(table_lme, 'exploration/output/table_lme.csv')

ggsave('exploration/visualization/lme_randomslope_vs_rainfall.png',
       plot=plot_lme_slope,
       width=12,
       height=4.5,
       units='in')

# -----------------------------------------------------------------------------
# End 34_lme_final