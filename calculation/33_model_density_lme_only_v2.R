# model_density_lme_only_v2
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
# fit lme to data
# -----------------------------------------------------------------------------
res_dens_minq_rslope <- workflow() %>%
  add_variables(
    outcomes = z_ln_density,
    predictors = any_of(c('site_code', 'z_ln_wk4_min_to_annual_med'))) %>%
  add_model( linear_reg() %>% set_engine('lmer'), 
             formula = z_ln_density ~ 
              z_ln_wk4_min_to_annual_med +
              (1 + z_ln_wk4_min_to_annual_med|site_code)) %>%
  fit(d)



# -----------------------------------------------------------------------------
# split data
# -----------------------------------------------------------------------------
# split baseline into training and testing sets
d_split <- d %>% 
  initial_split(prop=0.80, strata = site_code)
d_train <- training(d_split)
d_test <- testing(d_split)

# cross-validation sets
d_folds <- vfold_cv(d_train, v= 10, repeats=5)
d_boots <- bootstraps(d, times=100)

# -----------------------------------------------------------------------------
# mixed effect must be pre-processed manually
# -----------------------------------------------------------------------------
get_model <- function(x) { extract_fit_parsnip(x) %>% tidy() }
keep_pred <- control_resamples(save_pred = T, save_workflow = T)
mod_lme <- linear_reg() %>% set_engine('lmer')
make_res <- function(my_predictors, my_response, my_data) {
  formula_x <- reformulate(str_replace_all(my_predictors, 'site_code', 
                                           '(1|site_code)'),
                           response = my_response)
  
  res_lme_x <- workflow() %>%
    add_variables(
      outcomes = {{my_response}},
      predictors = any_of(my_predictors)) %>%
    add_model(mod_lme, formula = formula_x) %>%
    fit_resamples(resamples=my_data, control = keep_pred)
  
  return(res_lme_x) }

res_dens_minq <- make_res(
  my_predictors = c('site_code', 'z_ln_wk4_min_to_annual_med'),
  my_response = 'z_ln_density',
  my_data = d_boots) # Warning help('iSingular')

res_dens_maxq <- make_res(
  my_predictors = c('site_code', 'z_ln_wk4_max_to_annual_med'),
  my_response = 'z_ln_density',
  my_data = d_boots) # Warning help('iSingular')

res_dens_minq_rslope <- workflow() %>%
  add_variables(
    outcomes = z_ln_density,
    predictors = any_of(c('site_code', 'z_ln_wk4_min_to_annual_med'))) %>%
  add_model(mod_lme, formula = z_ln_density ~ 
              z_ln_wk4_min_to_annual_med +
              (1 + z_ln_wk4_min_to_annual_med|site_code)) %>%
  fit_resamples(resamples=d_boots, control = keep_pred)

res_dens_maxq_rslope <- workflow() %>%
  add_variables(
    outcomes = z_ln_density,
    predictors = any_of(c('site_code', 'z_ln_wk4_max_to_annual_med'))) %>%
  add_model(mod_lme, formula = z_ln_density ~ 
              z_ln_wk4_max_to_annual_med +
              (1 + z_ln_wk4_max_to_annual_med|site_code)) %>%
  fit_resamples(resamples=d_boots, control = keep_pred)

# workflow set mixed effects models
all_workflows <- as_workflow_set(
   int_dens_minq = res_dens_minq,
   int_dens_maxq = res_dens_maxq,
   rsl_dens_minq = res_dens_minq_rslope,
   rsl_dens_maxq = res_dens_minq_rslope)

# -----------------------------------------------------------------------------
# Tune and Evaluate Models
# -----------------------------------------------------------------------------
race_ctrl <- 
  control_race(save_pred=T, parallel_over = 'everything', save_workflow=T)

race_results <-
  all_workflows %>% workflow_map(
    'tune_race_anova', 
    seed=1503, 
    resamples = d_boots, 
    grid = 25, 
    control = race_ctrl )

plot_rank_rmse <- 
  autoplot( race_results, rank_metric = 'rmse', metric = 'rmse', select_best=T) +
  geom_text(aes(y=mean-1/8, label=wflow_id), hjust=1) +
  theme(legend.position = 'none') +
  ylim(c(0.5,1.25)) +
  coord_flip()

plot_rank_rsq <- 
  autoplot( race_results, rank_metric = 'rsq', metric = 'rsq', select_best=T) +
  geom_text(aes(y=mean-1/12, label=wflow_id), hjust=1) +
  theme(legend.position = 'none') +
  ylim(c(0.0,0.99)) +
  coord_flip()

(plot_rank_both <- plot_rank_rmse + plot_rank_rsq)

race_results %>%
  rank_results(rank_metric = 'rmse', select_best=T) %>%
  filter(.metric == 'rmse') %>%
  select(wflow_id, model, .config, rmse=mean, rank) %>%
  print(n=20)

# -----------------------------------------------------------------------------
# Post Processing: functions
# -----------------------------------------------------------------------------
make_best <- function(wfID) {
  race_results %>% 
    extract_workflow_set_result(wfID) %>% 
    select_best(metric = 'rmse') }

make_best_fit <- function(wfID) {
  x_best <- make_best(wfID)
  
  race_results %>%
    extract_workflow(wfID) %>%
    finalize_workflow(x_best) %>%
    last_fit(split = d_split) }

make_best_fit_baseline <- function(wfID) {
  x_best <- make_best(wfID)
  
  race_results %>%
    extract_workflow(wfID) %>%
    finalize_workflow(x_best) %>%
    fit(d)
}

visualize_best_fit_test <- function(wfID) {
  x_best_fit <- make_best_fit(wfID)
  
  x_best_fit %>%
    collect_predictions() %>%
    ggplot(aes(x=z_ln_density, y=.pred)) +
    geom_abline(lty=2) +
    geom_point(alpha=0.5) +
    coord_obs_pred() +
    labs(x='Observed', y = 'Predicted') +
    ggtitle(paste('Validation Test:', wfID)) }

visualize_expected_v_observed <- function(wfID) {
  x_best_fit <- make_best_fit_baseline(wfID)
  
  d %>% 
    mutate(prediction = predict(x_best_fit, d) %>% as_vector(),
           is_baseline = ifelse(is_baseline == 'test', 'Hurricane', 'Baseline'),
           my_labels = ifelse(is_baseline=='Hurricane', month, NA)) %>%
    ggplot(aes(x=z_ln_density, y=prediction, color=is_baseline, fill=is_baseline)) +
    geom_abline(lty=2) +
    geom_point(size = 2, shape = 21, alpha =.3) +
    geom_point(size = 2, shape = 21, fill=NA) +
    coord_obs_pred() +
    labs(x='Observed', y = 'Predicted') +
    ggtitle(wfID)}

visualize_time_series <- function(wfID) { 
  x_best_fit <- make_best_fit_baseline(wfID)
  
  d %>% 
    mutate(prediction = predict(x_best_fit, d) %>% as_vector()) %>%
    add_hurricane() %>%
    mutate(ln_hu_flood_ratio = log(hu_flood_ratio)) %>%
    ggplot(aes(x=collection_period, y=(z_ln_density-prediction), 
               color = ln_hu_flood_ratio, 
               fill = ln_hu_flood_ratio)) +
    geom_hline(yintercept=0, lty=2, lwd=.2, color='white') +
    geom_vline(xintercept=ymd('2017-08-27'), color = 'red', lty=2, lwd=.2) +
    geom_smooth(aes(group=NA),
                method='loess', span=.3, color=NA, lty=1, lwd=.4,
                show.legend = F)+
    geom_point(size = 4, shape = 21, alpha =.3) +
    geom_point(size = 4, shape = 21, fill=NA) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction =-1) +
    scale_x_date(date_breaks = "4 month", date_labels = "%y:%m") +
    labs(y='Observed - Expected', x ='Year : Month') +
    ggtitle(wfID) }

# -----------------------------------------------------------------------------
# Post Processing: iterations
# -----------------------------------------------------------------------------
# Mixed Effect (int_dens_maxq)
best_lme <- make_best('int_dens_maxq')
best_lme_fit <- make_best_fit('int_dens_maxq')
plot_test_best_lme <- visualize_best_fit_test('int_dens_maxq')
plot_timeseries_density_lme <- visualize_time_series('int_dens_maxq')
extract_fit_parsnip(best_lme_fit) %>% tidy()

extract_fit_parsnip(best_lme_fit) %>% summary()

# Mixed Effect (rsl_dens_minq)
best_lme <- make_best('rsl_dens_minq')
best_lme_fit <- make_best_fit('rsl_dens_minq')
plot_test_best_lme <- visualize_best_fit_test('rsl_dens_minq')
plot_EvO_density_lme <- visualize_expected_v_observed('rsl_dens_minq')
plot_timeseries_density_lme <- visualize_time_series('rsl_dens_minq')
extract_fit_parsnip(best_lme_fit) %>% tidy()

# Mixed Effect (rsl_dens_maxq)
best_lme <- make_best('rsl_dens_maxq')
best_lme_fit <- make_best_fit('rsl_dens_maxq')
plot_test_best_lme <- visualize_best_fit_test('rsl_dens_maxq')
plot_EvO_density_lme <- visualize_expected_v_observed('rsl_dens_maxq')
plot_timeseries_density_lme <- visualize_time_series('rsl_dens_maxq')
extract_fit_parsnip(best_lme_fit) %>% tidy()

# Mixed Effect (int_dens_minq)
best_lme <- make_best('int_dens_minq')
best_lme_fit <- make_best_fit('int_dens_minq')
plot_test_best_lme <- visualize_best_fit_test('int_dens_minq')
plot_EvO_density_lme <- visualize_expected_v_observed('int_dens_minq')
plot_timeseries_density_lme <- visualize_time_series('int_dens_minq')
extract_fit_parsnip(best_lme_fit) %>% tidy()

#------------------------------------------------------------------------------
# Export figures
#------------------------------------------------------------------------------
# my_objects<-ls()
# my_figure_names <- my_objects[str_detect(my_objects,'plot_')]
# 
# my_figures<-list(
#   plot_EvO_density_glm, plot_EvO_density_lme, plot_EvO_density_rf, plot_EvO_density_svm, plot_EvO_density_xgb, plot_rank_both, 
#   plot_rank_rmse, plot_rank_rsq, plot_test_best_glm, plot_test_best_lme, plot_test_best_rf, plot_test_best_svm, plot_test_best_xgb, plot_timeseries_density_glm, plot_timeseries_density_lme, plot_timeseries_density_rf, plot_timeseries_density_svm, plot_timeseries_density_xgb, plot_xgboost_important_vars)
# 
# names(my_figures) <- str_replace_all(my_figure_names, 'plot_', '33_density_rapid_')
# 
# for (i in 1:length(my_figures)) {
#   my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
#   my_object <- my_figures[[i]]
#   ggsave(my_place,
#          plot = my_object,
#          width = 10,
#          height = 10,
#          units = c("in")) }

#------------------------------------------------------------------------------
# End 33_model_density_00mo