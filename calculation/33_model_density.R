# model_density
# Sean Kinard
# 2023-07-10
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

d_fish <- read_csv('exploration/output/model_combined_vars.csv') %>%
  mutate(ln_density = log(density))

# reserve hurricane data for predicting expected values
d_hurricane <- d_fish %>% filter(is_baseline == 'test') %>% select(-is_baseline)
d_baseline <- d_fish %>% filter(is_baseline == 'baseline') %>% select(-is_baseline)

# -----------------------------------------------------------------------------
# split data
# -----------------------------------------------------------------------------
# split baseline into training and testing sets
d_split <- d_baseline %>% 
  initial_split(prop=0.80, strata = site_code)
d_train <- training(d_split)
d_test <- testing(d_split)

# cross-validation sets
d_folds <- vfold_cv(d_train, v= 10, repeats=5)
d_boots <- bootstraps(d_fish %>% filter(is_baseline=='baseline'), times=100)

# -----------------------------------------------------------------------------
# stepwise regression: Identify good predictors
# -----------------------------------------------------------------------------
# rich_full <- lm(ln_density ~ annualrain + ammonia + total_algae + conductivity +
#                   d_oxygen + depth_mx + doc + gravel + ln_q_CDC + ln_q_FS +
#                   ln_q_HFPP15+ ln_q_HFPP7 + ln_q_WCDC + ln_q_max + ln_q_med +
#                   ln_q_min + ln_q_mu + ln_q_q25 + ln_q_q75 + ln_q_rsd +
#                   ln_q_sd + ln_rain + nitrate + ph + phosphate + q_DS +
#                   q_HFPP3 + q_LFPP + sand + silt + temperature + total_algae +
#                   turbidity + width,
#                 data= d_fish%>% filter(is_baseline=='baseline'))
# 
# step <- stepAIC(rich_full, direction='both')
# step$coefficients %>% colnames()

apriori <- c('annualrain', 'conductivity', 'nitrate', 'd_oxygen', 
             'ln_q_med', 'ln_q_rsd', 'temperature', 'ln_rain')

astep<-c("annualrain", "gravel", "ln_q_CDC", "ln_q_HFPP15", "ln_q_med", 
          "ln_q_sd", "ln_rain", "nitrate", "ph", "q_HFPP3", "temperature", 
          "width")

bstep<-c("annualrain", "total_algae", "conductivity", "depth_mx", 
               "doc", "ln_q_CDC", "ln_q_HFPP7", "ln_q_med", "ln_q_min", 
               "ln_q_mu", "ln_q_q75", "ln_q_rsd", "ln_rain", "nitrate", 
               "phosphate", "q_DS", "q_HFPP3", "sand", "silt", "temperature", 
               "width")

all_pred <- c("annualrain", "ammonia", "bluegreen_cyano", "canopy", "conductivity", "d_oxygen", "depth_mx", "diatoms", "doc", "gravel", "green_algae", "ln_q_CDC", "ln_q_FS", "ln_q_HFPP15", "ln_q_HFPP7", "ln_q_WCDC", "ln_q_max", "ln_q_med", "ln_q_min", "ln_q_mu", "ln_q_q25", "ln_q_q75", "ln_q_rsd", "ln_q_sd", "ln_rain", "nitrate", "ph", "phosphate", "q_DS", "q_HFPP3", "q_LFPP", "sand", "silt", 
"temperature", "total_algae", "turbidity", "width")

# -----------------------------------------------------------------------------
# Set outcomes and predictors
# -----------------------------------------------------------------------------
# no preprocessing
vars_apriori <- workflow_variables(outcomes=ln_density, 
                                   predictors = any_of(apriori))
vars_astep <- workflow_variables(outcomes=ln_density, 
                                  predictors = any_of(astep))
vars_bstep <- workflow_variables(outcomes=ln_density, 
                                       predictors = any_of(bstep))
vars_all <- workflow_variables(outcomes=ln_density, 
                               predictors = any_of(all_pred))

# normalized
norm_apriori <- vars_apriori %>% step_normalize(all_predictors())
norm_astep <- vars_astep %>% step_normalize(all_predictors())
norm_bstep <- vars_bstep %>% step_normalize(all_predictors())
norm_all <- vars_all %>% step_normalize(all_predictors())

# quadratic and 2-way interactions
poly_apriori <- norm_apriori %>% step_poly(all_predictors()) %>%
  step_interact(~ all_predictors():all_predictors())
poly_astep <- norm_astep %>% step_poly(all_predictors()) %>%
  step_interact(~ all_predictors():all_predictors())
poly_bstep <- norm_bstep %>% step_poly(all_predictors())%>%
  step_interact(~ all_predictors():all_predictors())
poly_all <- norm_all %>% step_poly(all_predictors())%>%
  step_interact(~ all_predictors():all_predictors())
# -----------------------------------------------------------------------------
# model engines
# -----------------------------------------------------------------------------
linear_reg_spec <- 
  linear_reg(penalty=tune(), mixture=tune()) %>%
  set_engine('glmnet')

nnet_spec <-
  mlp(hidden_units=tune(), penalty=tune(), epochs=tune()) %>%
  set_engine("nnet", MaxNwts=2600) %>%
  set_mode('regression')

mars_spec <- 
  mars(prod_degree=tune()) %>% 
  set_engine('earth') %>%
  set_mode('regression')

svm_r_spec <- 
  svm_rbf(cost = tune(), rbf_sigma=tune()) %>%
  set_engine('kernlab') %>%
  set_mode('regression')

svm_p_spec <- 
  svm_poly(cost = tune(), degree=tune()) %>%
  set_engine('kernlab') %>%
  set_mode('regression')

knn_spec <- 
  nearest_neighbor(neighbors=tune(), dist_power=tune(), weight_func=tune()) %>%
  set_engine('kknn') %>%
  set_mode('regression')

cart_spec <- 
  decision_tree(cost_complexity=tune(), min_n=tune()) %>%
  set_engine('rpart') %>%
  set_mode('regression')

bag_cart_spec <- 
  bag_tree() %>%
  set_engine('rpart', times=50L) %>%
  set_mode('regression')

rf_spec <-
  rand_forest(mtry=tune(), min_n=tune(), trees=1000) %>%
  set_engine('ranger') %>%
  set_mode('regression')

xgb_spec <-
  boost_tree(tree_depth=tune(), learn_rate=tune(), loss_reduction=tune(),
             min_n=tune(), sample_size=tune(), trees=tune()) %>%
  set_engine('xgboost') %>%
  set_mode('regression')

cubist_spec <- 
  cubist_rules(committees=tune(), neighbors=tune()) %>%
  set_engine('Cubist')

# -----------------------------------------------------------------------------
# workflow set
# -----------------------------------------------------------------------------
no_pre_proc <-
  workflow_set(
    preproc = list(apr = vars_apriori,
                   #astep = vars_astep,
                   bas = vars_bstep,
                   all = vars_all),
    models = list(#MARS = mars_spec,
                  #CART = cart_spec,
                  #CART_bagged = bag_cart_spec,
                  #Cubist = cubist_spec,
                  RF = rf_spec,
                  XGB = xgb_spec) )

normalized <-
  workflow_set(
    preproc = list(aprN = norm_apriori,
                   #astepN = norm_astep,
                   basN = norm_bstep,
                   allN = norm_all),
    models = list(#neural_network = nnet_spec,
                  SVM_radial = svm_r_spec,
                  SVM_poly = svm_p_spec,
                  KNN = knn_spec  ) )

with_features <- 
  workflow_set(
    preproc = list(aprP = poly_apriori,
                   #astepP = poly_astep,
                   basP = poly_bstep,
                   allP = poly_all),
    models = list(GLM = linear_reg_spec,
                  KNN = knn_spec) )

all_workflows <-
  bind_rows(no_pre_proc, normalized, with_features)


# -----------------------------------------------------------------------------
# mixed effect must be pre-processed manually
# -----------------------------------------------------------------------------
get_model <- function(x) { extract_fit_parsnip(x) %>% tidy() }
keep_pred <- control_resamples(save_pred = T, save_workflow = T)
mod_lme <- linear_reg() %>% set_engine('lmer')

res_lme_apriori <- workflow() %>%
  add_variables(
    outcomes = ln_density,
    predictors = any_of(apriori)) %>%
  add_model(mod_lme, formula = ln_density ~ (1|annualrain) +
              conductivity + nitrate + d_oxygen + ln_q_med + ln_q_rsd + 
              temperature + ln_rain) %>%
  fit_resamples(resamples=d_folds, control = keep_pred)

res_lme_astep <- workflow() %>%
  add_variables(
    outcomes = ln_density,
    predictors = any_of(astep)) %>%
  add_model(mod_lme, formula = ln_density ~ (1|annualrain) +
              gravel + ln_q_CDC + ln_q_HFPP15 + ln_q_med + ln_q_sd + 
              ln_rain + nitrate + ph + q_HFPP3 + temperature + width) %>%
  fit_resamples(resamples=d_folds, control = keep_pred)

res_lme_bstep <- workflow() %>%
  add_variables(
    outcomes = ln_density,
    predictors = any_of(bstep)) %>%
  add_model(mod_lme, formula = ln_density ~ (1|annualrain) +
              total_algae + conductivity + depth_mx + 
            doc + ln_q_CDC + ln_q_HFPP7 + ln_q_med + ln_q_min + 
            ln_q_mu + ln_q_q75 + ln_q_rsd + ln_rain + nitrate + 
            phosphate + q_DS + q_HFPP3 + sand + silt + temperature + 
            width) %>%
  fit_resamples(resamples=d_folds, control = keep_pred)

# add mixed effects model
all_workflows <- as_workflow_set(
  apr_LME = res_lme_apriori,
  bstep_LME = res_lme_bstep,
  astep_LME = res_lme_astep) %>%
  bind_rows(all_workflows)

# -----------------------------------------------------------------------------
# Tune and Evaluate Models
# -----------------------------------------------------------------------------
race_ctrl <- 
  control_race(save_pred=T, parallel_over = 'everything', save_workflow=T)

race_results <-
  all_workflows %>% workflow_map(
    'tune_race_anova', 
    seed=1503, 
    resamples = d_folds, 
    grid = 25, 
    control = race_ctrl )

plot_rank_rmse <- 
  autoplot( race_results, rank_metric = 'rmse', metric = 'rmse', select_best=T) +
  geom_text(aes(y=mean-1/8, label=wflow_id), hjust=1) +
  theme(legend.position = 'none') +
  ylim(c(0.5,1.5)) +
  coord_flip()

plot_rank_rsq <- 
  autoplot( race_results, rank_metric = 'rsq', metric = 'rsq', select_best=T) +
  geom_text(aes(y=mean-1/12, label=wflow_id), hjust=1) +
  theme(legend.position = 'none') +
  ylim(c(0.0,0.7)) +
  coord_flip()

(plot_rank_both <- plot_rank_rmse + plot_rank_rsq)

race_results %>%
  rank_results(rank_metric = 'rmse', select_best=T) %>%
  filter(.metric == 'rmse') %>%
  select(wflow_id, model, .config, rmse=mean, rank)

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
    fit(d_baseline)
}

visualize_best_fit_test <- function(wfID) {
  x_best_fit <- make_best_fit(wfID)
  
  x_best_fit %>%
    collect_predictions() %>%
    ggplot(aes(x=ln_density, y=.pred)) +
    geom_abline(lty=2) +
    geom_point(alpha=0.5) +
    coord_obs_pred() +
    labs(x='Observed', y = 'Predicted') +
    ggtitle(paste('Validation Test:', wfID)) }

visualize_expected_v_observed <- function(wfID) {
  x_best_fit <- make_best_fit_baseline(wfID)
  
  d_fish %>% 
    mutate(prediction = predict(x_best_fit, d_fish) %>% as_vector(),
           is_baseline = ifelse(is_baseline == 'test', 'Hurricane', 'Baseline'),
           my_labels = ifelse(is_baseline=='Hurricane', month, NA)) %>%
    ggplot(aes(x=ln_density, y=prediction, color=is_baseline, fill=is_baseline)) +
    geom_abline(lty=2) +
    geom_point(size = 2, shape = 21, alpha =.3) +
    geom_point(size = 2, shape = 21, fill=NA) +
    coord_obs_pred() +
    labs(x='Observed', y = 'Predicted') +
    ggtitle(wfID)}

visualize_time_series <- function(wfID) { 
  x_best_fit <- make_best_fit_baseline(wfID)
  
  d_fish %>% 
    mutate(prediction = predict(x_best_fit, d_fish) %>% as_vector()) %>%
    add_hurricane() %>%
    mutate(ln_hu_flood_ratio = log(hu_flood_ratio)) %>%
    ggplot(aes(x=collection_period, y=(ln_density-prediction), 
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
# Mixed Effect (lme)
best_lme <- make_best('bstep_LME')
best_lme_fit <- make_best_fit('bstep_LME')
best_baseline_fit_lme <- make_best_fit_baseline('bstep_LME')
plot_test_best_lme <- visualize_best_fit_test('bstep_LME')
plot_EvO_density_lme <- visualize_expected_v_observed('bstep_LME')
plot_timeseries_density_lme <- visualize_time_series('bstep_LME')

# LASSO Penalized glm (glm)
best_glm <- make_best('basP_GLM')
best_glm_fit <- make_best_fit('basP_GLM')
best_baseline_fit_glm <- make_best_fit_baseline('basP_GLM')
plot_test_best_glm <- visualize_best_fit_test('basP_GLM')
plot_EvO_density_glm <- visualize_expected_v_observed('basP_GLM')
plot_timeseries_density_glm <- visualize_time_series('basP_GLM')
extract_fit_parsnip(best_glm_fit) %>% tidy() %>% arrange(desc(abs(estimate))) %>% 
  print(n=20)

# Support Vector Machine (svm)
best_svm <- make_best('svm_rbf')
best_svm_fit <- make_best_fit('svm_rbf')
best_baseline_fit_svm <- make_best_fit_baseline('svm_rbf')
plot_test_best_svm <- visualize_best_fit_test('svm_rbf')
plot_EvO_density_svm <- visualize_expected_v_observed('svm_rbf')
plot_timeseries_density_svm <- visualize_time_series('svm_rbf')

# Random Forest (rf)
best_rf <- make_best('bas_RF')
best_rf_fit <- make_best_fit('bas_RF')
best_baseline_fit_rf <- make_best_fit_baseline('bas_RF')
plot_test_best_rf <- visualize_best_fit_test('bas_RF')
plot_EvO_density_rf <- visualize_expected_v_observed('bas_RF')
plot_timeseries_density_rf <- visualize_time_series('bas_RF')

# Boosted Trees Ensemble (xgb)
best_xgb <- make_best('all_XGB')
best_xgb_fit <- make_best_fit('all_XGB')
best_baseline_fit_xgb <- make_best_fit_baseline('all_XGB')
plot_test_best_xgb <- visualize_best_fit_test('all_XGB')
plot_EvO_density_xgb <- visualize_expected_v_observed('all_XGB')
plot_timeseries_density_xgb <- visualize_time_series('all_XGB')
plot_xgboost_important_vars <- best_xgb_fit %>% 
  extract_fit_parsnip() %>% 
  vip() +
  ggtitle("Variable Importance for Xgboost model")
# xgboost interpretting predictor effects:
# library(DALEXtra)
# explained_xgboost <-  explain_tidymodels(
#   model =    best_xgboost_fit,
#   data =  d_train %>% select(-ln_density) ,
#   y =  d_train$ln_density,
#   type = "regression",
#   predict_function_target_column = 53,
#   verbose = F)

#------------------------------------------------------------------------------
# Export figures
#------------------------------------------------------------------------------
my_objects<-ls()
my_figure_names <- my_objects[str_detect(my_objects,'plot_')]

my_figures<-list(
plot_EvO_density_glm, plot_EvO_density_lme, plot_EvO_density_rf, plot_EvO_density_svm, plot_EvO_density_xgb, plot_rank_both, 
plot_rank_rmse, plot_rank_rsq, plot_test_best_glm, plot_test_best_lme, plot_test_best_rf, plot_test_best_svm, plot_test_best_xgb, plot_timeseries_density_glm, plot_timeseries_density_lme, plot_timeseries_density_rf, plot_timeseries_density_svm, plot_timeseries_density_xgb, plot_xgboost_important_vars)

names(my_figures) <- str_replace_all(my_figure_names, 'plot_', '33_density_')

for (i in 1:length(my_figures)) {
  my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
  my_object <- my_figures[[i]]
  ggsave(my_place,
         plot = my_object,
         width = 10,
         height = 10,
         units = c("in")) }

#------------------------------------------------------------------------------
# End 33_model_density