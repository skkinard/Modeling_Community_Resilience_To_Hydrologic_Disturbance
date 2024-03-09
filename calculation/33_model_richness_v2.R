# model_data_preparation
# Sean Kinard
# 2023-07-04
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
tidymodels_prefer()
theme_set(dark_theme_grey(base_size = 14))

d_fish <- read_csv('exploration/output/model_combined_vars.csv')

# reserve hurricane data for predicting expected values
d_hurricane <- d_fish %>% filter(is_baseline == 'test') %>% select(-is_baseline)

# -----------------------------------------------------------------------------
# split data
# -----------------------------------------------------------------------------
# split baseline into training and testing sets
d_split <- d_fish %>% 
  filter(is_baseline=='baseline') %>% 
  select(-is_baseline) %>%
  initial_split(prop=0.90, strata = site_code)
d_train <- training(d_split)
d_test <- testing(d_split)

# cross-validation sets
d_folds <- vfold_cv(d_train, v= 10, repeats=5)
d_boots <- bootstraps(d_fish %>% filter(is_baseline=='baseline'), times=100)

# -----------------------------------------------------------------------------
# stepwise regression: Identify good predictors
# -----------------------------------------------------------------------------
# rich_full <- lm(richness ~ annualrain + ammonia + total_algae + conductivity +
#                   d_oxygen + depth_mx + doc + gravel + ln_q_CDC + ln_q_FS +
#                   ln_q_HFPP15+ ln_q_HFPP7 + ln_q_WCDC + ln_q_max + ln_q_med + 
#                   ln_q_min + ln_q_mu + ln_q_q25 + ln_q_q75 + ln_q_rsd + 
#                   ln_q_sd + ln_rain + nitrate + ph + phosphate + q_DS + 
#                   q_HFPP3 + q_LFPP + sand + silt + temperature + total_algae +
#                   turbidity + width,
#                 data= d_fish %>% filter(is_baseline=='baseline'))
# 
# step <- stepAIC(rich_full, direction='both')
# step$coefficients

apriori <- c('annualrain', 'conductivity', 'nitrate', 'd_oxygen', 
                'ln_q_med', 'ln_q_rsd', 'temperature', 'ln_rain')

alldat <- c('annualrain', 'gravel', 'ln_q_CDC', 'ln_q_WCDC',  
              'q_LFPP', 'turbidity', 'width')

baselinedat <- c('annualrain', 'gravel', 'ln_q_CDC', 'ln_q_HFPP15', 'ln_q_HFPP7', 
                   'ln_q_WCDC', 'ln_q_med', 'ln_q_min', 'ln_q_q25', 'nitrate', 
                   'phosphate', 'q_DS', 'q_LFPP', 'silt', 'width')

# -----------------------------------------------------------------------------
# variables
# -----------------------------------------------------------------------------
# no preprocessing
vars_apriori <- workflow_variables(outcomes=richness, 
                                   predictors = any_of(apriori))
vars_alldat <- workflow_variables(outcomes=richness, 
                                   predictors = any_of(alldat))
vars_baselinedat <- workflow_variables(outcomes=richness, 
                                   predictors = any_of(baselinedat))

# normalized
norm_apriori <- vars_apriori %>% step_normalize(all_predictors())
norm_alldat <- vars_alldat %>% step_normalize(all_predictors())
norm_baselinedat <- vars_baselinedat %>% step_normalize(all_predictors())

# quadratic and 2-way interactions
poly_apriori <- norm_apriori %>% step_poly(all_predictors()) %>%
  step_interact(~ all_predictors():all_predictors())
poly_alldat <- norm_alldat %>% step_poly(all_predictors()) %>%
  step_interact(~ all_predictors():all_predictors())
poly_baselinedat <- norm_baselinedat %>% step_poly(all_predictors())%>%
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
                   all = vars_alldat,
                   bas = vars_baselinedat),
    models = list(MARS = mars_spec,
                  CART = cart_spec,
                  CART_bagged = bag_cart_spec,
                  RF = rf_spec,
                  boosting = xgb_spec,
                  Cubist = cubist_spec) )

normalized <-
  workflow_set(
    preproc = list(aprN = norm_apriori,
                   allN = norm_alldat,
                   basN = norm_baselinedat),
    models = list(SVM_radial = svm_r_spec,
                  SVM_poly = svm_p_spec,
                  KNN = knn_spec,
                  neural_network = nnet_spec) )

with_features <- 
  workflow_set(
    preproc = list(aprP = poly_apriori,
                   allP = poly_alldat,
                   basP = poly_baselinedat),
    models = list(linear_reg = linear_reg_spec,
                  KNN = knn_spec) )

all_workflows <-
  bind_rows(no_pre_proc, normalized, with_features)
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
  geom_text(aes(y=mean-1/2, label=wflow_id), hjust=1) +
  theme(legend.position = 'none') +
  coord_flip()

plot_rank_rsq <- 
  autoplot( race_results, rank_metric = 'rsq', metric = 'rsq', select_best=T) +
  geom_text(aes(y=mean-1/12, label=wflow_id), hjust=1) +
  theme(legend.position = 'none') +
  coord_flip()

plot_rank_rmse + plot_rank_rsq

# -----------------------------------------------------------------------------
# Extract Top models
# -----------------------------------------------------------------------------
race_results %>%
  rank_results(rank_metric = 'rmse', select_best=T) %>%
  filter(.metric == 'rmse') %>%
  select(wflow_id, model, .config, rmse=mean, rank)

best_results <-
  race_results %>%
  extract_workflow_set_result('allN_KNN') %>%
  select_best(metric = 'rmse')

allN_KNN_test_results <-
  race_results %>%
  extract_workflow('allN_KNN') %>%
  finalize_workflow(best_results) %>%
  last_fit(split = d_split)

collect_metrics(allN_KNN_test_results)

allN_KNN_test_results %>%
  collect_predictions() %>%
  ggplot(aes(x=richness, y=.pred)) +
  geom_abline(lty=2) +
  geom_point(alpha=0.5) +
  coord_obs_pred() +
  labs(x='Observed', y = 'Predicted')
# -----------------------------------------------------------------------------
allN_KNN_results <-
  race_results %>%
  extract_workflow_set_result('all_RF') %>%
  select_best(metric = 'rmse')

all_RF_test_results <-
  race_results %>%
  extract_workflow('all_RF') %>%
  finalize_workflow(allN_KNN_results) %>%
  last_fit(split = d_split)

collect_metrics(all_RF_test_results)

all_RF_test_results %>%
  collect_predictions() %>%
  ggplot(aes(x=richness, y=.pred)) +
  geom_abline(lty=2) +
  geom_point(alpha=0.5) +
  coord_obs_pred() +
  labs(x='Observed', y = 'Predicted')
