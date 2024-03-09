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
tidymodels_prefer()
theme_set(dark_theme_grey(base_size = 14))

d_fish <- read_csv('exploration/output/model_combined_vars.csv')

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

# reserve hurricane data for predicting expected values
d_hurricane <- d_fish %>% filter(is_baseline == 'test') %>% select(-is_baseline)

# cross-validation sets
d_folds <- vfold_cv(d_train, v= 10, repeats=5)
d_boots <- bootstraps(d_fish %>% 
                        filter(is_baseline=='baseline'), times=100)

# -----------------------------------------------------------------------------
# stepwise regression: Identify good predictors
# -----------------------------------------------------------------------------
rich_full <- lm(richness ~ annualrain + ammonia + total_algae + conductivity +
                  d_oxygen + depth_mx + doc + gravel + ln_q_CDC + ln_q_FS +
                  ln_q_HFPP15+ ln_q_HFPP7 + ln_q_WCDC + ln_q_max + ln_q_med + 
                  ln_q_min + ln_q_mu + ln_q_q25 + ln_q_q75 + ln_q_rsd + 
                  ln_q_sd + ln_rain + nitrate + ph + phosphate + q_DS + 
                  q_HFPP3 + q_LFPP + sand + silt + temperature + total_algae +
                  turbidity + width,
                data= d_fish %>% filter(is_baseline=='baseline'))

step <- stepAIC(rich_full, direction='both')
step$coefficients

aprio_vars <- c('annualrain', 
                'conductivity', 'nitrate', 'd_oxygen', 
                'ln_q_med', 'ln_q_rsd', 'temperature', 'ln_rain')

all_vars <- c('annualrain', 'gravel', 'ln_q_CDC', 'ln_q_WCDC',  
              'q_LFPP', 'turbidity', 'width')

baseline_vars <- c('annualrain', 'gravel', 'ln_q_CDC', 'ln_q_HFPP15', 'ln_q_HFPP7', 
                   'ln_q_WCDC', 'ln_q_med', 'ln_q_min', 'ln_q_q25', 'nitrate', 
                   'phosphate', 'q_DS', 'q_LFPP', 'silt', 'width')

test_vars<-c('annualrain', 'conductivity', 'doc', 'ln_q_FS', 'ln_q_WCDC', 
             'ln_q_max', 'ln_rain', 'nitrate', 'ph', 'q_DS', 
             'q_HFPP3', 'silt', 'turbidity', 'width')

# -----------------------------------------------------------------------------
# Preprocessing
# -----------------------------------------------------------------------------
# Recipes
rec_slim <- recipe(richness ~ gravel + ln_q_CDC + 
                     ln_q_WCDC + q_LFPP + turbidity + width,
                    data = d_train) 

rec_baseline <- recipe(richness ~ annualrain + gravel + ln_q_CDC + 
                         ln_q_HFPP15 + ln_q_HFPP7 + 
                         ln_q_WCDC + ln_q_med + ln_q_min + ln_q_q25 + nitrate + 
                         phosphate + q_DS + q_LFPP + silt + width,
                       data = d_train)

rec_baseline_pca_q <- rec_baseline %>%
  step_pca(contains('q_'), num_comp=2)

# Parsnip Models
mod_lm <- linear_reg() %>% set_engine('lm')

mod_rf <- rand_forest(trees=1000) %>% 
  set_engine('ranger') %>% 
  set_mode('regression')

# -----------------------------------------------------------------------------
# model workflow
# -----------------------------------------------------------------------------
# create workflow set
preproc <- list(slim = rec_slim,
                baseline = rec_baseline,
                pca_q = rec_baseline_pca_q)

lm_models <- workflow_set(preproc, list(lm = mod_lm,
                                        rf = mod_rf), cross=T)

get_model <- function(x) { extract_fit_parsnip(x) %>% tidy() }

keep_pred <- control_resamples(save_pred = T, save_workflow = T)

lm_models <- lm_models %>% workflow_map("fit_resamples", seed=1001, verbose=T, 
                                        resamples=d_folds, control=keep_pred)

# -----------------------------------------------------------------------------
# mixed effect must be pre-processed manually 
# -----------------------------------------------------------------------------
# mod_lme <- linear_reg() %>% set_engine('lmer')
# 
# wrk_lme_baseline <- workflow() %>%
#   add_variables(
#     outcomes = richness,
#     predictors = c(annualrain, gravel, ln_q_CDC, ln_q_HFPP15, ln_q_HFPP7, 
#                      ln_q_WCDC, ln_q_med, ln_q_min, ln_q_q25, nitrate, 
#                      phosphate, q_DS, q_LFPP, silt, width)) %>%
#   add_model(mod_lme, formula = richness ~ (1|annualrain) + 
#               gravel + ln_q_CDC + ln_q_HFPP15 + ln_q_HFPP7 + 
#               ln_q_WCDC + ln_q_med + ln_q_min + ln_q_q25 + nitrate + 
#               phosphate + q_DS + q_LFPP + silt + width)
# 
# wrk_lme_slim <- workflow() %>%
#   add_variables(
#     outcomes = richness,
#     predictors = c(annualrain, gravel, ln_q_CDC, ln_q_WCDC, q_LFPP, 
#                    turbidity, width)) %>%
#   add_model(mod_lme, formula = richness ~ (1|annualrain) + 
#               gravel + ln_q_CDC + ln_q_WCDC +  q_LFPP)
# 
# # add mixed effects model
# res_lme_baseline <- wrk_lme_baseline %>% 
#   fit_resamples(resamples=d_folds, control = keep_pred)
# lm_models <- as_workflow_set(baseline_lme = res_lme_baseline) %>% 
#   bind_rows(lm_models)
# 
# res_lme_slim <- wrk_lme_slim %>% 
#   fit_resamples(resamples=d_folds, control = keep_pred)
# lm_models <- as_workflow_set(slim_lme = res_lme_slim) %>% 
#   bind_rows(lm_models)

# -----------------------------------------------------------------------------
# Performance Assessment
# -----------------------------------------------------------------------------
# Performance Statistic Comparison
table_PSC <- collect_metrics(lm_models)

plot_PSC_rsq <- autoplot(lm_models, metric = 'rsq') +
  geom_text_repel(aes(label=wflow_id), nudge_x=1/4, nudge_y=1/100)

# Performance Statistic Correlation within Resamples
rsq_indiv_estimates <- collect_metrics(lm_models, summarize = F) %>%
  filter(.metric == 'rsq') %>%
  mutate(id3 = paste(str_replace_all(id, 'epeat', ''),
                     str_replace_all(id2, 'old', ''),
                     sep='_'))
rsq_wider <- rsq_indiv_estimates %>% select(wflow_id, .estimate, id3) %>%
  pivot_wider(id_cols='id3', names_from='wflow_id', values_from='.estimate')

table_PS_rsq_corr <- corrr::correlate(rsq_wider%>% select(-id3), quiet=T)

plot_PS_rsq_corr <- rsq_indiv_estimates %>%
  mutate(wflow_id = reorder(wflow_id, .estimate)) %>%
  mutate(is_pca_q=ifelse(str_detect(wflow_id, 'baseline'), 
                         'baseline', 'pca_q')) %>%
  ggplot(aes(x=wflow_id, y=.estimate, group=id3, color = id2)) +
  facet_wrap(~is_pca_q, scales='free') +
  geom_line(alpha=.4, lwd=1.2) +
  theme(legend.position='none')

# -----------------------------------------------------------------------------
# Hypothesis Testing: formal model assessment
# -----------------------------------------------------------------------------
# Bayesian: Posterior distributions for the coefficient of determination
rsq_anova <- perf_mod(lm_models, metric='rsq', 
                      prior_intercept=rstanarm::student_t(df=1),
                      chains=4, iter = 10000, seed=1102)

model_post <- rsq_anova %>% tidy(seed=1103)

# glimpse(model_post)
plot_coef_posterior_rsq <- model_post %>%
  mutate(model = forcats::fct_inorder(model)) %>%
  ggplot(aes(x = posterior)) +
  geom_histogram(bins=50, color = 'white', fill='skyblue2', alpha=.5) +
  facet_wrap(~ model, ncol=1)

plot_PSC_bayes_rsq <- autoplot(rsq_anova) +
  geom_text_repel(aes(label=workflow), nudge_x=1/4, nudge_y=1/100) +
  theme(legend.position='none')

# # Bayesian: compare 2 models
# rqs_diff <- contrast_models(rsq_anova,
#                             list_1 = 'pca_q_lm',
#                             list_2 = 'baseline_lme',
#                             seed = 1104)
# 
# plot_lm_Vs_lme <- rqs_diff %>%
#   as_tibble() %>%
#   ggplot(aes(x=difference)) +
#   geom_vline(xintercept=0, lty=2) +
#   geom_histogram(bins=50, color = 'skyblue2', fill = 'skyblue3', alpha = 0.4)
# 
# table_lm_Vs_lme <- summary(rqs_diff) %>% select(-starts_with('pract'))

plot_prob_prac_equiv <- autoplot(rsq_anova, type='ROPE', size =0.02) +
  geom_text_repel(aes(label=workflow)) +
  theme(legend.position = 'none')

# -----------------------------------------------------------------------------
# extract predictions to visualize model fit and residuals
# -----------------------------------------------------------------------------
wflow_slim_lm <- lm_models %>%  extract_workflow(id = 'slim_lm')
ctrl <- control_resamples(extract=get_model)
res_slim_lm <- wflow_slim_lm %>% fit_resamples(resamples = d_folds, control=ctrl)
res_slim_lm$.extracts[[1]][[1]]
all_coef_slim_lm <- map_dfr(res_slim_lm$.extracts, ~ .x[[1]][[1]])

plot_box_coeff_lm_slim <- all_coef_slim_lm %>%
  ggplot(aes(term, estimate)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2, color = 'grey50') +
  xlab(element_blank())

(plot_PSC_bayes_rsq + plot_prob_prac_equiv) / plot_box_coeff_lm_slim

# -----------------------------------------------------------------------------
# model tuning: grid search
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# model tuning: iterative search
# -----------------------------------------------------------------------------
