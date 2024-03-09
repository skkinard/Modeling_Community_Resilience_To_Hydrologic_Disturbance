# 35_sem
# Sean Kinard
# 2023-08-10
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')
library(tidySEM)
library(lavaan)
theme_set(theme_bw(base_size = 14))

# Data
d_outcome <- read_csv('exploration/output/model_outcome_vars.csv') %>%
  select(-qtr)

d_predictor <- read_csv('exploration/output/model_predictor_vars.csv') %>%
  mutate(flood_days = exp(ln_flood_total_days),
         low_flow_days = exp(ln_lf_days))

d_maj_flood <- read_csv(
  'data/06_fill_data/site_flow_4week_flood_duration.csv') %>%
  select(site_code, collection_period, contains('is_'))

# vectors
my_index <- c('site_code', 'collection_period')
my_predictors <- c('annualrain', 'sqrt_rain', 'ln_wk4_max', 'ln_wk4_min',
                   'low_flow_days', 'flood_days')
# -----------------------------------------------------------------------------
# data prep
# -----------------------------------------------------------------------------
# predictor inspection
d_test <- d_predictor %>%
  select(any_of(my_index), any_of(my_predictors)) %>%
  select(-annualrain) %>%
  pivot_longer(cols=-any_of(my_index),
               names_to='predictor', values_to='xvalue') 
plot_hist_predictors <- d_test %>%
  ggplot(aes(xvalue)) +
  geom_histogram(bins=50) +
  facet_wrap(~predictor)

# Log Response Ratio
d_out_long <- d_outcome%>%select(-simpson) %>%
  pivot_longer(cols=-any_of(my_index),
               names_to='outcome', values_to='xvalue') 
d_LRR <- d_out_long %>%
  group_by(site_code, outcome) %>%
  dplyr::summarize(mu = mean(xvalue, na.rm=T)) %>%
  right_join(d_out_long) %>%
  mutate(LRR = log(xvalue/mu)) %>%
  select(-mu, -xvalue) %>%
  pivot_wider(names_from=outcome, values_from=LRR, names_prefix = 'LRR_')

# merge
dc <- d_predictor %>% 
  select(any_of(my_index), any_of(my_predictors)) %>%
  left_join(d_LRR)

# filter to flood and drought only
only_flood <- d_maj_flood %>%
  filter(is_maj_flood == 'yes') %>% unique() %>%
  select(site_code, collection_period)
only_drought <- d_maj_flood %>%
  filter(is_drought == 'yes') %>% unique() %>%
  select(site_code, collection_period)
only_drought_flood <- d_maj_flood %>%
  filter(is_drought=='yes' | is_maj_flood == 'yes') %>% unique() %>%
  select(site_code, collection_period)

# model data
dc_f <-left_join(only_flood, dc) # only floods
dc_d <-left_join(only_drought, dc) # only droughts
dc_fd <- left_join(only_drought_flood, dc) # floods and droughts

my_outcomes <- d_LRR %>%
  ungroup() %>%
  select(-site_code, -collection_period) %>%
  colnames()

make_fit_flood <-function(x_formula) {sem(x_formula, data = dc_f)}
make_fit_drought <-function(x_formula) {sem(x_formula, data = dc_d)}
make_smry <- function(x_fit) {summary(x_fit, standardized = TRUE)}

sem_flood <- tibble(outcome=my_outcomes) %>%
  mutate(formula = paste(
    "\nln_wk4_max ~ sqrt_rain + annualrain\nflood_days ~ ln_wk4_max + annualrain\n",
    outcome, 
    " ~ ln_wk4_max + flood_days + annualrain\n", sep=''),
    sem_fit = map(formula, make_fit_flood),
    sem_smry = map(sem_fit, make_smry),
    sem_graph = map(sem_fit, graph_sem),
    dataset = 'flood') %>%
  select(dataset, outcome, contains('sem'), formula)

sem_drought <- tibble(outcome=my_outcomes) %>%
  mutate(formula = paste(
    "\nln_wk4_min ~ sqrt_rain + annualrain\nlow_flow_days ~ ln_wk4_min + annualrain\n",
    outcome, 
    " ~ ln_wk4_min + low_flow_days + annualrain\n", sep=''),
    sem_fit = map(formula, make_fit_drought),
    sem_smry = map(sem_fit, make_smry),
    sem_graph = map(sem_fit, graph_sem),
    dataset = 'drought') %>%
  select(dataset, outcome, contains('sem'), formula)

# biomass

sem_flood$sem_smry[[1]]
sem_drought$sem_smry[[1]]
sem_flood$sem_graph[[1]] + 
  ggtitle(paste(
    'Flood Events (p = ', sem_flood$sem_smry[[1]]$test$standard$pvalue%>%round(3), ")",
    sep='')) + 
  sem_drought$sem_graph[[1]] + ggtitle(paste(
    'Drought Events (p = ', sem_drought$sem_smry[[1]]$test$standard$pvalue%>%round(3), ")",
    sep=''))

# centroid
sem_flood$sem_smry[[2]]
sem_drought$sem_smry[[2]]
sem_flood$sem_graph[[2]] + 
  ggtitle(paste(
    'Flood Events (p = ', sem_flood$sem_smry[[2]]$test$standard$pvalue%>%round(3), ")",
    sep='')) + 
  sem_drought$sem_graph[[2]] + ggtitle(paste(
    'Drought Events (p = ', sem_drought$sem_smry[[2]]$test$standard$pvalue%>%round(3), ")",
    sep=''))

# density
sem_flood$sem_smry[[3]]%>%str()
sem_drought$sem_smry[[3]]
sem_flood$sem_graph[[3]] + 
  ggtitle(paste(
    'Flood Events (p = ', sem_flood$sem_smry[[3]]$test$standard$pvalue%>%round(3), ")",
    sep='')) + 
  sem_drought$sem_graph[[3]] + ggtitle(paste(
    'Drought Events (p = ', sem_drought$sem_smry[[3]]$test$standard$pvalue%>%round(3), ")",
    sep=''))
# richness
sem_flood$sem_smry[[4]]
sem_drought$sem_smry[[4]]
sem_flood$sem_graph[[4]] + 
  ggtitle(paste(
    'Flood Events (p = ', sem_flood$sem_smry[[4]]$test$standard$pvalue%>%round(3), ")",
    sep='')) + 
  sem_drought$sem_graph[[4]] + ggtitle(paste(
    'Drought Events (p = ', sem_drought$sem_smry[[4]]$test$standard$pvalue%>%round(3), ")",
    sep=''))

# shannon
sem_flood$sem_smry[[5]]
sem_drought$sem_smry[[5]]
sem_flood$sem_graph[[5]] + 
  ggtitle(paste(
    'Flood Events (p = ', sem_flood$sem_smry[[5]]$test$standard$pvalue%>%round(3), ")",
    sep='')) + 
  sem_drought$sem_graph[[5]] + ggtitle(paste(
    'Drought Events (p = ', sem_drought$sem_smry[[5]]$test$standard$pvalue%>%round(3), ")",
    sep=''))

# -----------------------------------------------------------------------------
# Unified Model attempt #1
# -----------------------------------------------------------------------------
mod_uni <- '
Flood.Stress =~ ln_wk4_max + flood_days + sqrt_rain + annualrain
Drought.Stress =~ ln_wk4_min + low_flow_days + sqrt_rain + annualrain
LRRcentroid ~  Flood.Stress + Drought.Stress + annualrain 
'

make_fit_uni <-function(x_formula) {sem(x_formula, data = dc)}
make_sem_graph <- function(x_fit) {graph_sem(model=x_fit, angle=130)}

sem_uni <- tibble(outcome=my_outcomes) %>%
  mutate(formula = paste(
    "\nFlood.Stress =~ ln_wk4_max + flood_days + sqrt_rain + annualrain\nDrought.Stress =~ ln_wk4_min + low_flow_days + sqrt_rain + annualrain\n",
    outcome, 
    " ~  Flood.Stress + Drought.Stress + annualrain\n", sep=''),
    sem_fit = map(formula, make_fit_uni),
    sem_smry = map(sem_fit, make_smry),
    sem_graph = map(sem_fit, make_sem_graph),
    dataset = 'flood') %>%
  select(dataset, outcome, contains('sem'), formula)

sem_uni$sem_smry[[1]]
sem_uni$sem_graph[[1]] + 
  ggtitle(paste(
    'All Data (p = ', sem_uni$sem_smry[[1]]$test$standard$pvalue%>%round(3), ")",
    sep=''))

sem_uni$sem_smry[[2]]
sem_uni$sem_graph[[2]] + 
  ggtitle(paste(
    'All Data (p = ', sem_uni$sem_smry[[2]]$test$standard$pvalue%>%round(3), ")",
    sep=''))


# -----------------------------------------------------------------------------
# Generate Dictionary
# -----------------------------------------------------------------------------
d_tidy <- dc_fd
colnames(d_tidy) = str_replace_all(colnames(d_tidy), 'LRR_', '')
d_tidy <- d_tidy %>%
  rename(rain = sqrt_rain,
         flood_magnitude= ln_wk4_max,
         drought_magnitude = ln_wk4_min,
         drought_duration = low_flow_days,
         flood_duration = flood_days,
         centroid = centroid_dist) %>%
  select(-site_code, -collection_period)


model <- tidy_sem(d_tidy%>%select(-any_of(c(
  'centroid', 'density', 'richness', 'shannon'))))
model
dictionary(model)

# -----------------------------------------------------------------------------
# Generate Syntax
# -----------------------------------------------------------------------------
model %>%
  measurement() -> model
model
syntax(model)

# -----------------------------------------------------------------------------
# Add Paths
# -----------------------------------------------------------------------------
mod_uni2 <- '
flood_magnitude ~ annualrain + rain
flood_duration ~ annualrain + rain
drought_magnitude ~ annualrain + rain
drought_duration ~ annualrain + rain
biomass ~ flood_magnitude + flood_duration + drought_magnitude + drought_duration + annualrain
'

fit <- sem(mod_uni2, data = d_tidy)
fit %>% summary(standardized = TRUE)
fit %>% graph_sem(angle=170)


fit_uni <- model %>%
  add_paths(mod_uni2) %>%
  estimate_lavaan()

fit_uni %>% summary()
fit_uni %>% graph_sem(angle=130)


# -----------------------------------------------------------------------------
# Run the model
# -----------------------------------------------------------------------------
model %>%
  estimate_lavaan()

# -----------------------------------------------------------------------------
# Access the Dictionary, data and syntax
# -----------------------------------------------------------------------------
dictionary(model)
syntax(model)

# -----------------------------------------------------------------------------
# Modify the dictionary and syntax
# -----------------------------------------------------------------------------
dictionary(model) %>%
  mutate(label = ifelse(label == "vis", "Visual", label))

syntax(model) %>%
  mutate(lhs = ifelse(lhs == "spe" & op == "=~", "tex", lhs)) %>%
  filter(!(lhs == "spe" | rhs == "spe")) -> syntax(model)

estimate_lavaan(model)

# -----------------------------------------------------------------------------
# Adding Paths
# -----------------------------------------------------------------------------
model %>%
  add_paths("vis ~ tex") %>%
  estimate_lavaan() %>%
  summary(estimates = TRUE)

model %>%
  add_paths("vis ~ tex", vis =~ spe_1) %>%
  estimate_lavaan()