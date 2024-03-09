# 35_sem_v2
# Sean Kinard
# 2023-08-10
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')
library(tidySEM)
library(lavaan)
theme_set(theme_bw(base_size = 20))

# Data
d_outcome <- read_csv('exploration/output/model_outcome_vars.csv') %>%
  select(-qtr)

d_flow <- read_csv('data/06_fill_data/site_flow_4week_flood_duration.csv')
d_env <- read_csv('exploration/output/model_predictor_vars.csv')

my_index <- c('site_code', 'collection_period')
# -----------------------------------------------------------------------------
# Data Prep
# -----------------------------------------------------------------------------
# flood and drought metrics
d_predictor <- d_flow %>%
  mutate(RR_f.mag = (wk4_max-q_med_yr)/q_med_yr,
         RR_d.mag = (wk4_min-q_med_yr)/q_med_yr,
         f.dur = flood_duration,
         d.dur = lf_drought_duration) %>%
  select(site_code, collection_period, 
         RR_f.mag, RR_d.mag, f.dur, d.dur, contains('is_')) %>%
  add_rain()

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
  pivot_wider(names_from=outcome, values_from=LRR, names_prefix = 'LRR_') %>% 
  ungroup()

# merge
dc <- left_join(d_predictor, d_LRR) %>% unique() 

# filter to flood and drought only
d_fl <- filter(dc, is_maj_flood == 'yes') %>% unique() 
d_dr <- filter(dc, is_drought == 'yes') %>% unique() 
d_fldr <- filter(dc, is_drought=='yes' | is_maj_flood == 'yes')  %>% unique() 

n_flood <- d_fl %>%
  fix_site_order() %>%
  group_by(site_code) %>%
  dplyr::summarize(n_flood=length(site_code))

n_drought <- d_dr %>%
  fix_site_order() %>%
  group_by(site_code) %>%
  dplyr::summarize(n_drought=length(site_code))

n_flood_drought <- left_join(n_flood, n_drought) %>% 
  add_rain() %>%
  rename(Site = site_code, Rain.30yr=annualrain, 
         Floods=n_flood, Droughts=n_drought)

raw <- left_join(d_predictor, d_LRR) %>% left_join(d_env) %>% unique() %>%
  rename(F.Magnitude = RR_f.mag,
         D.Magnitude = RR_d.mag,
         F.Duration = f.dur,
         D.Duration = d.dur,
         Rain.Month = sqrt_rain,
         Rain.30yr = annualrain,
         Sp.Richness = LRR_richness,
         SW.Index = LRR_shannon,
         Abundance = LRR_density,
         Biomass = LRR_biomass,
         Centroid.Dist = LRR_centroid_dist)
my_outcomes <- c('Sp.Richness', 'SW.Index', 'Abundance', 'Biomass', 
                 'Centroid.Dist')

raw_fl <- raw %>% filter(is_maj_flood == 'yes')
raw_dr <- raw %>% filter(is_drought == 'yes')

# -----------------------------------------------------------------------------
# SEM functions
# -----------------------------------------------------------------------------
make_fit_raw_fl <-function(x_formula) {sem(x_formula, data = raw_fl)}

make_fit_raw_dr <-function(x_formula) {sem(x_formula, data = raw_dr)}

make_smry <- function(x_fit) {summary(x_fit, standardized = TRUE)}

make_pest <- function(x_fit) {
  x_fit %>%
    parameterestimates(se=F, standardized=T) %>%
    as_tibble() %>%
    filter(op != "~~") %>%
    rename(Response = lhs, Predictor = rhs, Estimate=est, SE=std.all) %>%
    select(Response, Predictor, Estimate, SE)  }

make_performace <- function(x_fit) {x_fit %>%
    model_performance(c("Chi2", "RMSEA", "CFI", "SRMR")) %>% as_tibble() }

make_sem_graph <- function(x_fit) {
  graph_data <- prepare_graph(x_fit)
  edges(graph_data) <- edges(graph_data) %>%
    mutate( label = est_std) %>%
    filter( ! arrow %in% c('both', 'none')) %>%
    mutate(colour = ifelse(est_std <=0, "red", "white")) %>%
    mutate(size = abs(as.numeric(est_std))*4)
  nodes(graph_data) <- nodes(graph_data) %>%
    mutate(label_family = 'sans',
           label_size = 4,
           label_color = case_when(
             name %in% c('D.Magnitude', 'D.Duration') ~ 'yellow',
             TRUE ~ 'cyan'))
  graph_data %>% plot() }

# -----------------------------------------------------------------------------
# SEM fits
# -----------------------------------------------------------------------------
sem_flood <- tibble(outcome=my_outcomes) %>%
  mutate(formula = paste(
    "F.Magnitude ~ Rain.Month
    F.Duration ~ Rain.Month + F.Magnitude
    ", outcome,  " ~ F.Magnitude + F.Duration + Rain.30yr", 
    sep=''),
    sem_fit = map(formula, make_fit_raw_fl),
    sem_smry = map(sem_fit, make_smry),
    sem_pest = map(sem_fit, make_pest),
    sem_perf = map(sem_fit, make_performace),
    sem_graph = map(sem_fit, make_sem_graph),
    dataset = 'flood') %>%
  select(dataset, outcome, contains('sem'), formula)

sem_drought <- tibble(outcome=my_outcomes) %>%
  mutate(formula = paste(
    "D.Magnitude ~ Rain.Month
    D.Duration ~ Rain.Month + D.Magnitude
    ", outcome, " ~ D.Magnitude + D.Duration + Rain.30yr", 
    sep=''),
    sem_fit = map(formula, make_fit_raw_dr),
    sem_smry = map(sem_fit, make_smry),
    sem_pest = map(sem_fit, make_pest),
    sem_perf = map(sem_fit, make_performace),
    sem_graph = map(sem_fit, make_sem_graph),
    dataset = 'drought') %>%
  select(dataset, outcome, contains('sem'), formula)

# -----------------------------------------------------------------------------
# SEM Plots
# -----------------------------------------------------------------------------
pretty_sem <- function(x_sem, x_outcome) {
  x_filtered <- x_sem %>% filter(outcome == x_outcome)
  
  chi_label <- paste("chi^2~'='~", x_filtered$sem_perf[[1]]$Chi2%>%round(2), sep='')
  rmsea_label <- paste("RMSEA = ", x_filtered$sem_perf[[1]]$RMSEA%>%round(2), sep='')
  cfi_label <- paste("CFI = ", x_filtered$sem_perf[[1]]$CFI%>%round(2), sep='')
  srmr_label <- paste("RMSEA = ", x_filtered$sem_perf[[1]]$SRMR%>%round(2), sep='')
 
  if(x_filtered$sem_perf[[1]]$Chi2 < 0.05 |
     x_filtered$sem_perf[[1]]$RMSEA > 0.08 |
     x_filtered$sem_perf[[1]]$CFI < 0.9 ) {
    
    x_filtered$sem_graph[[1]] + 
      annotate("text", x = 5.5, y = 2.7, size = 4, hjust = 0, color = 'red2', parse=T,
               label = chi_label) +
      annotate("text", x = 5.5, y = 2, size = 4, hjust = 0, color = 'red2',
               label = paste(rmsea_label,cfi_label,
                             sep='\n'))
  } else {
  x_filtered$sem_graph[[1]] + 
    annotate("text", x = 5.5, y = 2.7, size = 4, hjust = 0, parse=T,
             label = chi_label) +
    annotate("text", x = 5.5, y = 2, size = 4, hjust = 0,
             label = paste(rmsea_label,cfi_label,
                           sep='\n')) } }

sem_post <- function(x_sem, x_outcome) {
  x_filtered <- x_sem %>% filter(outcome == x_outcome)
  x_filtered$sem_smry[[1]] }

# individual plots
sem_flood %>% sem_post('Abundance')
plot_sem_f_abu <- sem_flood %>% pretty_sem('Abundance')

sem_flood %>% sem_post('Biomass')
plot_sem_f_bio <- sem_flood %>% pretty_sem('Biomass')

sem_flood %>% sem_post('Sp.Richness')
plot_sem_f_ric <- sem_flood %>% pretty_sem('Sp.Richness')

sem_flood %>% sem_post('SW.Index')
plot_sem_f_sha <- sem_flood %>% pretty_sem('SW.Index')

sem_flood %>% sem_post('Centroid.Dist')
plot_sem_f_cen <- sem_flood %>% pretty_sem('Centroid.Dist')

sem_drought %>% sem_post('Abundance')
plot_sem_d_abu <- sem_drought %>% pretty_sem('Abundance')

sem_drought %>% sem_post('Biomass')
plot_sem_d_bio <- sem_drought %>% pretty_sem('Biomass')

sem_drought %>% sem_post('Sp.Richness')
plot_sem_d_ric <- sem_drought %>% pretty_sem('Sp.Richness')

sem_drought %>% sem_post('SW.Index')
plot_sem_d_sha <- sem_drought %>% pretty_sem('SW.Index')

sem_drought %>% sem_post('Centroid.Dist')
plot_sem_d_cen <- sem_drought %>% pretty_sem('Centroid.Dist')

# multi-panel plots
plot_multi_floods <- (plot_sem_f_abu + plot_sem_f_bio) /
  (plot_sem_f_sha + plot_sem_f_cen)

plot_multi_droughts <- (plot_sem_d_abu + plot_sem_d_bio) /
  (plot_sem_d_sha + plot_sem_d_cen)

plot_multi_goodfit <- 
  (plot_sem_f_bio + plot_sem_f_sha) /
  (plot_sem_d_cen + plot_sem_d_sha) +
  plot_annotation(tag_levels = 'A')

plot_multi_badfit <- 
  (plot_sem_f_abu + plot_sem_f_cen) /
  (plot_sem_d_abu + plot_sem_d_bio) +
  plot_annotation(tag_levels = 'A')

# -----------------------------------------------------------------------------
# SEM Table
# -----------------------------------------------------------------------------
# parameter estimate table
add_model_fit <- function(df) {
  df %>%
    select(outcome, sem_perf) %>%
    unnest(cols = sem_perf) %>%
    mutate(model_fit = ifelse(
      Chi2 > 0.05 & RMSEA < 0.08 & CFI > 0.9, 'Good', 'Poor'  )) %>%
    right_join(df)  }

extract_pest <- function(df) {
  df %>%
    add_model_fit() %>%
    select(model_fit, sem_pest) %>% 
    unnest(cols = c(sem_pest)) %>%
    mutate(Estimate = Estimate %>% round(10),
           SE = SE %>% round(5)) %>%
    unique() %>%
    mutate(Response = str_replace_all(Response, 'F.Mag', 'Mag'),
           Response = str_replace_all(Response, 'F.Dur', 'Dur'),
           Response = str_replace_all(Response, 'D.Mag', 'Mag'),
           Response = str_replace_all(Response, 'D.Dur', 'Dur'),
           Predictor = str_replace_all(Predictor, 'F.Mag', 'Mag'),
           Predictor = str_replace_all(Predictor, 'F.Dur', 'Dur'),
           Predictor = str_replace_all(Predictor, 'D.Mag', 'Mag'),
           Predictor = str_replace_all(Predictor, 'D.Dur', 'Dur'))}

table_sem_pest <- sem_flood %>%
  extract_pest() %>%
  rename(Estimate.F=Estimate, SE.F = SE) %>%
  left_join(
    sem_drought %>%
      extract_pest() %>%
      rename(Estimate.D=Estimate, SE.D = SE)) %>%
  filter(model_fit == 'Good') %>%
  unique() %>%
  select(-model_fit)

#-------------------------------------------------------------------------------
# End 35_sem_v2
