# environment_split
# Sean Kinard
# 2023-06-20

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

d_flo <- read_csv("data/02_clean_data/site_flow.csv")
d_4wf <- read_csv("data/06_fill_data/site_flow_4week_flood_duration.csv")
d_rai <- read_csv("data/04_match_data/site_rain_match.csv")
d_wat <- read_csv("data/06_fill_data/site_watershed_fill.csv")
d_tra <- read_csv("data/04_match_data/site_transect_match.csv")

#------------------------------------------------------------------------------
# isolate long-term predictors
#------------------------------------------------------------------------------
# prep watershed metrics
d_wat_slim <- d_wat %>% select(site_code, lat:otherland) %>% unique()

# prep 30 year q metrics
q30 <- d_flo %>% select(site_code, contains('al_')) %>% unique()
colnames(q30) <- str_replace_all(colnames(q30), 'al_', 'q_')

# prep transect study average
mu_tra <- d_tra %>% pivot_longer(
  cols=-c('site_code', 'collection_period', 'collection_date'), 
  names_to = 'xvar', values_to = 'xval') %>%
  group_by(site_code, xvar) %>%
  dplyr::summarize(xmu=mean(xval, na.rm=T)) %>%
  pivot_wider(names_from=xvar, values_from = xmu) %>% unique()

# merge long-term
lterm <- left_join(d_wat_slim, q30) %>%
  left_join(mu_tra)

# remove variables with zero variance across sites
zero_variance <- lterm %>% 
  pivot_longer(cols=!site_code, names_to = 'xvar', 
               values_to='xval') %>%
  group_by(xvar) %>%
  dplyr::summarize(variance = var(xval, na.rm=T)) %>% 
  filter(near(variance,0)) %>%
  pull(xvar)

lterm <- lterm %>% select( - any_of(zero_variance), -q_n)

#------------------------------------------------------------------------------
# isolate short-term predictors
#------------------------------------------------------------------------------
# monthly transect
d_tra <- d_tra %>% select(-collection_date)

# prep rain metrics
rain_slim <- d_rai %>% select(site_code, collection_period, rain_mL_m2) %>%
  rename(rain = rain_mL_m2)

# 4-week prior flow stats
qm <- d_4wf %>%
  select(-site_period, -contains('is_'), -contains('date'), - lf_min_return) %>% 
  unique() %>%
  mutate(flood_duration = ifelse(flood_duration<0, 0, flood_duration))

# merge short-term
sterm <- left_join(d_tra, rain_slim) %>%
  left_join(qm) 

# remove duplicates by averaging site x collection_period x variable
sterm <- sterm %>% pivot_longer(cols = -c('site_code', 'collection_period'),
                       names_to = 'xvar', values_to = 'xval') %>%
  group_by(site_code, collection_period, xvar) %>%
  dplyr::summarize(xval = mean(xval, na.rm=T) ) %>%
  pivot_wider(names_from = 'xvar', values_from = 'xval')

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(lterm, 'exploration/output/environment_longterm.csv')
write_csv(sterm, 'exploration/output/environment_shortterm.csv')

#------------------------------------------------------------------------------
# End environment_split