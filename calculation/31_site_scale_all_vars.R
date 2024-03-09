# 31_site_scale_all_vars
# Sean Kinard
# 2023-07-09
#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

d <- read_csv('exploration/output/model_combined_vars.csv')

trash_vars <- c('year', 'qtr', 'month', 'day', 'nday', 'year_week', 
               'year_month', 'is_baseline', 'annualrain')
id_vars <- c('site_code', 'collection_period')

#------------------------------------------------------------------------------
# site_scale
#------------------------------------------------------------------------------
# inspect biovars for normality
d %>%
  select(any_of(id_vars), density, biomass, richness, shannon, simpson, centroid_dist) %>%
  pivot_longer(cols=-any_of(id_vars), names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  geom_histogram() +
  facet_wrap(~xvar, scales='free')

# ln transform biomass, centroid_dist, density
d <- d %>%
  mutate(ln_density = log(density),
         ln_biomass = log(biomass),
         sqrt_centroid_dist = sqrt(centroid_dist)) %>%
  select(-density, -biomass, -centroid_dist)

d %>%
  select(any_of(id_vars), ln_density, ln_biomass, sqrt_centroid_dist, 
         richness, shannon, simpson) %>%
  pivot_longer(cols=-any_of(id_vars), names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  geom_histogram() +
  facet_wrap(~xvar, scales='free')

d_scaled <- d %>%
  select(-any_of(trash_vars)) %>%
  pivot_longer(cols = -id_vars, names_to='xvar', values_to='xval') %>%
  group_by(site_code, xvar) %>%
  dplyr::mutate('sc_xval' = scale(xval, scale=FALSE)%>%as.vector())%>%
  ungroup() %>%
  select(-xval) %>%
  pivot_wider(names_from='xvar', values_from='sc_xval', names_prefix='z_')

d_scaled_bio <- d_scaled %>%
  select(site_code, collection_period, 
        z_richness, z_shannon, z_simpson, 
        z_ln_density, z_ln_biomass, z_sqrt_centroid_dist)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d_scaled, 'exploration/output/31_combined_vars_SCALED.csv')
write_csv(d_scaled_bio, 'exploration/output/31_bio_vars_SCALED.csv')

#------------------------------------------------------------------------------
# End 31_site_scale_all_vars