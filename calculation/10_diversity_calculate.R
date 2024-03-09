# 10_diversity_calculate
# Sean Kinard
# 2023-06-23

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')
library(iNEXT)

d_biomass <- read_csv('data/06_fill_data/fish_biomass.csv')
d_density <- read_csv('data/04_match_data/fish_density.csv') %>%
  combo_create_period()

#------------------------------------------------------------------------------
# Tidy Data
#------------------------------------------------------------------------------
# transform variable_of_interest to pseudo-abundance
trash <- c('AFDMg', 'AFDMg_m2', 'abun', 'f_density', 'stream_area',
           'collection_date', 'site_period')

bio <- d_biomass %>% 
  mutate(pseudo_abun = ceiling(AFDMg_m2 *1000)) %>%
  select(-any_of(trash))

den <- d_density %>% 
  mutate(pseudo_abun = ceiling(f_density *1000)) %>%
  select(-any_of(trash))

# format: rows = species,  column = 'site_code collection_date'
prep_inext <- function(df) {
  df %>%
    unite('UID', site_code, collection_period, sep=' ') %>%
    pivot_wider(names_from=UID, values_from = pseudo_abun, values_fill = 0) %>%
    column_to_rownames(var = "lowest_taxon") %>%
    as.data.frame() }

bio_prep <- bio %>% prep_inext()

den_prep <- den %>% prep_inext()

#------------------------------------------------------------------------------
# Calculate Hill Numbers
#------------------------------------------------------------------------------
# supported data formats (abundance, incidence_raw, incidence_frequency)
# 'abundance' -> Individual‐based abundance data (datatype="abundance"): Input data for each assemblage/site include species abundances in an empirical sample of n individuals (“reference sample”). When there are N assemblages, input data consist of an S by N abundance matrix, or N lists of species abundances.

# Create iNEXt objects: biomass
bio_inext <- iNEXT(bio_prep)

# Extract Hill numbers from iNEXT objects: biomass
bio_hill <- bio_inext$AsyEst %>%
  separate(Assemblage, c('site_code', 'collection_period'), sep = ' ') %>%
  as_tibble()
bio_div <- bio_hill %>%
  select(site_code:Observed) %>%
  pivot_wider(names_from = Diversity, values_from = Observed) %>%
  dplyr::rename(richness = `Species richness`,
                shannon = `Shannon diversity`,
                simpson = `Simpson diversity`)

# Create iNEXt objects: density
den_inext <- iNEXT(den_prep)

# Extract Hill numbers from iNEXT objects: density
den_hill <- den_inext$AsyEst %>%
  separate(Assemblage, c('site_code', 'collection_period'), sep = ' ') %>%
  as_tibble()

den_div <- den_hill %>%
  select(site_code:Observed) %>%
  pivot_wider(names_from = Diversity, values_from = Observed) %>%
  dplyr::rename(richness = `Species richness`,
                shannon = `Shannon diversity`,
                simpson = `Simpson diversity`)

#------------------------------------------------------------------------------
# Impute missing spring 2017 richness and shannon using PeerJ publication
#------------------------------------------------------------------------------
sp17_d <- read_csv(
  "data/00_source_data/Kinard_data/sp17_fish_diversity_estimates.csv") %>%
  as_tibble() %>% r_friendly_colnames() %>% select(-contains('_'))
sp17_e <- read_csv("data/00_source_data/Kinard_data/sp17_site_x_env.csv") %>%
  as_tibble() %>% r_friendly_colnames() %>% select(staid, ap)
sp_17_c <- left_join(sp17_d, sp17_e) %>% rename(annualrain=ap,
                                                richness=rich) %>%
  mutate(type='Kinardetal2020')

den_sp17 <- den_div %>% 
  filter(collection_period < ymd("2017-06-01")) %>%
  add_rain() %>%
  select(-collection_period, -site_code) %>%
  mutate(type='density-based')

# fit model using 13 sites from spring surveys
lm_richness <- lm(richness~annualrain, data = sp_17_c)
lm_shannon <- lm(shannon~annualrain, data = sp_17_c)
lm_simpson <- lm(simpson~annualrain, data = sp_17_c)

# fit model using 7 sites using density estimates
lm2_richness <- lm(richness~annualrain, data = den_sp17)
lm2_shannon <- lm(shannon~annualrain, data = den_sp17)
lm2_simpson <- lm(simpson~annualrain, data = den_sp17)

# use averaged models to predict sites monitored after hurricane harvey (n=9)
div_sim <- tibble(site_code = my_sites) %>% add_rain() %>% 
  mutate(r_1 = predict(lm_richness, cur_data()),
         s_1 = predict(lm_shannon, cur_data()),
         i_1 = predict(lm_simpson, cur_data()),
         r_2 = predict(lm2_richness, cur_data()),
         s_2 = predict(lm2_shannon, cur_data()),
         i_2 = predict(lm2_simpson, cur_data()),
         r_3 = predict(lm2_richness, cur_data()),
         s_3 = predict(lm2_shannon, cur_data()),
         i_3 = predict(lm2_simpson, cur_data()),
         type = 'simulated') %>%
  rowwise() %>%
  mutate(richness = ceiling(mean(c_across(c('r_1', 'r_2', 'r_3',)))),
         shannon = mean(c_across(c('s_1', 's_2', 's_3'))),
         simpson = mean(c_across(c('i_1', 'i_2', 'i_3')))) %>%
  ungroup() %>%
  select(-contains('_1'), - contains('_2'), - contains('_3'))

# conservative simulated vs real (weighted 2:1 density-based:kinardetal2020)
div_simulated_sp17 <- den_sp17 %>%
  full_join(div_sim%>%select(-site_code)) %>%
  full_join(sp_17_c%>%select(-staid)) %>%
  pivot_longer(cols=c(richness, shannon, simpson), 
               names_to='diversity', values_to = 'x') %>%
  ggplot(aes(annualrain, x, color = type, fill=type)) +
  facet_wrap(~diversity, scales='free', ncol=1) +
  geom_smooth(method='lm', se=F,show.legend = F, lty=2, linewidth=.4) +
  geom_point(shape=21, size = 4, alpha=.1) +
  geom_point(shape=21, size = 4, fill=NA) +
  dark_theme_grey(base_size=14) +
  scale_color_manual(values=c( 'yellow', 'skyblue', 'white')) +
  scale_fill_manual(values=c('yellow', 'skyblue', 'white'))

# add simulated div to density based diversity estimates:
sp17_fill <- full_join(den_div %>% 
            mutate(collection_period = ymd(collection_period)) %>%
            filter(collection_period < ymd("2017-06-01")),
          div_sim %>%
            filter(annualrain>86) %>%
            mutate(collection_period = ymd('2017-04-01')) %>%
            select(-type, - annualrain)) 
  
# add fill data to density diversity
den_div <- den_div %>% 
  mutate(collection_period = ymd(collection_period)) %>%
  filter(collection_period > ymd("2017-06-01")) %>%
  full_join(sp17_fill) %>%
  arrange(collection_period,site_code)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(bio_div, 'exploration/output/fish_diversity_biomass.csv')
write_csv(den_div, 'exploration/output/fish_diversity_density.csv')

#------------------------------------------------------------------------------
# End 10_diversity_calculate
