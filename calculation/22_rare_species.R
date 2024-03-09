# rare_species
# Sean Kinard
# 2023-06-28
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

d_density <- read_csv("data/06_fill_data/fish_density_fill.csv")

d_bio_categories <- read_csv("data/02_clean_data/fish_biological_scales_of_comparison.csv")

d_lte <- read_csv("exploration/output/environment_longterm_with_ln.csv")
d_ste <- read_csv("exploration/output/environment_shortterm_with_ln.csv")
d_hur <- read_csv("data/02_clean_data/hurricane_stats.csv")

#------------------------------------------------------------------------------
# id species unique to certain time periods
#------------------------------------------------------------------------------
rare_occurance <- d_density %>%
  mutate(year_quarter = yearquarter(collection_period)) %>%
  group_by(year_quarter, lowest_taxon) %>%
  dplyr::summarize(presence = ifelse(sum(density)>0, 1, 0)) %>%
  ungroup() %>%
  group_by(lowest_taxon) %>%
  dplyr::summarize(occurance = sum(presence)) %>%
  filter(occurance<5) %>%
  pull(lowest_taxon)

rare_occurance_table <- d_density %>%
  filter(lowest_taxon %in% rare_occurance) %>%
  mutate(year_quarter = yearquarter(collection_period)) %>%
  group_by(year_quarter, lowest_taxon) %>%
  dplyr::summarize(presence = ifelse(sum(density)>0, 1, 0)) %>%
  pivot_wider(names_from=year_quarter, values_from=presence)

#------------------------------------------------------------------------------
# Probability-based occurance
#------------------------------------------------------------------------------
# presence or absence
d_present <-  d_density %>%
  group_by(site_code, collection_period, lowest_taxon) %>%
  dplyr::summarize(present = ifelse(density>0, 1, NA)) %>%
  pivot_wider(names_from=lowest_taxon, 
              values_from=present,
              values_fill=0) %>%
  pivot_longer(cols=contains(' '), names_to='lowest_taxon', values_to='presence') 

# Probability of Occurance 
prob_occur <- d_present %>%
  group_by(lowest_taxon, site_code) %>%
  dplyr::summarize(prob_present = sum(presence / length(presence))) %>%
  mutate(prob_present = ifelse(is.na(prob_present), 0 , prob_present))

species_occur_prob_vs_site <- prob_occur %>%
  add_rain() %>%
  mutate(annualrain = round(annualrain,0)%>% as.factor()) %>%
  fix_site_order() %>%
  filter(prob_present > 0) %>%
  ggplot(aes(prob_present, fill=annualrain, color = annualrain)) +
  facet_wrap(~site_code, ncol=3) +
  geom_density(alpha=.1) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  dark_theme_grey(base_size=14) +
  ylab('Probability Density') +
  xlab('Species Occurance Probability')

d_rare_tot <- left_join(d_present, prob_occur) %>%
  filter(prob_present <.25 & prob_present > 0) %>%
  left_join(d_density) %>%
  group_by(site_code, collection_period) %>%
  dplyr::summarize(rare_tot = sum(density, na.rm=T)) %>%
  ungroup()

d_rare_tot_scaled <- d_rare_tot %>%
  pivot_wider(values_from=rare_tot, names_from = site_code) %>%
  mutate(AR = scale(AR),
         EM = scale(EM),
         GC = scale(GC),
         MR = scale(MR),
         PD = scale(PD),
         PL = scale(PL),
         SF = scale(SF),
         TR = scale(TR),
         WM = scale(WM)) %>%
  pivot_longer(cols=AR:WM, names_to='site_code', values_to='rare_tot_scaled') %>%
  na.omit() 

rare_tot_vs_site <- d_rare_tot_scaled %>%
  add_rain() %>%
  create_site_group() %>%
  mutate(Rain = round(annualrain,0)%>% as.factor()) %>%
  fix_site_order() %>%
  ggplot(aes(collection_period, rare_tot_scaled, fill=Rain, color = Rain)) +
  facet_wrap(~site_code, ncol=3) +
  geom_smooth(method='loess', se=F, span=.8) +
  geom_point(alpha=.1, size=3, shape =21) +
  geom_point(fill=NA, size=3, shape = 21) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  dark_theme_grey(base_size=14) +
  ylab('Total Rare Species (scaled)') +
  xlab('Time')

rare_tot_scaled_vs_time <- d_rare_tot_scaled %>%
  left_join(d_hur) %>%
  mutate(hu_max = round(hu_max,0)%>% as.factor()) %>%
  fix_site_order() %>%
  ggplot(aes(collection_period, rare_tot_scaled, fill=site_code, 
             color = site_code)) +
  facet_wrap(~hu_max, ncol=3) +
  geom_smooth(method='loess', se=F, span=.4) +
  geom_point(alpha=.1, size=3, shape =21) +
  geom_point(fill=NA, size=3, shape = 21) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  dark_theme_grey(base_size=14) +
  ylab('Total Rare Species (scaled)') +
  xlab('Time') +
  geom_vline(xintercept = ymd("2017-08-27"), linewidth=.3, color = 'white') +
  geom_vline(xintercept = ymd(c("2018-08-27", "2019-08-27", "2020-08-27")),
             linewidth=.1)


#------------------------------------------------------------------------------
# present in recovery vs baseline
#------------------------------------------------------------------------------

recovery_start <- ymd("2017-08-20")
recovery_end <- ymd("2018-10-01")

# present @ recovery
present_at_recovery <- d_density %>%
  filter(collection_period > recovery_start & 
           collection_period < recovery_end) %>%
  mutate(year_quarter = yearquarter(collection_period)) %>%
  group_by(year_quarter, lowest_taxon) %>%
  dplyr::summarize(present_at_recovery = ifelse(sum(density)>0, 1, 0)) %>%
  ungroup() %>%
  select(lowest_taxon, present_at_recovery) %>%
  unique()

# present @ baseline
present_at_baseline <- d_density %>%
  filter(collection_period < recovery_start | 
           collection_period > recovery_end) %>%
  mutate(year_quarter = yearquarter(collection_period)) %>%
  group_by(year_quarter, lowest_taxon) %>%
  dplyr::summarize(present_at_baseline = ifelse(sum(density)>0, 1, 0)) %>%
  ungroup() %>%
  select(lowest_taxon, present_at_baseline) %>%
  unique()

only_present_during_recovery <- left_join(present_at_recovery,
                                          present_at_baseline) %>%
  filter(is.na(present_at_baseline))

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
my_tables <- list(prob_occur, rare_occurance_table, 
                   only_present_during_recovery)

names(my_tables) <- c('rare_spe_prob_occur', 'rare_spe_occurance_table', 
                      'rare_spe_only_recovery')

for (i in 1:length(my_tables)) {
  my_place <- paste('exploration/output/', names(my_tables[i]), ".csv", sep='')
  my_object <- my_tables[[i]]
  write_csv(my_object, my_place) }


my_figures <- list(
  species_occur_prob_vs_site, rare_tot_vs_site, rare_tot_scaled_vs_time)

names(my_figures) <- c('rare_species_occur_prob_vs_site', 'rare_tot_vs_site',
                       'rare_tot_scaled_vs_time')

for (i in 1:length(my_figures)) {
  my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
  my_object <- my_figures[[i]]
  ggsave(my_place,
         plot = my_object,
         width = 6,
         height = 9,
         units = c("in")) }

#------------------------------------------------------------------------------
# End rare_species

