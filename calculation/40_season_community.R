# season_community
# Sean Kinard
# 2023-08-28
#------------------------------------------------------------------------------
# setup
#-----------------------------------------------------------------------------
source('exploration/toolkit.R')

d_spp <- read_csv('data/02_clean_data/fish_species.csv') %>%
  select(lowest_taxon, family, genus, species) %>%
  mutate(
    lowest_taxon = lowest_taxon %>% 
      str_replace_all("cyanoguttatum", "cyanoguttatus"),
    species = species %>%
      str_replace_all("cyanoguttatum", "cyanoguttatus"))

d_den <- read_csv('data/06_fill_data/fish_density_fill.csv') %>% left_join(d_spp)

d <- read_csv('exploration/output/model_outcome_vars.csv') %>%
  add_rain() %>%
  fix_site_order() %>%
  add_season() %>%
  mutate(ln_biomass=log(biomass),
         ln_density=log(density),
         sqrt_centroid_distance = sqrt(centroid_dist))

d_scaled <- read_csv('exploration/output/31_bio_vars_SCALED.csv') %>%
  add_rain() %>%
  fix_site_order() %>%
  add_season() %>%
  select(-z_simpson)

colnames(d_scaled) <- str_replace_all(colnames(d_scaled), 'z_', '')

d_flow <- read_csv('data/06_fill_data/site_flow_4week_flood_duration.csv')

my_outcomes <- c('ln_density', 'ln_biomass', 'richness', 'shannon')
my_index <- c('site_code', 'collection_period')

#-----------------------------------------------------------------------------
# functions
#-----------------------------------------------------------------------------
make_lm_stats <- function(xdata) {
lm(value~annualrain, xdata) %>%
  glance() }

extract_p <- function(xdata) {
  lm(value~annualrain, xdata) %>%
    glance() %>%
    pull(p.value) }

visualize_season <- function(xdata) {
  basic_plot <- xdata %>%
    ggplot(aes(x=annualrain, y=value)) +
    geom_point(size=3, fill = 'grey70', alpha=.5, shape=21) +
    stat_poly_eq(label.x = "middle",
                 label.y = "top",
                 use_label(c("R2", "p"))) +
    labs(y=pretty_x(xdata$xvar[1]),
         x='Annual Rainfall (cm)') +
    theme_bw(base_size=14)
  
  if(extract_p(xdata) < 0.05) {
    basic_plot + 
      geom_smooth(method='lm', se=F, lwd=.5, color='red', lty=2)
      } else {
        basic_plot } }

visualize_multi_outcome_season <- function(xseason) {
  dx <- d_lm %>% filter(season==xseason)
  
  p1 <- dx %>% filter(outcome == 'ln_density')
  p2 <- dx %>% filter(outcome == 'ln_biomass')
  p3 <- dx %>% filter(outcome == 'richness')
  p4 <- dx %>% filter(outcome == 'shannon') 
  
  p1 <- p1$lm_plot[[1]]
  p2 <- p2$lm_plot[[1]]
  p3 <- p3$lm_plot[[1]]
  p4 <- p4$lm_plot[[1]]
  p1+p2+p3+p4+plot_layout(ncol=2)+
    plot_annotation(title = xseason) }

make_lm_stats_2 <- function(xdata) {
  lm(prcnt_com~annualrain, xdata) %>% 
    tidy() %>%
    filter(term=='annualrain') %>%
    cbind(lm(prcnt_com~annualrain, xdata) %>%glance()%>%select(r.squared)) %>%
    select(estimate, r.squared, p.value) }

#-----------------------------------------------------------------------------
# Season panels
#-----------------------------------------------------------------------------
d_lm <- d %>%
  select(site_code, annualrain, collection_period, season, 
         any_of(my_outcomes)) %>%
  pivot_longer(cols=any_of(my_outcomes), 
               names_to='outcome', 
               values_to = 'value') %>%
  mutate(xvar=outcome) %>%
  group_by(outcome, season) %>%
  nest() %>%
  mutate(lm_p = map(data, make_lm_stats),
         lm_plot = map(data, visualize_season)) %>%
  unnest(lm_p)

plot_winter_lm <- visualize_multi_outcome_season('Winter')
plot_spring_lm <- visualize_multi_outcome_season('Spring')
plot_summer_lm <- visualize_multi_outcome_season('Summer')
plot_fall_lm <- visualize_multi_outcome_season('Fall')

d_lm %>%
  mutate(is_sig=ifelse(p.value<.05, 'Yes', 'No'),
         outcome=paste(pretty_x(outcome), "vs AR", sep=' ')) %>%
  ggplot(aes(y=statistic, x=season, fill=is_sig)) +
  facet_wrap(~outcome, scale='free') +
  geom_point(shape=21, size=3) +
  scale_fill_manual(values=c('salmon', 'black')) +
  labs(fill = expression(paste(italic("p"), " < 0.05" )),
       y="Slope", 
       x='Season')

table_season_lm <- d_lm %>% select(season, outcome, statistic, r.squared, p.value, nobs)

#-----------------------------------------------------------------------------
# Trend vs Season
#-----------------------------------------------------------------------------
make_lm_stats_3 <- function(xdata) {
  lm(value~annualrain, xdata) %>%
    tidy() %>%
    filter(term == 'annualrain')%>%
    select(-term) }

d_lm2 <- d %>%
  select(site_code, annualrain, collection_period, season, 
         any_of(my_outcomes)) %>%
  pivot_longer(cols=any_of(my_outcomes), 
               names_to='outcome', 
               values_to = 'value') %>%
  mutate(xvar=outcome) %>%
  group_by(outcome, season) %>%
  nest() %>%
  mutate(lm_p = map(data, make_lm_stats_3)) %>%
  unnest(lm_p)

plot_trend_v_season <- d_lm2 %>%
  mutate(is_sig=ifelse(p.value<.05, 'Yes', 'No'),
         outcome=pretty_x(outcome)) %>%
  ggplot(aes(y=estimate, x=season, fill=is_sig)) +
  facet_wrap(~outcome, scale='free') +
  geom_point(shape=21, size=3) +
  scale_fill_manual(values=c('salmon', 'black')) +
  labs(fill = expression(paste(italic("p"), " < 0.05" )),
       y="Slope (per cm Annual Rainfall)", 
       x='Season')

#-----------------------------------------------------------------------------
# Seasonal changes within sites
#-----------------------------------------------------------------------------
d_scaled_long <- d_scaled %>%
  pivot_longer(cols=any_of(my_outcomes), 
               names_to='xname', values_to='xvalue')

d_scaled_long_summary <- d_scaled_long %>%
  group_by(site_code, annualrain, season,xname) %>%
  dplyr::summarize(z_mu = mean(xvalue),
                  z_sd = sd(xvalue)) %>% 
  ungroup() %>%
  fix_site_order()

site_group_colors <- c('white', '#A1CAF6FF', '#4C6FA1FF')

plot_season_box <- d_scaled_long %>%
  create_site_group() %>%
  mutate(xname=pretty_x(xname)) %>%
  ggplot(aes(x=season, y=xvalue, fill=site_group)) +
  facet_wrap(~xname, scales='free') +
  geom_hline(yintercept=0, lty=2, lwd=.5, color='red') +
  geom_boxplot(alpha=.5) +
  scale_fill_manual(values=site_group_colors) +
  labs(x=element_blank(), y = 'Z-Score') +
  theme_bw(base_size=14) +
  theme(legend.title=element_blank())

table_season_box <- d_scaled_long %>%
  create_site_group() %>%
  group_by(site_group, xname) %>%
  dplyr::summarize(x_mean=mean(xvalue),
                   x_sd=sd(xvalue),
                   x_n = length(xvalue))

#-----------------------------------------------------------------------------
# LRR Heat Map
#-----------------------------------------------------------------------------
# Log Response Ratio
d_out_long <- d%>%select(any_of(my_index), density, biomass, richness, shannon) %>%
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

d_LRR_summary <- d_LRR %>%
  mutate(Month=month(collection_period)) %>%
  pivot_longer(cols=contains('LRR'), names_to='xname', values_to = 'xvalue') %>%
  group_by(site_code, Month, xname) %>%
  dplyr::summarize(x_mu=mean(xvalue, na.rm=T)) %>%
  pivot_wider(names_from='xname', values_from='x_mu') %>%
  add_rain() %>%
  mutate(Rainfall=round(annualrain,0)%>%as.factor()) %>%
  unique()

# fill missing linear interpolation
site_x_date_cross <- tibble(site_code = list(my_sites),
                            Month = list(1:12)) %>%
  unnest(site_code) %>%
  unnest(Month)

d_LRR_summary_fill <- left_join(site_x_date_cross, d_LRR_summary) %>%
  select(-contains('ain')) %>%
  rowwise() %>%
  as_tsibble(key=site_code,
             index=Month) %>%  
  na_ma(weighting="linear") %>%
  ungroup() %>%
  as_tibble() %>%
  pivot_longer(cols=-any_of(c('site_code', 'Month')),
               names_to='xname',
               values_to = 'x') %>%
  group_by(site_code, Month, xname) %>%
  dplyr::summarize(x=mean(x ,na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from = x) %>%
  ungroup() %>%
  add_rain() %>%
  mutate(Rainfall=round(annualrain,0)%>%as.factor())

visualise_heatmap <- function(xoutcome) {
  my_title <- str_replace_all(xoutcome, 'LRR_','')%>%str_to_title()
  df <- d_LRR_summary_fill
  colnames(df) <- str_replace_all(colnames(df), xoutcome, 'LRR')
  
  df %>%
    ggplot(aes(x=Rainfall, y=Month, fill=LRR)) +
    geom_tile(color='white', na.rm = T) +
    scale_fill_viridis_c(option="inferno", direction=-1, limits = c(-4, 4)) +
    theme_bw(base_size=14) +
    scale_y_continuous(breaks=1:12) +
    labs(title=my_title, x=element_blank(), y=element_blank()) +
    coord_flip() }

plot_season_heatmap <- ((visualise_heatmap('LRR_density') + xlab('Rainfall (cm/yr)')) +
(visualise_heatmap('LRR_biomass') )) /
((visualise_heatmap('LRR_richness') + ylab('Month') + xlab('Rainfall (cm/yr)')) +
(visualise_heatmap('LRR_shannon') + ylab('Month'))) +
  plot_layout(guides='collect')

#-----------------------------------------------------------------------------
# Flow Disturbances vs Season
#-----------------------------------------------------------------------------
disturbance <- d_flow %>% 
  select(contains('is_'), collection_period, site_code) %>%
  mutate(is_maj_flood=ifelse(is_maj_flood=='yes', 1,0),
         is_drought = ifelse(is_drought=='yes', 1,0)) %>%
  add_season() %>%
  group_by(season, site_code) %>%
  dplyr::summarize(n_floods = sum(is_maj_flood, na.rm=T),
                   n_droughts = sum(is_drought, na.rm=T),
                   n_obs = length(collection_period)) %>%
  arrange(site_code, season) %>%
  ungroup() %>%
  mutate(percent_flood = n_floods/n_obs*100,
         percent_drought=n_droughts/n_obs*100)
  
d_heat_flow <- disturbance %>%
  add_rain() %>%
  mutate(Rainfall=round(annualrain,0)%>%as.factor()) %>%
  pivot_longer(cols=contains('percent'), names_to = 'xname', values_to = 'Percent') %>%
  mutate(xname=str_replace_all(xname, 'percent_', '')%>%str_to_title())

 plot_heat_flow <- d_heat_flow %>%
  ggplot(aes(x=season, y=Rainfall, fill=Percent)) +
  facet_wrap(~xname) +
  geom_tile(color='white') +
  scale_fill_viridis_c(option="inferno", direction=-1) +
  theme_bw(base_size=14) +
  labs(x=element_blank())

 #-----------------------------------------------------------------------------
 # Composition vs Season
 #-----------------------------------------------------------------------------
d_den_fam <- d_den %>%
  mutate(family=ifelse(family %in% c('Centrarchidae', 'Poeciliidae'), family, 'Other')) %>%
  group_by(site_code, collection_period, family) %>%
  dplyr::summarize(density_m2=sum(density, na.rm=T)) 

 extract_regression_stats <- function(xdata) {
   lm(prcnt_com~annualrain, data = xdata) %>%
     glance() %>%
     select(r.squared, p.value) }
 
d_den_fam_prep <- d_den_fam %>%
  group_by(site_code, collection_period) %>%
  dplyr::summarize(total_density_m2=sum(density_m2, na.rm=T)) %>%
  right_join(d_den_fam) %>%
  mutate(prcnt_com=density_m2/total_density_m2*100) %>%
  filter(family!='Other') %>%
  add_season() %>%
  add_rain() %>%
  group_by(season, family) %>%
  nest() %>%
  mutate(lm_stats=map(data, extract_regression_stats)) %>%
  unnest(lm_stats) %>%
  unnest(data)

plot_season_composition <- d_den_fam_prep %>%
  ggplot(aes(x=annualrain, y=prcnt_com, color = family)) +
  facet_grid(season~family) +
  geom_hline(yintercept=0, lty=1, lwd=.3, color='grey') +
  geom_smooth(data = d_den_fam_prep%>% filter(p.value < 0.05),
              method='lm', se=F, lwd=.6, lty=2, show.legend=F, color='red3') +
  geom_point(shape=21, size=1, fill='grey35',color='black') +
  stat_poly_eq(color='black', 
               label.x = "middle",
               label.y = "top",
               use_label(c("R2", "p"))) +
  theme_bw(base_size=14) +
  labs(x='Rainfall (cm/yr)', y='% Total Abundance')

table_season_composition <- d_den_fam_prep %>%
  group_by(season, family) %>%
  nest() %>%
  mutate(lm_stats=map(data, make_lm_stats_2)) %>%
  unnest(lm_stats) %>%
  select(-data) %>%
  rename(Family=family, Season=season, Slope=estimate, r_squared=r.squared, 
         p_value=p.value) %>%
  arrange(Family, Season)

table_season_composition %>%
  ggplot(aes(y=abs(Slope), x=Season, color=Family, size=1/p_value)) +
  geom_point()

plot_composition_v_season <- table_season_composition %>%
  mutate(is_sig=ifelse(p_value<.05, 'Yes', 'No')) %>%
  ggplot(aes(y=Slope, x=Season, fill=is_sig, shape=Family)) +
  geom_point(size=3) +
  scale_fill_manual(values=c('salmon', 'black')) +
  scale_shape_manual(values=c(21,22)) +
  labs(fill = expression(paste(italic("p"), " < 0.05" )),
       y="Slope (% of Community per cm Annual Rainfall)", 
       x='Season')

#------------------------------------------------------------------------------
# Are Centrarchid Densities lower in Summer across the region?
#------------------------------------------------------------------------------

site_z_score <- function(x) {
  x %>%
    mutate(z_score = scale(sum_cent)) }

prepare_z_score <- function(my_family) {
  d_den %>%
    filter(family == my_family) %>%
    group_by(site_code, collection_period) %>%
    dplyr::summarize(sum_cent = sum(density)%>%log()) %>%
    ungroup() %>%
    group_by(site_code) %>%
    nest() %>%
    mutate(z_values = map(data, site_z_score)) %>%
    unnest(z_values) %>%
    add_season() %>%
    add_rain() }

visualize_prcnt_com_v_season <- function(my_data, my_title) {
  ggboxplot(my_data, x = "season", y = "z_score", 
            add = "mean", 
            ylab = "Abundance (z score)",
            xlab = "Season",
            title = my_title) +
    theme_bw(base_size=14) + 
    stat_compare_means(comparisons = list(c("Winter", "Spring"), 
                                          c("Winter", "Summer"), 
                                          c("Winter", "Fall")))
}

# Cyprinodontidae
d_cyp <- prepare_z_score('Cyprinodontidae')
cyp_anova <- aov(z_score~season, data = d_cyp) 
cyp_anova_table <- cyp_anova%>% tidy()
cyp_thsd <- TukeyHSD(cyp_anova)
cyp_thsd_table <- cyp_thsd %>% tidy()
cyp_anova_figure <- visualize_prcnt_com_v_season(d_cyp, "Cyprinodontidae")

# Centrarchidae
d_cen <- prepare_z_score('Centrarchidae')
cen_anova <- aov(z_score~season, data = d_cen) 
cen_anova_table <- cen_anova%>% tidy()
cen_thsd <- TukeyHSD(cen_anova)
cen_thsd_table <- cen_thsd %>% tidy()
cen_anova_figure <- visualize_prcnt_com_v_season(d_cen, "Centrarchidae")

# Poeciliidae
d_poe <- prepare_z_score('Poeciliidae')
poe_anova <- aov(z_score~season, data = d_poe) 
poe_anova_table <- poe_anova%>% tidy()
poe_thsd <- TukeyHSD(poe_anova)
poe_thsd_table <- poe_thsd %>% tidy()
poe_anova_figure <- visualize_prcnt_com_v_season(d_poe, "Poeciliidae")

# Leuciscidae
d_leu <- prepare_z_score('Leuciscidae')
leu_anova <- aov(z_score~season, data = d_leu) 
leu_anova_table <- leu_anova%>% tidy()
leu_thsd <- TukeyHSD(leu_anova)
leu_thsd_table <- leu_thsd %>% tidy()
leu_anova_figure <- visualize_prcnt_com_v_season(d_leu, "Leuciscidae")

#------------------------------------------------------------------------------
# Composition time series
#------------------------------------------------------------------------------

make_com_prcnt <- function(x) {
  total_density <- sum(x$density)
  
  x %>%
    group_by(family) %>%
    mutate(com_prcnt = density/total_density*100) %>%
    select(family, com_prcnt)
}

d_prcnt <- d_den %>%
  group_by( collection_period) %>%
  nest() %>%
  mutate(my_output = map(data, make_com_prcnt)) %>%
  unnest(my_output) %>%
  add_season() %>% ungroup()

plot_comp_month <- d_prcnt %>%
  filter(family %in% c('Centrarchidae', 'Cyprinodontidae', 'Leuciscidae', 'Poeciliidae')) %>%
  mutate(xmonth=month(collection_period)) %>%
  group_by(xmonth, family) %>%
  dplyr::summarize(com_prcnt_mu=mean(com_prcnt)) %>%
  ggplot(aes(x = xmonth, y = com_prcnt_mu, color=family, fill = family)) +
  geom_smooth(se=F, lwd=.5, method='loess', span=.6) +
  scale_fill_manual(values = c('cyan', 'blue', 'goldenrod', 'red')) +
  scale_color_manual(values = c('cyan', 'blue', 'goldenrod', 'red')) +
  theme_bw(base_size = 12)

#------------------------------------------------------------------------------
# Export Table
#------------------------------------------------------------------------------
# my_objects <- ls()
# my_csv_names <- my_objects[str_detect(my_objects, 'table_')]
# 
# my_tables <- list(
#   table_season_box, table_season_composition, table_season_lm)
# 
# names(my_tables) <- my_csv_names
# 
# for (i in 1:length(my_tables)) {
#   my_place <- paste('exploration/output/40_', names(my_tables[i]), ".csv", sep='')
#   my_object <- my_tables[[i]]
#   write_csv(my_object, my_place) }

#------------------------------------------------------------------------------
# Export figures
#------------------------------------------------------------------------------
# my_figure_names <- my_objects[str_detect(my_objects, 'plot_')]
# 
# my_figures <- list(
#   plot_fall_lm, plot_heat_flow, plot_season_box,
#   plot_season_composition, plot_season_heatmap, plot_spring_lm, plot_summer_lm,
#   plot_winter_lm)
# 
# names(my_figures) <- my_figure_names
# 
# for (i in 1:length(my_figures)) {
#   my_place <- paste('exploration/visualization/40_', names(my_figures[i]), ".png", sep='')
#   my_object <- my_figures[[i]]
#   ggsave(my_place,
#          plot = my_object,
#          width = 10,
#          height = 10,
#          units = c("in")) }

#------------------------------------------------------------------------------
# End 40_season_community