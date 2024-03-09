# 04_environment_long_term_preview
# Sean Kinard
# 2023-06-20

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

lterm <- read_csv('exploration/output/environment_longterm.csv')

# remove lat lon
lterm <- lterm %>% select(-any_of(location))

#------------------------------------------------------------------------------
# autocorrelation
#------------------------------------------------------------------------------
lte_correlation_table <- lterm %>% 
  correlate() %>%
  shave() %>% 
  fashion()
  
lte_correlation_plot <- lterm %>% 
  correlate() %>%
  shave() %>%
  rplot(colours = c("#F24405FF", "black", "#00FFFFFF")) +
  dark_theme_grey() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

#------------------------------------------------------------------------------
# histogram
#------------------------------------------------------------------------------
ex_hist <- function(my_data) {
  my_data %>% 
  pivot_longer(cols=-site_code, names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  geom_histogram(color = 'chartreuse3', fill = NA) +
  geom_histogram(fill = 'chartreuse3', color = NA, alpha=.2)+
  facet_wrap(~xvar, scales='free') +
  dark_theme_grey(base_size = 14) +
  xlab(element_blank()) +
  ylab(element_blank()) }

lte_histogram_landuse <- lterm %>% 
  select(site_code, any_of(landuse)) %>%
  ex_hist() # no log transforms necessary

lte_histogram_flow <- lterm %>% 
  select(site_code, any_of(flow)) %>%
  ex_hist() # need to log transform many q variables

lte_histogram_waterquality <- lterm %>% 
  select(site_code, any_of(water_quality)) %>%
  ex_hist() # need to log transform conductivity

lte_histogram_geomorph <- lterm %>% 
  select(site_code, any_of(geomorph)) %>%
  ex_hist()

#------------------------------------------------------------------------------
# natural log transformation
#------------------------------------------------------------------------------
ln_vars <- lterm %>%
  pivot_longer(cols=-site_code, names_to='xvar', values_to='xval') %>%
  group_by(xvar) %>%
  dplyr::summarize(range_over_mean = range(xval)/mean(xval)) %>%
  filter(range_over_mean > 4) %>%
  pull(xvar)

# visual verfication
plot_hist_lnvars_raw <- lterm %>% 
  select(site_code, any_of(ln_vars)) %>%
  pivot_longer(cols=-site_code, 
               names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  facet_wrap(~xvar, scales='free') +
  geom_histogram()

lterm_ln <- lterm %>% 
  select(site_code, any_of(ln_vars)) %>%
  pivot_longer(cols=-site_code, names_to='xvar', values_to='xval') %>%
  mutate(xval = log(xval)) %>%
  pivot_wider(names_from='xvar', values_from='xval')
colnames(lterm_ln) <- paste('ln_', colnames(lterm_ln), sep='')
lterm_ln <- lterm_ln %>% rename(site_code = ln_site_code)
lterm <- lterm %>% select(-any_of(ln_vars)) %>% left_join(lterm_ln)

ln_names <- paste('ln_', ln_vars, sep='')

lte_densiogram_ln_vars <- lterm %>% 
  select(site_code, any_of(ln_names)) %>%
  ex_hist()

sqrt_vars <- c('')

# restore ln_vars to categories
climate <-restore_category(climate)
landuse <- restore_category(landuse)
flow <- restore_category(flow)
water_quality <- restore_category(water_quality)
geomorph <- restore_category(geomorph)
algae <- restore_category(algae)

#------------------------------------------------------------------------------
# lte scatterplots
#------------------------------------------------------------------------------
scatter_plot <- function(my_data) {
  
  my_data %>%
    pivot_longer(cols = 2:length(colnames(my_data)),
                 names_to = "variable",
                 values_to = "value") %>%
    left_join(lterm %>% select(site_code, annualrain) %>% distinct()) %>%
    fix_site_order() %>%
    ggplot(aes(x=annualrain, y=value)) +
    facet_wrap(~variable, scale = 'free', ncol=2) +
    geom_smooth(method = "loess", se=FALSE, color="grey30", linetype=2,
                size=1.2, span=.8) +
    geom_point(aes(fill=site_code), color = 'black', size = 4, shape=21) +
    scale_fill_manual(values=my_colors) +
    labs(fill='Site') +
    xlab('Rain (cm/yr)') +
    ylab(element_blank()) +
    dark_theme_grey(base_size=14) }

lte_scatter_landuse <- lterm %>% 
  select(site_code, any_of(landuse)) %>%
  scatter_plot()

lte_scatter_flow <- lterm %>% 
  select(site_code, any_of(flow)) %>%
  scatter_plot() +
  facet_wrap(~variable, scale = 'free', ncol=3)

lte_scatter_waterquality <- lterm %>% 
  select(site_code, any_of(water_quality)) %>%
  scatter_plot() +
  facet_wrap(~variable, scale = 'free', ncol=3)

lte_scatter_geomorph <- lterm %>% 
  select(site_code, any_of(geomorph)) %>%
  scatter_plot()

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(lterm, 'exploration/output/environment_longterm_with_ln.csv')

# autocorrelation
write_csv(lte_correlation_table,
          'exploration/output/environment_longterm_autocorrelation.csv')
ggsave(filename = 'exploration/visualization/lte_correlation.png',
       lte_correlation_plot,
       width=9, height=9, units="in")

# scatter figures
ggsave(filename = 'exploration/visualization/lte_scatter_landuse.png',
       lte_scatter_landuse,
       width=7, height=7, units="in")
ggsave(filename = 'exploration/visualization/lte_scatter_flow.png',
       lte_scatter_flow,
       width=9, height=15, units="in")
ggsave(filename = 'exploration/visualization/lte_scatter_waterquality.png',
       lte_scatter_waterquality,
       width=10, height=9, units="in")
ggsave(filename = 'exploration/visualization/lte_scatter_geomorph.png',
       lte_scatter_geomorph,
       width=9, height=9, units="in")

#------------------------------------------------------------------------------
# End 04_environment_long_term_preview