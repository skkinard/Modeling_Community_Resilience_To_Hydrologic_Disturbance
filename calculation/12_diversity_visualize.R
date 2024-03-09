# diversity_visualize
# Sean Kinard
# 2023-06-23
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

d_bio <- read_csv('exploration/output/fish_diversity_biomass.csv')
d_den <- read_csv('exploration/output/fish_diversity_density.csv')
d_lte <- read_csv('exploration/output/environment_longterm_with_ln.csv')
d_ste <- read_csv('exploration/output/environment_shortterm_with_ln.csv')

#------------------------------------------------------------------------------
# Prep data
#------------------------------------------------------------------------------
# combine data
prep_data <- function(df, xname) {
  df %>% pivot_longer(
    cols=c('richness', 'shannon', 'simpson'),
    names_to='diversity', values_to= 'x') %>%
    mutate(type= xname) }

d_both <- full_join(prep_data(d_bio, 'biomass'),
                    prep_data(d_den, 'density')) %>%
  create_date_vars() %>% # add time variables
  left_join(d_lte %>% select(site_code, any_of(climate))) %>% # add climate
  left_join(d_ste) %>% # add environmental variables
  fix_site_order()

#------------------------------------------------------------------------------
# biomass_diversity vs density_diversity
#------------------------------------------------------------------------------
density_vs_biomass <- d_both %>%
  pivot_wider(names_from = type, values_from = x) %>%
  ggplot(aes(biomass, density, color = site_code)) +
  facet_wrap(~diversity, scales='free', nrow=3) +
  geom_abline(lty=2) +
  geom_point() +
  scale_color_manual(values = my_colors) +
  dark_theme_grey()

density_vs_biomass2 <- d_both %>%
  ggplot(aes(x = site_code, y = x, fill = type)) +
  facet_wrap(~year+diversity, ncol=3) +
  geom_boxplot(position=position_dodge(), color = 'grey50') +
  scale_fill_manual(values = c('red2', 'yellow3')) +
  dark_theme_grey(base_size=10) +
  ylab('Diversity') +
  labs(fill = 'Rain (cm/yr)') +
  xlab(element_blank())

# biomass appears underdispersed and prone to underestimation
d_div <- d_both %>% filter(type == 'density')

#------------------------------------------------------------------------------
# diversity vs site
#------------------------------------------------------------------------------
diversity_vs_site <- d_div %>%
  group_by(site_code, year, qtr, diversity) %>% # to even sample distribution
  dplyr::summarize(x=mean(x, na.rm=T)) %>%
  ggplot(aes(x = site_code, y = x, fill = site_code)) +
  facet_wrap(~diversity, scales='free') +
  geom_boxplot(position=position_dodge(), color = 'grey50') +
  scale_fill_manual(values = my_colors) +
  dark_theme_grey(base_size=18) +
  ylab('Diversity') +
  labs(fill = 'Rain (cm/yr)') +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position='bottom')

diversity_vs_site_year <- d_div %>%
  ggplot(aes(x = site_code, y = x, fill = site_code)) +
  facet_grid(year~diversity, scales='free') +
  geom_boxplot(position=position_dodge(), color = 'grey50') +
  scale_fill_manual(values = my_colors) +
  dark_theme_grey(base_size=18) +
  ylab('Diversity') +
  labs(fill = 'Rain (cm/yr)') +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position='bottom')

#------------------------------------------------------------------------------
# diversity vs time
#------------------------------------------------------------------------------
base_div <- function(my_data, xspan) {
  my_data %>%
    mutate(year_month = as_date(year_month)) %>%
    fix_site_order() %>%
    ggplot(aes(x = year_month, y = x, color = site_code, fill = site_code)) +
    scale_color_manual('Site', values = my_colors) +
    scale_fill_manual('Site', values=my_colors) +
    dark_theme_grey(base_size=14) +
    scale_x_date(date_breaks='3 months', date_labels = "%y-%m") +
    ylab(element_blank())+
    xlab('Year Month') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_vline(xintercept = ymd("2017-08-27"), linewidth=.3, color = 'white') +
    geom_vline(xintercept = ymd(c("2018-08-27", "2019-08-27", "2020-08-27")),
               linewidth=.1) +
    geom_jitter(shape = 21, size = 2, alpha = .4, show.legend=T,
                position=position_dodge(width=45)) +
    geom_jitter(shape = 21, size = 2, fill=NA, show.legend=T,
                position=position_dodge(width=45)) +
    geom_smooth(aes(group=NA), method = "loess", se=FALSE, color = 'purple',
                linewidth=1, span = xspan, show.legend=F) }

div_vs_time <- d_div %>%
  filter(diversity != 'simpson') %>%
  base_div(0.3) +
  facet_wrap(~pretty_x(diversity), scales = 'free_y', ncol=1)
  
div_rich_vs_time_site <- d_div %>%
  create_site_group() %>%
  filter(diversity == 'richness') %>%
  base_div(0.6) +
  facet_wrap(~site_group, scales = 'free_y', ncol=1)

div_sha_vs_time_site <- d_div %>%
  create_site_group() %>%
  filter(diversity == 'shannon') %>%
  base_div(0.6) +
  facet_wrap(~site_group, scales = 'free_y', ncol=1)

# MR 2019 (especially fall) is highly suspicious
# it looks like linear mixed effect model using site as a random effect will reveal a post-hurricane bump and steady rise in diversity over time at each site.

#------------------------------------------------------------------------------
# autocorrelation
#------------------------------------------------------------------------------
trash_time_vars <- c('collection_period', 'year', 'month', 'day', 
               'nday', 'qtr', 'year_week')

# fill missing predictors with linear interpolation
d_div_fill <- d_div %>%
  select(-any_of(trash_time_vars)) %>%
  arrange(year_month) %>%
  pivot_wider(names_from = diversity, values_from = x) %>%
  as_tsibble(key=site_code,
             index=year_month) %>%  
  na_ma(weighting="linear")

autocorrelate <- function(df) {
  df %>%
    select(where(is.numeric), -year_month) %>%
    correlate() %>%
    select(term, richness, shannon, simpson) %>%
    na.omit() }

div_cor <- d_div_fill %>%
  group_by(site_code) %>%
  nest() %>%
  mutate(data_cor = map(data, autocorrelate)) %>%
  unnest(data_cor) %>%
  select(-data)

div_correlation <- div_cor %>%
  select(-simpson) %>%
  pivot_longer(cols=c('richness', 'shannon'), 
               names_to='type', values_to='r') %>%
  rename(predictor=term) %>%
  ggplot(aes(x=r, y = predictor, color = site_code, 
             fill=site_code, alpha = abs(r), size= abs(r^2))) +
  facet_wrap(~type) +
  scale_color_manual('Site', values = my_colors) +
  scale_fill_manual('Site', values=my_colors) +
  dark_theme_grey(base_size=14) +
  geom_vline(xintercept = 0, color='grey90', lty=2, linewidth=.3) +
  geom_point(shape = 21, show.legend = F) 

# diversity at sites in the middle of the gradient are negatively related to flow variability
#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
# my_figures <- list(div_correlation,
#                    div_sha_vs_time_site,
#                    div_rich_vs_time_site,
#                    div_vs_time)
# 
# names(my_figures) <- c('div_correlation',
#                        'div_sha_vs_time_site',
#                        'div_rich_vs_time_site',
#                        'div_vs_time')
# 
# for (i in 1:length(my_figures)) {
#   my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
#   my_object <- my_figures[[i]]
#   ggsave(my_place,
#          plot = my_object,
#          width = 7,
#          height = 9,
#          units = c("in")) }

#------------------------------------------------------------------------------
# End diversity_visualize