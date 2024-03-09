# biomass_visualize
# Sean Kinard
# 2023-06-23
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

good_vars <- c('site_code', 'collection_period', 'lowest_taxon', 'biomass', 
               'density')

d_bio <- read_csv('data/06_fill_data/fish_biomass.csv') %>% 
  rename(biomass = AFDMg_m2) %>% 
  select(any_of(good_vars))

d_den <- read_csv('data/06_fill_data/fish_density_fill.csv') %>% 
  select(any_of(good_vars))

d_lte <- read_csv('exploration/output/environment_longterm_with_ln.csv')

d_ste <- read_csv('exploration/output/environment_shortterm_with_ln.csv')

#------------------------------------------------------------------------------
# Prep data
#------------------------------------------------------------------------------
# combine data
d_both <- full_join(d_bio, d_den) %>%
  pivot_longer(cols=c('biomass', 'density'), 
               names_to='abun_type', values_to='x') %>%
  create_date_vars() %>% # add time variables
  left_join(d_lte %>% select(site_code, any_of(climate))) %>% # add climate
  left_join(d_ste) %>% # add environmental variables
  fix_site_order()

#------------------------------------------------------------------------------
# biomass_abun_type vs density_abun_type
#------------------------------------------------------------------------------
density_vs_biomass <- full_join(d_bio%>%select(any_of(good_vars)),
                                d_den%>%select(any_of(good_vars))) %>% 
  fix_site_order() %>%
  group_by(site_code, collection_period) %>%
  dplyr::summarize(biomass = sum(biomass, na.rm=T),
                   density = sum(density, na.rm=T)) %>%
  ggplot(aes(biomass, density, color = site_code)) +
  geom_point() +
  scale_color_manual('Site', values = my_colors) +
  dark_theme_grey()

abn_density_vs_biomass2 <- d_both %>%
  ggplot(aes(x = site_code, y = x, fill = abun_type)) +
  geom_boxplot(position=position_dodge(), color = 'grey50') +
  scale_fill_manual('Type', values = c('red2', 'yellow3')) +
  dark_theme_grey(base_size=10) +
  ylab('abun_type') +
  xlab(element_blank()) +
  scale_y_log10()

# biomass varies more and is more conservative

#------------------------------------------------------------------------------
# abun_type vs site
#------------------------------------------------------------------------------
base_plot <- function(x) {
  x %>%
    ggplot() +
    scale_fill_manual(values = my_colors) +
    dark_theme_grey(base_size=18) +
    xlab(element_blank()) +
    ylab(element_blank()) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position='bottom') }

abun_type_vs_site <- d_both %>%
  group_by(site_code, year, qtr, abun_type) %>% # to even sample distribution
  dplyr::summarize(x=mean(x, na.rm=T)) %>%
  base_plot() +
  facet_wrap(~abun_type, scales='free') +
  geom_boxplot(aes(x = site_code, y = x, fill = site_code),
               position=position_dodge(), color = 'grey50')

abun_type_vs_site_year <- d_both %>%
  base_plot() +
  facet_grid(year~abun_type, scales='free') +
  geom_boxplot(aes(x = site_code, y = x, fill = site_code),
               position=position_dodge(), color = 'grey50') 
  
#------------------------------------------------------------------------------
# abun_type vs time
#------------------------------------------------------------------------------
base_abn <- function(my_data, my_x) {
  
  my_data$x_var <- my_data %>% pull(my_x)
  
  my_data %>%
    mutate(x_var = as_date(x_var)) %>%
    fix_site_order() %>%
    ggplot(aes(x = x_var, y = x, color = site_code, fill = site_code)) +
    scale_color_manual('Site', values = my_colors) +
    scale_fill_manual('Site', values=my_colors) +
    dark_theme_grey(base_size=14) +
    facet_wrap(~pretty_x(abun_type), scales = 'free_y', ncol=1) +
    scale_x_date(date_breaks='3 months', date_labels = "%y-%m") +
    ylab(element_blank())+
    xlab(pretty_x(my_x))+
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))}

abn_vs_time <- d_both %>%
  base_abn('year_month') +
  geom_vline(xintercept = ymd("2017-08-27"), linewidth=.3, color = 'white') +
  geom_vline(xintercept = ymd(c("2018-08-27", "2019-08-27", "2020-08-27")),
             linewidth=.1) +
  geom_jitter(shape = 21, size = 2, alpha = .4, show.legend=T,
              position=position_dodge(width=45)) +
  geom_jitter(shape = 21, size = 2, fill=NA, show.legend=T,
              position=position_dodge(width=45)) +
  geom_smooth(aes(group=NA), method = "loess", se=FALSE, color = 'purple',
              linewidth=1, span = 0.3, show.legend=F)

abn_den_vs_time_site <- d_both %>%
  create_site_group() %>%
  filter(abun_type == 'density') %>%
  base_abn('year_month') +
  facet_wrap(~site_group, ncol=1) +
  geom_vline(xintercept = ymd("2017-08-27"), linewidth=.3, color = 'white') +
  geom_vline(xintercept = ymd(c("2018-08-27", "2019-08-27", "2020-08-27")),
             linewidth=.3) +
  geom_jitter(shape = 21, size = 2, alpha = .4, show.legend=T,
              position=position_dodge(width=45)) +
  geom_jitter(shape = 21, size = 2, fill=NA, show.legend=T,
              position=position_dodge(width=45)) +
  geom_smooth(aes(group=site_group), color = 'purple',
              method = "loess", se=FALSE,
              linewidth=.4, span = .3, show.legend=F)

abn_bio_vs_time_site <- d_both %>%
  create_site_group() %>%
  filter(abun_type == 'biomass') %>%
  base_abn('year_month') +
  facet_wrap(~site_group, ncol=1) +
  geom_vline(xintercept = ymd("2017-08-27"), linewidth=.3, color = 'white') +
  geom_vline(xintercept = ymd(c("2018-08-27", "2019-08-27", "2020-08-27")),
             linewidth=.3) +
  geom_jitter(shape = 21, size = 2, alpha = .4, show.legend=T,
              position=position_dodge(width=45)) +
  geom_jitter(shape = 21, size = 2, fill=NA, show.legend=T,
              position=position_dodge(width=45)) +
  geom_smooth(aes(group=site_group), color = 'purple',
              method = "loess", se=FALSE,
              linewidth=.4, span = .3, show.legend=F)

#------------------------------------------------------------------------------
# autocorrelation
#------------------------------------------------------------------------------
trash_time_vars <- c('collection_period', 'year', 'month', 'day', 
                     'nday', 'qtr', 'year_week')

# fill missing predictors with linear interpolation
d_ste_fill <- d_ste %>%
  as_tsibble(key=site_code,
             index=collection_period) %>%  
  na_ma(weighting="linear") %>%
  ungroup() %>%
  as_tibble() %>%
  pivot_longer(cols=-any_of(c('site_code', 'collection_period')),
               names_to='xname',
               values_to = 'x') %>%
  group_by(site_code, collection_period, xname) %>%
  dplyr::summarize(x=mean(x ,na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from = x)

autocorrelate <- function(df) {
  df %>%
    select(where(is.numeric), -collection_period) %>%
    correlate() %>%
    select(term, biomass, density) %>%
    na.omit() }

abn_cor <- full_join(d_bio, d_den) %>%
  left_join(d_ste_fill) %>%
  group_by(site_code) %>%
  nest() %>%
  mutate(data_cor = map(data, autocorrelate)) %>%
  unnest(data_cor) %>%
  select(-data)

abn_correlation <- abn_cor %>%
  pivot_longer(cols=c('biomass', 'density'), 
               names_to='abun_type', values_to='r') %>%
  rename(predictor=term) %>%
  ggplot(aes(x=r, y = predictor, color = site_code, 
             fill=site_code, alpha = abs(r), size= abs(r^2))) +
  facet_wrap(~abun_type) +
  scale_color_manual('Site', values = my_colors) +
  scale_fill_manual('Site', values=my_colors) +
  dark_theme_grey(base_size=14) +
  geom_vline(xintercept = 0, color='grey90', lty=2, linewidth=.3) +
  geom_point(shape = 21, show.legend = F) 

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
# my_figures <- list(abn_correlation,
#                    abn_bio_vs_time_site,
#                    abn_den_vs_time_site,
#                    abn_vs_time,
#                    abn_density_vs_biomass2)
# 
# names(my_figures) <- c('abn_correlation',
#                        'abn_bio_vs_time_site',
#                        'abn_den_vs_time_site',
#                        'abn_vs_time',
#                        'abn_density_vs_biomass2')
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
# End biomass_visualize