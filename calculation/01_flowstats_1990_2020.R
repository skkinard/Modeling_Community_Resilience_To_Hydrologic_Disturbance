# daily_discharge_and_storm_stats_1990_2020
# Sean Kinard
# 2023-02-24

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')

d_mega <- read_csv('data/00_source_data/Kinard_data/daily_discharge_and_storm_stats_1990_2020.csv') %>%
  add_rain() %>%
  mutate(site_code=factor(site_code)) %>%
  fix_site_order() 

hurricane_start <- ymd('2017-08-15')
hurricane_end <- ymd('2017-09-15')
regionflood_start <- ymd('2018-06-15')
regionflood_end <- ymd('2018-07-25')

# -----------------------------------------------------------------------------
# Visualization: Time Series
# -----------------------------------------------------------------------------
ts_base <- function(my_data) {
  my_data %>%
    ggplot(aes(collection_date, y=q_region)) +
    geom_point(alpha=.2, color='skyblue') +
    dark_theme_gray(base_size=14) +
    scale_x_date(date_breaks = '1 year', date_labels = "%y") +
    scale_color_manual(values=my_colors) +
    ylab('Regional Discharge') +
    scale_y_log10() +
    xlab('Year')
}

regional_discharge_time_series <- d_mega %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_med, na.rm=T)) %>%
  ts_base() +
  geom_smooth(method='loess', se=F, span=.2, color='red3')

regional_discharge_time_series_facet <- d_mega %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_mu, na.rm=T)) %>%
  mutate(time_period = ifelse(collection_date <= ymd('1997-01-01'), '1990-1997',
                       ifelse(collection_date > ymd('1997-01-01') &
                              collection_date <= ymd('2004-01-01'), '1997-2004',
                      ifelse(collection_date >= ymd('2004-01-01') &
                             collection_date < ymd('2011-01-01'), '2004-2011', '2011-2021')))) %>%
  ts_base() +
  facet_wrap(~time_period, scales='free_x', ncol=1) +
  geom_smooth(method='loess', se=F, span=.001, color='red3')
  
# -----------------------------------------------------------------------------
# Visualization: Histogram
# -----------------------------------------------------------------------------
hist_base <- function(my_data) {
  my_data %>%
    ggplot(aes(log(q), fill=site_code, color = site_code)) +
    geom_density(alpha=.1) +
    scale_color_manual(values=my_colors) +
    scale_fill_manual(values=my_colors) +
    dark_theme_grey() +
    ylab(element_blank()) +
    xlab('Log Discharge')
}

flow_density <- d_mega %>%
  hist_base()

flow_density_sites <- d_mega %>%
  hist_base() +
  facet_wrap(~site_code, ncol=3)

# -----------------------------------------------------------------------------
# Visualization: Base Flows
# -----------------------------------------------------------------------------
box_yr_base <- function(my_data) {
  my_data %>%
    ggplot(aes(site_code, val, color=annualrain, fill=annualrain)) +
    geom_boxplot(alpha=.1) +
    geom_boxplot(fill=NA) +
    #geom_point(size=5, shape = 21, alpha=.1) +
    #geom_point(size=5, shape = 21, fill=NA) +
    paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
    paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
    dark_theme_grey(base_size=12) +
    scale_y_log10() +
    xlab(element_blank()) +
    ylab('Discharge (l/s)') +
    theme(legend.position = 'none')
  
}

base_flows <- d_mega %>%
  select(annualrain, ye_med, site_code) %>%
  unique() %>%
  dplyr::rename(val = ye_med) %>%
  box_yr_base()

# -----------------------------------------------------------------------------
# Visualization: Season
# -----------------------------------------------------------------------------

# weekly
seas_weekly <- d_mega %>%
  select(year, week, we_mu, al_mu, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(week, y=we_mu/al_mu, color=site_code)) +
  theme_bw(base_size=12) +
  ylab('weekly_avg / yr_avg') +
  scale_x_continuous(breaks = seq(from=1, to=53, by=4)) +
  scale_color_manual('Site', values=my_colors) +
  scale_y_log10() +
  geom_vline(xintercept=week(ymd('2017-08-27')),
             lty=2, lwd=.5, color = 'black') +
  geom_label(aes(x=week(ymd('2017-08-27')), y = 0.18, label='Hurricane'),
             fill='grey90', color = 'black', show.legend = F) +
  geom_smooth(aes(group=site_code), method = "loess", se=F, span=.4, 
              color = 'black', lwd=1.5) +
  geom_smooth(method = "loess", se=F, span=.4, alpha=.3) + 
  coord_cartesian(xlim = c(15, 40)) +
  theme(legend.position = 'bottom',
    legend.background = element_rect(color='black', fill="white", size=.3))
  

# monthly
seas_monthly <- d_mega %>%
  filter(year>=2000) %>%
  select(year, ymonth, month, ym_mu, ye_med, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(month, y=ym_mu/ye_med, color=site_code)) +
  geom_smooth(method = "loess", se=F, span=.1, alpha=.3) +
  geom_point(alpha=.1) +
  dark_theme_gray() +
  scale_x_continuous(breaks = seq(from=1, to=12, by=1)) +
  ylab('monthly_mu / yr_med') +
  scale_color_manual(values=my_colors) +
  scale_y_log10()

# monthly_faceted
annual_variation_in_month <- d_mega %>%
  filter(year>=2000) %>%
  select(year, ymonth, month, ym_mu, ye_med, site_code, annualrain) %>%
  unique() %>%
  ggplot(aes(year, y=ym_mu/ye_med, color=site_code)) +
  geom_smooth(method = "loess", se=F, color='purple', span=.1) +
  geom_point(alpha=.3) +
  dark_theme_gray() +
  ylab('weekly_avg / yr_avg') +
  #scale_x_continuous(breaks = seq(from=1, to=60, by=3)) +
  scale_color_manual(values=my_colors) +
  scale_y_log10() +
  facet_wrap(~month, scales='free_x', ncol=2)

# Consistent troughs in February and August
# consistent peaks in May and September

# -----------------------------------------------------------------------------
# Visualization: Variability
# -----------------------------------------------------------------------------
# 1990-2020: measures of intra-annual variation
q_variability <- d_mega %>%
  select(annualrain, site_code, ye_FS, ye_DS, mo_rsd, mo_WCDC, ye_LFPP,
         ye_HFPP3, ye_HFPP7, ye_HFPP15) %>%
  pivot_longer(cols=ye_FS:ye_HFPP15, names_to='var', values_to = 'val') %>%
  box_yr_base() +
  ylab(element_blank()) +
  facet_wrap(~var, scales='free')
  
# -----------------------------------------------------------------------------
# Visualization Hurricane Harvey
# -----------------------------------------------------------------------------
# 2016-2019 regional hydrograph
hu_context <- d_mega %>%
  filter(year %in% 2016:2019) %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_med, na.rm=T)) %>%
  mutate(collection_period = collection_date) %>%
  create_date_vars() %>%
  unique() %>%
  ggplot(aes(nday, y=q_region)) +
  geom_vline(aes(xintercept=271), alpha = .5, linetype=2) +
  geom_point(alpha=.4, color='skyblue') +
  dark_theme_gray(base_size=14) +
  scale_color_manual(values=my_colors) +
  ylab('Regional Discharge Ratio (daily mean / annual mean)') +
  scale_y_log10() +
  xlab('Day') +
  facet_wrap(~year) +
  geom_smooth(method='loess', se=F, span=.1, color='red3')
  
# 2017 hydrographs faceted by site
hu_sites <- d_mega %>%
  filter(year == 2017) %>%
  filter(nday %% 7 == 0) %>%
  select(site_code, annualrain, yw_mu, yw_med, collection_date) %>%
  unique() %>%
  ggplot(aes(collection_date, y=yw_mu/yw_med, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=ymd('2017-08-27')), alpha=.5, linetype=2) +
  #geom_point(size=3, shape=21, alpha=.3) +
  #geom_point(size=3, shape=21, fill=NA) +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge Ratio (Weekly Mean / Annual Median)') +
  scale_y_log10() +
  xlab('Year') +
  geom_smooth(method='loess', se=F, span=.1, linewidth=.4, alpha=.7) +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# Visualize hurricane duration
hu_duration <- d_mega %>%
  filter(year == 2017) %>%
  filter(collection_date > hurricane_start) %>%
  ggplot(aes(collection_date, y=q, 
             fill=site_code, color = site_code)) +
  facet_wrap(~site_code) +
  scale_color_manual('Site', values=my_colors) +
  scale_fill_manual('Site', values=my_colors) +
  theme_bw(base_size=12) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  ylab('Discharge (l/s)') +
  scale_y_log10() +
  xlab('Month') +
  theme(legend.position='none') +
  geom_vline(aes(xintercept=hu_flood_max_date), 
             alpha=.8, linetype=2, color='black') +
  geom_vline(aes(xintercept=hu_return_date), 
             alpha=.8, linetype=2, color='black') +
  geom_hline(data = filter(d_mega, collection_date==hurricane_start),
             aes(yintercept=al_med), 
             alpha=.8, linetype=2, color='black') +
  geom_line(color='black', lwd=1.5) +
  geom_line(lwd=1.1)
  

# -----------------------------------------------------------------------------
# Visualization June 2018 Regional Flood Event
# -----------------------------------------------------------------------------
# 2017-2020 regional hydrograph
rf_context <- d_mega %>%
  filter(year %in% 2017:2020) %>%
  group_by(collection_date) %>%
  dplyr::summarize(q_region = mean(q/al_med, na.rm=T)) %>%
  mutate(collection_period=collection_date) %>%
  create_date_vars() %>%
  ggplot(aes(nday, y=q_region)) +
  geom_vline(aes(xintercept=200), alpha=.5, linetype=2) +
geom_point(alpha=.4, color='skyblue') +
  dark_theme_gray(base_size=14) +
  scale_color_manual(values=my_colors) +
  ylab('Regional Discharge Ratio (daily mean / annual mean)') +
  scale_y_log10() +
  xlab('Day') +
  facet_wrap(~year) +
  geom_smooth(method='loess', se=F, span=.1, color='red3')

# 2018 sites faceted
rf_sites <- d_mega %>%
  filter(year == 2018) %>%
  ggplot(aes(collection_date, y=yw_mu/yw_med, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=rf_flood_max_date), alpha=.5, linetype=2) +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge Ratio (Weekly Mean / Annual Median)') +
  scale_y_log10() +
  xlab('Month') +
  geom_smooth(method='loess', se=F, span=.1, linewidth=.4, alpha=.7) +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# Visualize storm duration
rf_duration <- d_mega %>%
  filter(year == 2018) %>%
  filter(collection_date > regionflood_start) %>%
  select(site_code, annualrain, q, mo_med, 
         collection_date, rf_return_date,rf_flood_max_date) %>%
  unique() %>%
  ggplot(aes(collection_date, y=q, 
             fill=annualrain, color = annualrain)) +
  geom_vline(aes(xintercept=rf_flood_max_date), alpha=.5, linetype=2) +
  geom_vline(aes(xintercept=rf_return_date), alpha=.5, linetype=2) +
  geom_hline(data = filter(d_mega, collection_date==regionflood_start),
             aes(yintercept=al_med), alpha=.5, linetype=2) +
  geom_line() +
  dark_theme_gray(base_size=14) +
  scale_x_date(date_breaks = '1 month', date_labels = "%m") +
  paletteer::scale_fill_paletteer_c("grDevices::Zissou 1", direction=-1) +
  paletteer::scale_color_paletteer_c("grDevices::Zissou 1", direction=-1) +
  ylab('Discharge (l/s)') +
  scale_y_log10() +
  xlab('Month') +
  #geom_smooth(method='loess', se=F, span=.1, linewidth=.4, alpha=.7) +
  facet_wrap(~site_code) +
  theme(legend.position='none')

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
# # export megaframe
# write_csv(d_mega, 'Data/daily_discharge_and_storm_stats_1990_2020.csv')
# 
# # Site Summary Table
# site_summary_table <- fstats_alltime %>% 
#   left_join(RI_hur%>%select(-hu_mo_med)%>%unique()) %>%
#   left_join(RI_rfl%>%select(-rf_mo_med)%>%unique())
# # export summarytable
# write_csv(site_summary_table, 'Data/storm_stats_summary.csv')
# 
# ggsave('Figures/base_flows.pdf',
#        plot = base_flows,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/flow_density.pdf',
#        plot = flow_density,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/flow_density_sites.pdf',
#        plot = flow_density_sites,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/hu_duration.pdf',
#        plot = hu_duration,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/hu_sites.pdf',
#        plot = hu_sites,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/q_variability.pdf',
#        plot = q_variability,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/regional_discharge_time_series.pdf',
#        plot = regional_discharge_time_series,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/regional_discharge_time_series_facet.pdf',
#        plot = regional_discharge_time_series_facet,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/rf_context.pdf',
#        plot = rf_context,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/rf_duration.pdf',
#        plot = rf_duration,
#        width = 9,
#        height = 9,
#        units = c("in"))
# 
# ggsave('Figures/rf_sites.pdf',
#        plot = rf_sites,
#        width = 9,
#        height = 9,
#        units = c("in"))

# -----------------------------------------------------------------------------

# End daily_discharge_and_storm_stats_1990_2020











