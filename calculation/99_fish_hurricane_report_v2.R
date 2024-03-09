# 99_fish_hurricane_report_v2
# Sean Kinard
# 2023-07-10
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
source('exploration/toolkit.R')

d_fish <- read_csv('exploration/output/model_combined_vars_SCALED.csv') %>%
  add_rain() %>%
  is_baseline(end_test = "2018-12-31") # set hurricane test period limit

d_rain <- read_csv('data/02_clean_data/site_rain.csv') %>%
  mutate(rain_L_m2=rain_m3/watershed_km2/1000)
# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
# Broad Trends in Biological Response Time Series
plot_site_scaled_biovars_vs_time <- d_fish %>%
  fix_site_order() %>%
  select(collection_period, site_code, 
         ln_density, ln_biomass, richness, shannon) %>%
  pivot_longer(cols=-any_of(c('site_code', 'collection_period')),
               names_to='xname', values_to='xvalue') %>%
  ggplot(aes(collection_period, xvalue, color=site_code, fill = site_code)) +
  theme_bw(base_size=14) +
  theme(legend.position = 'bottom',
        legend.background = element_rect(color='black', fill="white", size=.3)) +
  facet_wrap(~xname, scales='free_y') +
  scale_color_manual('Site', values=my_colors) +
  scale_fill_manual('Site', values=my_colors) +
  geom_vline(xintercept=ymd('2017-08-27'), color='red', lwd=.3, lty=2, 
             color='white') +
  geom_point(shape = 21, size = 3, alpha = .6) +
  geom_point(shape = 21, size = 3, fill=NA) +
  geom_point(shape = 21, size = 3, color='black', fill=NA) +
  geom_smooth(aes(group=NA),
              method='loess', se=F, span=.3, lwd=1.5, color = 'black',
              show.legend = F) +
  geom_smooth(aes(group=NA),
              method='loess', se=F, span=.3, lwd=.7, color = 'chartreuse',
              show.legend = F) +
  labs(x=element_blank(), y='Site-Scaled Value')

# Considerable Month-Month Variation
plot_density_vs_month <- d_fish %>%
  fix_site_order() %>%
  filter(collection_period < ymd("2019-01-01") &
           collection_period > ymd("2017-08-01")) %>%
  ggplot(aes(collection_period, ln_density, color=site_code, fill=site_code )) +
  theme_bw(base_size=14) +
  scale_color_manual('Site', values=my_colors) +
  scale_fill_manual('Site', values=my_colors) +
  labs(x='Month', y='Density') +
  facet_wrap(~site_code, scales='free_y') +
  scale_x_date(breaks = '2 month', date_labels = "%m") +
  geom_line(lty=2, lwd=.3, color='grey25') +
  geom_point(shape = 21, size = 3, alpha = .6, show.legend = F) +
  geom_point(shape = 21, size = 3, fill=NA, show.legend = F) +
  geom_point(shape = 21, size = 3, color='black', fill=NA, show.legend = F)

# Seasonal Rains
d_rain_mu <- d_rain %>% 
  group_by(year, month) %>% 
  dplyr::summarize(xmu=mean(rain_L_m2))

year_colors <- c('yellow', 'orange', 'white', 'lightgreen', 'skyblue')

plot_rain_vs_month <- d_rain %>%
  ggplot(aes(x=as.numeric(month), y = rain_L_m2)) +
  theme_bw(base_size=14) +
  scale_x_continuous(breaks=seq(1:12)) +
  scale_color_manual('Site', values=my_colors) +
  labs(x='Month', y='Rain (L/m2)') +
  geom_vline(xintercept=8.75, lty=2, lwd=.4, color='black') +
  geom_label(aes(x=8.75, y = 450, label='Hurricane'),
             fill='grey90', color = 'black', show.legend = F) +
  geom_violin(aes(group=month), fill='grey40') +
  
  geom_label(data = d_rain_mu, aes(as.numeric(month), xmu, 
                                   label=substr(year, 3, 4), fill=as.factor(year)),
             color='black') +
  geom_smooth(aes(group=NA), method='loess', span=.3, se=F,
                                         lty=1, lwd=1.5, color='black', show.legend = F)+
  geom_smooth(aes(group=NA), method='loess', span=.3, se=F,
              lty=1, lwd=1.1, color='purple', show.legend = F) +
  scale_fill_manual('Year', values=year_colors)+
  theme(legend.position='none')
  

# -----------------------------------------------------------------------------
# Tables
# -----------------------------------------------------------------------------

# Uneven Time Series Concerns
table_fish_dates <- d_fish %>%
  group_by(collection_period) %>%
  summarize(n=length(site_code)) %>%
  ungroup() %>%
  mutate(year=year(collection_period),
         month=month(collection_period)) %>%
  select(-collection_period) %>%
  arrange(month) %>%
  pivot_wider(names_from=year, values_from=n) %>%
  select(month, `2017`, `2018`, `2019`, `2020`)

plot_fish_dates <- d_fish %>%
  group_by(collection_period) %>%
  summarize(n=length(site_code)) %>%
  ungroup() %>%
  mutate(year=year(collection_period),
         month=month(collection_period)) %>%
  select(-collection_period) %>%
  arrange(month) %>%
  pivot_wider(names_from=year, values_from=n, values_fill = 0) %>%
  pivot_longer(cols=-month, names_to='year', values_to='xvalue') %>%
  ggplot(aes(month, xvalue)) +
  facet_wrap(~year, ncol=1) +
  theme_bw(base_size=14) +
  scale_x_continuous(breaks=seq(1:12)) +
  scale_y_continuous(breaks=seq(0,10, by=2)) +
  labs(x='Month', y='Sites Sampled') +
  scale_color_manual('Year', values=year_colors) +
  scale_fill_manual('Year', values=year_colors) +
  geom_vline(xintercept=8.75, lty=2, lwd=.4, color='black') +
  geom_label(aes(x=8.75, y = 8, label='HH'),
             fill='grey90', color = 'black', show.legend = F) +
  # geom_jitter(size=3, shape = 24, color='black')
  geom_line(lty=2, lwd=.5, color = 'black') +
  geom_point(shape=21, size=4, color='black', fill='palegreen3')
    
# -----------------------------------------------------------------------------
# Broad Trends in Biological Response Time Series
# -----------------------------------------------------------------------------
# Abundance based metrics indicate a region-wide bump in the first quarter of 2018. There is no way to attribute the bump to the hurricane since there is evidence of a similar bump occuring in the fourth quarter of 2020. This may indicate that there is a seasonal bump in fish count and/or mass in the winter quarters.

# Biomass drops from spring 2017 values after the hurricane and remained low for one month. However, after 2018, biomass remained low throughout 2019. Biomass recovered to pre-hurricane levels in the second quarter of 2020. 

# The trends indicated by the smoothed moving average are deceiving, since the localized average captures trajectories in 2018 (because of frequent sampling), but not other years (because of infrequent sampling). However, it is impossible to see any month-to-month or seasonal variation in the other years. Comparing storm-effects to the singular set of 7 site, pre-storm data is tempting. But the only monthly data we have (2018) indicates a drastic decline in fish abundance and richness between March and August (likely due to annual hydrologic drought in July). So the fact that the pre-storm data was collected in a different season most likely does not represent actual pre-storm communities in mid August. 

# plot_site_scaled_biovars_vs_time

# -----------------------------------------------------------------------------
# Considerable Month-Month Variation
# -----------------------------------------------------------------------------
# Month-to-month variation is extraordinary in some of the sampled streams, as assessed in 2018. Wild drops and recoveries of fish densities occur throughout the year at 6/9 sites.

# plot_density_vs_month

# -----------------------------------------------------------------------------
# Uneven Time Series Concerns
# -----------------------------------------------------------------------------
# Post-2018 data also does not provide obvious baseline conditions. There is no way to accurately define seasonal biological trends, whose timing can vary weeks to months. Splitting years in quarters can bias comparisons if one seasonal grouping is offset by a month. Also, TERRG sampling intervals were uneven and failed to capture the 'spring' or second quarter: samples in 2020 were collected in January (m=01), late May(m=06), late September (m=10), and November (m=11). So spring comparisons are limited to 2017 (with only 7 sites), 2018, and 2019. Fall samples include 2017, 2019 and 2020 . But we are missing fall 2018 which was impeded by technical, weather, and manpower issues. 2019 is missing quarter one and three. In summary, seasonal comparisons are data deficient across the board.

# table_fish_dates

# -----------------------------------------------------------------------------
# Insufficient Baseline Data For Disturbance Comparison
# -----------------------------------------------------------------------------
# There are few options for defining 'baseline' conditions to compare with the alleged storm-effects. Pre-storm data lacks 22% of the sites, and is 4 months removed from the event, during which time we expect seasonal drying to drastically lower abundance and diversity within streams. The expected seasonal effect would produce the same observed effect of the storm. So the temporal distance and inadequate coverage renders the pre-storm data insufficient for disturbance assessment. Non-testing period data (data after 2018) is so infrequent, that it's impossible to ascertain seasons and associated baselines, causing a similar issue to using the pre-storm data. So the final option is to aggregrate all data outside of the testing period and create a baseline akin to an 'annual average'. But the annual averages for 2017-prestorm, 2019, and 2020 would contain 1, 2, and 4 time periods respectively. So, one of the annual averages is not an average at all, the second hardly qualifies as an average, and the third year barely qualifies. Lastly, Most of the data prospected for establishing 'baseline' conditions occur after the storm, and should therefore be influenced by the storm.

# -----------------------------------------------------------------------------
# T-tests identify 'abnormal' Collection Periods After Hurricane Harvey
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Hypothetical Baseline Calculations
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Hypothetical Disturbance Calculations
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Hypothetical Disturbance Hypothesis tests
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Hypothetical Disturbance Findings
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
#  Hinderances to Disturbance Inferences
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Descriptive Models for Sub-Tropical Coastal Stream Fish Communities
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Descriptive Model Inferences
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Mitigating Hindrances to Descriptive Inferences
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Conclusions
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Supplemental Figures
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Supplemental Tables
# -----------------------------------------------------------------------------