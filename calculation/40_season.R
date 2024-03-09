# Season
# Sean Kinard
# 2023-08-28
#------------------------------------------------------------------------------
# setup
#-----------------------------------------------------------------------------
source('exploration/toolkit.R')

d_rain <- read_csv("data/06_fill_data/site_rain_match.csv") %>%
  mutate(collection_period = case_when( # undo lag-1 from other models
    month(collection_period) == 1 ~ paste(
      (year(collection_period)-1), '-12-01', sep='') %>% ymd(),
    T ~ paste(year(collection_period), 
              (month(collection_period)-1), 
              '01', sep='-')%>% ymd())) %>%
  rename(rain=rain_mL_m2) %>%
  select(rain, site_code, collection_period)
  
d_temp <- read_csv("data/04_match_data/site_transect_match.csv") %>%
  select(site_code, collection_period, temperature)

dc <- left_join(d_temp, d_rain) %>% filter(collection_period > ymd('2017-10-01')) %>%
  scale_by_site(xvar='rain') %>%
  scale_by_site(xvar='temperature')

colnames(dc) <- str_replace_all(colnames(dc), '_scaled', '')

#-----------------------------------------------------------------------------
# Visualize rain and temperature
#-----------------------------------------------------------------------------
my_season <- tibble(
  label = c('Winter', 'Spring', 'Summer', 'Fall'),
  x = c(11, 4, 7, 9)%>%month(),
  y=2)

dc_long <- dc %>%
  add_rain() %>%
  mutate(Month=month(collection_period)) %>%
  pivot_longer(cols=c(rain, temperature), names_to='xvar', values_to='xval') %>%
  mutate(xvar=str_to_title(xvar))

dc_summary <- dc_long %>%
  group_by(Month, xvar) %>%
  dplyr::summarize(Z_mu = mean(xval),
                   Z_sd = sd(xval)) %>%
  ungroup()

table_season <-  dc_summary %>%
  mutate(label = paste(round(Z_mu,2), " (", round(Z_sd,2), ")", sep='')) %>%
  select(-contains('Z_')) %>%
  pivot_wider(names_from=xvar, values_from=label) 

plot_season <- dc_long %>%
  ggplot(aes(Month, y=xval)) +
  geom_hline(yintercept=0, lty=2, lwd=.5, color='grey70') +
  geom_vline(xintercept = my_season$x, color='grey70') +
  geom_label(data=my_season,
             aes(x=x, y=y, label=label), fill='white', nudge_x=.6) +
  geom_smooth(aes(color=xvar), method='loess', span=.3, lty=2, lwd=.4, se=F) +
  geom_point(data=dc_summary, aes(x=Month, y=Z_mu, fill = xvar),
             shape=21, color='black', size=4) +
  scale_color_manual(values=c('blue', 'red')) +
  scale_fill_manual(values=c('blue', 'red')) +
  scale_x_continuous(breaks=1:12) +
  labs(x='Month', y='Z-Score') +
  theme_bw(base_size=14) +
  theme(legend.position=c(.1, .9),
        legend.title = element_blank())

# Seasons
# Cool and Dry 12-4 ("Winter")
# Hot and Wet 5-7 ("Spring")
# Hot and Dry 8-9 ("Summer")
# Hot and Wet 10-11 ("Fall")

#-----------------------------------------------------------------------------
# Export
#-----------------------------------------------------------------------------
# ggsave('exploration/visualization/plot_season.png',
#        plot = plot_season,
#        width = 10,
#        height = 5,
#        units = c("in"))
# 
# write_csv(table_season, 'exploration/output/table_season.csv')

#-----------------------------------------------------------------------------
# End Season