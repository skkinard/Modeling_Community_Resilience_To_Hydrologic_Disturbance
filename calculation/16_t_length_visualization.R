# t_length_visualization
# Sean Kinard
# 2023-06-26
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')
library(tsibble)
library(ggpubr)

d_len <- read_csv('data/06_fill_data/fish_length_fill.csv')

d_lte <- read_csv('exploration/output/environment_longterm_with_ln.csv')

d_ste <- read_csv('exploration/output/environment_shortterm_with_ln.csv')
#------------------------------------------------------------------------------
# Calculations: L_mu (average length)
#------------------------------------------------------------------------------
# make L_mu for community , taxonomic family, taxonomic species
my_avg <- function(my_data) {
  
  L_com <- my_data %>%
    group_by(site_code, collection_period) %>%
    dplyr::summarize(L_mu_com = mean(lengthmm),
                     L_me_com = median(lengthmm),
                     L_sd_com = sd(lengthmm),
                     L_n_com = length(lengthmm)) %>% ungroup()
  
  L_fam <- my_data %>%
    group_by(site_code, collection_period, family) %>%
    dplyr::summarize(L_mu_fam = mean(lengthmm),
                     L_me_fam = median(lengthmm),
                     L_sd_fam = sd(lengthmm),
                     L_n_fam = length(lengthmm)) %>% ungroup()
  
  L_spe <- my_data %>%
    group_by(site_code, collection_period, lowest_taxon) %>%
    dplyr::summarize(L_mu_spe = mean(lengthmm),
                     L_me_spe = median(lengthmm),
                     L_sd_spe = sd(lengthmm),
                     L_n_spe = length(lengthmm)) %>%
    filter(L_mu_spe>0) %>%
    left_join(select(my_data, family, lowest_taxon) %>% unique())
  
  my_output <- L_spe %>%
    left_join(L_fam) %>%
    left_join(L_com)
  
  return(my_output) }

d <- my_avg(d_len)

#------------------------------------------------------------------------------
# Prep data for figures
#------------------------------------------------------------------------------
d <- d %>% create_date_vars()

# Add some environmental predictors
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

d <- d %>% add_rain() %>% left_join(d_ste_fill) %>% fix_site_order()
  
#------------------------------------------------------------------------------
# t_length vs site
#------------------------------------------------------------------------------
t_length_com_vs_site <-  d %>%
  ggplot(aes(x = site_code, y = L_mu_com, fill = site_code)) +
  facet_wrap(~year, ncol=2) +
  geom_boxplot(color='grey75', show.legend = F) +
  scale_fill_manual(values=my_colors) +
  dark_theme_grey(base_size=14) +
  ylab('Total Length (mm)') +
  xlab(element_blank())

#------------------------------------------------------------------------------
# t_length vs time
#------------------------------------------------------------------------------
base_t_length <- function(my_data, my_x, xspan) {
  df <- my_data
  colnames(df) <- str_replace_all(colnames(df), my_x, 'XXX')
  
  df %>%
    fix_site_order() %>%
    ggplot(aes(x = collection_period, y = XXX,
               color = site_code, fill = site_code)) +
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

d_com <- d %>% 
  select(site_code, collection_period, L_mu_com) %>%
  unique() %>%
  left_join(d_ste_fill) %>% 
  fix_site_order() %>%
  create_site_group() %>%
  create_date_vars() %>%
  add_rain()

t_length_com_vs_time <-  d_com %>% 
  base_t_length('L_mu_com', .3) +
  facet_wrap(~site_group, ncol=1)

autocorrelate <- function(df) {
  df %>%
    select(where(is.numeric), -year_month) %>%
    correlate() %>%
    select(term, L_mu_com) %>%
    na.omit() %>%
    arrange(desc(abs(L_mu_com))) }

t_length_cor <- d_com %>%
  group_by(site_code) %>%
  nest() %>%
  mutate(data_cor = map(data, autocorrelate)) %>%
  unnest(data_cor) %>%
  select(-data)

t_length_correlation <- t_length_cor %>%
  rename(r=L_mu_com, predictor=term) %>%
  ggplot(aes(x=r, y = predictor, color = site_code, 
             fill=site_code, alpha = abs(r), size= abs(r^2))) +
  scale_color_manual('Site', values = my_colors) +
  scale_fill_manual('Site', values=my_colors) +
  dark_theme_grey(base_size=14) +
  geom_vline(xintercept = 0, color='grey90', lty=2, linewidth=.3) +
  geom_point(shape = 21, show.legend = F) 

#------------------------------------------------------------------------------
# Taxa-specific densities vs time
#------------------------------------------------------------------------------
my_families <- c('Centrarchidae', 'Chichlidae', 'Lepisosteidae', 'Leuciscidae', 'Poeciliidae')

my_centrarchids <- c('L. auritus',  "L. macrochirus", "L. megalotis", 
                     "L. cyanellus", "L. gulosus")

small_gape <- c('H. cyanoguttatum', 'L. auritus',  "L. macrochirus", "L. megalotis")

large_gape <- c("L. cyanellus", "L. gulosus", "M. salmoides")

t_length_fam_vs_time <- d %>% 
  filter(family %in% my_families) %>%
  base_t_length('L_mu_fam', .3) +
  facet_wrap(~family, ncol=1, scales='free_y')

t_length_poec_vs_time <- d %>% 
  filter(family == 'Poeciliidae') %>%
  mutate(lowest_taxon = fct_relevel(
    lowest_taxon, 'G. affinis', 'P. latipinna', 'P. formosa')) %>%
  base_t_length('L_mu_spe', 0.4) +
  facet_wrap(~lowest_taxon, ncol=1)

t_length_sgap_vs_time <- d %>% 
  filter(lowest_taxon %in% small_gape) %>%
  base_t_length('L_mu_spe', 0.4) +
  facet_wrap(~lowest_taxon, ncol=1)

t_length_lgap_vs_time <- d %>% 
  filter(lowest_taxon %in% large_gape) %>%
  base_t_length('L_mu_spe', 0.4) +
  facet_wrap(~lowest_taxon, ncol=1)
 
#------------------------------------------------------------------------------
# Annual t_length (AD)
#------------------------------------------------------------------------------

# base plot function
L_base <- function(mydata) {
  
  mydata %>%
    ggplot(aes(x=annualrain, y = y_var)) +
    geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=2,
                formula =  y ~ poly(x, 2)) +
    geom_smooth(method = "lm", se=FALSE, color="grey20", linetype=1,
                formula =  y ~ x) +
    geom_point(shape = 21, size = 4, fill = 'white', alpha = .5)+
    stat_cor(label.y = 115, size = 6)+ 
    stat_regline_equation(label.y = 125, size = 6) + 
    ylab('t_length (mm)') +
    xlab('Annual Rainfall (cm)') +
    dark_theme_grey(base_size = 14) }

t_length_com_vs_rainfall_annual <- d %>%
  select(! c(L_mu_fam, L_mu_spe, family, lowest_taxon)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_com) %>%
  L_base() + facet_wrap(~year)

t_length_cent_vs_rainfall_annual <- d %>%
  filter(lowest_taxon %in% my_centrarchids) %>%
  select(! c(L_mu_fam,family)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_spe) %>%
  L_base() + facet_wrap(~year, ncol=2)

t_length_poec_vs_rainfall_annual <- d %>%
  filter(family == 'Poeciliidae') %>%
  select(! c(L_mu_fam,family)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_spe) %>%
  L_base() + facet_wrap(~year, ncol=2)

caption_t_length <- 'Fish total length (mm) plotted against Rainfall (cm/yr). The vertical axis is log-transformed. Data comes from monthly sampling in 2017 and 2018, and quarterly sampling in 2020.'

#------------------------------------------------------------------------------
# Quarterly t_length (QD)
#------------------------------------------------------------------------------
t_length_com_vs_rainfall_quarterly <- d %>%
  select(! c(L_mu_fam, L_mu_spe, family, lowest_taxon)) %>%
  unique() %>%
  dplyr::rename(y_var = L_mu_com) %>%
  L_base() +
  facet_wrap(~year+qtr)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
# write_csv(d, 'exploration/output/fish_Length_stats.csv')
# write_csv(t_length_cor, 'exploration/output/t_length_cor.csv')
# 
# 
# my_figures <- list(t_length_com_vs_rainfall_quarterly,
#                    t_length_poec_vs_rainfall_annual,
#                    t_length_cent_vs_rainfall_annual,
#                    t_length_com_vs_rainfall_annual,
#                    t_length_lgap_vs_time,
#                    t_length_sgap_vs_time,
#                    t_length_poec_vs_time,
#                    t_length_fam_vs_time,
#                    t_length_com_vs_time,
#                    t_length_correlation,
#                    t_length_com_vs_site)
# 
# names(my_figures) <- c('t_length_com_vs_rainfall_quarterly',
#                        't_length_poec_vs_rainfall_annual',
#                        't_length_cent_vs_rainfall_annual',
#                        't_length_com_vs_rainfall_annual',
#                        't_length_lgap_vs_time',
#                        't_length_sgap_vs_time',
#                        't_length_poec_vs_time',
#                        't_length_fam_vs_time',
#                        't_length_com_vs_time',
#                        't_length_correlation',
#                        't_length_com_vs_site')
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
# End t_length_visualization