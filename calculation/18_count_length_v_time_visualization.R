# count_length_v_time_visualization
# Sean Kinard
# 2023-06-27

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')

d_length <- read_csv("data/06_fill_data/fish_length_fill.csv")

d_bio_categories <- read_csv(
  "data/02_clean_data/fish_biological_scales_of_comparison.csv")

# -----------------------------------------------------------------------------
# Count and Length vs Time by groups
# -----------------------------------------------------------------------------
# Attempt multiple scaled metrics in the same plot
group_scale <- function(df, grp, grp_x, my_var)
{ # filters by a group and group_x and then scales the aggregated value by site, including instance of zero catch
  df <- df %>% 
    dplyr::rename(g = {{grp}})
  df <- df %>%
    dplyr::rename(v = {{my_var}})
  
  my_output <- df %>%
    filter(g == grp_x) %>%
    group_by(site_code, collection_period) %>%
    dplyr::summarize(value = sum(v)) %>%
    ungroup() %>%
    right_join(d_length%>%select(site_code, collection_period)%>% unique()) %>% 
    mutate(value = ifelse(is.na(value), 0 ,value)) %>%
    arrange(site_code, collection_period) %>%
    pivot_wider(values_from=value, names_from = site_code) %>%
    mutate(AR = scale(AR),
           EM = scale(EM),
           GC = scale(GC),
           MR = scale(MR),
           PD = scale(PD),
           PL = scale(PL),
           SF = scale(SF),
           TR = scale(TR),
           WM = scale(WM)) %>%
    pivot_longer(cols=AR:WM, names_to='site_code', values_to='value') %>%
    na.omit() 
  
  return(my_output)
}

group_scale_plot <- function(grp, grp_x) {
  
  d_count_s <- d_length %>%
    group_by(site_code, collection_period, lowest_taxon) %>%
    dplyr::summarize(cnt=length(lengthmm)) %>%
    ungroup() %>%
    left_join(d_bio_categories) %>%
    group_scale(grp=grp, grp_x = grp_x, my_var='cnt') %>%
    dplyr::rename(count_s=value)
  
  d_size_s <- d_length %>%
    group_by(site_code, collection_period, lowest_taxon) %>%
    dplyr::summarize(mu=mean(lengthmm)) %>%
    ungroup() %>%
    left_join(d_bio_categories) %>%
    group_scale(grp=grp, grp_x = grp_x, my_var='mu') %>%
    dplyr::rename(size_s=value)
  
  d_count_s %>%
    left_join(d_size_s) %>%
    fix_site_order() %>%
    pivot_longer(cols=contains('_s'), 
                 names_to='xvar', values_to='x') %>%
    mutate(xvar = str_replace_all(xvar, '_s','')) %>%
    ggplot(aes(collection_period, x, fill=xvar, color=xvar)) +
    geom_line(linewidth=.25, linetype=2) +
    geom_point(shape=21, size = 3, alpha = .1) +
    geom_point(shape=21, size = 3, fill=NA) +
    dark_theme_grey(base_size=12) +
    ylab(element_blank()) +
    theme(legend.position='top') +
    scale_fill_manual(grp_x, values=c('darkgoldenrod1', 'skyblue')) +
    scale_color_manual(grp_x, values=c('darkgoldenrod1', 'skyblue')) +
    facet_wrap(~site_code, ncol=1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

# Family
n_length_Poeciliidae <- group_scale_plot('family', 'Poeciliidae')
n_length_Centrarchidae <- group_scale_plot('family', 'Centrarchidae')
n_length_Ictaluridae <- group_scale_plot('family', 'Ictaluridae')
n_length_Lepisosteidae <- group_scale_plot('family', 'Lepisosteidae')
n_length_Leuciscidae <- group_scale_plot('family', 'Leuciscidae')
n_length_Cichlidae <- group_scale_plot('family', 'Cichlidae')
n_length_Cyprinodontidae <- group_scale_plot('family', 'Cyprinodontidae')

# Species
n_length_macrochirus <- group_scale_plot('lowest_taxon', 'L. macrochirus')
n_length_megalotis <- group_scale_plot('lowest_taxon', 'L. megalotis')
n_length_gulosus <- group_scale_plot('lowest_taxon', 'L. gulosus')
n_length_cyanellus <- group_scale_plot('lowest_taxon', 'L. cyanellus')
n_length_cyanoguttatum <- group_scale_plot('lowest_taxon', 'H. cyanoguttatus')
n_length_latipinna <- group_scale_plot('lowest_taxon', 'P. latipinna')
n_length_variegatus <- group_scale_plot('lowest_taxon', 'C. variegatus')
n_length_vigilax <- group_scale_plot('lowest_taxon', 'P. vigilax')

# Trophic Category
n_length_Herbivore <- group_scale_plot('diet_simple', 'Herbivore')
n_length_Omnivore <- group_scale_plot('diet_simple', 'Mixed')
n_length_Piscivore <- group_scale_plot('diet_simple', 'Carnivore')

# Oxygen Sensitivity
n_length_sensitive <- group_scale_plot('hypoxia_tolerance', 'Sensitive')
n_length_mesotolerant <- group_scale_plot('hypoxia_tolerance', 'Mesotolerant')
n_length_tolerant <- group_scale_plot('hypoxia_tolerance', 'Tolerant')

# Reproductive Strategy
n_length_SimpleNest <- group_scale_plot('reproduction', 'SimpleNest')
n_length_Bearer <- group_scale_plot('reproduction', 'Bearer')
n_length_ComplexNest <- group_scale_plot('reproduction', 'ComplexNest')
n_length_Migratory <- group_scale_plot('reproduction', 'Migratory')
n_length_Broadcaster <- group_scale_plot('reproduction', 'Broadcaster')

# Endemism
n_length_Local <- group_scale_plot('geographic_range', 'Local')
n_length_Restricted <- group_scale_plot('geographic_range', 'Restricted')
n_length_Unlimited <- group_scale_plot('geographic_range', 'Unlimited')

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
# my_objects <- ls()
# my_objects[str_detect(my_objects, 'n_length')]
# 
# my_figures <- list(
#   n_length_Bearer, n_length_Broadcaster, n_length_Centrarchidae, 
#   n_length_Cichlidae, n_length_ComplexNest, n_length_cyanellus, 
#   n_length_cyanoguttatum, n_length_Cyprinodontidae, n_length_gulosus,
#   n_length_Herbivore, n_length_Ictaluridae, n_length_latipinna, 
#   n_length_Lepisosteidae, n_length_Leuciscidae, n_length_Local, 
#   n_length_macrochirus, n_length_megalotis, n_length_mesotolerant, 
#   n_length_Migratory, n_length_Omnivore, n_length_Piscivore, 
#   n_length_Poeciliidae, n_length_Restricted, n_length_sensitive, 
#   n_length_SimpleNest, n_length_tolerant, n_length_Unlimited, 
#   n_length_variegatus, n_length_vigilax)
# 
# names(my_figures) <- my_objects[str_detect(my_objects, 'n_length')]
# 
# for (i in 1:length(my_figures)) {
#   my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
#   my_object <- my_figures[[i]]
#   ggsave(my_place,
#          plot = my_object,
#          width = 10,
#          height = 15,
#          units = c("in")) }

# -----------------------------------------------------------------------------
# End count_length_v_time_visualization