# trait_density_vs_time
# Sean Kinard
# 2023-06-28
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')
library(BiodiversityR)

d_density <- read_csv("data/06_fill_data/fish_density_fill.csv")
#d_biomass <- read_csv("data/06_fill_data/fish_biomass.csv")
d_bio_categories <- read_csv("data/02_clean_data/fish_biological_scales_of_comparison.csv")

#------------------------------------------------------------------------------
# trait density functions
#------------------------------------------------------------------------------
create_d_trait <- function(xtrait) {
  my_data <- d_density %>% left_join(d_bio_categories)
  colnames(my_data) <- str_replace_all(colnames(my_data), xtrait, 'XXX')
  
  output <- my_data %>%
    group_by(site_code, collection_period, XXX) %>%
    dplyr::summarize(density = sum(density, na.rm=T)) %>%
    ungroup()
  
  colnames(output) <- str_replace_all(colnames(output), 'XXX', xtrait)
  return(output) }

trait_density_plot <- function(df, xgroup, xspan=.5) {
  my_data <- df
  colnames(my_data) <- str_replace_all(colnames(my_data), xgroup, 'XTRAIT')
  
  my_data %>%
    scale_by_site('density') %>%
    fix_site_order() %>%
    ggplot(aes(x=collection_period, y=density_scaled, 
               color = site_code, fill=site_code)) +
    facet_wrap(~XTRAIT, scales='free', ncol=1) +
    geom_point(size = 4, alpha=.1, shape = 21) +
    geom_point(size = 4, fill=NA, shape = 21) +
    geom_smooth(aes(group=NA), method='loess', se=F, span=xspan,
                color='purple', linewidth=1, show.legend=F) +
    dark_theme_grey(base_size = 14) +
    scale_color_manual(values=my_colors) +
    scale_fill_manual(values=my_colors) }

#------------------------------------------------------------------------------
# trait density visualization
#------------------------------------------------------------------------------
# reproduction
trait_density_reproduction_table <- create_d_trait('reproduction')
trait_density_reproduction_vs_time <- trait_density_reproduction_table %>% 
  filter(reproduction != 'Migratory') %>%
  trait_density_plot('reproduction')

# diet_simple
trait_density_diet_simple_table <- create_d_trait('diet_simple')
trait_density_diet_simple_vs_time <- trait_density_diet_simple_table %>% 
  trait_density_plot('diet_simple')

# substrate_simple
trait_density_substrate_simple_table <- create_d_trait('substrate_simple')
trait_density_substrate_simple_vs_time <- trait_density_substrate_simple_table %>% 
  trait_density_plot('substrate_simple')

# macrohabitat_simple
trait_density_macrohabitat_simple_table <- create_d_trait('macrohabitat_simple')
trait_density_macrohabitat_simple_vs_time <- trait_density_macrohabitat_simple_table %>% 
  trait_density_plot('macrohabitat_simple')

# mesohabitat_simple
trait_density_mesohabitat_simple_table <- create_d_trait('mesohabitat_simple')
trait_density_mesohabitat_simple_vs_time <- trait_density_mesohabitat_simple_table %>% 
  trait_density_plot('mesohabitat_simple')

# hypoxia_tolerance
trait_density_hypoxia_tolerance_table <- create_d_trait('hypoxia_tolerance')
trait_density_hypoxia_tolerance_vs_time <- trait_density_hypoxia_tolerance_table %>% 
  trait_density_plot('hypoxia_tolerance')

# geographic_range
trait_density_geographic_range_table <- create_d_trait('geographic_range')
trait_density_geographic_range_vs_time <- trait_density_geographic_range_table %>% 
  trait_density_plot('geographic_range')

#------------------------------------------------------------------------------
# Export tables
#------------------------------------------------------------------------------
my_objects <- ls()
my_csv_names <- my_objects[str_detect(my_objects, '_table')]

my_csvs <- list( 
  trait_density_diet_simple_table, trait_density_geographic_range_table,
  trait_density_hypoxia_tolerance_table, trait_density_macrohabitat_simple_table,
  trait_density_mesohabitat_simple_table, trait_density_reproduction_table, 
  trait_density_substrate_simple_table)

names(my_csvs) <- my_csv_names

for (i in 1:length(my_csvs)) {
  my_place <- paste('exploration/output/', names(my_csvs[i]), ".csv", sep='')
  my_object <- my_csvs[[i]]
  write_csv(my_object, my_place) }

#------------------------------------------------------------------------------
# Export figures
#------------------------------------------------------------------------------
# my_figure_names <- my_objects[str_detect(my_objects, '_vs_time')]
# 
# my_figures <- list(
#   trait_density_diet_simple_vs_time, trait_density_geographic_range_vs_time,
#   trait_density_hypoxia_tolerance_vs_time, trait_density_macrohabitat_simple_vs_time,
#   trait_density_mesohabitat_simple_vs_time, trait_density_reproduction_vs_time,
#   trait_density_substrate_simple_vs_time)
# 
# names(my_figures) <- my_figure_names
# 
# for (i in 1:length(my_figures)) {
#   my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
#   my_object <- my_figures[[i]]
#   ggsave(my_place,
#          plot = my_object,
#          width = 10,
#          height = 15,
#          units = c("in")) }

#------------------------------------------------------------------------------
# End trait_density_vs_time