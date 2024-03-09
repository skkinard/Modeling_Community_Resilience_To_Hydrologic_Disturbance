# site_map
# Sean Kinard
# 1-25-2023

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidygeocoder)
library(maps)
library(ggrepel)
library(paletteer)
library(ozmaps) 
library(grid)
library(gt)

# state data
state_map_data <- ne_states(country = 'united states of america', 
                            returnclass = "sf")

# rivers
rivers_tx <- read_sf(
  dsn = "data/00_source_data/Map_Data/Surface_Water/Surface_Water.shp", 
  layer = "Surface_Water")

# precip
precip <- read_sf(dsn = "data/00_source_data/Map_Data/precip1981_2010_a_tx.shp", 
                  layer = "precip1981_2010_a_tx") %>%
  rename('Rainfall (cm/yr)' = 'Precip_cm')

# oceans
oceans <- read_sf(
  dsn = 'data/00_source_data/Map_Data/natural_earth_vector/10m_physical/ne_10m_ocean.shp', 
  layer = "ne_10m_ocean")

# cities
cities <- read_sf(
  dsn = 'data/00_source_data/Map_Data/natural_earth_vector/10m_cultural/ne_10m_populated_places_simple.shp', 
  layer = "ne_10m_populated_places_simple") %>%
  filter(latitude > 27.25) %>%
  filter(latitude < 29.25) %>%
  filter(longitude > -98.25) %>%
  filter(longitude < -95.75) 

# urban
urban <- read_sf(
  dsn = 'data/00_source_data/Map_Data/natural_earth_vector/10m_cultural/ne_10m_urban_areas.shp', 
  layer = "ne_10m_urban_areas")

# Site labels
d_site <- read_csv('data/02_clean_data/watershed.csv') %>%
  select(site_code, lat, lon)

#------------------------------------------------------------------------------
# Site Map
#------------------------------------------------------------------------------

site_map <- ggplot() +
  geom_sf(data = precip, aes(fill=`Rainfall (cm/yr)`), color = NA) + # precipitation layer
  geom_sf(data = oceans, fill = 'lightblue3') +                      # oceans layer
  geom_sf(data = state_map_data, fill = 'grey90', size = .4, alpha = .1) + # state map layer
  geom_sf(data = urban, fill = 'grey60') +                                 # urban area layer
  geom_sf(data = rivers_tx, color = 'deepskyblue1', alpha = .35) +         # surface waters layer
  geom_sf(data = cities, fill = 'black', shape =22, size = 3) +            # cities layer
  geom_label_repel(data = cities, aes(x=cities$longitude, y = cities$latitude, 
                                      label=cities$name), fill='grey80', alpha = .7) + # city labels
  geom_label_repel(data = d_site, aes(x = lon, y = lat, label = site_code),
                   fill = 'chartreuse', box.padding=.1, force=.2) +        # site labels
  coord_sf(xlim = c(-98.25, -95.75), ylim = c(27.25,29.25)) +              # limit lat and lon
  paletteer::scale_fill_paletteer_c("grDevices::Purple-Brown", 
                                    direction = -1,
                                    limits = c(50,125),
                                    breaks = c(50,75,100,125),
                                    name = "Rainfall \n(cm/yr)",
                                    guide = guide_colorbar(
                                      direction = "vertical",
                                      title.position = "top"))  +          # palette and legend settings
  theme_classic(base_size = 18)+                                           # theme
  ggsn::north(location = "topleft", scale = 0.8, symbol = 12,              
              x.min = -95.8, x.max = -95.6, y.min = 27.15, y.max = 27.3) + # add north arrow
  ggsn::scalebar(location = "bottomright", dist = 50,
                 dist_unit = "km", transform = TRUE, 
                 x.min=-97, x.max=-96, y.min=27.2, y.max=27.5,
                 st.bottom = FALSE, height = 0.025,
                 st.dist = 0.2, st.size = 5) +                             # add scale bar
  ylab(element_blank()) +                                                  # remove axes labels
  xlab(element_blank()) +                                                  # remove axes labels
  theme(legend.position = c(.9,.35),
        legend.background = element_rect(fill="lightblue2",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.key.size = unit(1, 'cm'),
        legend.title.align = .5,
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))                              # legend settings

#------------------------------------------------------------------------------
# Inset
#------------------------------------------------------------------------------

inset <- ggplot() +
  geom_sf(data = oceans, fill = 'lightblue3') +
  geom_sf(data = state_map_data, fill = 'white', size = 3, linewidth=1) +
  geom_rect(aes(xmin = -98.25, xmax = -95.75, ymin = 27.25, ymax = 29.25), 
            color = "red", fill = NA, linewidth=1.5) +
  geom_text(label='Texas', aes(x=-99, y = 31.5), size = 6) +
  xlim(c(-107, -93)) +
  ylim(c(26,37)) +
  labs(x = NULL, y = NULL) +
  theme_test() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "white"))


site_map
print(inset, vp = viewport(0.255, 0.85, width = 0.25, height = 0.25))

ggsave('exploration/visualization/site_map.png',
       plot = site_map,
       width = 9,
       height = 9,
       units = c("in"))

ggsave('exploration/visualization/site_map_inset.png',
       plot = inset,
       width = 2,
       height = 2,
       units = c("in"))

caption_site_map <- 'Study sites (in green) where fish, invertebrates, and environmental data were collected monthly following hurricane Harvey (August 2017) for 12 months. An overlay indicates the average annual precipitation (brown-purple) from USGS PRISM data (1981-2010). Cities (black squares) and urban areas (grey) were included for geographic reference. This map was made with Natural Earth.'

#------------------------------------------------------------------------------
# End site_map