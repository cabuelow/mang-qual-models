# miscellaneous plotting/mapping

library(tidyverse)
library(sf)
library(tmap)
sf_use_s2(FALSE)

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('outputs/master-dat.csv')

typ_points <- typ_points %>% 
  left_join(filter(spatial_dat, pressure_def == 1)) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

# historical observations of net change

map <- tm_shape(world_mang) +
  tm_fill(col = 'gray85') +
  tm_shape(typ_points) + 
  tm_dots('land_net_change', title = '',  palette = c('Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4')) +
  tm_layout(frame = T, legend.position = c(0.2, 0))
map
tmap_save(map, 'outputs/maps/historical-landward-net-change.png', width = 8, height = 2)

map <- tm_shape(world_mang) +
  tm_fill(col = 'gray85') +
  tm_shape(typ_points) + 
  tm_dots('sea_net_change', title = '',  palette = c('Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4')) +
  tm_layout(frame = T, legend.position = c(0.2, 0))
map
tmap_save(map, 'outputs/maps/historical-seaward-net-change.png', width = 8, height = 2)
