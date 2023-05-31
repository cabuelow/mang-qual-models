# map spatial model outcomes

library(sf)
library(tmap)
library(tidyverse)
sf_use_s2(FALSE)

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")

# which model outcomes to map?

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
allout <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv'))

# join outcomes to spatial typologies

landward <- typ_points %>% 
  left_join(filter(allout, var == 'LandwardMang'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

seaward <- typ_points %>% 
  left_join(filter(allout, var == 'SeawardMang'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

world_mang <- st_crop(World, xmin = -150, ymin = -40, xmax = 180, ymax = 33)

# map

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(landward, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(landward, !is.na(Prob_change))) +
  tm_dots('Prob_change', 
          palette = 'Spectral',
          breaks = c(-100,-50,0,50,100),
          midpoint = 0,
          title = '',
          legend.is.portrait = T) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'A) Landward mangrove',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
lmap

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(seaward, !is.na(Prob_change))) +
  tm_dots('Prob_change', 
          palette = 'Spectral',
          breaks = c(-100,-50,0,50,100),
          midpoint = 0,
          title = '',
          legend.is.portrait = T) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'B) Seaward mangrove',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
smap

maps <- tmap_arrange(lmap, smap, ncol = 1)

tmap_save(maps, 'outputs/map_change.png', width = 10, height = 5)





