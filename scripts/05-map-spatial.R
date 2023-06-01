# map spatial model outcomes

library(sf)
library(tmap)
library(tidyverse)
library(scales)
sf_use_s2(FALSE)
tmap_mode('plot')

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")

# which model outcomes to map?

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv'))

# join outcomes to spatial typologies

landward_forecast <- typ_points %>% 
  left_join(filter(dat, var == 'LandwardMang' & cast == 'forecast'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

seaward_forecast <- typ_points %>% 
  left_join(filter(dat, var == 'SeawardMang' & cast == 'forecast'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

landward_hindcast <- typ_points %>% 
  left_join(filter(dat, var == 'LandwardMang' & cast == 'hindcast'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

seaward_hindcast <- typ_points %>% 
  left_join(filter(dat, var == 'SeawardMang' & cast == 'hindcast'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

world_mang <- st_crop(World, xmin = -150, ymin = -40, xmax = 180, ymax = 33)

# map forecasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(landward_forecast, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(landward_forecast, !is.na(Prob_change))) +
  tm_dots('Prob_gain_neutral', 
          palette = 'Spectral',
          breaks = c(0,25,50,75,100),
          labels = c("-100 to -75", "-75 to -50", "50 to 75", '75 to 100'),
          title = '',
          legend.is.portrait = T, 
          alpha = 0.5) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'A) Landward mangrove forecast',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
lmap

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward_forecast, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(seaward_forecast, !is.na(Prob_change))) +
  tm_dots('Prob_gain_neutral', 
          palette = 'Spectral',
          breaks = c(0,25,50,75,100),
          labels = c("-100 to -75", "-75 to -50", "50 to 75", '75 to 100'),
          title = '',
          legend.is.portrait = T,
          alpha = 0.5) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'B) Seaward mangrove forecast',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
smap

maps <- tmap_arrange(lmap, smap, ncol = 1)
maps

tmap_save(maps, paste0('outputs/maps/forecast_map_', chosen_model_name, '.png'), width = 10, height = 5)

# map hindcasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  #tm_shape(filter(landward_hindcast, is.na(Prob_change))) +
  #tm_dots('darkgrey') +
  tm_shape(filter(landward_hindcast, !is.na(Prob_change))) +
  tm_dots('Prob_gain_neutral', 
          palette = 'Spectral',
          breaks = c(0,25,50,75,100),
          labels = c("-100 to -75", "-75 to -50", "50 to 75", '75 to 100'),
          title = '',
          legend.is.portrait = T, 
          alpha = 0.5) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'A) Landward mangrove hindcast',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
lmap

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward_forecast, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(seaward_forecast, !is.na(Prob_change))) +
  tm_dots('Prob_gain_neutral', 
          palette = 'Spectral',
          breaks = c(0,25,50,75,100),
          labels = c("-100 to -75", "-75 to -50", "50 to 75", '75 to 100'),
          title = '',
          legend.is.portrait = T,
          alpha = 0.5) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'B) Seaward mangrove hindcast',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
smap

maps <- tmap_arrange(lmap, smap, ncol = 1)
maps

tmap_save(maps, paste0('outputs/maps/hindcast_map_', chosen_model_name, '.png'), width = 10, height = 5)




