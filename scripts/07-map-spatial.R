# map spatial model outcomes

library(sf)
library(tmap)
library(tidyverse)
library(scales)
library(cowplot)
library(ggplotify)
library(tmaptools)
library(patchwork)
source('scripts/helpers/models_v2.R')
sf_use_s2(FALSE)

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  filter(pressure_def == thresh_press)

# set pressure and ambiguity thresholds for mapping
accuracy <- read.csv('outputs/validation/accuracy-threshold-vary.csv')
thresh_sea <- filter(accuracy, mangrove == 'Seaward', Overall_accuracy == max(filter(accuracy, mangrove == 'Seaward')$Overall_accuracy))[1,'ambig_threshold']
thresh_land <- filter(accuracy, mangrove == 'Landward', Overall_accuracy == max(filter(accuracy, mangrove == 'Landward')$Overall_accuracy))[1,'ambig_threshold']
thresh_press <- filter(accuracy, mangrove == 'Landward', Overall_accuracy == max(filter(accuracy, mangrove == 'Landward')$Overall_accuracy))[1,'pressure_def']

# which model outcomes to map?
names(models) # names of available models
chosen_model_name <- 'mangrove_model'

# uncalibrated hind/forecasts
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  filter(pressure_def == thresh_press) %>% 
  mutate(Prob_gain_neutrality = Prob_gain + Prob_neutral)

# calibrated forecasts
datcal <- read.csv(paste0('outputs/validation/calibrated_forecast', chosen_model_name, '_spatial.csv')) %>% 
  mutate(Prob_gain_neutrality = Prob_gain + Prob_neutral)
  
# classify outcomes according to ambiguity threshold
  
land <- dat %>% 
  filter(var == 'LandwardMang') %>% 
  mutate(Change = case_when(Prob_gain_neutrality > thresh_land ~ 'Gain_neutrality',
                            Prob_loss < -thresh_land ~ 'Loss',
                            .default = 'Ambiguous')) %>% 
  inner_join(select(spatial_dat, Type, pressure_def, land_net_change), by = c('Type', 'pressure_def'))

sea <- dat %>% 
  filter(var == 'SeawardMang') %>% 
  mutate(Change = case_when(Prob_gain_neutrality > thresh_sea ~ 'Gain_neutrality',
                            Prob_loss < -thresh_sea ~ 'Loss',
                            .default = 'Ambiguous')) %>% 
  inner_join(select(spatial_dat, Type, pressure_def, sea_net_change), by = c('Type', 'pressure_def'))

landcal <- datcal %>% 
  filter(var == 'LandwardMang') %>% 
  mutate(Change = case_when(Prob_gain_neutrality > thresh_land ~ 'Gain_neutrality',
                            Prob_loss < -thresh_land ~ 'Loss',
                            .default = 'Ambiguous')) %>% 
  inner_join(select(spatial_dat, Type, pressure_def, land_net_change), by = c('Type'))

seacal <- datcal %>% 
  filter(var == 'SeawardMang') %>% 
  mutate(Change = case_when(Prob_gain_neutrality > thresh_sea ~ 'Gain_neutrality',
                            Prob_loss < -thresh_sea ~ 'Loss',
                            .default = 'Ambiguous')) %>% 
  inner_join(select(spatial_dat, Type, pressure_def, sea_net_change), by = c('Type'))

# summarise the number of units hindcast/forecast to be each change class

landsummary <- land %>%
  right_join(select(spatial_dat, Type)) %>% 
  mutate(Change = ifelse(is.na(Change), 'No pressures', Change)) %>% 
  mutate(num = 1) %>% 
  group_by(cast, Change) %>% 
  summarise(percent_total = (sum(num)/length(unique(spatial_dat$Type)))*100,
            number_units = sum(num))
sum(landsummary$percent_total) # should equal 100
write.csv(landsummary, 'outputs/summary-stats/land.csv', row.names = F)

seasummary <- sea %>%
  right_join(select(spatial_dat, Type)) %>% 
  mutate(Change = ifelse(is.na(Change), 'No pressures', Change)) %>% 
  mutate(num = 1) %>% 
  group_by(cast, Change) %>% 
  summarise(percent_total = (sum(num)/length(unique(spatial_dat$Type)))*100,
            number_units = sum(num))
sum(seasummary$percent_total) # should equal 100
write.csv(seasummary, 'outputs/summary-stats/sea.csv', row.names = F)

# join to spatial typologies for plotting

landward_hindcast <- typ_points %>% 
  left_join(filter(land, cast == 'hindcast'), by = 'Type') %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

seaward_hindcast <- typ_points %>% 
  left_join(filter(sea, cast == 'hindcast'), by = 'Type') %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map hindcasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(landward_hindcast, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(landward_hindcast, Change == 'Gain_neutrality' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
  tm_shape(filter(landward_hindcast, Change == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(landward_hindcast, Change == 'Loss' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'B) Landward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_', chosen_model_name, '.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward_hindcast, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(seaward_hindcast, Change == 'Loss' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gai_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(seaward_hindcast, Change == 'Gain_neutrality' & !is.na(Change))) +
  tm_dots('Change', 
         palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
        alpha = 0.5, 
       title = '',
      legend.show = F, 
     size = 0.025) +
  tm_shape(filter(seaward_hindcast, Change == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Change', 
         palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
        alpha = 0.5, 
     title = '',
    legend.show = F, 
   size = 0.025) +
tm_layout(legend.outside = F,
          #legend.outside.position = 'bottom',
          legend.position = c(0.13, 0.01),
          title.position = c(0.01,0.45),
          legend.title.size = 0.45,
          legend.text.size = 0.35,
          main.title = 'A) Seaward hindcast',
          main.title.size = 0.45,
          frame = T,
          legend.bg.color = 'white',
          legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', chosen_model_name, '.png'), width = 5, height = 3)

maps <- tmap_arrange(smap, lmap, ncol = 1)
maps

tmap_save(maps, paste0('outputs/maps/hindcast_map_', chosen_model_name, '.png'), width = 6, height = 3)

# map forecasts

landward_forecast <- typ_points %>% 
  left_join(filter(land, cast == 'forecast'), by = 'Type') %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

seaward_forecast <- typ_points %>% 
  left_join(filter(sea, cast == 'forecast'), by = 'Type') %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

landward_forecast_cal <- typ_points %>% 
  left_join(landcal, by = 'Type') %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

seaward_forecast_cal <- typ_points %>% 
  left_join(seacal, by = 'Type') %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map forecasts - uncalibrated

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(landward_forecast, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(landward_forecast, Change == 'Gain_neutrality' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
  tm_shape(filter(landward_forecast, Change == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(landward_forecast, Change == 'Loss' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'B) Landward forecast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', chosen_model_name, '.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward_forecast, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(seaward_forecast, Change == 'Loss' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  #tm_shape(filter(seaward_forecast, Sea_Change == 'Gain' & !is.na(Sea_Change))) +
  #tm_dots('Sea_Change', 
   #       palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
    #      alpha = 0.5, 
     #     title = '',
      #    legend.show = F, 
       #   size = 0.025) +
  tm_shape(filter(seaward_forecast, Change == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'A) Seaward forecast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', chosen_model_name, '.png'), width = 5, height = 3)

maps <- tmap_arrange(smap, lmap, ncol = 1)
maps

tmap_save(maps, paste0('outputs/maps/forecast_map_', chosen_model_name, '.png'), width = 6, height = 3)

# map forecasts - calibrated

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(landward_forecast_cal, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(landward_forecast_cal, Change == 'Gain_neutrality' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
  tm_shape(filter(landward_forecast_cal, Change == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(landward_forecast_cal, Change == 'Loss' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'B) Landward forecast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', chosen_model_name, '_calibrated.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward_forecast_cal, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(seaward_forecast, Change == 'Loss' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  #tm_shape(filter(seaward_forecast_cal, Sea_Change == 'Gain' & !is.na(Sea_Change))) +
  #tm_dots('Sea_Change', 
  #       palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
  #      alpha = 0.5, 
  #     title = '',
  #    legend.show = F, 
  #   size = 0.025) +
  tm_shape(filter(seaward_forecast_cal, Change == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Change', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'A) Seaward forecast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', chosen_model_name, '_calibrated.png'), width = 5, height = 3)

maps <- tmap_arrange(smap, lmap, ncol = 1)
maps

tmap_save(maps, paste0('outputs/maps/forecast_map_', chosen_model_name, '_calibrated.png'), width = 6, height = 3)

# summarise characteristics of typologies with gains, losses, or ambiguity

# land
land_sum <- spatial_dat %>% 
  filter(pressure_def == thresh_press) %>% 
  select(Type, fut_csqueeze, sed_supp, fut_slr, fut_dams, fut_gwsub, fut_storms, fut_drought, fut_ext_rain, 
         Tidal_Class, prop_estab) %>% 
  rowwise() %>% dplyr::mutate(geomorph = strsplit(Type, split="_")[[1]][1]) %>% 
  pivot_wider(names_from = 'geomorph', values_from = 'geomorph') %>% 
  mutate_at(vars(Delta:OpenCoast), ~ifelse(is.na(.), 0, 1)) %>% 
  inner_join(filter(land, cast == 'forecast' & !is.na(Change)), by = 'Type') %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'Low', 1, .)) %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'Medium', 2, .)) %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'High', 3, .)) %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'None', 0, .)) %>% 
  mutate(Tidal_Class = ifelse(Tidal_Class == 'Micro', 1, Tidal_Class)) %>% 
  mutate(Tidal_Class = ifelse(Tidal_Class == 'Meso', 2, Tidal_Class)) %>% 
  mutate(Tidal_Class = ifelse(Tidal_Class == 'Macro', 3, Tidal_Class)) %>% 
  mutate_at(vars(fut_csqueeze:OpenCoast), ~rescale(as.integer(.), c(0,1))) %>% 
  pivot_longer(fut_csqueeze:OpenCoast, names_to = 'variable', values_to = 'val') %>% 
  group_by(Change, variable) %>% 
  summarise(val = mean(val),
            total_change_cat = n(),
            percent_change_cat = (n()/nrow(filter(land, !is.na(Change))))*100) %>% 
  mutate(variable = recode(variable, 'Tidal_Class' = 'Tidal range',
                           'sed_supp' = 'Sediment supply',
                           'prop_estab' = 'Propagule establishment',
                           'OpenCoast' = 'Open coast',
                           'fut_storms' = 'Intense storms',
                           'fut_slr' = 'Sea-level rise',
                           'fut_dams' = 'Dams',
                           'fut_gwsub' = 'Subsidence',
                           'fut_ext_rain' = 'Extreme rainfall',
                           'fut_drought' = 'Drought', 
                           'fut_csqueeze' = 'Coastal squeeze')) %>% 
  mutate(group = ifelse(variable %in% c('Sea-level rise', 'Intense storms', 'Extreme rainfall',
                                        'Drought'), 'C) Climate', NA),
         group = ifelse(variable %in% c('Subsidence',  'Dams',
                                        'Coastal squeeze','Sediment supply'), 'D) Anthropogenic', group),
         group = ifelse(variable %in% c('Tidal range', 'Propagule establishment',
                                        'Delta', 'Lagoon', 'Open coast', 'Estuary'), 'E) Biophysical', group)) %>% 
  mutate(variable = factor(variable, levels = c('Sea-level rise', 'Intense storms', 'Extreme rainfall',
                                                'Drought',
                                                'Subsidence',  'Dams',
                                                'Coastal squeeze','Sediment supply', 'Tidal range', 'Propagule establishment',
                                                'Delta', 'Lagoon', 'Open coast', 'Estuary'))) %>% 
  mutate(group = factor(group, levels = c('C) Climate', 'D) Anthropogenic', 'E) Biophysical')))

a <- ggplot(land_sum) +
  geom_tile(aes(x = Change, y = variable, fill = val) ) +
  scale_fill_distiller(palette = 'PuBuGn', direction = 1,  name = '') +
  facet_wrap(~group, scales = 'free') +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 1)) +
  xlab('') +
  ylab('') +
  #ggtitle('B) Landward') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, hjust = 0.9))
a
ggsave('outputs/maps/landward-map-characteristics.png', width = 8, height = 2)

#sea
sea_sum <- spatial_dat %>% 
  filter(sensitivity == 2) %>% 
  select(Type, fut_csqueeze, sed_supp, fut_slr, fut_dams, fut_gwsub, fut_storms, fut_drought, fut_ext_rain, 
         Tidal_Class, prop_estab) %>% 
  rowwise() %>% dplyr::mutate(geomorph = strsplit(Type, split="_")[[1]][1]) %>% 
  pivot_wider(names_from = 'geomorph', values_from = 'geomorph') %>% 
  mutate_at(vars(Delta:OpenCoast), ~ifelse(is.na(.), 0, 1)) %>% 
  inner_join(filter(sea, !is.na(Sea_Change)), by = 'Type') %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'Low', 1, .)) %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'Medium', 2, .)) %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'High', 3, .)) %>% 
  mutate_at(vars(fut_csqueeze:sed_supp, prop_estab), ~ifelse(. == 'None', 0, .)) %>% 
  mutate(Tidal_Class = ifelse(Tidal_Class == 'Micro', 1, Tidal_Class)) %>% 
  mutate(Tidal_Class = ifelse(Tidal_Class == 'Meso', 2, Tidal_Class)) %>% 
  mutate(Tidal_Class = ifelse(Tidal_Class == 'Macro', 3, Tidal_Class)) %>%  
  mutate_at(vars(fut_csqueeze:OpenCoast), ~rescale(as.integer(.), c(0,1))) %>% 
  pivot_longer(fut_csqueeze:OpenCoast, names_to = 'variable', values_to = 'val') %>% 
  group_by(Sea_Change, variable) %>% 
  summarise(val = mean(val),
            total_change_cat = n(),
            percent_change_cat = (n()/nrow(filter(sea, !is.na(Sea_Change))))*100) %>% 
  mutate(variable = recode(variable, 'Tidal_Class' = 'Tidal range',
                           'sed_supp' = 'Sediment supply',
                           'prop_estab' = 'Propagule establishment',
                           'OpenCoast' = 'Open coast',
                           'fut_storms' = 'Intense storms',
                           'fut_slr' = 'Sea-level rise',
                           'fut_dams' = 'Dams',
                           'fut_gwsub' = 'Subsidence',
                           'fut_ext_rain' = 'Extreme rainfall',
                           'fut_drought' = 'Drought', 
                           'fut_csqueeze' = 'Coastal squeeze')) %>% 
  mutate(group = ifelse(variable %in% c('Sea-level rise', 'Intense storms', 'Extreme rainfall',
                                        'Drought'), 'Climate', NA),
         group = ifelse(variable %in% c('Subsidence',  'Dams',
                                        'Coastal squeeze','Sediment supply'), 'Anthropogenic', group),
         group = ifelse(variable %in% c('Tidal range', 'Propagule establishment',
                                        'Delta', 'Lagoon', 'Open coast', 'Estuary'), 'Biophysical', group)) %>% 
  mutate(variable = factor(variable, levels = c('Sea-level rise', 'Intense storms', 'Extreme rainfall',
                                                'Drought',
                                                'Subsidence',  'Dams',
                                                'Coastal squeeze','Sediment supply', 'Tidal range', 'Propagule establishment',
                                                'Delta', 'Lagoon', 'Open coast', 'Estuary'))) %>% 
  mutate(group = factor(group, levels = c('Climate', 'Anthropogenic', 'Biophysical')))

b <- ggplot(sea_sum) +
  geom_tile(aes(x = Sea_Change, y = variable, fill = val)) +
  scale_fill_distiller(palette = 'Greens', direction = 1,  name = '') +
  facet_wrap(~group, scales = 'free') +
  #facet_grid(rows = vars(group), scales = 'free', space = 'free_y') +
  xlab('') +
  ylab('') +
  #ggtitle('A) Seaward') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.9))
b
ggsave('outputs/maps/seaward-map-characteristics.png', width = 6, height = 2)

c <- b/a
c

ggsave('outputs/maps/map-characteristics.png', width = 5.3, height = 7.5)


