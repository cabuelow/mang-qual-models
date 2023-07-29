# make forecasts of mangrove loss or gain using valid and highly likely matrices identified 
# during hindcast calibration and validation under different biophysical/pressure scenarios
# then ask what happens if there is a sustained increase in mangrove propagules via
# management or restoration? solve the valid and highly likely matrices given pressures + restoration

library(QPress)
library(tidyverse)
library(sf)
library(tmap)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
sf_use_s2(FALSE)

# read in spatial data (mangrove typological units)
typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  mutate(land_net_change_obs = ifelse(land_net_change == 'Gain', 1, -1),
         sea_net_change_obs = ifelse(sea_net_change == 'Gain', 1, -1)) %>% 
  mutate(csqueeze = recode(csqueeze, 'Medium' = 'M', 'High' = 'L', 'Low' = 'H'), # note counterintuitive notation here
         csqueeze_1 = ifelse(csqueeze == 'None', 0, 1), 
         fut_csqueeze = recode(fut_csqueeze, 'Medium' = 'M', 'High' = 'L', 'Low' = 'H'), # note counterintuitive notation here
         fut_csqueeze_1 = ifelse(fut_csqueeze == 'None', 0, 1), 
         sed_supp = recode(sed_supp, 'Medium' = 'M', 'Low' = 'L', 'High' = 'H'), 
         fut_dams = ifelse(fut_dams == 1, 'L', sed_supp), 
         prop_estab = recode(prop_estab, 'Medium' = 'M', 'High' = 'H', 'Low' = 'L'),
         Tidal_Class = recode(Tidal_Class, 'Micro' = 'H', 'Meso' = 'M', 'Macro' = 'L')) %>% 
  mutate(land_net_change = case_when(land_net_change == 'Gain' ~ 'Gain_neutrality',
                                     land_net_change == 'Neutral' ~ 'Gain_neutrality',
                                     .default = land_net_change)) %>% 
  mutate(sea_net_change = case_when(sea_net_change == 'Gain' ~ 'Gain_neutrality',
                                    sea_net_change == 'Neutral' ~ 'Gain_neutrality',
                                    .default = sea_net_change)) %>% 
  mutate(Geomorphology = sub("\\_.*", "", Type))

# import the final set of calibrated posterior predictions for each biophysical/pressure scenario 
# using optimal pressure definition and calibrated ambiguity threshold
press <- 4 # optimal pressure definition threshold
thresh <- 85 # optimal ambiguity threshold
preds <- read.csv(paste0('outputs/predictions/final-calibrated-predictions_', press, '_', thresh, '.csv'))

spatial_pred <- spatial_dat %>% 
  dplyr::select(pressure_def, Type, fut_csqueeze, fut_csqueeze_1, sed_supp, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, hist_gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T) %>% 
  filter(pressure_def == press) %>% 
  left_join(preds, by = 'scenario')

# map final 'all data' hindcasts

preds <- typ_points %>% 
  left_join(spatial_pred) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map hindcasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Landward == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Landward == 'Loss' & !is.na(Change))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  #tm_shape(filter(preds, Landward == 'Gain' & !is.na(Change))) +
  #tm_dots('Landward', 
   #       palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
    #      alpha = 0.5, 
     #     title = '',
      #    legend.show = F, 
       #   size = 0.025) +
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
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', chosen_model_name, '_all-data.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  #tm_shape(filter(preds, Seaward == 'Ambiguous' & !is.na(Change))) +
  #tm_dots('Seaward', 
   #       palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
    #      alpha = 0.5, 
     #     title = '',
      #    legend.show = F, 
       #   size = 0.0015) +
  tm_shape(filter(preds, Seaward == 'Loss' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  #tm_shape(filter(preds, Seaward == 'Gain' & !is.na(Change))) +
  #tm_dots('Seaward', 
   #       palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
    #      alpha = 0.5, 
     #     title = '',
      #    legend.show = F, 
       #   size = 0.025) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'B) Seaward forecast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', chosen_model_name, '_all-data.png'), width = 5, height = 3)

