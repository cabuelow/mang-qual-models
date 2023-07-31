# make forecasts of mangrove loss or gain using calibrated posterior hindcasts
# then ask what happens if there is a sustained increase in mangrove propagules via
# management or restoration? solve the matrices under each scenario with increased propagules

library(QPress)
library(tidyverse)
library(sf)
library(tmap)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
sf_use_s2(FALSE)
chosen_model_name <- 'mangrove_model'

# read in spatial data (mangrove typological units)
typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('outputs/master-dat.csv')

# import the final set of calibrated posterior hindcasts for each biophysical/pressure scenario 
# using optimal pressure definition and calibrated ambiguity threshold
press_thresh <- 4 # optimal pressure definition threshold
thresh <- 70 # optimal ambiguity threshold
naive_outcomes <- read.csv('outputs/validation/naive_outcomes.csv')
post_prob <- read.csv('outputs/validation/matrix-posterior-prob.csv')

# make posterior forecasts using naive hindcasts and posterior probabilities
#TODO: here rename sed_supp as 'fut_dams' to allow for future dams
spatial_pred <- spatial_dat %>% # here renaming future pressures as historical pressures so can join to posterior hindcasts
  dplyr::select(pressure_def, Type, fut_csqueeze, fut_csqueeze_1, sed_supp, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  relocate(press_gwsub, .after = press_csqueeze_1) %>% # relocate columns so in same order as historical pressures
  relocate(press_hist_drought, .after = press_gwsub) %>% 
  relocate(press_hist_ext_rain, .after = press_hist_drought) %>% 
  relocate(press_storms, .after = press_hist_ext_rain) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.') %>% 
  filter(pressure_def == press_thresh) %>% 
  left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
  left_join(post_prob, by = c('scenario', 'nsim')) %>% 
  mutate(LandwardMang = ifelse(LandwardMang == -1, 0, LandwardMang), # here turn losses into a 0 so just calculating the probability of gain/neutrality
         SeawardMang = ifelse(SeawardMang == -1, 0, SeawardMang)) %>% 
  mutate(LandwardMang = LandwardMang*matrix_post_prob,
         SeawardMang = SeawardMang*matrix_post_prob) %>% 
  group_by(pressure_def, Type, land_net_change_obs, sea_net_change_obs) %>% 
  summarise(LandwardMang = (sum(LandwardMang)/sum(matrix_post_prob))*100,
            SeawardMang = (sum(SeawardMang)/sum(matrix_post_prob))*100) %>% 
  mutate(Landward = case_when(is.na(LandwardMang) ~ NA,
                              LandwardMang >= thresh ~ 'Gain_neutrality',
                              LandwardMang < 100-thresh ~ 'Loss',
                              .default = 'Ambiguous'),
         Seaward = case_when(is.na(SeawardMang) ~ NA,
                             SeawardMang >= thresh ~ 'Gain_neutrality',
                             SeawardMang < 100-thresh ~ 'Loss',
                             .default = 'Ambiguous')) %>% 
  mutate(ambig_threshold = thresh)

# map final 'all data' forecasts

preds <- typ_points %>% 
  left_join(spatial_pred) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

# map forecasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Landward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Landward == 'Ambiguous' & !is.na(Landward))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Landward == 'Loss' & !is.na(Landward))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Landward == 'Gain_neutrality' & !is.na(Landward))) +
  tm_dots('Landward', 
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
  tm_shape(filter(preds, is.na(Seaward))) +
  tm_dots('darkgrey', size = 0.001) +
  #tm_shape(filter(preds, Seaward == 'Ambiguous' & !is.na(Seaward))) +
  #tm_dots('Seaward', 
   #       palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
    #      alpha = 0.5, 
     #     title = '',
      #    legend.show = F, 
       #   size = 0.0015) +
  tm_shape(filter(preds, Seaward == 'Loss' & !is.na(Seaward))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Gain_neutrality' & !is.na(Seaward))) +
  tm_dots('Seaward', 
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
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', chosen_model_name, '_all-data.png'), width = 5, height = 3)

# now take biophysical matrices and solve if mangrove propagules increase, plus future pressures
model <- models$mangrove_model # choose correct network model
chosen_model_name <- 'mangrove_model'

# find future biophysical and pressure scenarios

dat <- spatial_dat %>% # get unique biophysical/pressure combinations historically, given 5 different pressure definitions
  dplyr::select(pressure_def, Type, fut_csqueeze, fut_csqueeze_1, sed_supp, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>%
  rename(csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  relocate(press_gwsub, .after = press_csqueeze_1) %>% # relocate columns so in same order as historical pressures
  relocate(press_hist_drought, .after = press_gwsub) %>% 
  relocate(press_hist_ext_rain, .after = press_hist_drought) %>% 
  relocate(press_storms, .after = press_hist_ext_rain) %>% 
  select(-c(pressure_def,Type)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.') %>% 
  distinct()

# read in matrices  
tmp <- readRDS('outputs/simulation-outcomes/scenario_matrices.RDS')
matrices <- lapply(tmp, function(x){x[[2]]})
matrix_index <- data.frame(index = 1:length(matrices), scenario = unlist(lapply(tmp, function(x){names(x)[1]})))

# loop through each unique biophysical-pressure scenario and solve matrices in the relevant biophysical model
# to obtain naive forecasts of loss/gain

tmp2 <- list()
system.time( # takes 1.6 hours
  for(i in 1:nrow(dat)){
    bio_model <- matrices[[filter(matrix_index, scenario == as.character(dat[i,'scenario']))$index]]
    pressures <- data.frame(press = unlist(strsplit(as.character(dat[i,'press']), '\\.'))) %>% 
      mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                            'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
    tmp <- vector("list", dim(bio_model)[3])
    for(j in 1:dim(bio_model)[3]){
      tmp[[j]] <- data.frame(solver(bio_model[,,j], chosen_model, c('LandwardAvailableProp', 'SeawardAvailableProp', pressures$press)), nsim = j)  
    }
    tmp2[[i]] <- do.call(rbind, tmp) %>% 
      mutate(scenario = as.character(dat[i,'scenario']),
             press = as.character(dat[i,'press']))
  })
naive_outcomes_restore <- do.call(rbind, tmp2) %>% pivot_wider(names_from = 'node', values_from = 'outcome')
write.csv(naive_outcomes_restore, 'outputs/validation/naive_outcomes_restore.csv', row.names = F)
naive_outcomes_restore <- read.csv('outputs/validation/naive_outcomes_restore.csv')

# make posterior forecast for each future scenario, using matrix calibrated posterior probabilities

spatial_pred <- spatial_dat %>% # here renaming future pressures as historical pressures so can join to posterior hindcasts
  dplyr::select(pressure_def, Type, fut_csqueeze, fut_csqueeze_1, sed_supp, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  relocate(press_gwsub, .after = press_csqueeze_1) %>% # relocate columns so in same order as historical pressures
  relocate(press_hist_drought, .after = press_gwsub) %>% 
  relocate(press_hist_ext_rain, .after = press_hist_drought) %>% 
  relocate(press_storms, .after = press_hist_ext_rain) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.') %>% 
  filter(pressure_def == press_thresh) %>% 
  left_join(naive_outcomes_restore, by = c('scenario', 'press')) %>% 
  left_join(post_prob, by = c('scenario', 'nsim')) %>% 
  mutate(LandwardMang = ifelse(LandwardMang == -1, 0, LandwardMang), # here turn losses into a 0 so just calculating the probability of gain/neutrality
         SeawardMang = ifelse(SeawardMang == -1, 0, SeawardMang)) %>% 
  mutate(LandwardMang = LandwardMang*matrix_post_prob,
         SeawardMang = SeawardMang*matrix_post_prob) %>% 
  group_by(pressure_def, Type, land_net_change_obs, sea_net_change_obs) %>% 
  summarise(LandwardMang = (sum(LandwardMang)/sum(matrix_post_prob))*100,
            SeawardMang = (sum(SeawardMang)/sum(matrix_post_prob))*100) %>% 
  mutate(Landward = case_when(is.na(LandwardMang) ~ NA,
                              LandwardMang >= thresh ~ 'Gain_neutrality',
                              LandwardMang < 100-thresh ~ 'Loss',
                              .default = 'Ambiguous'),
         Seaward = case_when(is.na(SeawardMang) ~ NA,
                             SeawardMang >= thresh ~ 'Gain_neutrality',
                             SeawardMang < 100-thresh ~ 'Loss',
                             .default = 'Ambiguous')) %>% 
  mutate(ambig_threshold = thresh)

# map final 'all data' forecasts

preds <- typ_points %>% 
  left_join(spatial_pred) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map forecasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Landward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Landward == 'Ambiguous' & !is.na(Landward))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Landward == 'Loss' & !is.na(Landward))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Landward == 'Gain_neutrality' & !is.na(Landward))) +
  tm_dots('Landward', 
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
            main.title = 'B) Landward forecast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', chosen_model_name, '_all-data_restore.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Seaward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Ambiguous' & !is.na(Seaward))) +
  tm_dots('Seaward', 
         palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
        alpha = 0.5, 
       title = '',
      legend.show = F, 
     size = 0.0015) +
  tm_shape(filter(preds, Seaward == 'Loss' & !is.na(Seaward))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Gain_neutrality' & !is.na(Seaward))) +
  tm_dots('Seaward', 
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
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', chosen_model_name, '_all-data_restore.png'), width = 5, height = 3)


