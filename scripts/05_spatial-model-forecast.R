# make forecasts of mangrove loss or gain using calibrated posterior hindcasts
# then ask what happens if there is a sustained increase in mangrove propagules via
# management or restoration? solve the matrices under each scenario with increased propagules

library(QPress)
library(tidyverse)
library(sf)
library(tmap)
library(foreach)
library(doParallel)
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
go <- 1 # which coastal dev threshold?
press <- 3 # which pressure definition threshold?
thresh <- 75 # which ambiguity threshold?
rm_e <- 'N' # remove erosion from validation? Y or N
naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes_', go,'_', rm_e,'.csv'))
post_prob <- read.csv(paste0('outputs/validation/matrix-posterior-prob', go, '_', rm_e, '_', press, '_', thresh, '.csv'))

# make posterior forecasts using naive hindcasts and posterior probabilities

spatial_pred <- spatial_dat %>% # here renaming future pressures as historical pressures so can join to posterior hindcasts
  filter(Cdev_thresh == go & pressure_def == press) %>% 
  dplyr::select(pressure_def, Type, fut_dams, fut_csqueeze, fut_csqueeze_1, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(sed_supp = fut_dams, csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab))
# reorder column names so always in same order regardless of filtering
spatial_pred <- spatial_pred[,c('pressure_def', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 
                          'land_net_change_obs', 'sea_net_change_obs',  'vals_gwsub', 'vals_hist_drought',
                          'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
                          'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                          'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
spatial_pred <- spatial_pred %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  mutate(press = ifelse(!is.na(press_no_press), 'none', press)) %>% 
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
write.csv(spatial_pred, paste0('outputs/predictions/forecast-predictions', go, '_', rm_e, '_', press, '_', thresh, '.csv'), row.names = F)

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
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', go, '_', rm_e, '_', press, '_', thresh,'_all-data.png'), width = 5, height = 1)

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
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', go, '_', rm_e, '_', press, '_', thresh, '_all-data.png'), width = 5, height = 1)

# now take biophysical matrices and solve if mangrove propagules increase, plus future pressures
chosen_model <- models$mangrove_model # choose correct network model

# find future biophysical and pressure scenarios

dat <- spatial_dat %>% # get unique biophysical/pressure combinations historically, given 5 different pressure definitions
  filter(Cdev_thresh == go & pressure_def == press) %>% 
  dplyr::select(pressure_def, Type, fut_csqueeze, fut_csqueeze_1, fut_dams, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>%
  rename(sed_supp = fut_dams, csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  select(-c(pressure_def, Type))
# reorder column names so always in same order regardless of filtering
dat <- dat[,c('csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'vals_gwsub', 'vals_hist_drought',
                                'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
                                'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                                'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
dat <- dat %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  distinct()

# read in matrices  
tmp <- readRDS(paste0('outputs/simulation-outcomes/scenario_matrices_', go,'_', rm_e, '.RDS'))
matrices <- lapply(tmp, function(x){x[[2]]})
matrix_index <- data.frame(index = 1:length(matrices), scenario = unlist(lapply(tmp, function(x){names(x)[1]})))

# loop through each unique biophysical-pressure scenario and solve matrices in the relevant biophysical model
# to obtain naive forecasts of loss/gain given management/conservation scenarios

dat2 <- dat %>% mutate(split = rep(1:5, (nrow(.)+4)/5)[1:nrow(.)]) # add a column to split so can parallelise
scenarios <- list(c('LandwardAvailableProp', 'SeawardAvailableProp'), 'Sediment', 'Coastalsqueeze')
system.time(
for(b in seq_along(scenarios)){
  scn <- scenarios[[b]]
  if(scn[1] != 'Coastalsqueeze'){
cl <- makeCluster(5)
registerDoParallel(cl)
results <- foreach(h = seq_along(unique(dat2$split)), .packages = c('tidyverse', 'QPress')) %dopar% {
datsub <- dat2 %>% filter(split == h)
tmp2 <- list()
  for(i in 1:nrow(datsub)){
    bio_model <- matrices[[filter(matrix_index, scenario == as.character(datsub[i,'scenario']))$index]]
    pressures <- data.frame(press = unlist(strsplit(as.character(datsub[i,'press']), '\\.'))) %>% 
      mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                            'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
    tmp <- vector("list", dim(bio_model)[3])
    for(j in 1:dim(bio_model)[3]){
      tmp[[j]] <- data.frame(solver(bio_model[,,j], chosen_model, c(scn, pressures$press)), nsim = j)  
    }
    tmp2[[i]] <- do.call(rbind, tmp) %>% mutate(scenario = as.character(datsub[i,'scenario']), press = as.character(datsub[i,'press']))
  }
do.call(rbind, tmp2)
}
stopCluster(cl)
results <- do.call(rbind, results) %>% pivot_wider(names_from = 'node', values_from = 'outcome')
write.csv(results, paste0('outputs/predictions/naive_outcomes_scenario_', scn[1], '_', go, '_', rm_e, '_', press, '_', thresh, '.csv'), row.names = F)
  }else{
    cl <- makeCluster(5)
    registerDoParallel(cl)
    results <- foreach(h = seq_along(unique(dat2$split)), .packages = c('tidyverse', 'QPress')) %dopar% {
      datsub <- dat2 %>% filter(split == h)
      tmp2 <- list()
      for(i in 1:nrow(datsub)){
        bio_model <- matrices[[filter(matrix_index, scenario == as.character(datsub[i,'scenario']))$index]]
        pressures <- data.frame(press = unlist(strsplit(as.character(datsub[i,'press']), '\\.'))) %>% 
          mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                                'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
        tmp <- vector("list", dim(bio_model)[3])
        for(j in 1:dim(bio_model)[3]){
          new_mat <- bio_model[,,j]
          new_mat[8,9] <- 1 # update so interaction between sea level rise and landward mangrove is always high (i.e., 1)
          tmp[[j]] <- data.frame(solver(bio_model[,,j], chosen_model, pressures$press), nsim = j)  
        }
        tmp2[[i]] <- do.call(rbind, tmp) %>% mutate(scenario = as.character(datsub[i,'scenario']), press = as.character(datsub[i,'press']))
      }
      do.call(rbind, tmp2)
    }
    stopCluster(cl)
    results <- do.call(rbind, results) %>% pivot_wider(names_from = 'node', values_from = 'outcome')
    write.csv(results, paste0('outputs/predictions/naive_outcomes_scenario_', scn[1], '_', go, '_', rm_e, '_', press, '_', thresh, '.csv'), row.names = F)
}})

# make posterior forecast for each future scenario, using matrix calibrated posterior probabilities

spatial_pred <- spatial_dat %>% # here renaming future pressures as historical pressures so can join to posterior hindcasts
  filter(Cdev_thresh == go & pressure_def == press) %>% 
  dplyr::select(pressure_def, Type, fut_dams, fut_csqueeze, fut_csqueeze_1, Tidal_Class, prop_estab, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(sed_supp = fut_dams, csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, ant_slr = fut_slr, gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>%
  mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab))
# reorder column names so always in same order regardless of filtering
spatial_pred <- spatial_pred[,c('pressure_def', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 
                                  'land_net_change_obs', 'sea_net_change_obs',  'vals_gwsub', 'vals_hist_drought',
                                  'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 'vals_csqueeze_1','press_no_press', 'press_gwsub',
                                  'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                                  'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
spatial_pred <- spatial_pred %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  mutate(press = ifelse(!is.na(press_no_press), 'none', press))

# identify the unique 'no pressure' bio scenario combinations and add in 1000 sims as neutral outcomes to naive hindcasts
nsim <- 1000
d <- select(filter(spatial_pred, press == 'none'), press, scenario) %>% distinct()
naive_outcomes_restore <- rbind(naive_outcomes_restore, data.frame(nsim = rep(1:nsim, nrow(d)), 
                                                   scenario = rep(d$scenario, each = nsim), 
                                                   press = rep(d$press, each = nsim),
                                                   LandwardMang = 1, # assuming no pressure means neutral
                                                   SeawardMang = 1))
write.csv(naive_outcomes_restore, paste0('outputs/validation/naive_outcomes_restore_', go, '_', rm_e, '_', press, '_', thresh, '.csv'), row.names = F)
naive_outcomes_restore <- read.csv(paste0('outputs/validation/naive_outcomes_restore_', go, '_', rm_e, '_', press, '_', thresh, '.csv'))

spatial_pred <- spatial_pred %>% 
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
write.csv(spatial_pred, paste0('outputs/predictions/forecast-predictions', go, '_', rm_e, '_', press, '_', thresh, '_restore.csv'), row.names = F)

# map final 'all data' forecasts

preds <- typ_points %>% 
  left_join(spatial_pred) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map forecasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  #tm_shape(filter(preds, is.na(Landward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Landward == 'Gain_neutrality' & !is.na(Landward))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
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
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', go, '_', rm_e, '_', press, '_', thresh, '_all-data_restore.png'), width = 5, height = 1)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  #tm_shape(filter(preds, is.na(Seaward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Gain_neutrality' & !is.na(Seaward))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.025) +
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
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', go, '_', rm_e, '_', press, '_', thresh, '_all-data_restore.png'), width = 5, height = 1)


