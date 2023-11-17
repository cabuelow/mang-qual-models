# make forecasts of mangrove loss or gain using matrix posterior probabilities for fitting

library(QPress)
library(tidyverse)
library(sf)
library(tmap)
library(foreach)
library(doParallel)
library(RColorBrewer)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
set.seed(123) # set random number generator so reproducible
sf_use_s2(FALSE)
chosen_model_name <- 'mangrove_model'

# read in spatial data (mangrove typological units)
typ_points <- st_read('data/typologies/Mangrove_Typology_v3.14_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('data/master-dat.csv')
pal <- brewer.pal(11, 'Spectral') # color palette

# import the final set of calibrated posterior hindcasts for each biophysical/pressure scenario 
# using optimal pressure definition and calibrated ambiguity threshold
press <- 4 # which pressure definition threshold?
thresh <- 75 # which ambiguity threshold?
naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes.csv'))
post_prob <- read.csv(paste0('outputs/validation/matrix-posterior-prob_', press, '_', thresh, '.csv'))

# make posterior forecasts using naive hindcasts and posterior probabilities, without future sea level rise

spatial_pred <- spatial_dat %>% # here renaming future pressures as historical pressures so can join to posterior hindcasts
  filter(pressure_def == press) %>% 
  dplyr::select(pressure_def, Type, fut_dams, fut_csqueeze, fut_csqueeze_1, Tidal_Class, prop_estab, climate, fut_cdev, #fut_slr, 
                fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(sed_supp = fut_dams, csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, cdev = fut_cdev, #ant_slr = fut_slr, 
         gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  mutate(no_press = gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  #mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  pivot_longer(cols = c(csqueeze_1,gwsub:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  #pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab),
         climate_2 = paste0('climate_', .$climate),
         cdev_2 = paste0('cdev_', .$cdev))
# reorder column names so always in same order regardless of filtering
spatial_pred <- spatial_pred[,c('pressure_def', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'climate', 'cdev',
                          'land_net_change_obs', 'sea_net_change_obs',  'vals_gwsub', 'vals_hist_drought',
                          'vals_hist_ext_rain', 'vals_storms', #'vals_ant_slr', 
                          'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
                          'press_hist_drought', 'press_hist_ext_rain', 'press_storms', #'press_ant_slr',
                          'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2', 'climate_2', 'cdev_2')]
spatial_pred <- spatial_pred %>% 
  unite('scenario', csqueeze_2:cdev_2, na.rm = T, sep = '.') %>% 
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
write.csv(spatial_pred, paste0('outputs/predictions/forecast-predictions_', press, '_', thresh, '.csv'), row.names = F)
#spatial_pred <- read.csv(paste0('outputs/predictions/forecast-predictions_', press, '_', thresh, '.csv'))

# summarise predictions
datsum <- ungroup(spatial_pred) %>% 
  mutate(Landward = paste0('Landward_', .$Landward),
         Seaward = paste0('Seaward_', .$Seaward)) %>% 
  mutate(Landward_seaward = paste0(.$Landward, '.', .$Seaward)) %>% 
  group_by(Landward_seaward) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
write.csv(datsum, paste0('outputs/predictions/forecast-predictions_', press, '_', thresh, '_summary.csv'), row.names = F)    

# map final 'all data' forecasts

preds <- typ_points %>% 
  left_join(spatial_pred) %>%
  #mutate_at(vars(LandwardMang, SeawardMang), ~ifelse(. >=100-thresh & . < thresh, 100-thresh, .)) %>% 
  mutate_at(vars(LandwardMang, SeawardMang), ~./100) %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

# map forecasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(preds, is.na(Landward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, !is.na(Landward))) +
  tm_bubbles('LandwardMang', 
          palette = pal[1:10],
          midpoint = 0.5,
          breaks = seq(0,1,0.1),
          size = 'LandwardMang',
          scale = 0.25,
          alpha = 0.5, 
          border.alpha = 0,
          legend.size.show = F,
          legend.col.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = 'B) Landward baseline forecast',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', press, '_', thresh,'_all-data_NoSLR.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(preds, is.na(Seaward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, !is.na(Seaward))) +
  tm_bubbles('SeawardMang', 
             palette = pal[1:10],
             midpoint = 0.5,
             breaks = seq(0,1,0.1),
             size = 'SeawardMang',
             scale = 0.25,
             alpha = 0.5, 
             border.alpha = 0,
             legend.size.show = F,
             legend.col.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = 'A) Seaward baseline forecast',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', press, '_', thresh,'_all-data_NoSLR.png'), width = 5, height = 1, dpi = 5000)

# now take biophysical matrices and solve with future SLR and 
chosen_model <- models$mangrove_model # choose correct network model

# find future biophysical and pressure scenarios
slr_lenient <- spatial_dat %>% filter(pressure_def == 1) %>% select(Type, fut_slr) # subset data for 'very lenient' future slr

# get unique biophysical/pressure combinations historically, given 5 different pressure definitions
dat <- spatial_dat %>% 
  filter(pressure_def == press) %>% 
  select(-fut_slr) %>% # remove future slr and replace with 'very lenient' fut slr
  left_join(slr_lenient, by = 'Type') %>% 
  dplyr::select(pressure_def, Type, fut_csqueeze, fut_csqueeze_1, fut_dams, Tidal_Class, prop_estab, climate, fut_cdev, fut_slr, 
                fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>%
  rename(sed_supp = fut_dams, csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, cdev = fut_cdev, ant_slr = fut_slr, 
         gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>% 
  mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  pivot_longer(cols = c(csqueeze_1, ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab),
         climate_2 = paste0('climate_', .$climate),
         cdev_2 = paste0('cdev_', .$cdev)) %>% 
  select(-c(pressure_def, Type))
# reorder column names so always in same order regardless of filtering
dat <- dat[,c('csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'climate', 'cdev', 'vals_gwsub', 'vals_hist_drought',
                                'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 
                                 'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
                                'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                                'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2', 'climate_2', 'cdev_2')]
dat <- dat %>% 
  unite('scenario', csqueeze_2:cdev_2, na.rm = T, sep = '.') %>% 
  unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  distinct()

# get hindcast matrices

tmp <- readRDS(paste0('outputs/simulation-outcomes/scenario_matrices.RDS'))
matrices <- lapply(tmp, function(x){x[[2]]})
matrix_index <- data.frame(index = 1:length(matrices), scenario = unlist(lapply(tmp, function(x){names(x)[1]})))

# get forecast matrices with no historical analogue

bio_dat <- dat %>% select(csqueeze:cdev, scenario) %>% distinct() %>% filter(!scenario %in% matrix_index$scenario) # unique biophysical contexts with no historical analogue

# simulate matrices for each biophysical context with no historical analogue

nsim <- 1000 # number of sims
tmp <- list()
system.time(
  for(k in 1:nrow(bio_dat)){
    if(bio_dat[k,'climate'] == 'arid'){
      aridity <- 'Y'
    }else{aridity <- 'N'}
    if(bio_dat[k,'cdev'] == 'None'){
      cdevel <- c(0,0)
    }else if(bio_dat[k,'cdev'] == 'Low'){
      cdevel <- c(0,0.33)
    }else if(bio_dat[k,'cdev'] == 'Medium'){
      cdevel <- c(0.34,0.66)
    }else if(bio_dat[k,'cdev'] == 'High'){
      cdevel <- c(0.67,1)
    }
    tmp[[k]] <- sim_mod(bio_dat[k,], nsim, arid = aridity, prob = c(0.1, 0.5), cdev = cdevel)
    names(tmp[[k]]) <- bio_dat[k,]$scenario
  }
)
saveRDS(tmp, paste0('outputs/simulation-outcomes/scenario_matrices_forecast.RDS'))
#tmp <- readRDS(paste0('outputs/simulation-outcomes/scenario_matrices_forecast.RDS'))
matrices_forecast <- lapply(tmp, function(x){x[[2]]})
matrix_index_forecast <- data.frame(index = (max(matrix_index$index)+1):(max(matrix_index$index)+length(matrices_forecast)), scenario = unlist(lapply(tmp, function(x){names(x)[1]})))

# bind hindcast and foreacst matrices

matrices <- c(matrices, matrices_forecast)
matrix_index <- rbind(matrix_index, matrix_index_forecast)

# loop through each unique biophysical-pressure scenario and solve matrices in the relevant biophysical model
# to obtain naive forecasts of loss/gain given management/conservation scenarios

dat2 <- dat %>% mutate(split = rep(1:5, (nrow(.)+4)/5)[1:nrow(.)]) # add a column to split so can parallelise
scenarios <- list('LandwardAvailableProp', 'Hydrology', 'SubVol', 'Coastalsqueeze', 'SeaLevelRise')
system.time( # takes ~1 hour
for(b in seq_along(scenarios)){
  scn <- scenarios[[b]]
  if(scn[1] == 'Hydrology'){
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
          new_mat[8,7] <- if(new_mat[8,7] != 0 & new_mat[8,7] < 0.67){runif(1, 0.67, 1)}else if(new_mat[8,7] == 0){0}else{new_mat[8,7]} # update so interaction between landward propagules and landward mangrove is always high (i.e., between 0.67-1), increasing hydro connectivity
          new_mat[11,10] <- if(new_mat[11,10] != 0 & new_mat[11,10] < 0.67){runif(1, 0.67, 1)}else if(new_mat[11,10] == 0){0}else{new_mat[11,10]}  # update so interaction between seaward propagules and seaward mangrove is always high (i.e., between 0.67-1), increasing hydro connectivity
          tmp[[j]] <- data.frame(solver(new_mat, chosen_model, c(pressures$press, 'LandwardAvailableProp', 'SeawardAvailableProp')), nsim = j)  
        }
        tmp2[[i]] <- do.call(rbind, tmp) %>% mutate(scenario = as.character(datsub[i,'scenario']), press = as.character(datsub[i,'press']))
      }
      do.call(rbind, tmp2)
    }
    stopCluster(cl)
    results <- do.call(rbind, results) #%>% pivot_wider(names_from = 'node', values_from = 'outcome')
    write.csv(results, paste0('outputs/predictions/naive_outcomes_scenario_', scn[1], '_', press, '_', thresh, '.csv'), row.names = F)
  }else if(scn[1] %in% c('LandwardAvailableProp','SubVol')){
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
    results <- do.call(rbind, results) #%>% pivot_wider(names_from = 'node', values_from = 'outcome')
    write.csv(results, paste0('outputs/predictions/naive_outcomes_scenario_', scn[1], '_', press, '_', thresh, '.csv'), row.names = F)
  }else if(scn[1] == 'Coastalsqueeze'){
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
          new_mat[8,9] <- runif(1, 0.67, 1) # update so interaction between sea level rise and landward mangrove is always high (i.e., 0.67 and 1)
          tmp[[j]] <- data.frame(solver(new_mat, chosen_model, pressures$press), nsim = j)  
        }
        tmp2[[i]] <- do.call(rbind, tmp) %>% mutate(scenario = as.character(datsub[i,'scenario']), press = as.character(datsub[i,'press']))
      }
      do.call(rbind, tmp2)
    }
    stopCluster(cl)
    results <- do.call(rbind, results) #%>% pivot_wider(names_from = 'node', values_from = 'outcome')
    write.csv(results, paste0('outputs/predictions/naive_outcomes_scenario_', scn[1], '_', press, '_', thresh, '.csv'), row.names = F)
}else if(scn[1] == 'SeaLevelRise'){
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
        tmp[[j]] <- data.frame(solver(bio_model[,,j], chosen_model, pressures$press), nsim = j)  
      }
      tmp2[[i]] <- do.call(rbind, tmp) %>% mutate(scenario = as.character(datsub[i,'scenario']), press = as.character(datsub[i,'press']))
    }
    do.call(rbind, tmp2)
  }
  stopCluster(cl)
  results <- do.call(rbind, results) #%>% pivot_wider(names_from = 'node', values_from = 'outcome')
  write.csv(results, paste0('outputs/predictions/naive_outcomes_scenario_', scn[1], '_', press, '_', thresh, '.csv'), row.names = F)
}
  })

# choose management/conservation scenario and map outcomes

for(i in seq_along(scenarios)){
  
scenario <- scenarios[[i]][1]

naive_outcomes_restore <- read.csv(paste0('outputs/predictions/naive_outcomes_scenario_', scenario, '_', press, '_', thresh, '.csv')) %>% 
  pivot_wider(names_from = 'node', values_from = 'outcome')

# make posterior forecast for each future scenario, using matrix calibrated posterior probabilities

spatial_pred <- spatial_dat %>% # here renaming future pressures as historical pressures so can join to posterior hindcasts
  filter(pressure_def == press) %>% 
  select(-fut_slr) %>% # remove future slr and replace with 'very lenient' fut slr
  left_join(slr_lenient, by = 'Type') %>% 
  dplyr::select(pressure_def, Type, fut_dams, fut_csqueeze, fut_csqueeze_1, Tidal_Class, prop_estab, climate, fut_cdev, fut_slr, 
                fut_gwsub, fut_drought, fut_ext_rain, fut_storms, land_net_change_obs, sea_net_change_obs) %>%
  rename(sed_supp = fut_dams, csqueeze = fut_csqueeze, csqueeze_1 = fut_csqueeze_1, cdev = fut_cdev, ant_slr = fut_slr, 
         gwsub = fut_gwsub, hist_drought = fut_drought, hist_ext_rain = fut_ext_rain, storms = fut_storms) %>%
  mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab),
         climate_2 = paste0('climate_', .$climate),
         cdev_2 = paste0('cdev_', .$cdev))
# reorder column names so always in same order regardless of filtering
spatial_pred <- spatial_pred[,c('pressure_def', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'climate', 'cdev',
                                  'land_net_change_obs', 'sea_net_change_obs',  'vals_gwsub', 'vals_hist_drought',
                                  'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 
                                'vals_csqueeze_1','press_no_press', 'press_gwsub',
                                  'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                                  'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2', 'climate_2', 'cdev_2')]
spatial_pred <- spatial_pred %>% 
  unite('scenario', csqueeze_2:cdev_2, na.rm = T, sep = '.') %>% 
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
write.csv(naive_outcomes_restore, paste0('outputs/validation/naive_outcomes_restore_', press, '_', thresh, '_', scenario, '.csv'), row.names = F)
naive_outcomes_restore <- read.csv(paste0('outputs/validation/naive_outcomes_restore_', press, '_', thresh, '_', scenario, '.csv'))

spatial_pred_fit <- spatial_pred %>% 
  left_join(naive_outcomes_restore, by = c('scenario', 'press')) %>% 
  left_join(post_prob, by = c('scenario', 'nsim')) %>% 
  mutate(LandwardMang = ifelse(LandwardMang == -1, 0, LandwardMang), # here turn losses into a 0 so just calculating the probability of gain/neutrality
         SeawardMang = ifelse(SeawardMang == -1, 0, SeawardMang),
         matrix_post_prob = ifelse(is.na(matrix_post_prob), 1, matrix_post_prob)) %>% # here make matrix probabilities 1 for those that are missing a historical analogue (so essentially are unfitted predictions)
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
if(scenario == 'SubVol'){spatial_pred_fit <- spatial_pred_fit %>% mutate(Class = sub('\\_.*', '', Type)) %>% mutate(Seaward = ifelse(Class == 'OpenCoast', 'Loss', Seaward))} # in subvol scenario, only allow reduced risk of loss or gain/neutrality in geomorphologies other than open coast
write.csv(spatial_pred_fit, paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenario, '_fit.csv'), row.names = F)
spatial_pred_fit <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[i]][1], '_fit.csv'))

spatial_pred_unfit <- spatial_pred %>% 
  left_join(naive_outcomes_restore, by = c('scenario', 'press')) %>% 
  left_join(post_prob, by = c('scenario', 'nsim')) %>% 
  mutate(LandwardMang = ifelse(LandwardMang == -1, 0, LandwardMang), # here turn losses into a 0 so just calculating the probability of gain/neutrality
         SeawardMang = ifelse(SeawardMang == -1, 0, SeawardMang),
         matrix_post_prob = 1) %>% 
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
if(scenario == 'SubVol'){spatial_pred_unfit <- spatial_pred_unfit %>% mutate(Class = sub('\\_.*', '', Type)) %>% mutate(Seaward = ifelse(Class == 'OpenCoast', 'Loss', Seaward))} # in subvol scenario, only allow reduced risk of loss or gain/neutrality in geomorphologies other than open coast
write.csv(spatial_pred_unfit, paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenario, '_unfit.csv'), row.names = F)
spatial_pred_unfit <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[i]][1], '_unfit.csv'))
titlea_fit <- c('C) Seaward forecast with increased propagules - fit', 'D) Seaward forecast with improved hydrology - fit','E) Seaward forecast with sediment addition - fit', 'F) Seaward forecast with removal of coastal barriers - fit', 'A) Seaward baseline forecast')
titleb_fit <- c('D) Landward forecast with increased propagules - fit', 'D) Landward forecast with improved hydrology - fit', 'F) Landward forecast with sediment addition - fit', 'H) Landward forecast with removal of coastal barriers - fit', 'C) Landward baseline forecast')
titlea_unfit <- c('C) Seaward forecast with increased propagules - unfit', 'D) Seaward forecast with improved hydrology - unfit', 'E) Seaward forecast with sediment addition - unfit', 'G) Seaward forecast with removal of coastal barriers - unfit', 'A) Seaward baseline forecast - unfit')
titleb_unfit <- c('D) Landward forecast with increased propagules - unfit', 'D) Landward forecast with improved hydrology - unfit', 'F) Landward forecast with sediment addition - unfit', 'H) Landward forecast with removal of coastal barriers - unfit', 'C) Landward baseline forecast - unfit')

# summarise predictions
datsum <- spatial_pred_fit %>% 
  mutate(Landward = paste0('Landward_', .$Landward),
         Seaward = paste0('Seaward_', .$Seaward)) %>% 
  mutate(Landward_seaward = paste0(.$Landward, '.', .$Seaward)) %>% 
  group_by(Landward_seaward) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
write.csv(datsum, paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[i]][1], '_summary_fit.csv'), row.names = F)    

datsum <- spatial_pred_unfit %>% 
mutate(Landward = paste0('Landward_', .$Landward),
       Seaward = paste0('Seaward_', .$Seaward)) %>% 
  mutate(Landward_seaward = paste0(.$Landward, '.', .$Seaward)) %>% 
  group_by(Landward_seaward) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
write.csv(datsum, paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[i]][1], '_summary_unfit.csv'), row.names = F) 

# map final 'all data' forecasts

preds_fit <- typ_points %>% 
  left_join(spatial_pred_fit) %>%
  #mutate_at(vars(LandwardMang, SeawardMang), ~ifelse(. >100-thresh & . <thresh, 100-thresh, .)) %>% 
  mutate_at(vars(LandwardMang, SeawardMang), ~./100) %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

preds_unfit <- typ_points %>% 
  left_join(spatial_pred_unfit) %>%
  #mutate_at(vars(LandwardMang, SeawardMang), ~ifelse(. >100-thresh & . <thresh, 100-thresh, .)) %>% 
  mutate_at(vars(LandwardMang, SeawardMang), ~./100) %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map forecasts
#fit
lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Landward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds_fit) +
  tm_bubbles('LandwardMang', 
           palette = pal[1:10], 
           midpoint = 0.5,
           breaks = seq(0,1,0.1),
           size = 'LandwardMang',
           scale = 0.25,
           alpha = 0.5, 
           border.alpha = 0,
           legend.size.show = F,
           legend.col.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.001),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = titleb_fit[i],
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]), title = 'Probability',
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', press, '_', thresh, '_all-data', '_', scenarios[[i]][1], '_fit.png'), width = 5, height = 1, dpi = 5000)

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Landward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds_fit) +
  tm_dots('Landward', 
             palette = c('Loss' = pal[1], 'Ambiguous' = pal[5], 'Gain_neutrality' = pal[10]), 
             alpha = 0.5, 
             size = 0.0015,
             border.alpha = 0,
             legend.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.001),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = 'B) Landward baseline forecast',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[c(1,5,10)]),
                labels =  c('Gain/Neutrality', 'Ambiguous', 'Loss'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', press, '_', thresh, '_all-data', '_', scenarios[[i]][1], '_fit_class.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Seaward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds_fit) +
  tm_bubbles('SeawardMang', 
             palette = pal[1:10], 
             midpoint = 0.5,
             breaks = seq(0,1,0.1),
             size = 'SeawardMang',
             scale = 0.25,
             alpha = 0.5, 
             border.alpha = 0,
             legend.size.show = F,
             legend.col.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = titlea_fit[i],
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),title = 'Probability',
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', press, '_', thresh, '_all-data', '_', scenarios[[i]][1], '_fit.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Landward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds_fit) +
  tm_dots('Seaward', 
          palette = c('Loss' = pal[1], 'Ambiguous' = pal[5], 'Gain_neutrality' = pal[10]), 
          alpha = 0.5, 
          size = 0.0015,
          border.alpha = 0,
          legend.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.001),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = 'A) Seaward baseline forecast',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[c(1,5,10)]),
                labels =  c('Gain/Neutrality', 'Ambiguous', 'Loss'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', press, '_', thresh, '_all-data', '_', scenarios[[i]][1], '_fit_class.png'), width = 5, height = 1, dpi = 5000)

# unfit
lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Landward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds_unfit) +
  tm_bubbles('LandwardMang', 
             palette = pal[1:10],
             midpoint = 0.5,
             breaks = seq(0,1,0.1),
             size = 'LandwardMang',
             scale = 0.25,
             alpha = 0.5, 
             border.alpha = 0,
             legend.size.show = F,
             legend.col.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = titleb_unfit[i],
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', press, '_', thresh, '_all-data', '_', scenarios[[i]][1], '_unfit.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Seaward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds_unfit) +
  tm_bubbles('SeawardMang', 
             palette = pal[1:10],
             midpoint = 0.5,
             breaks = seq(0,1,0.1),
             size = 'SeawardMang',
             scale = 0.25,
             alpha = 0.5, 
             border.alpha = 0,
             legend.size.show = F,
             legend.col.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = titlea_unfit[i],
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', press, '_', thresh, '_all-data', '_', scenarios[[i]][1], '_unfit.png'), width = 5, height = 1, dpi = 5000)
}

# get the landward and seaward forecasts for each scenario, and show where the additional gains or reduced risk of loss would be

scenarios <- list('LandwardAvailableProp', 'Hydrology', 'SubVol', 'Coastalsqueeze')
baseline <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', 'SeaLevelRise', '_fit.csv')) %>%
  filter(!is.na(Landward) & !is.na(Seaward)) %>% 
  rename('Landward_base' = 'LandwardMang', 'Seaward_base' = 'SeawardMang')
transplant <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[1]][1], '_fit.csv')) %>% 
  filter(!is.na(Landward) & !is.na(Seaward)) %>% 
  inner_join(select(baseline,Landward_base, Seaward_base, Type)) %>% 
  mutate(Transplant_Landward_gain = ifelse(Landward_base < thresh & LandwardMang > thresh, 'Plant', NA),
         Transplant_Landward_risk_reduced = ifelse(Landward_base < 100-thresh & LandwardMang > 50 & LandwardMang < thresh, 'Plant', NA),
         Transplant_Seaward_gain = ifelse(Seaward_base < thresh & SeawardMang > thresh, 'Plant', NA),
         Transplant_Seaward_risk_reduced = ifelse(Seaward_base < 100-thresh & SeawardMang > 50 & SeawardMang < thresh, 'Plant', NA))
hydrology <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[2]][1], '_fit.csv')) %>% 
  filter(!is.na(Landward) & !is.na(Seaward)) %>% 
  inner_join(select(baseline,Landward_base, Seaward_base, Type)) %>% 
  mutate(Hydrology_Landward_gain = ifelse(Landward_base < thresh & LandwardMang > thresh, 'Hydrology', NA),
         Hydrology_Landward_risk_reduced = ifelse(Landward_base < 100-thresh & LandwardMang > 50 & LandwardMang < thresh,  'Hydrology', NA),
         Hydrology_Seaward_gain = ifelse(Seaward_base < thresh & SeawardMang > thresh, 'Hydrology', NA),
         Hydrology_Seaward_risk_reduced = ifelse(Seaward_base < 100-thresh & SeawardMang > 50 & SeawardMang < thresh, 'Hydrology', NA))
sediment <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[3]][1], '_fit.csv')) %>% 
  filter(!is.na(Landward) & !is.na(Seaward)) %>% 
  inner_join(select(baseline,Landward_base, Seaward_base, Type)) %>% 
  mutate(Sediment_Landward_gain = ifelse(Landward_base < thresh & LandwardMang > thresh, 'Sediment', NA),
         Sediment_Landward_risk_reduced = ifelse(Landward_base < 100-thresh & LandwardMang > 50 & LandwardMang < thresh, 'Sediment', NA),
         Sediment_Seaward_gain = ifelse(Seaward_base < thresh & SeawardMang > thresh, 'Sediment', NA),
         Sediment_Seaward_risk_reduced = ifelse(Seaward_base < 100-thresh & SeawardMang > 50 & SeawardMang < thresh, 'Sediment', NA))
barriers <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', scenarios[[4]][1], '_fit.csv')) %>% 
  filter(!is.na(Landward) & !is.na(Seaward)) %>% 
  inner_join(select(baseline,Landward_base, Seaward_base, Type)) %>% 
  mutate(Barriers_Landward_gain = ifelse(Landward_base < thresh & LandwardMang > thresh, 'Barriers', NA),
         Barriers_Landward_risk_reduced = ifelse(Landward_base < 100-thresh & LandwardMang > 50 & LandwardMang < thresh,  'Barriers', NA),
         Barriers_Seaward_gain = ifelse(Seaward_base < thresh & SeawardMang > thresh, 'Barriers', NA),
         Barriers_Seaward_risk_reduced = ifelse(Seaward_base < 100-thresh & SeawardMang > 50 & SeawardMang < thresh,  'Barriers', NA))

all_scen <- transplant %>% left_join(hydrology, by = 'Type') %>% left_join(sediment, by = 'Type') %>% left_join(barriers, by = 'Type') %>% 
  unite('Landward_scenario_gain', Transplant_Landward_gain, Hydrology_Landward_gain, Sediment_Landward_gain, Barriers_Landward_gain, na.rm = T, sep = '_') %>% 
  unite('Landward_scenario_reduced_risk', Transplant_Landward_risk_reduced, Hydrology_Landward_risk_reduced, Sediment_Landward_risk_reduced, Barriers_Landward_risk_reduced, na.rm = T, sep = '_') %>% 
  unite('Seaward_scenario_gain', Transplant_Seaward_gain, Hydrology_Seaward_gain, Sediment_Seaward_gain, Barriers_Seaward_gain, na.rm = T, sep = '_') %>% 
  unite('Seaward_scenario_reduced_risk', Transplant_Seaward_risk_reduced, Hydrology_Seaward_risk_reduced, Sediment_Seaward_risk_reduced, Barriers_Seaward_risk_reduced, na.rm = T, sep = '_') %>% 
  mutate_at(vars(Landward_scenario_gain, Landward_scenario_reduced_risk,
                 Seaward_scenario_gain, Seaward_scenario_reduced_risk), ~ifelse(. == '', NA, .)) %>% 
  select(Type, Landward_scenario_gain, Landward_scenario_reduced_risk, Seaward_scenario_gain, Seaward_scenario_reduced_risk)
write.csv(all_scen, 'outputs/predictions/scenario-forecasts-all.csv', row.names = F)
all_scen <- read.csv('outputs/predictions/scenario-forecasts-all.csv')

# summarise predictions

# landward
datsum <- all_scen %>% 
  group_by(Landward_scenario_gain) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
datsum
write.csv(datsum, paste0('outputs/summary-stats/forecast-predictions_', press, '_', thresh, 'all_scenario_gain_summary-landward.csv'), row.names = F)   

datsum <- all_scen %>% 
  group_by(Landward_scenario_reduced_risk) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
datsum
write.csv(datsum, paste0('outputs/summary-stats/forecast-predictions_', press, '_', thresh, 'all_scenario_reduced_risk_summary-landward.csv'), row.names = F)   

# seaward
datsum <- all_scen %>% 
  group_by(Seaward_scenario_gain) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
datsum
write.csv(datsum, paste0('outputs/summary-stats/forecast-predictions_', press, '_', thresh, 'all_scenario_gain_summary-seaward.csv'), row.names = F)   

datsum <- all_scen %>% 
  group_by(Seaward_scenario_reduced_risk) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
datsum
write.csv(datsum, paste0('outputs/summary-stats/forecast-predictions_', press, '_', thresh, 'all_scenario_reduced_risk_summary-seaward.csv'), row.names = F)   

# landward and seaward

datsum <- all_scen %>% 
  mutate_at(c('Landward_scenario_gain', 'Landward_scenario_reduced_risk', 
              'Seaward_scenario_gain','Seaward_scenario_reduced_risk'), ~ifelse(!is.na(.), 1, 0)) %>% 
  mutate(gain_sum = Landward_scenario_gain + Seaward_scenario_gain,
         reduced_risk_sum = Landward_scenario_reduced_risk + Seaward_scenario_reduced_risk) %>% 
  mutate(gain = ifelse(gain_sum > 0, 1, 0),
         reduced_risk = ifelse(reduced_risk_sum > 0, 1, 0),
         reduced_risk_only = ifelse(gain_sum == 0 & reduced_risk_sum > 0, 1, 0)) %>% 
  pivot_longer(c(gain:reduced_risk_only), names_to = 'type', values_to = 'change') %>% 
  group_by(type) %>% 
  summarise(n = sum(change), percent = 100*(sum(change)/((nrow(.)/3))))
datsum
write.csv(datsum, paste0('outputs/summary-stats/forecast-predictions_', press, '_', thresh, 'all_scenario_summary.csv'), row.names = F)   

# map

scenario_change <- typ_points %>% 
  inner_join(all_scen) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

# landward 
lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(scenario_change, !is.na(Landward_scenario_gain))) +
  tm_dots('Landward_scenario_gain', 
          palette = c('Barriers' = 'darkcyan', 'Plant' = 'yellowgreen', 'Hydrology' = 'darkgoldenrod2', 'Plant_Barriers' = 'darkorange3', 'Hydrology_Barriers' = 'darkorchid4', 'Plant_Hydrology' = 'deeppink4', 'Plant_Hydrology_Barriers' = 'black', 'Plant_Hydrology_Sediment_Barriers' = 'navyblue'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.015) +
  tm_shape(filter(scenario_change, !is.na(Landward_scenario_reduced_risk))) +
  tm_dots('Landward_scenario_reduced_risk', 
          palette = c('Barriers' = 'darkcyan', 'Plant' = 'yellowgreen', 'Hydrology' = 'darkgoldenrod2', 'Plant_Barriers' = 'darkorange3', 'Hydrology_Barriers' = 'darkorchid4', 'Plant_Hydrology' = 'deeppink4', 'Plant_Hydrology_Barriers' = 'black', 'Plant_Hydrology_Sediment_Barriers' = 'navyblue'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_layout(legend.outside = F,
            legend.position = c(0.04, 0.01),
            legend.width = 1,
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title =  "D) Landward scenario-based forecast of net gain/neutrality or reduced certainty of loss relative to baseline",
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col =  c('white', 'white', 'white'), alpha = 0.1, 
                title = 'Action',
                labels =  c('B = Removal of barriers to landward migration', 'L = Increased landward propagules (assisted dispersal or enrichment)',
                            'EC = Improved ecological connectivity (assuming natural recruitment of propagules)'), border.alpha = 0, size = 0.3) +
  tm_add_legend('symbol', col =  c('darkcyan', 'yellowgreen','darkgoldenrod2', 'darkorchid4',  'deeppink4', 'black'), alpha = 0.8, is.portrait = F,
                labels =  c('B', 'L', 'EC', 'B or EC', 'L or EC', 'L or B or EC'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', press, '_', thresh, '_all-data', '_gain_reduced_risk_scenario.png'), width = 5, height = 1, dpi = 5000)

# seaward
smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(scenario_change, !is.na(Seaward_scenario_gain))) +
  tm_dots('Seaward_scenario_gain', 
          palette = c('Sediment' = 'plum4',  'Hydrology' = 'darkgoldenrod2', 'Hydrology_Sediment' = 'hotpink3', 'Plant_Hydrology' = 'black', 'Plant_Hydrology_Sediment' = 'midnightblue'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.015) +
  tm_shape(filter(scenario_change, !is.na(Seaward_scenario_reduced_risk))) +
  tm_dots('Seaward_scenario_reduced_risk', 
          palette = c('Sediment' = 'plum4',  'Hydrology' = 'darkgoldenrod2', 'Hydrology_Sediment' = 'hotpink3','Plant_Hydrology' = 'black', 'Plant_Hydrology_Sediment' = 'midnightblue'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_layout(legend.outside = F,
            legend.position = c(0.04, 0.01),
            title.position = c(0.01,0.45),
            legend.width = 1,
            legend.title.size = 0.35,
            legend.text.size = 0.25,
            main.title = "B) Seaward scenario-based forecast of net gain/neutrality or reduced certainty of loss relative to baseline",
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col =  c('white', 'white'), alpha = 0.8, title = 'Action',
                labels =  c('S = Sediment addition/trapping', 'L = Increased landward propagules (assisted dispersal or enrichment)',  'EC = Improved ecological connectivity (assuming natural recruitment of propagules)'), border.alpha = 0, size = 0.3) +
  tm_add_legend('symbol', col =  c('plum4', 'darkgoldenrod2',  'hotpink3', 'black', 'midnightblue'), alpha = 0.8,
                labels =  c('S', 'EC', 'S or EC', 'L or EC', 'S or EC or L'), border.alpha = 0, size = 0.25, is.portrait = F)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', press, '_', thresh, '_all-data', '_gain_reduced_risk_scenario.png'), width = 5, height = 1, dpi = 5000)

# end here

