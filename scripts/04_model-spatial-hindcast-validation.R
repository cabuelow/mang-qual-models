# spatial model hindcast calibration and cross validation

library(QPress)
library(tidyverse)
library(foreach)
library(doParallel)
library(caret)
library(ggh4x)
library(sf)
library(tmap)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
set.seed(123) # set random number generator to make results reproducible
sf_use_s2(FALSE)

# read in spatial data (mangrove typological units)

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('outputs/master-dat.csv')
drivers <- read.csv('data/typologies/SLR_Data.csv')

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model

# remove erosion from validation?
rm_e <- 'N'

system.time( # takes ~ 2 hours
for(go in 1:3){ # loop over coastal development threshold

  if(rm_e == 'N'){  
# simulate a set of matrices for each mangrove network model (i.e., biophysical setting) 
# and store the outcome for landward/seaward mangroves, i.e., gain/neutrality or loss

dat <- spatial_dat %>% # get unique biophysical/pressure combinations historically, given 5 different pressure definitions
  filter(Cdev_thresh == go) %>% 
  dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, #ant_slr, 
                gwsub, hist_drought, hist_ext_rain, storms) %>% 
  #pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  pivot_longer(cols = c(csqueeze_1,gwsub:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  select(-c(pressure_def,Type))
# reorder column names so always in same order regardless of filtering
dat <- dat[,c('csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'vals_gwsub', 'vals_hist_drought',
              'vals_hist_ext_rain', 'vals_storms', #'vals_ant_slr', 
              'vals_csqueeze_1', 'press_gwsub',
              'press_hist_drought', 'press_hist_ext_rain', 'press_storms', #'press_ant_slr',
              'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
dat <- dat %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  distinct()

bio_dat <- dat %>% select(csqueeze:prop_estab, scenario) %>% distinct() # unique biophysical contexts

# simulate matrices for each biophysical context

nsim <- 1000 # number of sims
tmp <- list()
system.time(
  for(k in 1:nrow(bio_dat)){
    tmp[[k]] <- sim_mod(bio_dat[k,], nsim)
    names(tmp[[k]]) <- bio_dat[k,]$scenario
  }
)
saveRDS(tmp, paste0('outputs/simulation-outcomes/scenario_matrices_', go,'_', rm_e, '.RDS'))
#tmp <- readRDS('outputs/simulation-outcomes/scenario_matrices.RDS')
matrices <- lapply(tmp, function(x){x[[2]]})
names(matrices) <- bio_dat$scenario
matrix_index <- data.frame(index = 1:length(matrices), scenario = unlist(lapply(tmp, function(x){names(x)[1]})))

# loop through each unique biophysical-pressure scenario and solve matrices in the relevant biophysical model
# to obtain naive hindcasts of loss/gain

dat2 <- dat %>% mutate(split = rep(1:5, (nrow(.)+4)/5)[1:nrow(.)]) # add a column to split so can parallelise
cl <- makeCluster(5)
registerDoParallel(cl)
system.time( # takes 22 mins
results <- foreach(h = seq_along(unique(dat2$split)), .packages = c('tidyverse', 'QPress')) %dopar% {
  datsub <- dat2 %>% filter(split == h)
  tmp2 <- list() # list for storing results
for(i in 1:nrow(datsub)){
  bio_model <- matrices[[filter(matrix_index, scenario == as.character(datsub[i,'scenario']))$index]]
  pressures <- data.frame(press = unlist(strsplit(as.character(datsub[i,'press']), '\\.'))) %>% 
    mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                          'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
  tmp <- vector("list", dim(bio_model)[3])
  for(j in 1:dim(bio_model)[3]){
    tmp[[j]] <- data.frame(solver(bio_model[,,j], chosen_model, pressures$press), nsim = j)  
  }
  tmp2[[i]] <- data.frame(do.call(rbind, tmp) %>% 
    mutate(scenario = as.character(datsub[i,'scenario']),
           press = as.character(datsub[i,'press'])))
}
  do.call(rbind, tmp2) %>% pivot_wider(names_from = 'node', values_from = 'outcome') %>% data.frame()
})
stopCluster(cl)
naive_outcomes <- do.call(rbind, results)

# split data into training and test folds

kfold <- 5 # number of folds

shuffled_dat <- spatial_dat %>% # shuffle the data and split
  filter(Cdev_thresh == go) %>% 
  left_join(data.frame(Type = .[1:length(unique(.$Type)),]$Type[sample(1:length(unique(.$Type)))],
                       k = c(rep(1:kfold, each = round(length(unique(.$Type))/kfold)),1:kfold)[1:length(unique(.$Type))]), by = 'Type') %>% 
  dplyr::select(pressure_def, k, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, #ant_slr, 
                gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  #mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  #pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  pivot_longer(cols = c(csqueeze_1,gwsub:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab))
# reorder column names so always in same order regardless of filtering
shuffled_dat <- shuffled_dat[,c('pressure_def', 'k', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'land_net_change_obs', 'sea_net_change_obs',
              'vals_gwsub', 'vals_hist_drought', 'vals_hist_ext_rain', 'vals_storms', #'vals_ant_slr', 
              'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
              'press_hist_drought', 'press_hist_ext_rain', 'press_storms', #'press_ant_slr',
              'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
shuffled_dat <- shuffled_dat %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  mutate(press = ifelse(!is.na(press_no_press), 'none', press))

# identify the unique 'no pressure' bio scenario combinations and add in 1000 sims as neutral outcomes to naive hindcasts
d <- select(filter(shuffled_dat, press == 'none'), press, scenario) %>% distinct()
naive_outcomes <- rbind(naive_outcomes, data.frame(nsim = rep(1:nsim, nrow(d)), 
                                                   scenario = rep(d$scenario, each = nsim), 
                                                   press = rep(d$press, each = nsim),
                                                   LandwardMang = 1,
                                                   SeawardMang = 1))
write.csv(naive_outcomes, paste0('outputs/validation/naive_outcomes_', go, '_', rm_e,'.csv'), row.names = F)
#naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes_', go, '_', rm_e,'.csv'))

# map the folds

#m.dat <- typ_points %>% left_join(filter(shuffled_dat, pressure_def == 1)) %>% mutate(k = factor(k)) #%>% st_transform(crs = 'ESRI:54030')
#world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33) 
#map <- tm_shape(world_mang) +
 # tm_fill(col = 'gray85') +
  #tm_shape(filter(m.dat, !is.na(k))) + 
  #tm_dots('k', title = '', legend.is.portrait = F, palette = 'Set2', size = 0.01) +
  #tm_layout(frame = T, legend.position = c(0.2, 0))
#map
#tmap_save(map, 'outputs/maps/k-fold-map.png', width = 8, height = 2)

# join naive outcomes/hindcasts to each unit in a training fold, validate against observed mangrove loss or gain
# calculate likelihood/posterior probability of each matrix in a biophysical model based on observed mangrove loss or gain
# use training posterior probabilities to make posterior hindcasts in test units and quantify accuracy
# do this looping over kfolds and pressure definition and ambiguity thresholds

threshold <- seq(60, 90, by = 5) # range of thresholds for defining when a prediction is ambiguous or not
cl <- makeCluster(5)
registerDoParallel(cl)
system.time( # ~30 mins
  results <- foreach(i = 1:kfold, .packages = c('tidyverse', 'caret'), .errorhandling = 'remove') %dopar% {
    acc <- list() # list to store accuracy outcomes
    preds <- list() # list to store predictions
    for(j in seq_along(threshold)){
      thresh <- threshold[j] # set ambiguity threshold
      acc2 <- list() # list to store accuracy outcomes
      preds2 <- list() # list to store predictions
      for(h in seq_along(unique(shuffled_dat$pressure_def))){
        
        # get training units for a fold and pressure definition
        train_units <- shuffled_dat %>% 
          filter(k != i & pressure_def == h)
        
        # calculate likelihood/posterior probability of each matrix in each bio model given observed loss or gain
        train_post_prob <- train_units %>% 
          left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
          mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
          group_by(scenario, nsim) %>% 
          summarise(matrix_post_prob = mean(valid)) 
        
        # get test units and make posterior predictions/hindcasts using posterior probability of each biomodel matrix
        test_units <- shuffled_dat %>% 
          filter(k == i & pressure_def == h)
        
        # join naive hindcasts and posterior probabilities, calculate posterior hindcasts for test units
        test_post_hindcasts <- test_units %>% 
          left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
          left_join(train_post_prob, by = c('scenario', 'nsim')) %>% 
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
                                     .default = 'Ambiguous'),
                 Change = case_when(LandwardMang >= thresh & SeawardMang >= thresh ~ "Gain",
                                    LandwardMang < 100-thresh & SeawardMang < 100-thresh ~ "Loss",
                                    Landward == 'Ambiguous' ~ "Ambiguous",
                                    Seaward == 'Ambiguous' ~ "Ambiguous",
                                    LandwardMang >= thresh & SeawardMang < 100-thresh ~ "Landward Gain & Seaward Loss",
                                    LandwardMang < 100-thresh & SeawardMang >= thresh ~ "Landward Loss & Seaward Gain",
                                    .default = NA),
                 Change_obs = case_when(land_net_change_obs == 1 & sea_net_change_obs == 1 ~ "Gain",
                                        land_net_change_obs == -1 & sea_net_change_obs == -1 ~ "Loss",
                                        land_net_change_obs == 1 & sea_net_change_obs == -1 ~ "Landward Gain & Seaward Loss",
                                        land_net_change_obs == -1 & sea_net_change_obs == 1 ~ "Landward Loss & Seaward Gain"),
                 land_net_change = ifelse(land_net_change_obs == -1, 'Loss', 'Gain_neutrality'),
                 sea_net_change = ifelse(sea_net_change_obs == -1, 'Loss', 'Gain_neutrality')) %>% 
          mutate(ambig_threshold = thresh)
        preds2[[h]] <- test_post_hindcasts
        
        test_set <- test_post_hindcasts %>% 
          filter(Landward != 'Ambiguous' & Seaward != 'Ambiguous')
        
        results <- data.frame(mangrove = 'Landward',
                              k = i,
                              pressure_def = h,
                              ambig_threshold = thresh, 
                              calc_accuracy(test_set$Landward, test_set$land_net_change)$accuracy.results)
        results2 <- data.frame(mangrove = 'Seaward',
                               k = i,
                               pressure_def = h,
                               ambig_threshold = thresh, 
                               calc_accuracy(test_set$Seaward, test_set$sea_net_change)$accuracy.results)
        results3 <- data.frame(mangrove = 'Seaward & Landward',
                               k = i,
                               pressure_def = h,
                               ambig_threshold = thresh, 
                               calc_accuracy(test_set$Change, test_set$Change_obs)$accuracy.results)
        acc2[[h]] <- rbind(results, results2, results3)
      }
      acc[[j]] <- do.call(rbind, acc2)
      preds[[j]] <- do.call(rbind, preds2)
    }
    list(do.call(rbind, acc), do.call(rbind, preds))
  })
stopCluster(cl)
saveRDS(results, paste0('outputs/validation/accuracy_', go,'_', rm_e, '.RDS'))
#results <- readRDS(paste0('outputs/validation/accuracy_', go,'_', rm_e, '.RDS'))
accuracy <- do.call(rbind, lapply(results, function(x)x[[1]])) %>% 
  mutate(pressure_def = recode(pressure_def, '1' = 'Very lenient',
                               '2' = 'Lenient',
                               '3' = 'Moderate',
                               '4' = 'Strict', 
                               '5' = 'Very strict')) %>% 
  mutate(pressure_def = factor(pressure_def, levels = c('Very lenient', 'Lenient', 'Moderate','Strict', 'Very strict')))
test_hindcasts <- do.call(rbind, lapply(results, function(x)x[[2]]))

# heatmap of accuracy metrics for combinations of pressure and ambiguity thresholds

# indidvidual kfolds
accuracy %>% 
  filter(mangrove != 'Seaward & Landward') %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  ggplot() +
  aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy', limits = c(0,100)) +
  facet_nested_wrap(~factor(k) + factor(mangrove) + factor(class) + factor(metric), nrow = 5) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold (%)') +
  theme_classic()

ggsave(paste0('outputs/validation/accuracy-heatmap_kfold_', go, '_', rm_e,'.png'), width = 15, height = 11)

# averaged across kfolds
accuracy_sum <- accuracy %>% 
  filter(mangrove != 'Seaward & Landward') %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'A) Seaward', mangrove == 'Landward' ~ 'B) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('A) Seaward', 'B) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  group_by(mangrove, pressure_def, ambig_threshold, class, metric) %>% 
  summarise(accuracy = mean(accuracy, na.rm = T)) 

accuracy_sum %>% 
  ggplot() +
  aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy', limits = c(0,100)) +
  facet_nested_wrap(~factor(mangrove) + factor(class) + factor(metric), nrow = 2) +
  ylab('Pressure definition threshold') +
  xlab('Ambiguity probability threshold (%)') +
  theme_classic() +
  theme(axis.title = element_text(size = 10))

ggsave(paste0('outputs/validation/accuracy-heatmap_kfold_averaged_', go, '_', rm_e,'.png'), width = 10, height = 4.5)
  }else{
    # simulate a set of matrices for each mangrove network model (i.e., biophysical setting) 
    # and store the outcome for landward/seaward mangroves, i.e., gain/neutrality or loss
    
    dat <- spatial_dat %>% # get unique biophysical/pressure combinations historically, given 5 different pressure definitions
      filter(Cdev_thresh == go) %>% 
      dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
      pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
      filter(vals == 1) %>% 
      pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
      mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
             sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
             Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
             prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
      select(-c(pressure_def,Type))
    # reorder column names so always in same order regardless of filtering
    dat <- dat[,c('csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'vals_gwsub', 'vals_hist_drought',
                  'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 'vals_csqueeze_1', 'press_gwsub',
                  'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                  'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
    dat <- dat %>% 
      unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
      unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
      distinct()
    
    bio_dat <- dat %>% select(csqueeze:prop_estab, scenario) %>% distinct() # unique biophysical contexts
    
    # simulate matrices for each biophysical context
    
    nsim <- 1000 # number of sims
    tmp <- list()
    system.time(
      for(k in 1:nrow(bio_dat)){
        tmp[[k]] <- sim_mod(bio_dat[k,], nsim)
        names(tmp[[k]]) <- bio_dat[k,]$scenario
      }
    )
    saveRDS(tmp, paste0('outputs/simulation-outcomes/scenario_matrices_', go,'_', rm_e, '.RDS'))
    #tmp <- readRDS('outputs/simulation-outcomes/scenario_matrices.RDS')
    matrices <- lapply(tmp, function(x){x[[2]]})
    names(matrices) <- bio_dat$scenario
    matrix_index <- data.frame(index = 1:length(matrices), scenario = unlist(lapply(tmp, function(x){names(x)[1]})))
    
    # loop through each unique biophysical-pressure scenario and solve matrices in the relevant biophysical model
    # to obtain naive hindcasts of loss/gain
    
    dat2 <- dat %>% mutate(split = rep(1:5, (nrow(.)+4)/5)[1:nrow(.)]) # add a column to split so can parallelise
    cl <- makeCluster(5)
    registerDoParallel(cl)
    system.time( # takes 22 mins
      results <- foreach(h = seq_along(unique(dat2$split)), .packages = c('tidyverse', 'QPress')) %dopar% {
        datsub <- dat2 %>% filter(split == h)
        tmp2 <- list() # list for storing results
        for(i in 1:nrow(datsub)){
          bio_model <- matrices[[filter(matrix_index, scenario == as.character(datsub[i,'scenario']))$index]]
          pressures <- data.frame(press = unlist(strsplit(as.character(datsub[i,'press']), '\\.'))) %>% 
            mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                                  'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
          tmp <- vector("list", dim(bio_model)[3])
          for(j in 1:dim(bio_model)[3]){
            tmp[[j]] <- data.frame(solver(bio_model[,,j], chosen_model, pressures$press), nsim = j)  
          }
          tmp2[[i]] <- data.frame(do.call(rbind, tmp) %>% 
                                    mutate(scenario = as.character(datsub[i,'scenario']),
                                           press = as.character(datsub[i,'press'])))
        }
        do.call(rbind, tmp2) %>% pivot_wider(names_from = 'node', values_from = 'outcome') %>% data.frame()
      })
    stopCluster(cl)
    naive_outcomes <- do.call(rbind, results)
    
    # split data into training and test folds
    
    kfold <- 5 # number of folds
    
    shuffled_dat <- spatial_dat %>% # shuffle the data and split
      left_join(select(drivers, Type, Erosion), by = 'Type') %>% 
      filter(Erosion < 0.1) %>% 
      filter(Cdev_thresh == go) %>% 
      left_join(data.frame(Type = .[1:length(unique(.$Type)),]$Type[sample(1:length(unique(.$Type)))],
                           k = c(rep(1:kfold, each = round(length(unique(.$Type))/kfold)),1:kfold)[1:length(unique(.$Type))]), by = 'Type') %>% 
      dplyr::select(pressure_def, k, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
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
    shuffled_dat <- shuffled_dat[,c('pressure_def', 'k', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'land_net_change_obs', 'sea_net_change_obs',
                                    'vals_gwsub', 'vals_hist_drought', 'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
                                    'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                                    'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
    shuffled_dat <- shuffled_dat %>% 
      unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
      unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
      mutate(press = ifelse(!is.na(press_no_press), 'none', press))
    
    # identify the unique 'no pressure' bio scenario combinations and add in 1000 sims as neutral outcomes to naive hindcasts
    d <- select(filter(shuffled_dat, press == 'none'), press, scenario) %>% distinct()
    naive_outcomes <- rbind(naive_outcomes, data.frame(nsim = rep(1:nsim, nrow(d)), 
                                                       scenario = rep(d$scenario, each = nsim), 
                                                       press = rep(d$press, each = nsim),
                                                       LandwardMang = 1,
                                                       SeawardMang = 1))
    write.csv(naive_outcomes, paste0('outputs/validation/naive_outcomes_', go, '_', rm_e,'.csv'), row.names = F)
    #naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes_', go, '_', rm_e,'.csv'))
    
    # map the folds
    
    #m.dat <- typ_points %>% left_join(filter(shuffled_dat, pressure_def == 1)) %>% mutate(k = factor(k)) #%>% st_transform(crs = 'ESRI:54030')
    #world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33) 
    #map <- tm_shape(world_mang) +
    # tm_fill(col = 'gray85') +
    #tm_shape(filter(m.dat, !is.na(k))) + 
    #tm_dots('k', title = '', legend.is.portrait = F, palette = 'Set2', size = 0.01) +
    #tm_layout(frame = T, legend.position = c(0.2, 0))
    #map
    #tmap_save(map, 'outputs/maps/k-fold-map.png', width = 8, height = 2)
    
    # join naive outcomes/hindcasts to each unit in a training fold, validate against observed mangrove loss or gain
    # calculate likelihood/posterior probability of each matrix in a biophysical model based on observed mangrove loss or gain
    # use training posterior probabilities to make posterior hindcasts in test units and quantify accuracy
    # do this looping over kfolds and pressure definition and ambiguity thresholds
    
    threshold <- seq(60, 90, by = 5) # range of thresholds for defining when a prediction is ambiguous or not
    cl <- makeCluster(5)
    registerDoParallel(cl)
    system.time( # ~30 mins
      results <- foreach(i = 1:kfold, .packages = c('tidyverse', 'caret'), .errorhandling = 'remove') %dopar% {
        acc <- list() # list to store accuracy outcomes
        preds <- list() # list to store predictions
        for(j in seq_along(threshold)){
          thresh <- threshold[j] # set ambiguity threshold
          acc2 <- list() # list to store accuracy outcomes
          preds2 <- list() # list to store predictions
          for(h in seq_along(unique(shuffled_dat$pressure_def))){
            
            # get training units for a fold and pressure definition
            train_units <- shuffled_dat %>% 
              filter(k != i & pressure_def == h)
            
            # calculate likelihood/posterior probability of each matrix in each bio model given observed loss or gain
            train_post_prob <- train_units %>% 
              left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
              mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
              group_by(scenario, nsim) %>% 
              summarise(matrix_post_prob = mean(valid)) 
            
            # get test units and make posterior predictions/hindcasts using posterior probability of each biomodel matrix
            test_units <- shuffled_dat %>% 
              filter(k == i & pressure_def == h)
            
            # join naive hindcasts and posterior probabilities, calculate posterior hindcasts for test units
            test_post_hindcasts <- test_units %>% 
              left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
              left_join(train_post_prob, by = c('scenario', 'nsim')) %>% 
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
                                         .default = 'Ambiguous'),
                     Change = case_when(LandwardMang >= thresh & SeawardMang >= thresh ~ "Gain",
                                        LandwardMang < 100-thresh & SeawardMang < 100-thresh ~ "Loss",
                                        Landward == 'Ambiguous' ~ "Ambiguous",
                                        Seaward == 'Ambiguous' ~ "Ambiguous",
                                        LandwardMang >= thresh & SeawardMang < 100-thresh ~ "Landward Gain & Seaward Loss",
                                        LandwardMang < 100-thresh & SeawardMang >= thresh ~ "Landward Loss & Seaward Gain",
                                        .default = NA),
                     Change_obs = case_when(land_net_change_obs == 1 & sea_net_change_obs == 1 ~ "Gain",
                                            land_net_change_obs == -1 & sea_net_change_obs == -1 ~ "Loss",
                                            land_net_change_obs == 1 & sea_net_change_obs == -1 ~ "Landward Gain & Seaward Loss",
                                            land_net_change_obs == -1 & sea_net_change_obs == 1 ~ "Landward Loss & Seaward Gain"),
                     land_net_change = ifelse(land_net_change_obs == -1, 'Loss', 'Gain_neutrality'),
                     sea_net_change = ifelse(sea_net_change_obs == -1, 'Loss', 'Gain_neutrality')) %>% 
              mutate(ambig_threshold = thresh)
            preds2[[h]] <- test_post_hindcasts
            
            test_set <- test_post_hindcasts %>% 
              filter(Landward != 'Ambiguous' & Seaward != 'Ambiguous')
            
            results <- data.frame(mangrove = 'Landward',
                                  k = i,
                                  pressure_def = h,
                                  ambig_threshold = thresh, 
                                  calc_accuracy(test_set$Landward, test_set$land_net_change)$accuracy.results)
            results2 <- data.frame(mangrove = 'Seaward',
                                   k = i,
                                   pressure_def = h,
                                   ambig_threshold = thresh, 
                                   calc_accuracy(test_set$Seaward, test_set$sea_net_change)$accuracy.results)
            results3 <- data.frame(mangrove = 'Seaward & Landward',
                                   k = i,
                                   pressure_def = h,
                                   ambig_threshold = thresh, 
                                   calc_accuracy(test_set$Change, test_set$Change_obs)$accuracy.results)
            acc2[[h]] <- rbind(results, results2, results3)
          }
          acc[[j]] <- do.call(rbind, acc2)
          preds[[j]] <- do.call(rbind, preds2)
        }
        list(do.call(rbind, acc), do.call(rbind, preds))
      })
    stopCluster(cl)
    saveRDS(results, paste0('outputs/validation/accuracy_', go,'_', rm_e, '.RDS'))
    #results <- readRDS('outputs/validation/accuracy.RDS')
    accuracy <- do.call(rbind, lapply(results, function(x)x[[1]])) %>% 
      mutate(pressure_def = recode(pressure_def, '1' = 'Very lenient',
                                   '2' = 'Lenient',
                                   '3' = 'Moderate',
                                   '4' = 'Strict', 
                                   '5' = 'Very strict')) %>% 
      mutate(pressure_def = factor(pressure_def, levels = c('Very lenient', 'Lenient', 'Moderate','Strict', 'Very strict')))
    test_hindcasts <- do.call(rbind, lapply(results, function(x)x[[2]]))
    
    # heatmap of accuracy metrics for combinations of pressure and ambiguity thresholds
    
    # indidvidual kfolds
    accuracy %>% 
      filter(mangrove != 'Seaward & Landward') %>% 
      mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
      mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
      pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
      mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
      distinct() %>% 
      ggplot() +
      aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
      geom_tile() +
      scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy', limits=c(0,100)) +
      facet_nested_wrap(~factor(k) + factor(mangrove) + factor(class) + factor(metric), nrow = 5) +
      ylab('Pressure definition') +
      xlab('Ambiguity probability threshold (%)') +
      theme_classic()
    
    ggsave(paste0('outputs/validation/accuracy-heatmap_kfold_', go, '_', rm_e,'.png'), width = 15, height = 11)
    
    # averaged across kfolds
    accuracy_sum <- accuracy %>% 
      filter(mangrove != 'Seaward & Landward') %>% 
      mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
      mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
      pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
      mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
      distinct() %>% 
      group_by(mangrove, pressure_def, ambig_threshold, class, metric) %>% 
      summarise(accuracy = mean(accuracy, na.rm = T)) 
    
    accuracy_sum %>% 
      ggplot() +
      aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
      geom_tile() +
      scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy', limits=c(0,100)) +
      facet_nested_wrap(~factor(mangrove) + factor(class) + factor(metric), nrow = 2) +
      ylab('Pressure definition') +
      xlab('Ambiguity probability threshold (%)') +
      theme_classic()
    
    ggsave(paste0('outputs/validation/accuracy-heatmap_kfold_averaged_', go, '_', rm_e,'.png'), width = 10, height = 4.5)
}
}
)

# End here

# accuracy for simultaneous seaward and landward 

accuracy_sum2 <- accuracy %>% 
  filter(mangrove == 'Seaward & Landward') %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  group_by(mangrove, pressure_def, ambig_threshold, class, metric) %>% 
  summarise(accuracy = mean(accuracy)) 

accuracy_sum2 %>% 
  ggplot() +
  aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(class) + factor(metric), nrow = 2) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/accuracy-heatmap_kfold_averaged_overall.png',  width = 10, height = 2)

filter(accuracy_sum2, accuracy == max(filter(accuracy_sum2, metric == 'Overall_accuracy' & class == 'Gain_neutrality & Loss')$accuracy))
