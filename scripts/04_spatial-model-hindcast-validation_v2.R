# spatial model hindcast and cross validation

library(QPress)
library(tidyverse)
#library(foreach)
#library(doParallel)
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

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model
chosen_model_name <- 'mangrove_model'

# simulate a set of matrices for each mangrove network model (i.e., biophysical setting and pressure combinations), 
# and store the outcome for landward/seaward mangroves, i.e., gain/neutrality or loss

# identify unique biophysical setting/pressure scenarios

press_dat <- spatial_dat %>% 
  dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  select(-c(pressure_def,Type)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T) %>% 
  distinct()

# simulate matrices for each scenario

nsim <- 100 # number of sims
tmp <- list()
system.time(
for(k in 1:nrow(press_dat)){ #TODO: not sure why apply won't work here over rows instead of forloop
    tmp[[k]] <- sim_mod(press_dat[k,], numsims = nsim)
    names(tmp[[k]]) <- press_dat[k,]$scenario
}
)
saveRDS(tmp, 'outputs/simulation-outcomes/scenario_matrices.RDS')
tmp <- readRDS('outputs/simulation-outcomes/scenario_matrices.RDS')

# extract outcomes and matrices for each scenario

outcomes <- do.call(rbind, lapply(tmp, function(x){data.frame(x[[3]], scenario = names(x)[1])})) %>% 
  filter(var %in% c('LandwardMang', 'SeawardMang')) %>% 
  mutate(outcome = case_when(outcome >= 0 ~ 1,
                             outcome < 0 ~ -1)) %>% 
  pivot_wider(names_from = 'var', values_from = 'outcome')
matrices <- lapply(tmp, function(x){x[[4]]})
names(matrices) <- press_dat$scenario

# calculate likelihood/posterior probability of each matrix based on observed mangrove loss or gain

kfold <- 5 # number of folds

shuffled_dat <- spatial_dat %>% 
  left_join(data.frame(Type = spatial_dat[1:length(unique(spatial_dat$Type)),]$Type[sample(1:length(unique(spatial_dat$Type)))],
                       k = rep(1:kfold, each = length(unique(spatial_dat$Type))/kfold)), by = 'Type') %>% 
  dplyr::select(pressure_def, k, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T)

matrix_likelihood <- shuffled_dat %>% 
  left_join(outcomes) %>% 
  mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
  group_by(pressure_def, k, scenario, nsim, LandwardMang, SeawardMang) %>% 
  summarise(matrix_post_prob = mean(valid)) 

# plot for just one training

ggplot(filter(matrix_likelihood, k != 2, scenario %in% unique(matrix_likelihood$scenario)[1:6])) +
  geom_point(aes(x = nsim, y = matrix_post_prob, col = factor(pressure_def)), alpha = 0.5) +
  facet_wrap(~scenario) +
  theme_classic()
  
# make posterior predictions accounting for matrix likelihood as a weighted sum using training data
# go through each k fold unit, make a new predictions using get relevant model outcomes and likelihoods from training set

preds <- list()
acc <- list()
for(i in 1:kfold){
  # make predictions using training data
  train_preds <- matrix_likelihood %>% 
    filter(k != i) %>% 
    mutate(LandwardMang = LandwardMang*matrix_post_prob,
           SeawardMang = SeawardMang*matrix_post_prob) %>% 
    group_by(pressure_def, scenario) %>% 
    summarise(LandwardMang = sum(LandwardMang)/sum(matrix_post_prob),
              SeawardMang = sum(SeawardMang)/sum(matrix_post_prob)) %>% 
    mutate(LandwardMang = case_when(LandwardMang >= 0.8 ~ 1,
                                    LandwardMang < -0.8 ~ -1),
           SeawardMang = case_when(SeawardMang >= 0.8 ~ 1,
                                   SeawardMang < -0.8 ~ -1))
  
  # join to test data based by pressure def and scenario
  test_set <- shuffled_dat %>% 
    filter(k == i) %>% 
    left_join(train_preds) %>% 
    mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0))
  
  acc2 <- list()
  for(j in seq_along(unique(test_set$pressure_def))){
    test_set2 <- test_set %>% filter(pressure_def == j)
    results <- data.frame(mangrove = 'Landward',
                          k = i,
                          pressure_def = j,
                          #ambig_threshold = thresh, 
                          calc_accuracy(test_set2$LandwardMang, test_set2$land_net_change_obs)$accuracy.results)
    results2 <- data.frame(mangrove = 'Seaward',
                           k = i,
                           pressure_def = j,
                           #ambig_threshold = thresh, 
                           calc_accuracy(test_set2$SeawardMang, test_set2$sea_net_change_obs)$accuracy.results)
    acc2[[j]] <- rbind(results, results2)
  }
  acc[[i]] <- do.call(rbind, acc2)
  preds[[i]] <- test_set
}

accuracy <- do.call(rbind, acc)
test_preds <- do.call(rbind, preds)
write.csv(accuracy, 'outputs/validation/hindcast-accuracy_pressure-ambiguity_kfold.csv', row.names = F)

# heatmap of accuracy metrics for combinations of pressure and ambiguity thresholds

# indidvidual kfolds
accuracy %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  ggplot() +
  #aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  aes(x = 1, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(k) + factor(mangrove) + factor(class) + factor(metric), nrow = 5) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/accuracy-heatmap_kfold.png', width = 12, height = 9)

# averaged across kfolds
accuracy %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  #group_by(mangrove, pressure_def, ambig_threshold, class, metric) %>% 
  group_by(mangrove, pressure_def, class, metric) %>% 
  summarise(accuracy = mean(accuracy)) %>% 
  ggplot() +
  #aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  aes(x = 1, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(mangrove) + factor(class) + factor(metric), nrow = 2) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/accuracy-heatmap_kfold_averaged.png',  width = 10, height = 4.5)

# map prediction matches and mismatches

preds <- typ_points %>% 
  left_join(filter(test_preds, pressure_def == 4)) %>% 
  mutate(Seaward_match = case_when(is.na(SeawardMang) ~'No Prediction',
                                   SeawardMang == sea_net_change_obs ~'Match', 
                                   SeawardMang != sea_net_change_obs ~ 'MisMatch'),
         Landward_match = case_when(is.na(LandwardMang) ~'No Prediction',
                                   LandwardMang == land_net_change_obs ~'Match', 
                                   LandwardMang != land_net_change_obs ~ 'MisMatch')) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(preds) +
  tm_dots('Landward_match', 
          palette = c('No Prediction' = 'lightgoldenrod', 'Mismatch' = 'black', 'Match' = 'palegreen4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'B) Landward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('black', 'lightgoldenrod', 'palegreen4'), 
                labels =  c('Mismatch','No Prediction', 'Match'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_match_', chosen_model_name, '.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(preds) +
  tm_dots('Seaward_match', 
          palette = c('No Prediction' = 'lightgoldenrod', 'Mismatch' = 'black', 'Match' = 'palegreen4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.35,
            main.title = 'A) Seaward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('black', 'lightgoldenrod', 'palegreen4'), 
                labels =  c('Mismatch','No Prediction', 'Match'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_match_', chosen_model_name, '.png'), width = 5, height = 3)
