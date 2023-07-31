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
source('scripts/helpers/spatial-helpers.R')
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
  #dplyr::select(pressure_def, Type, csqueeze, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
  #pivot_longer(cols = c(ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  select(-c(pressure_def,Type)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.') %>% 
  #unite('scenario', csqueeze_2:prop_estab_2, press_gwsub:press_ant_slr, na.rm = T, sep = '.') %>% 
  distinct()

# simulate matrices for each scenario

nsim <- 1000 # number of sims
tmp <- list()
system.time( # takes 4 mins
for(k in 1:nrow(press_dat)){ #TODO: not sure why apply won't work here over rows instead of forloop
    tmp[[k]] <- sim_mod(press_dat[k,], nsim)
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

# calculate likelihood/posterior probability of each matrix in a scenario model based on observed mangrove loss or gain

kfold <- 5 # number of folds

# shuffle the data and split
shuffled_dat <- spatial_dat %>% 
  left_join(data.frame(Type = .[1:length(unique(.$Type)),]$Type[sample(1:length(unique(.$Type)))],
                       k = rep(1:kfold, each = length(unique(.$Type))/kfold)), by = 'Type') %>% 
  dplyr::select(pressure_def, k, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  #dplyr::select(pressure_def, k, Type, csqueeze, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  #pivot_longer(cols = c(ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  #unite('scenario', csqueeze_2:prop_estab_2, press_gwsub:press_ant_slr, na.rm = T, sep = '.')
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.')

# map the folds

m.dat <- typ_points %>% left_join(filter(shuffled_dat, pressure_def == 1)) %>% mutate(k = factor(k)) #%>% st_transform(crs = 'ESRI:54030')
map <- tm_shape(world_mang) +
  tm_fill(col = 'gray85') +
  tm_shape(filter(m.dat, !is.na(k))) + 
  tm_dots('k', title = '', legend.is.portrait = F, palette = 'Set2', size = 0.01) +
  tm_layout(frame = T, legend.position = c(0.2, 0))
map
tmap_save(map, 'outputs/maps/k-fold-map.png', width = 8, height = 2)

# plot to check

#ggplot(filter(matrix_likelihood, k != 2, scenario %in% unique(matrix_likelihood$scenario)[1:6])) +
 # geom_point(aes(x = nsim, y = matrix_post_prob, col = factor(pressure_def)), alpha = 0.5) +
  #facet_wrap(~scenario) +
  #theme_classic()
  
# go through each k fold unit, make training posterior predictions using relevant scenario matrices 
# matrix likelihoods are calibrated against observed data
# use training posterior predictions as predictions in test set, and quantify accuracy

threshold <- seq(60, 90, by = 5) # range of thresholds for defining when a prediction is ambiguous or not
cl <- makeCluster(5)
registerDoParallel(cl)
system.time( # < 10mins
results <- foreach(i = 1:kfold, .packages = c('tidyverse', 'caret')) %dopar% {
  preds <- list()
  acc <- list()
  for(j in seq_along(threshold)){
    thresh <- threshold[j]
    preds2 <- list()
    acc2 <- list()
    for(h in seq_along(unique(shuffled_dat$pressure_def))){
 
       # make posterior predictions using training data
  train_preds <- shuffled_dat %>% 
      filter(k != i & pressure_def == h) %>% 
      left_join(outcomes, by = 'scenario') %>% 
      mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
      group_by(pressure_def, scenario, nsim, LandwardMang, SeawardMang) %>% 
      summarise(matrix_post_prob = mean(valid)) %>% 
    mutate(LandwardMang = ifelse(LandwardMang == -1, 0, LandwardMang), # here turn losses into a 0 so just calculating the probability of gain/neutrality
           SeawardMang = ifelse(SeawardMang == -1, 0, SeawardMang)) %>% 
    mutate(LandwardMang = LandwardMang*matrix_post_prob,
           SeawardMang = SeawardMang*matrix_post_prob) %>% 
    group_by(pressure_def, scenario) %>% 
    summarise(LandwardMang = (sum(LandwardMang)/sum(matrix_post_prob))*100,
              SeawardMang = (sum(SeawardMang)/sum(matrix_post_prob))*100) %>% 
    mutate(Landward = case_when(is.na(LandwardMang) ~ NA,
                                 LandwardMang >= thresh ~ 'Gain',
                                    LandwardMang < 100-thresh ~ 'Loss',
                                    .default = 'Ambiguous'),
           Seaward = case_when(is.na(SeawardMang) ~ NA,
                                  SeawardMang >= thresh ~ 'Gain',
                                   SeawardMang < 100-thresh ~ 'Loss',
                                   .default = 'Ambiguous'),
           Change = case_when(LandwardMang >= thresh & SeawardMang >= thresh ~ "Gain",
                              LandwardMang < 100-thresh & SeawardMang < 100-thresh ~ "Loss",
                              Landward == 'Ambiguous' | SeawardMang == 'Ambiguous' ~ "Ambiguous",
                              LandwardMang >= thresh & SeawardMang < thresh ~ "Landward Gain & Seaward Loss",
                              LandwardMang < 100-thresh & SeawardMang >= 100-thresh ~ "Landward Loss & Seaward Gain",
                              .default = NA))
  
  # join to test data by pressure def and scenario
  test_set <- shuffled_dat %>% 
    filter(k == i & pressure_def == h) %>% 
    mutate(Change_obs = case_when(land_net_change_obs == 1 & sea_net_change_obs == 1 ~ "Gain",
                                  land_net_change_obs == -1 & sea_net_change_obs == -1 ~ "Loss",
                                  land_net_change_obs == 1 & sea_net_change_obs == -1 ~ "Landward Gain & Seaward Loss",
                                  land_net_change_obs == -1 & sea_net_change_obs == 1 ~ "Landward Loss & Seaward Gain"),
           land_net_change_obs = ifelse(land_net_change_obs == 1, 'Gain', 'Loss'),
           sea_net_change_obs = ifelse(sea_net_change_obs == 1, 'Gain', 'Loss')) %>%
    left_join(train_preds) %>% 
    mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
    mutate(ambig_threshold = thresh)
  preds2[[h]] <- test_set
  
  test_set <- test_set %>% 
    filter(Landward != 'Ambiguous' & Seaward != 'Ambiguous')

    results <- data.frame(mangrove = 'Landward',
                          k = i,
                          pressure_def = h,
                          ambig_threshold = thresh, 
                          calc_accuracy(test_set$Landward, test_set$land_net_change_obs)$accuracy.results)
    results2 <- data.frame(mangrove = 'Seaward',
                           k = i,
                           pressure_def = h,
                           ambig_threshold = thresh, 
                           calc_accuracy(test_set$Seaward, test_set$sea_net_change_obs)$accuracy.results)
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
saveRDS(results, 'outputs/validation/accuracy.RDS')
results <- readRDS('outputs/validation/accuracy.RDS')
accuracy <- do.call(rbind, lapply(results, function(x)x[[1]]))
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
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(k) + factor(mangrove) + factor(class) + factor(metric), nrow = 5) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/accuracy-heatmap_kfold.png', width = 15, height = 11)

# averaged across kfolds
accuracy_sum <- accuracy %>% 
  filter(mangrove != 'Seaward & Landward') %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  group_by(mangrove, pressure_def, ambig_threshold, class, metric) %>% 
  summarise(accuracy = mean(accuracy)) 

accuracy_sum %>% 
  ggplot() +
  aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(mangrove) + factor(class) + factor(metric), nrow = 2) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/accuracy-heatmap_kfold_averaged.png',  width = 10, height = 4.5)

# both seaward and landward 

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

ggsave('outputs/validation/accuracy-heatmap_kfold_averaged_overall.png',  width = 10, height = 4.5)

filter(accuracy_sum2, accuracy == max(filter(accuracy_sum2, metric == 'Overall_accuracy' & class == 'Gain_neutrality & Loss')$accuracy))

# map prediction matches and mismatches

filter(accuracy_sum, pressure_def == 4, ambig_threshold == 85)
filter(accuracy_sum2, pressure_def == 4, ambig_threshold == 85)

preds <- typ_points %>% 
  left_join(filter(test_hindcasts, pressure_def == 4, ambig_threshold == 85)) %>% 
  mutate(Seaward_match = case_when(Seaward == 'Ambiguous' ~ 'Ambiguous',
                                   is.na(SeawardMang) ~'No Prediction',
                                   Seaward == sea_net_change_obs ~'Match', 
                                   Seaward != sea_net_change_obs ~ 'MisMatch'),
         Landward_match = case_when(Landward == 'Ambiguous' ~ 'Ambiguous',
                                   is.na(LandwardMang) ~'No Prediction',
                                   Landward == land_net_change_obs ~'Match', 
                                   Landward != land_net_change_obs ~ 'MisMatch')) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(preds) +
  tm_dots('Landward_match', 
          palette = c('No Prediction' = 'red', 'Ambiguous' = 'lightgoldenrod', 'MisMatch' = 'black', 'Match' = 'palegreen4'), 
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
  tm_add_legend('symbol', col =  c('red','black', 'lightgoldenrod', 'palegreen4'), 
                labels =  c('No Prediction', 'MisMatch', 'Ambiguous', 'Match'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_match_', chosen_model_name, '.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(preds) +
  tm_dots('Seaward_match', 
          palette = c('No Prediction' = 'red', 'Ambiguous' = 'lightgoldenrod', 'MisMatch' = 'black', 'Match' = 'palegreen4'), 
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
  tm_add_legend('symbol', col =  c('red','black', 'lightgoldenrod', 'palegreen4'), 
                labels =  c('No Prediction', 'MisMatch', 'Ambiguous', 'Match'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_match_', chosen_model_name, '.png'), width = 5, height = 3)

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
  tm_shape(filter(preds, Landward == 'Gain' & !is.na(Change))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
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
  tm_shape(filter(preds, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Seaward == 'Loss' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Gain' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
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
            main.title = 'B) Seaward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', chosen_model_name, '.png'), width = 5, height = 3)

# now make posterior predictions for each biophysical setting/pressure model using all the data (i.e. not split by kfolds)
# to be used for making forecasts
# here am only getting the predictions for pressure biophysical/pressure combinations using pressure definition 4

thresh <- 85 # optimal ambiguity threshold
press <- 4 # optimal pressure threshold
final_preds <- spatial_dat %>% 
  dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  #dplyr::select(pressure_def, Type, csqueeze, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  #pivot_longer(cols = c(ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.') %>% 
  #unite('scenario', csqueeze_2:prop_estab_2, press_gwsub:press_ant_slr, na.rm = T, sep = '.') %>% 
  filter(pressure_def == press) %>% 
  left_join(outcomes, by = 'scenario') %>% 
  mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
  group_by(pressure_def, scenario, nsim, LandwardMang, SeawardMang) %>% 
  summarise(matrix_post_prob = mean(valid)) %>% 
  mutate(LandwardMang = ifelse(LandwardMang == -1, 0, LandwardMang), # here turn losses into a 0 so just calcuting the probability of gain/neutrality
         SeawardMang = ifelse(SeawardMang == -1, 0, SeawardMang)) %>% 
  mutate(LandwardMang = LandwardMang*matrix_post_prob,
         SeawardMang = SeawardMang*matrix_post_prob) %>% 
  group_by(pressure_def, scenario) %>% 
  summarise(LandwardMang = (sum(LandwardMang)/sum(matrix_post_prob))*100,
            SeawardMang = (sum(SeawardMang)/sum(matrix_post_prob))*100) %>% 
  mutate(Landward = case_when(is.na(LandwardMang) ~ NA,
                              LandwardMang >= thresh ~ 'Gain',
                              LandwardMang < 100-thresh ~ 'Loss',
                              .default = 'Ambiguous'),
         Seaward = case_when(is.na(SeawardMang) ~ NA,
                             SeawardMang >= thresh ~ 'Gain',
                             SeawardMang < 100-thresh ~ 'Loss',
                             .default = 'Ambiguous'),
         Change = case_when(LandwardMang >= thresh & SeawardMang >= thresh ~ "Gain",
                            LandwardMang < 100-thresh & SeawardMang < 100-thresh ~ "Loss",
                            Landward == 'Ambiguous' | SeawardMang == 'Ambiguous' ~ "Ambiguous",
                            LandwardMang >= thresh & SeawardMang < thresh ~ "Landward Gain & Seaward Loss",
                            LandwardMang < 100-thresh & SeawardMang >= 100-thresh ~ "Landward Loss & Seaward Gain",
                            .default = NA)) %>% 
  mutate(ambig_threshold = thresh)
write.csv(final_preds, paste0('outputs/predictions/final-calibrated-predictions_',press, '_', thresh, '.csv'), row.names = F)

spatial_pred <- spatial_dat %>% 
  #dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  #pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  pivot_longer(cols = c(csqueeze_1,ant_slr:storms), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab)) %>% 
  unite('scenario', csqueeze_2:prop_estab_2, press_csqueeze_1:press_ant_slr, na.rm = T, sep = '.') %>% 
  #unite('scenario', csqueeze_2:prop_estab_2, press_gwsub:press_ant_slr, na.rm = T, sep = '.') %>% 
  filter(pressure_def == press) %>% 
  left_join(final_preds)

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
  tm_shape(filter(preds, Landward == 'Gain' & !is.na(Change))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
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
            main.title = 'B) Landward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_', chosen_model_name, '_all-data.png'), width = 5, height = 3)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Seaward == 'Loss' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Gain' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain' = 'deepskyblue4'), 
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
            main.title = 'B) Seaward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', chosen_model_name, '_all-data.png'), width = 5, height = 3)

