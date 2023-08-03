# map hindcast matches and mismatches for a given pressure and ambiguity threshold

library(tidyverse)
library(ggh4x)
library(sf)
library(tmap)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
set.seed(123) # set random number generator to make results reproducible
sf_use_s2(FALSE)

go <- 1 # which coastal dev threshold?
press <- 3 # which pressure definition threshold?
thresh <- 75 # which ambiguity threshold?
rm_e <- 'N' # remove erosion from validation? Y or N

# read in spatial data

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
world <- data("World")
spatial_dat <- read.csv('outputs/master-dat.csv')
drivers <- read.csv('data/typologies/SLR_Data.csv')
naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes_', go,'_', rm_e,'.csv'))
results <- readRDS(paste0('outputs/validation/accuracy_', go,'_', rm_e, '.RDS'))
accuracy <- do.call(rbind, lapply(results, function(x)x[[1]]))
test_hindcasts <- do.call(rbind, lapply(results, function(x)x[[2]]))

accuracy_sum <- accuracy %>% # summarise accuracy across folds
  filter(mangrove != 'Seaward & Landward') %>% 
  mutate(mangrove = case_when(mangrove == 'Seaward' ~ 'C) Seaward', mangrove == 'Landward' ~ 'D) Landward')) %>% 
  mutate(mangrove = factor(mangrove, levels = c('C) Seaward', 'D) Landward'))) %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  group_by(mangrove, pressure_def, ambig_threshold, class, metric) %>% 
  summarise(accuracy = mean(accuracy, na.rm = T)) 

filter(accuracy_sum, pressure_def == press, ambig_threshold == thresh) # look at accuracy for a given pressure and ambiguity threshold

if(rm_e == 'N'){
preds <- typ_points %>% # join hindcasts to spatial data
  left_join(filter(test_hindcasts, pressure_def == press, ambig_threshold == thresh)) %>% 
  mutate(Seaward_match = case_when(Seaward == 'Ambiguous' ~ 'Ambiguous',
                                   is.na(SeawardMang) ~'No Hindcast',
                                   Seaward == sea_net_change ~'Match', 
                                   Seaward != sea_net_change ~ 'Mis-match'),
         Landward_match = case_when(Landward == 'Ambiguous' ~ 'Ambiguous',
                                    is.na(LandwardMang) ~'No Hindcast',
                                    Landward == land_net_change ~'Match', 
                                    Landward != land_net_change ~ 'Mis-match')) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
}else{
  preds <- typ_points %>% # join hindcasts to spatial data
    left_join(select(drivers, Type, Erosion), by = 'Type') %>% 
    filter(Erosion < 0.1) %>% 
    left_join(filter(test_hindcasts, pressure_def == press, ambig_threshold == thresh)) %>% 
    mutate(Seaward_match = case_when(Seaward == 'Ambiguous' ~ 'Ambiguous',
                                     is.na(SeawardMang) ~'No Prediction',
                                     Seaward == sea_net_change ~'Match', 
                                     Seaward != sea_net_change ~ 'MisMatch'),
           Landward_match = case_when(Landward == 'Ambiguous' ~ 'Ambiguous',
                                      is.na(LandwardMang) ~'No Prediction',
                                      Landward == land_net_change ~'Match', 
                                      Landward != land_net_change ~ 'MisMatch')) %>%
    st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
}
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

# map 

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(preds) +
  tm_dots('Landward_match', 
          palette = c('No Hindcast' = 'red', 'Ambiguous' = 'lightgoldenrod', 'Mis-match' = 'black', 'Match' = 'palegreen4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.3,
            main.title = 'D) Landward hindcast matches and mis-matches with optimal thresholds',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('palegreen4','lightgoldenrod','black' ,'red'), 
                labels =  c( 'Match', 'Ambiguous', 'Mis-match','No Hindcast'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_match_', go, '_', rm_e, '_', press, '_', thresh, '.png'), width = 5, height = 1)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(preds) +
  tm_dots('Seaward_match', 
          palette = c('No Hindcast' = 'red', 'Ambiguous' = 'lightgoldenrod', 'Mis-match' = 'black', 'Match' = 'palegreen4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.01),
            title.position = c(0.01,0.45),
            legend.title.size = 0.45,
            legend.text.size = 0.3,
            main.title = 'C) Seaward hindcast matches and mis-matches with optimal thresholds',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('palegreen4','lightgoldenrod','black' ,'red'), 
                labels =  c( 'Match', 'Ambiguous', 'Mis-match','No Hindcast'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_match_', go, '_', rm_e, '_', press, '_', thresh, '.png'), width = 5, height = 1)

# map hindcasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Landward == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Landward == 'Loss' & !is.na(Change))) +
  tm_dots('Landward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Landward == 'Gain_neutrality' & !is.na(Change))) +
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
            main.title = 'B) Landward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_', go, '_', rm_e, '_', press, '_', thresh, '.png'), width = 5, height = 1)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(preds, is.na(Change))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Ambiguous' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F, 
          size = 0.0015) +
  tm_shape(filter(preds, Seaward == 'Loss' & !is.na(Change))) +
  tm_dots('Seaward', 
          palette = c('Ambiguous' = 'lightgoldenrod', 'Loss' = 'firebrick4', 'Gain_neutrality' = 'deepskyblue4'), 
          alpha = 0.5, 
          title = '',
          legend.show = F,
          size = 0.001) +
  tm_shape(filter(preds, Seaward == 'Gain_neutrality' & !is.na(Change))) +
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
            main.title = 'B) Seaward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', go, '_', rm_e, '_', press, '_', thresh, '.png'), width = 5, height = 1)

# now make posterior predictions for each biophysical setting/pressure model using all the data (i.e. not split by kfolds)
# to be used for making forecasts

# prepare all data prediction dataset
pred_dat <- spatial_dat %>% 
  filter(Cdev_thresh == go & pressure_def == press) %>% 
  dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, ant_slr, gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
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
pred_dat <- pred_dat[,c('pressure_def', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 
                         'land_net_change_obs', 'sea_net_change_obs',  'vals_gwsub', 'vals_hist_drought',
                          'vals_hist_ext_rain', 'vals_storms', 'vals_ant_slr', 'vals_csqueeze_1','press_no_press', 'press_gwsub',
                          'press_hist_drought', 'press_hist_ext_rain', 'press_storms', 'press_ant_slr',
                          'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2')]
pred_dat <- pred_dat %>% 
  unite('scenario', csqueeze_2:prop_estab_2, na.rm = T, sep = '.') %>% 
  unite('press',  press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  mutate(press = ifelse(!is.na(press_no_press), 'none', press))

# all data posterior probs
post_prob <- pred_dat %>% 
  left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
  mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
  group_by(scenario, nsim) %>% 
  summarise(matrix_post_prob = mean(valid)) 
write.csv(post_prob, paste0('outputs/validation/matrix-posterior-prob', go, '_', rm_e, '_', press, '_', thresh, '.csv'), row.names = F)

# make predictions 
final_preds <- pred_dat %>% 
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
write.csv(final_preds, paste0('outputs/predictions/final-calibrated-predictions_', go, '_', rm_e, '_', press, '_', thresh,'.csv'), row.names = F)

# map final 'all data' hindcasts

preds <- typ_points %>% 
  left_join(final_preds) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map hindcasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  #tm_shape(filter(preds, is.na(Landward))) +
  #tm_dots('darkgrey', size = 0.001) +
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
            main.title = 'B) Landward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
lmap

tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_', go, '_', rm_e, '_', press, '_', thresh, '_all-data.png'), width = 5, height = 1)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  #tm_shape(filter(preds, is.na(Seaward))) +
  #tm_dots('darkgrey', size = 0.001) +
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
            main.title = 'B) Seaward hindcast',
            main.title.size = 0.45,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('firebrick4', 'lightgoldenrod', 'deepskyblue4'), 
                labels =  c('Loss','Ambiguous', 'Gain/Neutrality'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', go, '_', rm_e, '_', press, '_', thresh, '_all-data.png'), width = 5, height = 1)
