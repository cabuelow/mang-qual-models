# map hindcast matches and mismatches for a given pressure and ambiguity threshold

library(tidyverse)
library(ggh4x)
library(sf)
library(tmap)
library(patchwork)
library(RColorBrewer)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
set.seed(123) # set random number generator to make results reproducible
sf_use_s2(FALSE)

press <- 4 # which pressure definition threshold?
thresh <- 75 # which ambiguity threshold?
pal <- brewer.pal(11, 'Spectral') # colour palette

# read in spatial data

typ_points <- st_read('data/typologies/Mangrove_Typology_v3.14_Composite_valid_centroids.gpkg')
world <- data("World")
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  
meow <- st_read('data/MEOW/meow_ecos.shp')
spatial_dat <- read.csv('data/master-dat.csv')
drivers <- read.csv('data/typologies/SLR_Data.csv')
naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes.csv'))
results <- readRDS(paste0('outputs/validation/accuracy.RDS'))
resamp_accuracy <- read.csv('outputs/validation/resampled_accuracy_summary.csv') %>% 
  mutate(metric = recode(metric, 'Producers_accuracy' = 'Producer accuracy',
                         'Users_accuracy' = 'User accuracy',
                         'Overall_accuracy' = 'Overall accuracy'))

# plot accuracy estimators

a <- ggplot(filter(resamp_accuracy, class == 'Gain_neutrality & Loss'), aes(x = metric, y = median, fill = mangrove)) +
  geom_bar(stat = 'identity', position='dodge') +
  geom_errorbar(aes(ymin=perc_0.025, ymax=perc_0.975), width=0.01, colour="black", position=position_dodge(.9)) +
  xlab('') +
  ylab('') +
  ylim(c(0,100)) +
  scale_fill_manual(values = c('gray30', 'gray')) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 10))+
  coord_flip() +
  ggtitle('A) Gain/neutrality & Loss')  +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 9),
        plot.margin = margin(c(0,0,0,0)))
a
b <- ggplot(filter(resamp_accuracy, class == 'Gain_neutrality'), aes(x = metric, y = median, fill = mangrove)) +
  geom_bar(stat = 'identity', position='dodge') +
  geom_errorbar(aes(ymin=perc_0.025, ymax=perc_0.975), width=0.01, colour="black", position=position_dodge(.9)) +
  xlab('') +
  ylab('') +
  ylim(c(0,100)) +
  scale_fill_manual(values = c('gray30', 'gray')) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 10))+
  coord_flip() +
  ggtitle('B) Gain/neutrality')  +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 9),
        plot.margin = margin(c(0,0,0,0)))
b
c <- ggplot(filter(resamp_accuracy, class == 'Loss'), aes(x = metric, y = median, fill = mangrove)) +
  geom_bar(stat = 'identity', position='dodge') +
  geom_errorbar(aes(ymin=perc_0.025, ymax=perc_0.975), width=0.01, colour="black", position=position_dodge(.9)) +
  xlab('') +
  ylab('Accuracy (%)') +
  ylim(c(0,100)) +
  scale_fill_manual(values = c('gray30', 'gray')) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 10))+
  coord_flip() +
  ggtitle('C) Loss')  +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        plot.title = element_text(size = 9),
        plot.margin = margin(c(0,0,0,0)))
c

d <- a/b/c + plot_layout(heights = c(0.5,1,1))
d

ggsave(paste0('outputs/validation/optimal-accuracy_', press, '_', thresh, '.png'), width = 2.5, height = 4)

# now make posterior predictions for each biophysical setting/pressure model using all the data (i.e. not split by kfolds)

# prepare all data prediction dataset
pred_dat <- spatial_dat %>% 
  filter(pressure_def == press) %>% 
  dplyr::select(pressure_def, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, climate, cdev,#ant_slr, 
                gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
  #mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% # here am removing antecedent sea level rise
  mutate(no_press = gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
  mutate(no_press = ifelse(no_press == 0, 1, 0)) %>% 
  #pivot_longer(cols = c(csqueeze_1,ant_slr:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  pivot_longer(cols = c(csqueeze_1,gwsub:storms, no_press), names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  pivot_wider(names_from = 'press', values_from = c('vals', 'press')) %>% 
  mutate(csqueeze_2 = paste0('Csqueeze_', .$csqueeze),
         sed_supp_2 = paste0('Sedsupp_', .$sed_supp),
         Tidal_Class_2 = paste0('TidalClass_', .$Tidal_Class),
         prop_estab_2 = paste0('Propestab_', .$prop_estab),
         climate_2 = paste0('climate_', .$climate),
         cdev_2 = paste0('cdev_', .$cdev))
# reorder column names so always in same order regardless of filtering
pred_dat <- pred_dat[,c('pressure_def', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'climate', 'cdev',
                         'land_net_change_obs', 'sea_net_change_obs',  'vals_gwsub', 'vals_hist_drought',
                          'vals_hist_ext_rain', 'vals_storms', #'vals_ant_slr', 
                        'vals_csqueeze_1','press_no_press', 'press_gwsub',
                          'press_hist_drought', 'press_hist_ext_rain', 'press_storms', #'press_ant_slr',
                          'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2', 'climate_2','cdev_2')]
pred_dat <- pred_dat %>% 
  unite('scenario', csqueeze_2:cdev_2, na.rm = T, sep = '.') %>% 
  unite('press',  press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
  mutate(press = ifelse(!is.na(press_no_press), 'none', press))

# all data posterior probs
post_prob <- pred_dat %>% 
  left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
  mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
  group_by(scenario, nsim) %>% 
  summarise(matrix_post_prob = mean(valid)) 
write.csv(post_prob, paste0('outputs/validation/matrix-posterior-prob_', press, '_', thresh, '.csv'), row.names = F)

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
write.csv(final_preds, paste0('outputs/predictions/final-calibrated-predictions_', press, '_', thresh,'.csv'), row.names = F)
#final_preds <- read.csv(paste0('outputs/predictions/final-calibrated-predictions_', press, '_', thresh,'.csv'))

final_preds_unfit <- pred_dat %>% 
  left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
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
write.csv(final_preds_unfit, paste0('outputs/predictions/final-calibrated-predictions_', press, '_', thresh,'_unfit.csv'), row.names = F)
#final_preds_unfit <- read.csv(paste0('outputs/predictions/final-calibrated-predictions_', press, '_', thresh,'_unfit.csv'))

# map final 'all data' hindcasts

preds <- typ_points %>% 
  left_join(final_preds) %>%
  #mutate_at(vars(LandwardMang, SeawardMang), ~ifelse(. >=100-thresh & . < thresh, 100-thresh, .)) %>% 
  mutate_at(vars(LandwardMang, SeawardMang), ~./100) %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map fit hindcasts

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(preds, is.na(Landward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds) +
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
            main.title = 'F) Landward hindcast - fit',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_', press, '_', thresh, '_all-data.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(preds, is.na(Seaward))) +
  tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds) +
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
            main.title = 'D) Seaward hindcast - fit',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', press, '_', thresh, '_all-data.png'), width = 5, height = 1, dpi = 5000)

# map unfit hindcasts

preds <- typ_points %>% 
  left_join(final_preds_unfit) %>%
  #mutate_at(vars(LandwardMang, SeawardMang), ~ifelse(. >=100-thresh & . < thresh, 100-thresh, .)) %>% 
  mutate_at(vars(LandwardMang, SeawardMang), ~./100) %>% 
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(LandwardMang))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds) +
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
            main.title = 'E) Landward hindcast - unfit',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
lmap

tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_', press, '_', thresh, '_all-data_unfit.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  #tm_shape(filter(preds, is.na(Seaward))) +
  #tm_dots('darkgrey', size = 0.001) +
  tm_shape(preds) +
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
            main.title = 'C) Seaward hindcast - unfit',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = rev(pal[1:10]),
                labels =  c('100-90% Gain/Neutrality', '90-80% Gain/Neutrality','80-70% Gain/Neutrality', '70-60% Gain/Neutrality', '60-50% Gain/Neutrality', '50-60% Loss', '60-70% Loss', '70-80% Loss', '80-90% Loss', '90-100% Loss'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_', press, '_', thresh, '_all-data_unfit.png'), width = 5, height = 1, dpi = 5000)

# wrangle hindcast matches and mismatches spatially

preds <- typ_points %>% # join hindcasts to spatial data
  left_join(final_preds) %>% 
  mutate(land_net_change = ifelse(land_net_change_obs == -1, 'Loss', 'Gain_neutrality'),
         sea_net_change = ifelse(sea_net_change_obs == -1, 'Loss', 'Gain_neutrality'))%>% 
  mutate(Seaward_match = case_when(Seaward == 'Ambiguous' ~ 'Ambiguous',
                                     is.na(SeawardMang) ~'No Hindcast',
                                     Seaward == sea_net_change ~'Match', 
                                     Seaward != sea_net_change ~ 'Mis-match'),
          Landward_match = case_when(Landward == 'Ambiguous' ~ 'Ambiguous',
                                      is.na(LandwardMang) ~'No Hindcast',
                                      Landward == land_net_change ~'Match', 
                                      Landward != land_net_change ~ 'Mis-match')) %>%
    st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)

# map 

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
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
            main.title = 'E) Landward hindcast matches and mis-matches with optimal thresholds',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('palegreen4','lightgoldenrod','black' ,'red'), 
                labels =  c( 'Match', 'Ambiguous', 'Mis-match','No Hindcast'), border.alpha = 0, size = 0.3)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-hindcast_map_match_', press, '_', thresh, '.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
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
            main.title = 'D) Seaward hindcast matches and mis-matches with optimal thresholds',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('palegreen4','lightgoldenrod','black' ,'red'), 
                labels =  c( 'Match', 'Ambiguous', 'Mis-match','No Hindcast'), border.alpha = 0, size = 0.3)
smap
tmap_save(smap, paste0('outputs/maps/seaward-hindcast_map_match_', press, '_', thresh, '.png'), width = 5, height = 1, dpi = 5000)

# characterise mis-matches according to geomorphology, driver of loss, and marine ecoregion

preds_df <- preds %>% 
  st_join(meow) %>% 
  st_drop_geometry() %>% 
  left_join(drivers)

sea_eco <- preds_df %>% 
  filter(Seaward_match != 'Ambiguous') %>% 
  mutate(mismatch = ifelse(Seaward_match == 'Mis-match', 1, 0)) %>% 
  group_by(ECOREGION) %>% 
  summarise(percent_mismatch = sum(mismatch)/n()*100) %>% 
  mutate(label = ifelse(percent_mismatch > 80, ECOREGION, NA))
sea_sum <- preds_df %>% 
  filter(ECOREGION %in% unique(sea_eco$label)) %>% 
  group_by(ECOREGION, sea_net_change) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = 'sea_net_change', values_from = 'n') %>% 
  mutate(Loss = ifelse(is.na(Loss), 0, Loss)) %>% 
  group_by(ECOREGION) %>% 
  summarise(percent_gain = Gain_neutrality/(Gain_neutrality + Loss))
sea_sum
sea_eco_sf <- meow %>% inner_join(sea_eco)

sea_m <- tm_shape(st_crop(World, st_bbox(sea_eco_sf))) +
  tm_fill('grey88') +
  tm_shape(sea_eco_sf) +
  tm_fill('percent_mismatch', title = 'Percent mis-match', legend.is.portrait = F) +
  tm_text('label', col = 'black', size = 0.15, just = 'left') +
  tm_layout(legend.position = c(0.13, 0.01),
            main.title.size = 0.4,
            legend.title.size = 0.45,
            legend.text.size = 0.3,
            legend.width = 0.2,
            legend.height = 0.2,
            main.title = 'A) Seaward')
sea_m
tmap_save(sea_m, 'outputs/maps/seaward-mismatch.png', width = 5, height = 1, dpi = 1000)

land_eco <- preds_df %>% 
  filter(Landward_match != 'Ambiguous') %>% 
  mutate(mismatch = ifelse(Landward_match == 'Mis-match', 1, 0)) %>% 
  group_by(ECOREGION) %>% 
  summarise(percent_mismatch = sum(mismatch)/n()*100) %>% 
  mutate(label = ifelse(percent_mismatch > 80, ECOREGION, NA))
land_sum <- preds_df %>% 
  filter(ECOREGION %in% unique(land_eco$label)) %>% 
  group_by(ECOREGION, land_net_change) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = 'land_net_change', values_from = 'n') %>% 
  mutate(Loss = ifelse(is.na(Loss), 0, Loss)) %>% 
  group_by(ECOREGION) %>% 
  summarise(percent_gain = Gain_neutrality/(Gain_neutrality + Loss))
land_sum
land_eco_sf <- meow %>% inner_join(land_eco)

land_m <- tm_shape(st_crop(World, st_bbox(land_eco_sf))) +
  tm_fill('grey88') +
  tm_shape(land_eco_sf) +
  tm_fill('percent_mismatch', title = 'Percent mis-match', legend.is.portrait = F) +
  tm_text('label', col = 'black', size = 0.15, just = 'left') +
  tm_layout(legend.position = c(0.13, 0.01),
            main.title.size = 0.4,
            legend.title.size = 0.45,
            legend.text.size = 0.3,
            legend.width = 0.2,
            legend.height = 0.2,
            main.title = 'B) Landward')
land_m
tmap_save(land_m, 'outputs/maps/landward-mismatch.png', width = 5, height = 1, dpi = 5000)

sea_geo <- preds_df %>% 
  group_by(Seaward_match, Class) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n/nrow(preds_df)*100)
sea_geo

land_geo <- preds_df %>% 
  group_by(Landward_match, Class) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n/nrow(preds_df)*100)
land_geo

sea_driver <- preds_df %>% 
  pivot_longer(cols = Erosion:Settlement) %>% 
  group_by(Seaward_match, name) %>% 
  summarise(mean_prop = median(value))
sea_driver

land_driver <- preds_df %>% 
  pivot_longer(cols = Erosion:Settlement) %>% 
  group_by(Landward_match, name) %>% 
  summarise(mean_prop = median(value))
land_driver

# map ecoregions with highe percentages of mismatches below

sea_m <- tm_shape(st_crop(World, st_bbox(sea_eco_sf))) +
  tm_fill('grey88') +
  tm_shape(filter(sea_eco_sf, percent_mismatch > 80)) +
  tm_fill('black', legend.show = F, alpha = 0.3) +
  #tm_text('label', col = 'black', size = 0.15, just = 'left') +
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
            main.title = 'D) Seaward hindcast matches and mis-matches with optimal thresholds',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('palegreen4','lightgoldenrod','black' ,'red'), 
                labels =  c( 'Match', 'Ambiguous', 'Mis-match','No Hindcast'), border.alpha = 0, size = 0.3)
sea_m
tmap_save(sea_m, 'outputs/maps/seaward-mismatch_ecoregion.png', width = 5, height = 1, dpi = 5000)

land_m <- tm_shape(st_crop(World, st_bbox(land_eco_sf))) +
  tm_fill('grey88') +
  tm_shape(filter(land_eco_sf, percent_mismatch > 80)) +
  tm_fill('black', legend.show = F, alpha = 0.3) +
  #tm_text('label', col = 'black', size = 0.15, just = 'left') +
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
            main.title = 'E) Landward hindcast matches and mis-matches with optimal thresholds',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8) +
  tm_add_legend('symbol', col =  c('palegreen4','lightgoldenrod','black' ,'red'), 
                labels =  c( 'Match', 'Ambiguous', 'Mis-match','No Hindcast'), border.alpha = 0, size = 0.3)

land_m
tmap_save(land_m, 'outputs/maps/landward-mismatch_ecoregion.png', width = 5, height = 1, dpi = 5000)

