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

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model
chosen_model_name <- 'mangrove_model'

# import the final set of hindcast predictions, posterior probabilities, and posterior predictions for each
# biophysical/pressure scenario using optimal pressure definition and calibrated ambiguity threshold



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