# randomly shuffle and split mangrove typological units into training and test folds
# loop through splits and use training data to calibrate model
# calibration is 2 steps, a. threshold calibration, b. parameter calibration
# then make hindcasts using test data and calibrated thresholds/parameters

library(QPress)
library(tidyverse)
library(scales)
library(foreach)
library(doParallel)
source('scripts/helpers/models_v2.R')
source('scripts/helpers/helpers_v2.R')
source('scripts/helpers/validation.R')
set.seed(123) # set random number generator to make results reproducible

# read in spatial data (mangrove typological units)

spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  mutate(land_net_change_obs = ifelse(land_net_change == 'Gain', 1, -1),
         sea_net_change_obs = ifelse(sea_net_change == 'Gain', 1, -1)) %>% 
  filter(pressure_def == 4) %>% 
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
drivers <- read.csv('data/typologies/SLR_Data.csv')

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model
chosen_model_name <- 'mangrove_model'

# define model relative edge constraints - which edge interaction strengths are greater than other
# in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
# under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
# the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 

# set up training vs. test folds

kfold <- 5 # number of folds
shuffled_dat <- spatial_dat[sample(1:nrow(spatial_dat)),] %>% 
  mutate(k = rep(1:kfold, each = nrow(spatial_dat)/kfold))

# loop through each fold

# use training folds to obtain hindcasts for range of pressure and ambiguity threshold definitions

# quantify accuracy of hindcasts and identify optimal thresholds

# use optimal thresholds and training set to make hindcasts and obtain valid interaction coefficients by comparing to historical observations

# use calibrated interaction coefficients and test fold to make independent hindcasts and quantify accuracy



