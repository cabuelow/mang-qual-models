# training testing cross validation

library(QPress)
library(tidyverse)
library(scales)
library(foreach)
library(doParallel)
source('scripts/helpers/models_v2.R')
source('scripts/helpers/helpers_v2.R')

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model
chosen_model_name <- 'mangrove_model'
spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  mutate(land_net_change_obs = ifelse(land_net_change == 'Gain', 1, -1),
         sea_net_change_obs = ifelse(sea_net_change == 'Gain', 1, -1)) %>% 
  filter(pressure_def == 3) %>% 
  mutate(csqueeze = recode(csqueeze, 'Medium' = 'M', 'High' = 'L', 'Low' = 'H'), # note counterintuitive notation here
         csqueeze_1 = ifelse(csqueeze == 'None', 0, 1), 
         fut_csqueeze = recode(fut_csqueeze, 'Medium' = 'M', 'High' = 'L', 'Low' = 'H'), # note counterintuitive notation here
         fut_csqueeze_1 = ifelse(fut_csqueeze == 'None', 0, 1), 
         sed_supp = recode(sed_supp, 'Medium' = 'M', 'Low' = 'L', 'High' = 'H'), 
         fut_dams = ifelse(fut_dams == 1, 'L', sed_supp), 
         prop_estab = recode(prop_estab, 'Medium' = 'M', 'High' = 'H', 'Low' = 'L'),
         Tidal_Class = recode(Tidal_Class, 'Micro' = 'H', 'Meso' = 'M', 'Macro' = 'L'))

# set simulation parameters

set.seed(123) # set random number generator to make results reproducible
numsims <- 100 # number of model simulations to run

# define relative edge constraints - which edge interaction strengths are greater than other
# in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
# under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
# the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 

# hindcast cross-validation
# randomly shuffle the units and split into 5 groups
# each group is used as a test set, and model is trained on remaining 75% of data
# store pressures and interaction strengths to summarise as range of interaction coefficients 
# make hindcasts using information on pressure-interaction strengths from the trained model and quantify accuracy

kfold <- 5 # number of folds
shuffled_dat <- spatial_dat[sample(1:nrow(spatial_dat)),] %>% 
  mutate(k = rep(1:5, each = nrow(spatial_dat)/kfold))

# loop through each fold and train/test model

cl <- makeCluster(5)
registerDoParallel(cl)

#tmpout <- list()
#tmpweights <- list()
system.time(
results <- foreach(i = 1:kfold, .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar% {
  
  train_dat <- shuffled_dat %>% filter(k != i)
  
  # loop through training sites, simulate and validate against observations, only keeping predictions that are valid
  tmpout2 <- list()
  tmpweights2 <- list()
  for(j in 1:nrow(train_dat)){
    
    datselect <- select(train_dat[j,], csqueeze_1, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
      pivot_longer(csqueeze_1:storms, names_to = 'press', values_to = 'vals') %>% 
      filter(vals == 1) %>% 
      mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                            'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
    datmonitor <- select(train_dat[j,], land_net_change_obs, sea_net_change_obs) %>% 
      pivot_longer(land_net_change_obs:sea_net_change_obs, names_to = 'monitor', values_to = 'vals') %>% 
      mutate(monitor = recode(monitor, 'land_net_change_obs' = 'LandwardMang', 'sea_net_change_obs' = "SeawardMang"))
    if(nrow(datselect) == 0){ # if there are no perturbations, go to next typology
      next
    }
    press.scenario <- rep(1, nrow(datselect))
    names(press.scenario) <- datselect$press
    monitor.scenario <- datmonitor$vals
    names(monitor.scenario) <- datmonitor$monitor
    
    # edge constraint scenarios
    
    if(train_dat[j,]$csqueeze == 'None'){
      datselect2 <- select(train_dat[j,], Tidal_Class, prop_estab) %>% 
        pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
      from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
      to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
      con.scenario <- c(datselect2$vals, datselect2$vals[2])
    }else{
      datselect2 <- select(train_dat[j,], Tidal_Class, prop_estab, csqueeze) %>% 
        pivot_longer(Tidal_Class:csqueeze, names_to = 'press', values_to = 'vals') 
      from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
      to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
      con.scenario <- c(datselect2$vals, datselect2$vals[2], datselect2$vals[3])
      con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
    }
    
    # select model for sediment supply
    
    if(train_dat[j,]$sed_supp == 'H'){
      model <- rel.edge.cons.scenarios[[1]]
    }else{
      model <- rel.edge.cons.scenarios[[2]]
    }
    
    # simulate outcomes
    
    sim <- system.sim_press_valid(numsims, constrainedigraph = model, 
                            from = from_vec,
                            to = to_vec,
                            class = con.scenario,
                            perturb = press.scenario,
                            monitor = monitor.scenario,
                            spatial = 'Y')
    
    tmpout2[[j]] <- sim$stableoutcome %>% 
      filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
      group_by(valid, pressures, var) %>% 
      summarise(Prob_gain = (sum(outcome>0)/n())*100,
                Prob_neutral = (sum(outcome==0)/n())*100,
                Prob_loss = (sum(outcome<0)/n())*-100) %>% 
      mutate(Type = train_dat[j, 'Type'],
             cast = 'training',
             kfold = i)
    
    tmpweights2[[j]] <- sim$stableweights %>% 
      mutate(Type = train_dat[j, 'Type'],
             cast = 'training',
             kfold = i)
  }
  #tmpout[[i]] <- do.call(rbind, tmpout2)
  data.frame(do.call(rbind, tmpweights2))
}
)

stopCluster(cl)
write.csv(results, paste0('outputs/validation/training_weights_', chosen_model_name, '_spatial.csv'), row.names = F)

# identify scenarios and locations that are invalid - i.e., model cannot predict for

# loop through the test folds and make hindcasts using weight ranges from relevant valid pressure scenarios
# quantify accuracy of outcomes




