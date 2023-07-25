# simulate models spatially, i.e., in each mangrove typological unit

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
dat <- read.csv('outputs/master-dat.csv')

# set simulation parameters

set.seed(123) # set random number generator to make results reproducible
numsims <- 1000 # number of model simulations to run

# loop over pressure definition thresholds
# takes 2.25 hours to run with 1000 sims

cl <- makeCluster(6)
registerDoParallel(cl)
system.time(
results <- foreach(k = seq_along(unique(dat$pressure_def)), 
                   .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar% {
  
  # read in and wrangle spatial data
  
  spatial_dat <- dat %>% 
    filter(pressure_def == k) %>% 
    mutate(csqueeze = recode(csqueeze, 'Medium' = 'M', 'High' = 'L', 'Low' = 'H'), # note counterintuitive notation here
           csqueeze_1 = ifelse(csqueeze == 'None', 0, 1), 
           fut_csqueeze = recode(fut_csqueeze, 'Medium' = 'M', 'High' = 'L', 'Low' = 'H'), # note counterintuitive notation here
           fut_csqueeze_1 = ifelse(fut_csqueeze == 'None', 0, 1), 
           sed_supp = recode(sed_supp, 'Medium' = 'M', 'Low' = 'L', 'High' = 'H'), 
           fut_dams = ifelse(fut_dams == 1, 'L', sed_supp), 
           prop_estab = recode(prop_estab, 'Medium' = 'M', 'High' = 'H', 'Low' = 'L'),
           Tidal_Class = recode(Tidal_Class, 'Micro' = 'H', 'Meso' = 'M', 'Macro' = 'L'))
  
  # run a model for each typological unit
  # define relative edge constraints - which edge interaction strengths are greater than other
  # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
  # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
  # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
  rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                  parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
  names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
  
  # iterate forecasts and hindcasts
  runs <- c('forecast', 'hindcast') # run both a forecast and a hindcast for each model
  
  tmp <- list() # tmp list to store results
  
  for(j in seq_along(runs)){
    run <- runs[[j]]
    tmp2 <- list() # tmp list to store results
    for(i in 1:nrow(spatial_dat)){
      
      if(run == 'forecast'){
        # perturbations for forecasting
        datselect <- select(spatial_dat[i,], fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
          pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
          filter(vals == 1) %>% 
          mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                                'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones'))
        
        if(nrow(datselect) == 0){ # if there are no perturbations, go to next typology
          next
        }
        
        press.scenario <- rep(1, nrow(datselect))
        names(press.scenario) <- datselect$press
        
        # edge constraint scenarios
        
        if(spatial_dat[i,]$fut_csqueeze == 'None'){
          datselect2 <- select(spatial_dat[i,], Tidal_Class, prop_estab) %>% 
            pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
          from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
          to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
          con.scenario <- c(datselect2$vals, datselect2$vals[2])
        }else{
          datselect2 <- select(spatial_dat[i,], Tidal_Class, prop_estab, fut_csqueeze) %>% 
            pivot_longer(Tidal_Class:fut_csqueeze, names_to = 'press', values_to = 'vals') 
          from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
          to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
          con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
        }
        
        # select model for sediment supply
        
        if(spatial_dat[i,]$fut_dams == 'L'){
          model <- rel.edge.cons.scenarios[[2]]
        }else if(spatial_dat[i,]$sed_supp == 'H'){
          model <- rel.edge.cons.scenarios[[1]]
        }else{
          model <- rel.edge.cons.scenarios[[2]]
        }
        
        # simulate outcomes
        
        sim <- system.sim_press(numsims, constrainedigraph = model, 
                                from = from_vec,
                                to = to_vec,
                                class = con.scenario,
                                perturb = press.scenario,
                                spatial = 'Y')
        
        out <- sim$stableoutcome %>% 
          filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
          group_by(var) %>% 
          summarise(Prob_gain = (sum(outcome>0)/n())*100,
                    Prob_neutral = (sum(outcome==0)/n())*100,
                    Prob_loss = (sum(outcome<0)/n())*-100)
        
        out$Type <- rep(spatial_dat[i, 'Type'], nrow(out))
        out$cast <- 'forecast'
        tmp2[[i]] <- out
      }else{
        
        # perturbations for hindcasting
        datselect <- select(spatial_dat[i,], csqueeze_1, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
          pivot_longer(csqueeze_1:storms, names_to = 'press', values_to = 'vals') %>% 
          filter(vals == 1) %>% 
          mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                                'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))

        if(nrow(datselect) == 0){ # if there are no perturbations, go to next typology
          next
        }
        
        press.scenario <- rep(1, nrow(datselect))
        names(press.scenario) <- datselect$press
        
        # edge constraint scenarios
        
        if(spatial_dat[i,]$csqueeze == 'None'){
          datselect2 <- select(spatial_dat[i,], Tidal_Class, prop_estab) %>% 
            pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
          from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
          to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
          con.scenario <- c(datselect2$vals, datselect2$vals[2])
        }else{
          datselect2 <- select(spatial_dat[i,], Tidal_Class, prop_estab, csqueeze) %>% 
            pivot_longer(Tidal_Class:csqueeze, names_to = 'press', values_to = 'vals') 
          from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
          to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
          con.scenario <- c(datselect2$vals, datselect2$vals[2], datselect2$vals[3])
          con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
        }
        
        # select model for sediment supply
        
        if(spatial_dat[i,]$sed_supp == 'H'){
          model <- rel.edge.cons.scenarios[[1]]
        }else{
          model <- rel.edge.cons.scenarios[[2]]
        }
        
        # simulate outcomes
        
        sim <- system.sim_press(numsims, constrainedigraph = model, 
                                from = from_vec,
                                to = to_vec,
                                class = con.scenario,
                                perturb = press.scenario,
                                spatial = 'Y')
        
        out <- sim$stableoutcome %>% 
          filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
          group_by(var) %>% 
          summarise(Prob_gain = (sum(outcome>0)/n())*100,
                    Prob_neutral = (sum(outcome==0)/n())*100,
                    Prob_loss = (sum(outcome<0)/n())*-100)
        
        out$Type <- rep(spatial_dat[i, 'Type'], nrow(out))
        out$cast <- 'hindcast'
        tmp2[[i]] <- out
      }
    }
    allout <- do.call(rbind, tmp2)
    tmp[[j]] <- allout
  }
  data.frame(pressure_def = k, do.call(rbind, tmp))
}
)
stopCluster(cl)

write.csv(results, paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv'), row.names = F)

