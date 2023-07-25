# training testing cross validation

library(QPress)
library(tidyverse)
library(scales)
library(foreach)
library(doParallel)
library(caret)
library(ggh4x)
source('scripts/helpers/models_v2.R')
source('scripts/helpers/helpers_v2.R')

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model
chosen_model_name <- 'mangrove_model'
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

# set simulation parameters

set.seed(123) # set random number generator to make results reproducible
numsims <- 100 # number of model simulations to run - takes 43 mins to do 100 in parallel
#TODO: to run 1000 sims takes up too much memory on laptop - any way to reduce results size?
# 7.4 GB results file from just 100 sims, if we bump up to 1000 wil by >70 GBs
# 1000 will take ~7.5 hrs to run in parallel

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
  mutate(k = rep(1:kfold, each = nrow(spatial_dat)/kfold))

# loop through each fold and train/test model

cl <- makeCluster(5)
registerDoParallel(cl)

system.time(
train_results <- foreach(i = 1:kfold, .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar% {
  
  train_dat <- shuffled_dat %>% filter(k != i)
  
  # loop through training sites, simulate and validate against observations, only keeping predictions that are valid
  tmpweights <- list() # list for storing weights
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
    
    outcomes <- sim$stableoutcome %>% 
      filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
      group_by(valid, pressures, var) %>% 
      summarise(Prob_gain = (sum(outcome>0)/n())*100,
                Prob_neutral = (sum(outcome==0)/n())*100,
                Prob_loss = (sum(outcome<0)/n())*-100) %>% 
      pivot_wider(id_cols = c('pressures'), names_from = 'var', values_from = c('Prob_gain', 'Prob_neutral', 'Prob_loss'))
    
    tmpweights[[j]] <- data.frame(sim$stableweights %>% 
      mutate(Type = train_dat[j, 'Type'],
             cast = 'training',
             kfold = i, 
             ), outcomes[,-1], 
      pivot_wider(datmonitor, names_from = 'monitor', values_from = 'vals', names_prefix = 'Observed_'))
  }
  data.frame(do.call(rbind, tmpweights))
}
)

stopCluster(cl)
write.csv(train_results, paste0('outputs/validation/training_weights_', chosen_model_name, '_spatial.csv'), row.names = F)
#train_results <- read.csv(paste0('outputs/validation/training_weights_', chosen_model_name, '_spatial.csv'))

# identify pressure combinations that are invalid - i.e., model cannot predict for - make these ambiguous?

nrow(filter(train_results, valid == 'valid'))/nrow(train_results) # what percentage of simulations are valid
invalid <- train_results %>% filter(valid == 'invalid') %>% 
  select(Type, pressures, Prob_gain_LandwardMang:Observed_SeawardMang) %>% 
  distinct()
length(unique(invalid$pressures)) # number of unique invalid pressure combinations
length(unique(invalid$Type)) # number of unique units with invalid predictions

# summarise interaction strength mean and 95% confidence interval under valid pressure and observed outcome combinations

param_sum <- train_results %>% 
  filter(valid == 'valid') %>%
  mutate(Geomorphology = sub("\\_.*", "", Type)) %>% 
  group_by(pressures, Geomorphology, param, Observed_LandwardMang, Observed_SeawardMang) %>% 
  summarise(weight_mean = mean(weight),
            weight_sd = sd(weight),
            weight_se = (sd(weight)/sqrt(n()))*1.96) %>% 
  mutate(Observed = case_when(Observed_LandwardMang == -1 & Observed_SeawardMang == -1 ~ 'Loss',
                   Observed_LandwardMang == 1 & Observed_SeawardMang == 1 ~ 'Gain',
                   Observed_LandwardMang == -1 & Observed_SeawardMang == 1 ~ 'Landward Loss & Seaward Gain',
                   Observed_LandwardMang == 1 & Observed_SeawardMang == -1 ~ 'Landward Gain & Seaward Loss'))

ggplot(param_sum) +
  geom_point(aes(x = weight_mean, y = pressures, col = Observed)) +
  geom_segment(aes(xend = weight_mean + weight_sd, x = weight_mean - weight_sd, y = pressures, yend = pressures, col = Observed), alpha = 0.5) +
  #facet_wrap(~param, nrow = 3) +
  facet_nested_wrap(~factor(Geomorphology) + factor(param), nrow = 4) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  ylab('') +
  xlab('') +
  theme_classic()

ggsave('outputs/validation/valid-parameter-weights.png', width = 30, height = 10)

# loop through the test sites and make hindcasts using weight ranges from relevant valid pressure scenarios and geomorphology

cl <- makeCluster(5)
registerDoParallel(cl)

system.time(
test_results <- foreach(i = 1:kfold, .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar% {
    
test_dat <- shuffled_dat %>% filter(k == i)
params <- train_results %>% # here am summarising over all valid simulations, regardless of observed outcome
  filter(valid == 'valid' & kfold != i) %>% # only use weights from training data 
  mutate(Geomorphology = sub("\\_.*", "", Type)) %>% 
  group_by(Geomorphology, pressures, param) %>% 
  summarise(weight_mean = mean(weight),
            weight_upp = mean(weight) + sd(weight),
            weight_low = mean(weight) - sd(weight)) %>% 
  mutate(weight_upp = ifelse(weight_mean < 0 & weight_upp > 0, 0, weight_upp),
         weight_low = ifelse(weight_mean > 0 & weight_low < 0, 0, weight_low))

# loop through test sites
tmp <- list() # list for storing outcomes
for(j in 1:nrow(test_dat)){
  
  datselect <- dplyr::select(test_dat[j,], csqueeze_1, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
    pivot_longer(csqueeze_1:storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                          'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
  if(nrow(datselect) == 0){ # if there are no perturbations, go to next
    next
  }
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario), collapse = '_') & Geomorphology == test_dat[j,]$Geomorphology)
  if(nrow(validweights) == 0){ # if there are no valid weights for this combination of pressures, go to next
    next
  }
  
  # edge constraint scenarios
  
  if(test_dat[j,]$csqueeze == 'None'){
    datselect2 <- dplyr::select(test_dat[j,], Tidal_Class, prop_estab) %>% 
      pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
    from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
    to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
    con.scenario <- c(datselect2$vals, datselect2$vals[2])
  }else{
    datselect2 <- dplyr::select(test_dat[j,], Tidal_Class, prop_estab, csqueeze) %>% 
      pivot_longer(Tidal_Class:csqueeze, names_to = 'press', values_to = 'vals') 
    from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
    to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
    con.scenario <- c(datselect2$vals, datselect2$vals[2], datselect2$vals[3])
    con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
  }
  
  # select model for sediment supply
  
  if(test_dat[j,]$sed_supp == 'H'){
    model <- rel.edge.cons.scenarios[[1]]
  }else{
    model <- rel.edge.cons.scenarios[[2]]
  }
  
  # simulate outcomes
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                                from = from_vec,
                                to = to_vec,
                                class = con.scenario,
                                perturb = press.scenario,
                                weights = validweights,
                                spatial = 'Y')
  
  tmp[[j]] <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = test_dat[j, 'Type'],
           kfold = i)
}
data.frame(do.call(rbind, tmp))
  }
)

stopCluster(cl)
write.csv(test_results, paste0('outputs/validation/test_outcomes_', chosen_model_name, '_spatial.csv'), row.names = F)

# classify outcomes and quanitfy accuracy

threshold <- seq(60, 90, by = 5) # range of thresholds for defining when a prediction is ambiguous or not

tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  tmp[[i]] <- test_results %>% 
    mutate(Prob_gain_neutrality = Prob_gain + Prob_loss) %>% 
    mutate(Change = case_when(Prob_gain_neutrality > thresh ~ 'Gain_neutrality',
                              Prob_loss < -thresh ~ 'Loss',
                              .default = 'Ambiguous')) %>%
    pivot_wider(id_cols = c('Type', 'kfold'), names_from = 'var', values_from = 'Change', names_prefix = 'Change_') %>% 
    inner_join(dplyr::select(spatial_dat, Type, pressure_def, land_net_change, sea_net_change), by = c('Type')) %>% 
    mutate(ambig_threshold = thresh)
}
test_class <- do.call(rbind, tmp)
write.csv(test_class, 'outputs/validation/test-crossvalidation-results.csv', row.names = F)

# split into land vs. sea and remove ambiguous responses and units where commodities and erosion are dominant drivers of loss, can't validate with our model
land_validate <- test_class %>% 
  left_join(drivers, by = 'Type') %>% 
  filter(Change_LandwardMang != 'Ambiguous' & Commodities < 0.1 & Erosion < 0.1) 
sea_validate <- test_class %>%
  left_join(drivers, by = 'Type') %>% 
  filter(Change_SeawardMang != 'Ambiguous' & Commodities < 0.1 & Erosion < 0.1) 

# net change hindcast accuracy
tmp2 <- list()
for(j in seq_along(unique(class$kfold))){
  land <- land_validate %>% filter(kfold == j)
  sea <- sea_validate %>% filter(kfold == j)
tmp <- list()
for(i in seq_along(threshold)){
    thresh <- threshold[i]
    landval <- land %>% filter(ambig_threshold == thresh)
    seaval <- sea %>% filter(ambig_threshold == thresh)
    results <- data.frame(validation = 'net',
                          mangrove = 'Landward',
                          ambig_threshold = thresh, 
                          calc_accuracy(landval$Change_LandwardMang, landval$land_net_change)$accuracy.results)
    results2 <- data.frame(validation = 'net',
                           mangrove = 'Seaward',
                           ambig_threshold = thresh, 
                           calc_accuracy(seaval$Change_SeawardMang, seaval$sea_net_change)$accuracy.results)
    tmp[[i]] <- rbind(results, results2)
}
tmp2[[j]] <- data.frame(kfold = j, do.call(rbind, tmp))
}

accuracy <- do.call(rbind, tmp2) %>%  mutate(mangrove = factor(mangrove, levels = c('Landward', 'Seaward')))
write.csv(accuracy, 'outputs/validation/cross-val-accuracy.csv', row.names = F)

# heatmap of accuracy metrics for combinations of pressure and ambiguity thresholds

accuracy %>% 
  filter(validation == 'net') %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
  distinct() %>% 
  ggplot() +
  aes(x = ambig_threshold, y = mangrove, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(kfold) + factor(class) + factor(metric), nrow = 5) +
  ylab('') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/cross-validation-accuracy-heatmap.png', width = 10, height = 8)

# forecast with calibrated paramater weights
    
params <- train_results %>% # here am summarising over all valid simulations, regardless of observed outcome
      filter(valid == 'valid') %>% # only use weights from training data 
      mutate(Geomorphology = sub("\\_.*", "", Type)) %>% 
      group_by(Geomorphology, pressures, param) %>% 
      summarise(weight_mean = mean(weight),
                weight_upp = mean(weight) + sd(weight),
                weight_low = mean(weight) - sd(weight)) %>% 
      mutate(weight_upp = ifelse(weight_mean < 0 & weight_upp > 0, 0, weight_upp),
             weight_low = ifelse(weight_mean > 0 & weight_low < 0, 0, weight_low))
    
# forecast

tmp <- list() # list for storing outcomes
for(i in 1:nrow(spatial_dat)){
  
  datselect <- select(spatial_dat[i,], fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
    pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                          'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones'))
  if(nrow(datselect) == 0){ # if there are no perturbations, go to next
    next
  }
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario), collapse = '_') & Geomorphology == spatial_dat[i,]$Geomorphology)
  if(nrow(validweights) == 0){ # if there are no valid weights for this combination of pressures, go to next
    next
  }
  
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
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                           from = from_vec,
                           to = to_vec,
                           class = con.scenario,
                           perturb = press.scenario,
                           weights = validweights,
                           spatial = 'Y')
  
  tmp[[i]] <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = spatial_dat[i, 'Type'])
}

forecast <- data.frame(do.call(rbind, tmp))
write.csv(forecast, paste0('outputs/validation/calibrated_forecast', chosen_model_name, '_spatial.csv'), row.names = F)

# forecast with calibrated paramater weights and restoration (increased sediment or increased propagules)

params <- train_results %>% # here am summarising over all valid simulations, regardless of observed outcome
  filter(valid == 'valid') %>% # only use weights from training data 
  mutate(Geomorphology = sub("\\_.*", "", Type)) %>% 
  group_by(Geomorphology, pressures, param) %>% 
  summarise(weight_mean = mean(weight),
            weight_upp = mean(weight) + sd(weight),
            weight_low = mean(weight) - sd(weight)) %>% 
  mutate(weight_upp = ifelse(weight_mean < 0 & weight_upp > 0, 0, weight_upp),
         weight_low = ifelse(weight_mean > 0 & weight_low < 0, 0, weight_low))

# forecast with increased sediment

tmp <- list() # list for storing outcomes
for(i in 1:nrow(spatial_dat)){
  
  datselect <- select(spatial_dat[i,], fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
    pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                          'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones')) %>% 
    rbind(data.frame(press = 'Sediment', vals = 1))
  if(nrow(datselect) == 0){ # if there are no perturbations, go to next
    next
  }
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario[-length(press.scenario)]), collapse = '_') & Geomorphology == spatial_dat[i,]$Geomorphology)
  if(nrow(validweights) == 0){ # if there are no valid weights for this combination of pressures, go to next
    next
  }
  
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
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                           from = from_vec,
                           to = to_vec,
                           class = con.scenario,
                           perturb = press.scenario,
                           weights = validweights,
                           spatial = 'Y')
  
  tmp[[i]] <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = spatial_dat[i, 'Type'])
}

forecast_restore <- data.frame(do.call(rbind, tmp))
write.csv(forecast_restore, paste0('outputs/validation/calibrated_forecast', chosen_model_name, '_spatial_sedrestore.csv'), row.names = F)

# forecast with increased propagules

tmp <- list() # list for storing outcomes
for(i in 1:nrow(spatial_dat)){
  
  datselect <- select(spatial_dat[i,], fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
    pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                          'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones')) %>% 
    rbind(data.frame(press = c('LandwardAvailableProp', 'SeawardAvailableProp'), vals = 1))
  if(nrow(datselect) == 0){ # if there are no perturbations, go to next
    next
  }
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario[-c(length(press.scenario)-1, length(press.scenario))]), collapse = '_') & Geomorphology == spatial_dat[i,]$Geomorphology)
  if(nrow(validweights) == 0){ # if there are no valid weights for this combination of pressures, go to next
    next
  }
  
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
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                           from = from_vec,
                           to = to_vec,
                           class = con.scenario,
                           perturb = press.scenario,
                           weights = validweights,
                           spatial = 'Y')
  
  tmp[[i]] <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = spatial_dat[i, 'Type'])
}

forecast_restoreprop <- data.frame(do.call(rbind, tmp))
write.csv(forecast_restoreprop, paste0('outputs/validation/calibrated_forecast', chosen_model_name, '_spatial_proprestore.csv'), row.names = F)
