# hindcast cross validation
# randomly shuffle and split mangrove typological units into training and test folds
# loop through splits and use training data to calibrate model
# calibration is 2 steps, a. threshold calibration, b. parameter calibration
# then make hindcasts using test data and calibrated thresholds/parameters

library(QPress)
library(tidyverse)
library(scales)
library(caret)
library(ggh4x)
library(foreach)
library(doParallel)
source('scripts/helpers/models.R')
source('scripts/helpers/helpers.R')
source('scripts/helpers/spatial-helpers.R')
set.seed(123) # set random number generator to make results reproducible

# read in spatial data (mangrove typological units)

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

nsim <- 2 # number of sims
kfold <- 5 # number of folds
shuffled_dat <- spatial_dat %>% 
  left_join(data.frame(Type = spatial_dat[1:length(unique(spatial_dat$Type)),]$Type[sample(1:length(unique(spatial_dat$Type)))],
                           k = rep(1:kfold, each = length(unique(spatial_dat$Type))/kfold)), by = 'Type') %>% 
  mutate(k_press = paste0(.$k, '_', .$pressure_def))

# use training folds to obtain hindcasts for range of pressure and ambiguity threshold definitions
cl <- makeCluster(6)
registerDoParallel(cl)

system.time(
train_press <- foreach(i = 1:kfold, .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar%{
  tmp <- list() # temp list for storing results
  for(j in seq_along(unique(shuffled_dat$pressure_def))){ # loop over pressure definitions
    train_dat <- shuffled_dat %>% filter(k != i & pressure_def == j)
    tmp2 <- list()
    for(k in 1:nrow(train_dat)){ #TODO: not sure why apply won't work here over rows instead of forloop
      tmp2[[k]] <- cast(train_dat[k,], numsims = nsim, type = 'hindcast')
    }
    tmp[[j]] <- data.frame(k = i, pressure_def = j, do.call(rbind, tmp2))
  }
  do.call(rbind, tmp)
}
)
stopCluster(cl)
write.csv(train_press, paste0('outputs/validation/hindcast-pressure-outcomes_', chosen_model_name, '.csv'), row.names = F)

# quantify accuracy of hindcasts across range of pressure and ambiguity defnitions and identify optimal thresholds

train_press <- train_press %>%
  mutate(Prob_gain_neutrality = Prob_gain + Prob_neutral)

threshold <- seq(60, 90, by = 5) # range of thresholds for defining when a prediction is ambiguous or not

# classify outcomes - gain, loss, ambiguous - for each ambiguity threshold

tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  tmp[[i]] <- train_press %>% 
    mutate(k_press = paste0(.$k, '_', .$pressure_def)) %>% 
    mutate(Change = case_when(Prob_gain_neutrality > thresh ~ 'Gain_neutrality',
                              Prob_loss < -thresh ~ 'Loss',
                              .default = 'Ambiguous')) %>%
    pivot_wider(id_cols = c('Type', 'k', 'pressure_def'), names_from = 'var', values_from = 'Change', names_prefix = 'Change_') %>% 
    left_join(select(shuffled_dat, Type, k, pressure_def, land_net_change, sea_net_change)) %>% 
    mutate(ambig_threshold = thresh)
}
class <- do.call(rbind, tmp)

# split into land vs. sea and remove ambiguous responses and units where commodities and erosion are dominant drivers of loss, can't validate with our model

land_validate <- class %>% 
  left_join(drivers, by = 'Type') %>% 
  filter(Change_LandwardMang != 'Ambiguous' & Commodities < 0.1 & Erosion < 0.1) 
sea_validate <- class %>%
  left_join(drivers, by = 'Type') %>% 
  filter(Change_SeawardMang != 'Ambiguous' & Commodities < 0.1 & Erosion < 0.1) 

# net change hindcast accuracy across different pressure and ambiguity thresholds
iter <- unique(shuffled_dat$k_press)
tmp2 <- list()
for(j in seq_along(unique(shuffled_dat$k_press))){
  land <- land_validate %>% filter(k_press == iter[j])
  sea <- sea_validate %>% filter(k_press == iter[j])
  tmp <- list()
  for(i in seq_along(threshold)){
    thresh <- threshold[i]
    landval <- land %>% filter(ambig_threshold == thresh)
    seaval <- sea %>% filter(ambig_threshold == thresh)
    results <- data.frame(validation = 'net',
                          mangrove = 'Landward',
                          k = landval[1,'k'],
                          pressure_def = landval[1,'pressure_def'],
                          ambig_threshold = thresh, 
                          calc_accuracy(landval$Change_LandwardMang, landval$land_net_change)$accuracy.results)
    results2 <- data.frame(validation = 'net',
                           mangrove = 'Seaward',
                           k = seaval[1,'k'],
                           pressure_def = seaval[1,'pressure_def'],
                           ambig_threshold = thresh, 
                           calc_accuracy(seaval$Change_SeawardMang, seaval$sea_net_change)$accuracy.results)
    tmp[[i]] <- rbind(results, results2)
  }
  tmp2[[j]] <- do.call(rbind, tmp)
}
accuracy <- do.call(rbind, tmp2)
write.csv(accuracy, 'outputs/validation/hindcast-accuracy_pressure-ambiguity_kfold.csv', row.names = F)

# heatmap of accuracy metrics for combinations of pressure and ambiguity thresholds

accuracy %>% 
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

ggsave('outputs/validation/accuracy-heatmap_kfold.png', width = 10, height = 4.5)

# use optimal thresholds and training set to make hindcasts and obtain valid interaction coefficients by comparing to historical observations
# pressure definition is based on optimal for landward above, taking average from across training fold and rounding down if needed

optimal_press <- round(mean(filter(accuracy, mangrove == 'Landward' & k != i, Overall_accuracy == max(filter(accuracy, mangrove == 'Landward' & k != i)$Overall_accuracy))$pressure_def))
cl <- makeCluster(6)
registerDoParallel(cl)
system.time(
  train_params <- foreach(i = 1:kfold, .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar%{
      train_dat <- shuffled_dat %>% filter(k != i & pressure_def == optimal_press)
      tmp <- list()
      for(k in 1:nrow(train_dat)){ #TODO: not sure why apply won't work here over rows instead of forloop
        tmp[[k]] <- train(train_dat[k,], numsims = nsim)
      }
      data.frame(k = i, pressure_def = optimal_press, do.call(rbind, tmp))
    }
)
stopCluster(cl)
write.csv(train_params, paste0('outputs/validation/training_weights_', chosen_model_name, '.csv'), row.names = F)

# use calibrated interaction coefficients and test fold to make independent hindcasts and quantify accuracy
# select valid interaction coefficients based on geomorphology and pressure presence

cl <- makeCluster(5)
registerDoParallel(cl)
system.time(
  test_results <- foreach(i = 1:kfold, .combine = rbind, .packages = c('QPress', 'tidyverse', 'scales')) %dopar% {
    
    test_dat <- shuffled_dat %>% filter(k == i, pressure_def == optimal_press)
    params <- train_params %>% # here am summarising over all valid simulations, regardless of observed outcome
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
    tmp[[j]] <- test(test_dat[j,], nsim, params)
    }
    data.frame(k = i, pressure_def = optimal_press, do.call(rbind, tmp))
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


