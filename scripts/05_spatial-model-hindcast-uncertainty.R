# quantify uncertainty in model hindcast accuracy using optimal thresholds
# randomly sample sites into 5 folds 200 times (which is 5*200 = 1000 test-training pairs, recommended as stable by Lyons et al. 2018)
# and obtain sampling distribution of each accuracy estimators (overall, producers, users)
# propagate uncertainty to estimates of number of units of loss, gain/neutrality, or ambiguity
# by taking the 95th percentiles of the commission and omission errors for each class, e.g.,
# num_gain_units(upper) = num_gain_units + (num_gain_units * omission_error(95 percentile)) # omission error is the rate of false negatives, so predictions of a class we are missing, so is the upper
# num_gain_units(lower) = num_gain_units - (num_gain_units * commission_error(95 percentile)) # commission error is the rate of false positives, so predictions of a class we shouldn't have, so is the lower

library(QPress)
library(tidyverse)
library(foreach)
library(doParallel)
library(sf)
source('scripts/helpers/models.R')
source('scripts/helpers/spatial-helpers_v2.R')
sf_use_s2(FALSE)
set.seed(123)

# set optimal thresholds
press <- 4 # which pressure definition threshold?
thresh <- 75 # which ambiguity threshold?

# read in and wrangle data
typ_points <- st_read('data/typologies/Mangrove_Typology_v3.14_Composite_valid_centroids.gpkg')
spatial_dat <- read.csv('data/master-dat.csv') %>% filter(pressure_def == press)

# randomly resample and split data into training and test folds
# loop through available models

for(i in seq_along(names(models))){
  
  chosen_model <- names(models[i]) # which model do you want to run?
  naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes_', chosen_model, '.csv'))
  nsamp <- 200 # number of random samples
  kfold <- 5 # number of folds
  tmp <- list()
  
  system.time( # approx 2 hours
    for(j in 1:nsamp){ 
      shuffled_dat <- spatial_dat %>% # shuffle the data and split
        left_join(data.frame(Type = .[1:length(unique(.$Type)),]$Type[sample(1:length(unique(.$Type)))],
                             k = c(rep(1:kfold, each = round(length(unique(.$Type))/kfold)),1:kfold)[1:length(unique(.$Type))]), by = 'Type') %>% 
        dplyr::select(pressure_def, k, Type, csqueeze, csqueeze_1, sed_supp, Tidal_Class, prop_estab, climate, cdev, #ant_slr, 
                      gwsub, hist_drought, hist_ext_rain, storms, land_net_change_obs, sea_net_change_obs) %>% 
        #mutate(no_press = ant_slr + gwsub + hist_drought + hist_ext_rain + storms + csqueeze_1) %>% 
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
      shuffled_dat <- shuffled_dat[,c('pressure_def', 'k', 'Type', 'csqueeze', 'sed_supp', 'Tidal_Class', 'prop_estab', 'climate', 'cdev', 'land_net_change_obs', 'sea_net_change_obs',
                                      'vals_gwsub', 'vals_hist_drought', 'vals_hist_ext_rain', 'vals_storms', #'vals_ant_slr', 
                                      'vals_csqueeze_1', 'press_no_press', 'press_gwsub',
                                      'press_hist_drought', 'press_hist_ext_rain', 'press_storms', #'press_ant_slr',
                                      'press_csqueeze_1', 'csqueeze_2', 'sed_supp_2', 'Tidal_Class_2', 'prop_estab_2', 'climate_2', 'cdev_2')]
      shuffled_dat <- shuffled_dat %>% 
        unite('scenario', csqueeze_2:cdev_2, na.rm = T, sep = '.') %>% 
        unite('press', press_gwsub:press_csqueeze_1, na.rm = T, sep = '.') %>% 
        mutate(press = ifelse(!is.na(press_no_press), 'none', press))
      
      # join naive outcomes/hindcasts to each unit in a training fold, validate against observed mangrove loss or gain
      # calculate likelihood/posterior probability of each matrix in a biophysical model based on observed mangrove loss or gain
      # use training posterior probabilities to make posterior hindcasts in test units and quantify accuracy
      # do this looping over kfolds
      
      cl <- makeCluster(5)
      registerDoParallel(cl)
      system.time(
        results <- foreach(i = 1:kfold, .packages = c('tidyverse', 'caret'), .errorhandling = 'remove') %dopar% {
          
          # get training units for a fold and pressure definition
          train_units <- shuffled_dat %>% filter(k != i)
          
          # calculate likelihood/posterior probability of each matrix in each bio model given observed loss or gain
          train_post_prob <- train_units %>% 
            left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
            mutate(valid = ifelse(land_net_change_obs == LandwardMang & sea_net_change_obs == SeawardMang, 1, 0)) %>% 
            group_by(scenario, nsim) %>% 
            summarise(matrix_post_prob = mean(valid)) 
          
          # get test units and make posterior predictions/hindcasts using posterior probability of each biomodel matrix
          test_units <- shuffled_dat %>% filter(k == i)
          
          # join naive hindcasts and posterior probabilities, calculate posterior hindcasts for test units
          test_post_hindcasts <- test_units %>% 
            left_join(naive_outcomes, by = c('scenario', 'press')) %>% 
            left_join(train_post_prob, by = c('scenario', 'nsim')) %>% 
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
                                       .default = 'Ambiguous'),
                   Change = case_when(LandwardMang >= thresh & SeawardMang >= thresh ~ "Gain",
                                      LandwardMang < 100-thresh & SeawardMang < 100-thresh ~ "Loss",
                                      Landward == 'Ambiguous' ~ "Ambiguous",
                                      Seaward == 'Ambiguous' ~ "Ambiguous",
                                      LandwardMang >= thresh & SeawardMang < 100-thresh ~ "Landward Gain & Seaward Loss",
                                      LandwardMang < 100-thresh & SeawardMang >= thresh ~ "Landward Loss & Seaward Gain",
                                      .default = NA),
                   Change_obs = case_when(land_net_change_obs == 1 & sea_net_change_obs == 1 ~ "Gain",
                                          land_net_change_obs == -1 & sea_net_change_obs == -1 ~ "Loss",
                                          land_net_change_obs == 1 & sea_net_change_obs == -1 ~ "Landward Gain & Seaward Loss",
                                          land_net_change_obs == -1 & sea_net_change_obs == 1 ~ "Landward Loss & Seaward Gain"),
                   land_net_change = ifelse(land_net_change_obs == -1, 'Loss', 'Gain_neutrality'),
                   sea_net_change = ifelse(sea_net_change_obs == -1, 'Loss', 'Gain_neutrality')) %>% 
            mutate(ambig_threshold = thresh)
          
          test_set <- test_post_hindcasts %>% 
            filter(Landward != 'Ambiguous' & Seaward != 'Ambiguous')
          
          results <- data.frame(mangrove = 'Landward',
                                resample = j,
                                k = i,
                                pressure_def = press,
                                ambig_threshold = thresh, 
                                calc_accuracy(test_set$Landward, test_set$land_net_change)$accuracy.results)
          results2 <- data.frame(mangrove = 'Seaward',
                                 resample = j,
                                 k = i,
                                 pressure_def = press,
                                 ambig_threshold = thresh, 
                                 calc_accuracy(test_set$Seaward, test_set$sea_net_change)$accuracy.results)
          results3 <- data.frame(mangrove = 'Seaward & Landward',
                                 resample = j,
                                 k = i,
                                 pressure_def = press,
                                 ambig_threshold = thresh, 
                                 calc_accuracy(test_set$Change, test_set$Change_obs)$accuracy.results)
          rbind(results, results2, results3)
        })
      stopCluster(cl)
      
      tmp[[j]] <- do.call(rbind, results)
    })
  accuracy <- do.call(rbind, tmp)
  write.csv(accuracy, paste0('outputs/validation/sampling_distribution_', chosen_model, '.csv'), row.names = F)
  
  # get summary stats for each metric
  
  accuracy_sum <- accuracy %>% 
    filter(mangrove != 'Seaward & Landward') %>% 
    pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
    mutate(class = ifelse(metric == 'Overall_accuracy', 'Gain_neutrality & Loss', class)) %>% 
    distinct() %>% 
    group_by(mangrove, class, metric) %>% 
    summarise(median = median(accuracy, na.rm = T),
              mean = mean(accuracy, na.rm = T),
              perc_0.025 = quantile(accuracy, 0.025),
              perc_0.975 = quantile(accuracy, 0.975),
              perc_0.05 =  quantile(accuracy, 0.05)) %>% 
    mutate(error_rate_upp = ifelse(metric == 'Producers_accuracy', (100-perc_0.05)/100, 0),
           error_rate_low = ifelse(metric == 'Users_accuracy', (100-perc_0.05)/100, 0))
  
  write.csv(accuracy_sum, paste0('outputs/validation/resampled_accuracy_summary_', chosen_model, '.csv'), row.names = F)
  
}
