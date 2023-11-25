# calculate 95% confidence intervals for the number of units in each forecast class

library(tidyverse)
source('scripts/helpers/models.R')
press <- 4 # which pressure definition threshold?
thresh <- 75 # which ambiguity threshold?

# loop through available models

for(i in seq_along(names(models))){
  
  chosen_model_name <- names(models[i]) # which model do you want to run?
  
  accuracy <- read.csv(paste0('outputs/validation/resampled_accuracy_summary_', chosen_model_name, '.csv'))
  spatial_pred_fit <- read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_SeaLevelRise_', chosen_model_name, '_fit.csv'))
  num_units <- nrow(filter(spatial_pred_fit, !is.na(Landward))) # total number of units for which we could make forecasts for (some unable to forecast due to lack of valid model, i.e., all matrix likelihoods sum to 0)
  
  # summarise number of units in each forecast class globally
  
  # baseline forecast
  
  datsum <- spatial_pred_fit %>% 
    pivot_longer(cols = c(Landward, Seaward), names_to = 'mangrove', values_to = 'forecast') %>% 
    mutate(n = 1) %>% 
    group_by(mangrove, forecast) %>% 
    summarise(n = sum(n)) %>% 
    filter(!is.na(forecast)) %>% # remove missing forecasts
    mutate(n_lower95 = ifelse(mangrove == 'Landward' & forecast == 'Gain_neutrality', 
                              round(n - n*filter(accuracy, mangrove == 'Landward' & class == 'Gain_neutrality' & metric == 'Users_accuracy')$error_rate_low), NA),
           n_upper95 = ifelse(mangrove == 'Landward' & forecast == 'Gain_neutrality', 
                              round(n + n*filter(accuracy, mangrove == 'Landward' & class == 'Gain_neutrality' & metric == 'Producers_accuracy')$error_rate_upp), NA),
           n_lower95 = ifelse(mangrove == 'Landward' & forecast == 'Loss', 
                              round(n - n*filter(accuracy, mangrove == 'Landward' & class == 'Loss' & metric == 'Users_accuracy')$error_rate_low), n_lower95),
           n_upper95 = ifelse(mangrove == 'Landward' & forecast == 'Loss', 
                              round(n + n*filter(accuracy, mangrove == 'Landward' & class == 'Loss' & metric == 'Producers_accuracy')$error_rate_upp), n_upper95),
           n_lower95 = ifelse(mangrove == 'Seaward' & forecast == 'Gain_neutrality', 
                              round(n - n*filter(accuracy, mangrove == 'Seaward' & class == 'Gain_neutrality' & metric == 'Users_accuracy')$error_rate_low), n_lower95),
           n_upper95 = ifelse(mangrove == 'Seaward' & forecast == 'Gain_neutrality', 
                              round(n + n*filter(accuracy, mangrove == 'Seaward' & class == 'Gain_neutrality' & metric == 'Producers_accuracy')$error_rate_upp), n_upper95),
           n_lower95 = ifelse(mangrove == 'Seaward' & forecast == 'Loss', 
                              round(n - n*filter(accuracy, mangrove == 'Seaward' & class == 'Loss' & metric == 'Users_accuracy')$error_rate_low), n_lower95),
           n_upper95 = ifelse(mangrove == 'Seaward' & forecast == 'Loss', 
                              round(n + n*filter(accuracy, mangrove == 'Seaward' & class == 'Loss' & metric == 'Producers_accuracy')$error_rate_upp), n_upper95)) %>% 
    mutate(n_upper95 = ifelse(mangrove == 'Landward' & forecast == 'Ambiguous', 
                              num_units - sum(filter(., mangrove == 'Landward')$n_lower95, na.rm = T), n_upper95),
           n_lower95 = ifelse(mangrove == 'Landward' & forecast == 'Ambiguous', 
                              num_units - sum(filter(., mangrove == 'Landward')$n_upper95, na.rm = T), n_lower95),
           n_upper95 = ifelse(mangrove == 'Seaward' & forecast == 'Ambiguous', 
                              num_units - sum(filter(., mangrove == 'Seaward')$n_lower95, na.rm = T), n_upper95),
           n_lower95 = ifelse(mangrove == 'Seaward' & forecast == 'Ambiguous', 
                              num_units - sum(filter(., mangrove == 'Seaward')$n_upper95, na.rm = T), n_lower95)) %>% 
    mutate(n_lower95 = ifelse(n_lower95 < 0, 0, n_lower95),
           n_upper95 = ifelse(n_upper95 > num_units, num_units, n_upper95)) %>% 
    mutate_at(c('n', 'n_lower95', 'n_upper95'), function(x){round((x/num_units)*100)})
  datsum
  write.csv(datsum, paste0('outputs/summary-stats/baseline-forecast_uncertainty_', chosen_model_name, '.csv'), row.names = F)
  
}
