# calculate 95% confidence intervals around the number of units in each forecast class

fit_summary <- read.csv(paste0('outputs/predictions/forecast-predictions', go, '_', rm_e, '_', press, '_', thresh, '_', 'SeaLevelRise_summary_fit.csv'))
accuracy <- read.csv('outputs/validation/resampled_accuracy_summary.csv')

fit_summary <- fit_summary %>% 
  mutate(n =)



spatial_pred_fit <- read.csv(paste0('outputs/predictions/forecast-predictions', go, '_', rm_e, '_', press, '_', thresh, '_', scenarios[[i]][1], '_fit.csv'))

# summarise predictions
datsum <- spatial_pred_fit %>% 
  mutate(Landward = paste0('Landward_', .$Landward),
         Seaward = paste0('Seaward_', .$Seaward)) %>% 
  mutate(Landward_seaward = paste0(.$Landward, '.', .$Seaward)) %>% 
  group_by(Landward_seaward) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
write.csv(datsum, paste0('outputs/predictions/forecast-predictions', go, '_', rm_e, '_', press, '_', thresh, '_', scenarios[[i]][1], '_summary_fit.csv'), row.names = F)    

datsum <- spatial_pred_unfit %>% 
  mutate(Landward = paste0('Landward_', .$Landward),
         Seaward = paste0('Seaward_', .$Seaward)) %>% 
  mutate(Landward_seaward = paste0(.$Landward, '.', .$Seaward)) %>% 
  group_by(Landward_seaward) %>% 
  summarise(n = n(), percent = 100*(n()/nrow(.)))
write.csv(datsum, paste0('outputs/predictions/forecast-predictions', go, '_', rm_e, '_', press, '_', thresh, '_', scenarios[[i]][1], '_summary_unfit.csv'), row.names = F) 
