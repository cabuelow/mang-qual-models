# compare hindcast probability of mangrove loss and gain with historical loss and gain to validate network model

library(tidyverse)
library(caret)
library(sf)
source('scripts/helpers/models_v2.R')

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')

# wrangle historical SRS observations of mangrove loss and gain into categories of change (no change, loss, gain, loss and gain)

spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  mutate(sea_gross_gain_loss = sea_gross_gain + sea_gross_loss,
        land_gross_gain_loss = land_gross_gain + land_gross_loss) %>% 
  mutate(sea_change_obs = ifelse(sea_gross_gain_loss == 2, 'Loss & Gain', NA),
         sea_change_obs = ifelse(sea_gross_gain_loss != 2 & sea_gross_gain ==1, 'Gain', sea_change_obs),
         sea_change_obs = ifelse(sea_gross_gain_loss != 2 & sea_gross_loss ==1, 'Loss', sea_change_obs),
         sea_change_obs = ifelse(sea_gross_gain_loss == 0, 'Gain', sea_change_obs), # here treating neutrality as gain, as in network model
         land_change_obs = ifelse(land_gross_gain_loss == 2, 'Loss & Gain', NA),
         land_change_obs = ifelse(land_gross_gain_loss != 2 & land_gross_gain ==1, 'Gain', land_change_obs),
         land_change_obs = ifelse(land_gross_gain_loss != 2 & land_gross_loss ==1, 'Loss', land_change_obs),
         land_change_obs = ifelse(land_gross_gain_loss == 0, 'Gain', land_change_obs)) %>%  # here treating neutrality as gain, as in network model 
  mutate(sea_gain_obs = ifelse(sea_gross_gain == 1, 'Gain', 'No Gain'),
         sea_loss_obs = ifelse(sea_gross_loss == 1, 'Loss', 'No Loss'),
         land_gain_obs = ifelse(land_gross_gain == 1, 'Gain', 'No Gain'),
         land_loss_obs = ifelse(land_gross_loss == 1, 'Loss', 'No Loss'))      
   #sea_gross_loss = ifelse(sea_gross_gain_loss == 2, 0, sea_gross_loss),
         #sea_gross_gain = ifelse(sea_gross_gain_loss == 2, 0, sea_gross_gain),
         #land_gross_loss = ifelse(land_gross_gain_loss == 2, 0, land_gross_loss),
         #land_gross_gain = ifelse(land_gross_gain_loss == 2, 0, land_gross_gain)) %>% 
  #mutate(sea_gross_gain_loss = ifelse(sea_gross_gain_loss == 2, 1, 0),
   #      land_gross_gain_loss = ifelse(land_gross_gain_loss == 2, 1, 0))

# which model outcomes to validate? Get outcomes for that model

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  filter(cast == 'hindcast') %>% 
  mutate(Prob_change = ifelse(Prob_gain_neutral > 50, Prob_gain_neutral, Prob_loss))

# join outcomes to typologies and compare probability of loss/gain with historical gross loss/gain (1996-2020)

threshold <- 75 # threshold for defining when a prediction is ambiguous or not

land <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Land_Gain = ifelse(LandwardMang > threshold, 'Gain', 'No Gain'),
         Land_Ambig = ifelse(LandwardMang <= threshold & LandwardMang >= -threshold, 'Ambiguous', 'Not Ambiguous'),
         Land_Loss = ifelse(LandwardMang < -threshold, 'Loss', 'No Loss')) %>% 
  mutate(Land_Change = ifelse(Land_Gain == 'Gain', 'Gain', NA),
         Land_Change = ifelse(Land_Ambig == 'Ambiguous', 'Ambiguous', Land_Change),
         Land_Change = ifelse(Land_Loss == 'Loss', 'Loss', Land_Change)) %>% 
  inner_join(select(spatial_dat, Type, sea_change_obs:land_loss_obs), by = 'Type') %>% 
  select(Type, Land_Gain:Land_Change, land_change_obs, land_gain_obs, land_loss_obs)
write.csv(land, 'outputs/land-validation-results.csv', row.names = F)

sea <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Sea_Gain = ifelse(SeawardMang > threshold, 'Gain', 'No Gain'),
         Sea_Ambig = ifelse(SeawardMang <= threshold & SeawardMang >= -threshold, 'Ambiguous', 'Not Ambiguous'),
         Sea_Loss = ifelse(SeawardMang < -threshold, 'Loss', 'No Loss')) %>% 
  mutate(Sea_Change = ifelse(Sea_Gain == 'Gain', 'Gain', NA),
         Sea_Change = ifelse(Sea_Ambig == 'Ambiguous', 'Ambiguous', Sea_Change),
         Sea_Change = ifelse(Sea_Loss == 'Loss', 'Loss', Sea_Change)) %>% 
  inner_join(select(spatial_dat, Type, sea_change_obs:land_loss_obs), by = 'Type') %>% 
  select(Type, Sea_Gain:Sea_Change, sea_change_obs, sea_gain_obs, sea_loss_obs)
write.csv(sea, 'outputs/sea-validation-results.csv', row.names = F)

# calculate overall prediction/classification accuracy
# and commission (users accuracy) and omission (producers accuracy) for each class

calc_accuracy <- function(x, x2){ # x is vector of predictions, x2 is reference vector
cont.table <- confusionMatrix(factor(x), factor(x2))$table # contingency table
commission <- diag(cont.table)/rowSums(cont.table)*100
omission <- diag(cont.table)/colSums(cont.table)*100
overall.accuracy <- sum(diag(cont.table))/sum(cont.table)*100
class.df <- data.frame(class = levels(factor(x2)), overall_accuracy = overall.accuracy, omission_accuracy = omission, commission_accuracy = commission)
accuracy_list <- list(class.df, cont.table)
names(accuracy_list) <- c('accuracy.results', 'contingency.table')
return(accuracy_list)
}

# individual classes
land_validate <- filter(land, Land_Ambig != 'Ambiguous') # get rid of ambiguous responses, can't validate
sea_validate <- filter(sea, Sea_Ambig != 'Ambiguous') # get rid of ambiguous responses, can't validate

calc_accuracy(land_validate$Land_Gain, land_validate$land_gain_obs)
calc_accuracy(land_validate$Land_Loss, land_validate$land_loss_obs)

calc_accuracy(sea_validate$Sea_Gain, sea_validate$sea_gain_obs)
calc_accuracy(sea_validate$Sea_Loss, sea_validate$sea_loss_obs)

# where are we not predicting well for seaward losses and seaward ambiguous responses? what are we predicting instead?

sea_poor_loss <- typ %>% 
  inner_join(filter(sea, Sea_Loss != sea_gross_loss))

sea_poor_gain_loss <- typ %>% 
  inner_join(filter(sea, Sea_Ambig != sea_gross_gain_loss))

qtm(sea_poor_loss, dots.col = 'sea_change_c')
qtm(sea_poor_gain_loss, dots.col = 'sea_change_c')


