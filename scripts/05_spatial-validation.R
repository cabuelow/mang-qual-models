# compare hindcast probability of mangrove loss and gain with historical loss and gain to validate network model

library(tidyverse)
library(caret)
source('scripts/helpers/models_v2.R')

# wrangle historical SRS observations of mangrove loss and gain into categories of change (no change, loss, gain, loss and gain)

spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  mutate(sea_change = sea_gain + sea_loss,
        land_change = land_gain + land_loss) %>% 
  mutate(sea_change_c = ifelse(sea_change == 2, 'Loss & Gain', NA),
         sea_change_c = ifelse(sea_change != 2 & sea_gain ==1, 'Gain', sea_change_c),
         sea_change_c = ifelse(sea_change != 2 & sea_loss ==1, 'Loss', sea_change_c),
         sea_change_c = ifelse(sea_change ==0, 'No change', sea_change_c),
         land_change_c = ifelse(land_change == 2, 'Loss & Gain', NA),
         land_change_c = ifelse(land_change != 2 & land_gain ==1, 'Gain', land_change_c),
         land_change_c = ifelse(land_change != 2 & land_loss ==1, 'Loss', land_change_c),
         land_change_c = ifelse(land_change ==0, 'No change', land_change_c))

# which model outcomes to validate? Get outcomes for that model

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  filter(cast == 'hindcast') %>% 
  mutate(Prob_change = ifelse(Prob_gain_neutral > 50, Prob_gain_neutral, Prob_loss))

# join outcomes to typologies and compare probability of loss/gain with historical gross loss/gain (1996-2020)

land <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Land_Gain = ifelse(LandwardMang > 75, 1, 0),
         Land_Ambig = ifelse(LandwardMang <= 50 & LandwardMang >= -75, 1, 0),
         Land_Loss = ifelse(LandwardMang < -75, 1, 0)) %>% 
  filter(Land_Ambig != 1) %>% # filter out ambiguous predictions
  inner_join(select(spatial_dat, Type, sea_gain:land_loss), by = 'Type') %>%
  as.data.frame()

sea <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Sea_Gain = ifelse(SeawardMang > 75, 1, 0),
         Sea_Ambig = ifelse(SeawardMang <= 50 & SeawardMang >= -75, 1, 0),
         Sea_Loss = ifelse(SeawardMang < -75, 1, 0)) %>% 
  filter(Sea_Ambig != 1) %>% # filter out ambiguous predictions
  inner_join(select(spatial_dat, Type, sea_gain:land_loss), by = 'Type') %>%
  as.data.frame()

# now do confusion matrix of non-ambiguous predictions of gain or loss, compared to net gain or loss

calc_accuracy <- function(x, x2){
  conMatrix <- confusionMatrix(factor(x, levels = c(0,1)), factor(x2, levels = c(0,1)))
  accuracy_percent <- conMatrix$overall[1]*100
  return(accuracy_percent)
  }

calc_accuracy(land$Land_Gain, land$land_gain)
calc_accuracy(land$Land_Loss, land$land_loss)
calc_accuracy(sea$Sea_Loss, sea$sea_loss)
calc_accuracy(sea$Sea_Gain, sea$sea_gain)