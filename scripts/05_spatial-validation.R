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
  mutate(sea_change_c = ifelse(sea_gross_gain_loss == 2, 'Loss & Gain', NA),
         sea_change_c = ifelse(sea_gross_gain_loss != 2 & sea_gross_gain ==1, 'Gain', sea_change_c),
         sea_change_c = ifelse(sea_gross_gain_loss != 2 & sea_gross_loss ==1, 'Loss', sea_change_c),
         sea_change_c = ifelse(sea_gross_gain_loss == 0, 'Gain', sea_change_c), # here treating neutrality as gain, as in network model
         land_change_c = ifelse(land_gross_gain_loss == 2, 'Loss & Gain', NA),
         land_change_c = ifelse(land_gross_gain_loss != 2 & land_gross_gain ==1, 'Gain', land_change_c),
         land_change_c = ifelse(land_gross_gain_loss != 2 & land_gross_loss ==1, 'Loss', land_change_c),
         land_change_c = ifelse(land_gross_gain_loss == 0, 'Gain', land_change_c), # here treating neutrality as gain, as in network model 
         sea_gross_loss = ifelse(sea_gross_gain_loss == 2, 0, sea_gross_loss),
         sea_gross_gain = ifelse(sea_gross_gain_loss == 2, 0, sea_gross_gain),
         land_gross_loss = ifelse(land_gross_gain_loss == 2, 0, land_gross_loss),
         land_gross_gain = ifelse(land_gross_gain_loss == 2, 0, land_gross_gain)) %>% 
  mutate(sea_gross_gain_loss = ifelse(sea_gross_gain_loss == 2, 1, 0),
         land_gross_gain_loss = ifelse(land_gross_gain_loss == 2, 1, 0))

# which model outcomes to validate? Get outcomes for that model

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  filter(cast == 'hindcast') %>% 
  mutate(Prob_change = ifelse(Prob_gain_neutral > 50, Prob_gain_neutral, Prob_loss))

# join outcomes to typologies and compare probability of loss/gain with historical gross loss/gain (1996-2020)

threshold <- 99 # threshold for defining when a prediction is ambiguous or not

land <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Land_Gain = ifelse(LandwardMang > threshold, 1, 0),
         Land_Ambig = ifelse(LandwardMang <= threshold & LandwardMang >= -threshold, 1, 0),
         Land_Loss = ifelse(LandwardMang < -threshold, 1, 0)) %>% 
  mutate(Land_Change = ifelse(Land_Gain == 1, 'Gain', NA),
         Land_Change = ifelse(Land_Ambig == 1, 'Loss & Gain', Land_Change),
         Land_Change = ifelse(Land_Loss == 1, 'Loss', Land_Change)) %>% 
  inner_join(select(spatial_dat, Type, sea_gross_gain:land_gross_gain_loss, sea_change_c, land_change_c), by = 'Type') %>% 
  select(Type, Land_Gain:Land_Change, land_gross_gain, land_gross_loss, land_gross_gain_loss, land_change_c)
write.csv(land, 'outputs/land-validation-results.csv', row.names = F)

sea <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Sea_Gain = ifelse(SeawardMang > threshold, 1, 0),
         Sea_Ambig = ifelse(SeawardMang <= threshold & SeawardMang >= -threshold, 1, 0),
         Sea_Loss = ifelse(SeawardMang < -threshold, 1, 0)) %>% 
  mutate(Sea_Change = ifelse(Sea_Gain == 1, 'Gain', NA),
         Sea_Change = ifelse(Sea_Ambig == 1, 'Loss & Gain', Sea_Change),
         Sea_Change = ifelse(Sea_Loss == 1, 'Loss', Sea_Change)) %>% 
  inner_join(select(spatial_dat, Type, sea_gross_gain:land_gross_gain_loss, sea_change_c, land_change_c), by = 'Type') %>% 
  select(Type, Sea_Gain:Sea_Change, sea_gross_gain, sea_gross_loss, sea_gross_gain_loss, sea_change_c)
write.csv(sea, 'outputs/sea-validation-results.csv', row.names = F)

# calculate overall prediction/classification accuracy
# and commission (users accuracy) and omission (producers accuracy) for each class

calc_accuracy <- function(x, x2){ # x is vector of predictions, x2 is reference vector
cont.table <- confusionMatrix(factor(x), factor(x2))$table # get a contingency table of reference values and predictions
commission <- diag(cont.table)/rowSums(cont.table)*100
omission <- diag(cont.table)/colSums(cont.table)*100
overall.accuracy <- sum(diag(cont.table))/sum(cont.table)*100
class.df <- data.frame(class = levels(factor(x2)), overall_accuracy = overall.accuracy, omission_accuracy = omission, commission_accuracy = commission)
accuracy_list <- list(class.df, cont.table)
names(accuracy_list) <- c('accuracy.results', 'contingency.table')
return(accuracy_list)
}

accuracy <- calc_accuracy(land$Land_Change, land$land_change_c)
accuracy

calc_accuracy(land$Land_Gain, land$land_gross_gain)
calc_accuracy(land$Land_Loss, land$land_gross_loss)
calc_accuracy(land$Land_Ambig, land$land_gross_gain_loss)
calc_accuracy(sea$Sea_Loss, sea$sea_gross_loss)
calc_accuracy(sea$Sea_Gain, sea$sea_gross_gain)
calc_accuracy(sea$Sea_Ambig, sea$sea_gross_gain_loss)
calc_accuracy(land$Land_Change, land$land_change_c)
calc_accuracy(sea$Sea_Change, sea$sea_change_c)

# where are we not predicting well for seaward losses and seaward ambiguous responses? what are we predicting instead?

sea_poor_loss <- typ %>% 
  inner_join(filter(sea, Sea_Loss != sea_gross_loss))

sea_poor_gain_loss <- typ %>% 
  inner_join(filter(sea, Sea_Ambig != sea_gross_gain_loss))

qtm(sea_poor_loss, dots.col = 'sea_change_c')
qtm(sea_poor_gain_loss, dots.col = 'sea_change_c')


