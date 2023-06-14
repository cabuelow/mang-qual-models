# compare hindcast probability of mangrove loss and gain with historical loss and gain to validate network model

library(tidyverse)
library(caret)
library(sf)
library(scales)
library(patchwork)
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
         land_loss_obs = ifelse(land_gross_loss == 1, 'Loss', 'No Loss'),
         sea_net_gain_obs = ifelse(sea_net_gain == 1, 'Gain', 'No Gain'),
         sea_net_loss_obs = ifelse(sea_net_loss == 1, 'Loss', 'No Loss'),
         land_net_gain_obs = ifelse(land_net_gain == 1, 'Gain', 'No Gain'),
         land_net_loss_obs = ifelse(land_net_loss == 1, 'Loss', 'No Loss'),
         land_net_change_obs = ifelse(land_net_gain == 1, 'Gain', 'Loss'),
         sea_net_change_obs = ifelse(sea_net_gain == 1, 'Gain', 'Loss'))

# which model outcomes to validate? Get outcomes for that model

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  filter(cast == 'hindcast') %>% 
  mutate(Prob_change = ifelse(Prob_gain_neutral > 50, Prob_gain_neutral, Prob_loss))

# join outcomes to typologies and compare probability of loss/gain with historical gross loss/gain (1996-2020)

threshold <- seq(60, 90, by = 5) # threshold for defining when a prediction is ambiguous or not

# land
tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  tmp[[i]] <- dat %>% 
    pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
    mutate(Land_Gain = ifelse(LandwardMang > thresh, 'Gain', 'No Gain'),
           Land_Ambig = ifelse(LandwardMang <= thresh & LandwardMang >= -thresh, 'Ambiguous', 'Not Ambiguous'),
           Land_Loss = ifelse(LandwardMang < -thresh, 'Loss', 'No Loss')) %>% 
    mutate(Land_Change = ifelse(Land_Gain == 'Gain', 'Gain', NA),
           Land_Change = ifelse(Land_Ambig == 'Ambiguous', 'Ambiguous', Land_Change),
           Land_Change = ifelse(Land_Loss == 'Loss', 'Loss', Land_Change)) %>% 
    inner_join(select(spatial_dat, Type, sea_change_obs:sea_net_change_obs), by = 'Type') %>% 
    select(Type, Land_Gain:Land_Change, land_change_obs, land_gain_obs, land_loss_obs, land_net_gain_obs, land_net_loss_obs, land_net_change_obs) %>% 
    mutate(ambig_threshold = thresh)
}
land <- do.call(rbind, tmp)
write.csv(land, 'outputs/validation/land-validation-results.csv', row.names = F)

# sea
tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  tmp[[i]] <- dat %>% 
    pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
    mutate(Sea_Gain = ifelse(SeawardMang > thresh, 'Gain', 'No Gain'),
           Sea_Ambig = ifelse(SeawardMang <= thresh & SeawardMang >= -thresh, 'Ambiguous', 'Not Ambiguous'),
           Sea_Loss = ifelse(SeawardMang < -thresh, 'Loss', 'No Loss')) %>% 
    mutate(Sea_Change = ifelse(Sea_Gain == 'Gain', 'Gain', NA),
           Sea_Change = ifelse(Sea_Ambig == 'Ambiguous', 'Ambiguous', Sea_Change),
           Sea_Change = ifelse(Sea_Loss == 'Loss', 'Loss', Sea_Change)) %>% 
    inner_join(select(spatial_dat, Type, sea_change_obs:sea_net_change_obs), by = 'Type') %>% 
    select(Type, Sea_Gain:Sea_Change, sea_change_obs, sea_gain_obs, sea_loss_obs, sea_net_gain_obs, sea_net_loss_obs, sea_net_change_obs) %>% 
    mutate(ambig_threshold = thresh)
}
sea <- do.call(rbind, tmp)
write.csv(sea, 'outputs/validation/sea-validation-results.csv', row.names = F)

# calculate overall prediction/classification accuracy
# and commission (users accuracy) and omission (producers accuracy) for each class

# function
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

# remove ambiguous responses
land_validate <- filter(land, Land_Ambig != 'Ambiguous') # get rid of ambiguous responses, can't validate
sea_validate <- filter(sea, Sea_Ambig != 'Ambiguous') # get rid of ambiguous responses, can't validate

# net change validation
tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  landval <- land_validate %>% filter(ambig_threshold == thresh)
  seaval <- sea_validate %>% filter(ambig_threshold == thresh)
  results <- data.frame(validation = 'net',
                        mangrove = 'Landward',
                        ambig_threshold = thresh, 
                        calc_accuracy(landval$Land_Change, landval$land_net_change_obs)$accuracy.results)
  results2 <- data.frame(validation = 'net',
                        mangrove = 'Seaward',
                        ambig_threshold = thresh, 
                        calc_accuracy(seaval$Sea_Change, seaval$sea_net_change_obs)$accuracy.results)
  tmp[[i]] <- rbind(results, results2)
}
net_val <- do.call(rbind, tmp)

# gross gains/losses validation
tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  landval <- land_validate %>% filter(ambig_threshold == thresh)
  seaval <- sea_validate %>% filter(ambig_threshold == thresh)
  results <- data.frame(validation = 'gross',
                        mangrove = 'Landward',
                        ambig_threshold = thresh, 
                        calc_accuracy(landval$Land_Gain, landval$land_gain_obs)$accuracy.results)
  results2 <- data.frame(validation = 'gross',
                        mangrove = 'Landward',
                        ambig_threshold = thresh, 
                        calc_accuracy(landval$Land_Loss, landval$land_loss_obs)$accuracy.results)
  results3 <- data.frame(validation = 'gross',
                         mangrove = 'Seaward',
                         ambig_threshold = thresh, 
                         calc_accuracy(seaval$Sea_Gain, seaval$sea_gain_obs)$accuracy.results)
  results4 <- data.frame(validation = 'gross',
                         mangrove = 'Seaward',
                         ambig_threshold = thresh, 
                         calc_accuracy(seaval$Sea_Loss, seaval$sea_loss_obs)$accuracy.results)
  tmp[[i]] <- rbind(results, results2, results3, results4)
}
gross_val <- do.call(rbind, tmp)
accuracy <- rbind(net_val, gross_val) %>% mutate(mangrove = factor(mangrove, levels = c('Seaward', 'Landward')))
write.csv(accuracy, 'outputs/validation/accuracy-threshold-vary.csv', row.names = F)

# plot accuracy results

# net
a <- accuracy %>% 
  filter(validation == 'net') %>% 
  pivot_longer(cols = c(overall_accuracy:commission_accuracy), names_to = 'accuracy_cat', values_to = 'accuracy') %>% 
  mutate(remove = ifelse(accuracy_cat == 'overall_accuracy' & class == 'Loss', 1, 0)) %>% 
  filter(remove == 0) %>% 
  mutate(cat = ifelse(accuracy_cat != 'overall_accuracy', paste0(class, '_', accuracy_cat), accuracy_cat)) %>% 
  mutate(cat = factor(cat, levels = c('overall_accuracy', 'Gain_commission_accuracy', 'Loss_commission_accuracy',
                                      'Gain_omission_accuracy', 'Loss_omission_accuracy'))) %>% 
  ggplot() +
  geom_line(aes(x = ambig_threshold, y = accuracy, col = cat), linewidth = 1, alpha = 0.5) +
  facet_wrap(~mangrove) +
  ylim(c(0,100)) +
  ylab('Accuracy') +
  xlab('Ambiguity probability threshold (%)') +
  scale_color_manual(values = c("black",'darkslategray3', 'darkolivegreen3',
                                'goldenrod3', 'lightcoral'),
                     labels = c('Overall accuracy', 'Gain commission accuracy',
                                'Loss commission accuracy', 'Gain omission accuracy',
                                'Loss omission accuracy')) +
  geom_vline(xintercept = 75, linetype = 'dashed', alpha = 0.4) +
  ggtitle('A) Net change') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(size = 11)) +
  guides(color=guide_legend(ncol=1))
a

# gross losses
b <- accuracy %>% 
  filter(validation == 'gross') %>% 
  mutate(type = ifelse(class == 'No Gain', 'Gain', NA),
         type = ifelse(class == 'Gain', 'Gain', type),
         type = ifelse(class == 'Loss', 'Loss', type),
         type = ifelse(class == 'No Loss', 'Loss', type)) %>% 
  filter(type == 'Loss') %>% 
  pivot_longer(cols = c(overall_accuracy:commission_accuracy), names_to = 'accuracy_cat', values_to = 'accuracy') %>% 
  mutate(remove = ifelse(accuracy_cat == 'overall_accuracy' & class == 'No Loss', 1, 0)) %>% 
  mutate(remove = ifelse(accuracy_cat == 'overall_accuracy' & class == 'No Gain', 1, remove)) %>% 
  filter(remove == 0) %>% 
  mutate(cat = ifelse(accuracy_cat != 'overall_accuracy', paste0(class, '_', accuracy_cat), accuracy_cat)) %>%
  mutate(cat = factor(cat, levels = c('overall_accuracy', 'Loss_commission_accuracy', 'No Loss_commission_accuracy',
                                      'Loss_omission_accuracy', 'No Loss_omission_accuracy'))) %>% 
  ggplot() +
  geom_line(aes(x = ambig_threshold, y = accuracy, col = cat), linewidth = 1, alpha = 0.5) +
  facet_wrap(~mangrove) +
  ylim(c(0,100)) +
  ylab('') +
  xlab('Ambiguity probability threshold (%)') +
  scale_color_manual(values = c("black",'darkslategray3', 'darkolivegreen3',
                                'goldenrod3', 'lightcoral'),
                     labels = c('Overall accuracy', 'Loss commission accuracy',
                                'No Loss commission accuracy', 'Loss omission accuracy',
                                'No Loss omission accuracy')) +
  geom_vline(xintercept = 75, linetype = 'dashed', alpha = 0.4) +
  ggtitle('B) Gross losses') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(size = 11)) +
  guides(color=guide_legend(ncol=1))
b

# gross gains
c <- accuracy %>% 
  filter(validation == 'gross') %>% 
  mutate(type = ifelse(class == 'No Gain', 'Gain', NA),
         type = ifelse(class == 'Gain', 'Gain', type),
         type = ifelse(class == 'Loss', 'Loss', type),
         type = ifelse(class == 'No Loss', 'Loss', type)) %>% 
  filter(type == 'Gain') %>% 
  pivot_longer(cols = c(overall_accuracy:commission_accuracy), names_to = 'accuracy_cat', values_to = 'accuracy') %>% 
  mutate(remove = ifelse(accuracy_cat == 'overall_accuracy' & class == 'No Loss', 1, 0)) %>% 
  mutate(remove = ifelse(accuracy_cat == 'overall_accuracy' & class == 'No Gain', 1, remove)) %>% 
  filter(remove == 0) %>% 
  mutate(cat = ifelse(accuracy_cat != 'overall_accuracy', paste0(class, '_', accuracy_cat), accuracy_cat)) %>% 
  mutate(cat = factor(cat, levels = c('overall_accuracy', 'Gain_commission_accuracy', 'No Gain_commission_accuracy',
                                      'Gain_omission_accuracy', 'No Gain_omission_accuracy'))) %>% 
  ggplot() +
  geom_line(aes(x = ambig_threshold, y = accuracy, col = cat), linewidth = 1, alpha = 0.5) +
  facet_wrap(~mangrove) +
  ylim(c(0,100)) +
  ylab('') +
  xlab('Ambiguity probability threshold (%)') +
  scale_color_manual(values = c("black",'darkslategray3', 'darkolivegreen3',
                                'goldenrod3', 'lightcoral'),
                     labels = c('Overall accuracy', 'Gain commission accuracy',
                                'No Gain commission accuracy', 'Gain omission accuracy',
                                'No Gain omission accuracy')) +
  geom_vline(xintercept = 75, linetype = 'dashed', alpha = 0.4) +
  ggtitle('C) Gross gains') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(size = 11)) +
  guides(color=guide_legend(ncol=1))
c

a+b+c
ggsave('outputs/validation/accuracy-plot.png', width = 11, height = 4)

