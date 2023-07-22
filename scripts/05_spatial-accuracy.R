# compare hindcast probability of mangrove loss and gain with historical loss and gain to validate network model

library(tidyverse)
library(caret)
library(sf)
library(scales)
library(patchwork)
library(ggh4x)
source('scripts/helpers/models_v2.R')

drivers <- read.csv('data/typologies/SLR_Data.csv')

# wrangle historical SRS observations of mangrove loss and gain into categories of change

spatial_dat <- read.csv('outputs/master-dat.csv')

# which model outcomes to validate? Get outcomes for that model

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  filter(cast == 'hindcast')# %>% 
  #mutate(Prob_change = ifelse(Prob_gain > 50, Prob_gain, Prob_loss))

# join outcomes to typologies and compare probability of loss/gain with historical gross loss/gain (1996-2020)

threshold <- seq(60, 90, by = 5) # range of thresholds for defining when a prediction is ambiguous or not

tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  tmp[[i]] <- dat %>% 
    mutate(Change = case_when(Prob_gain > thresh ~ 'Gain',
                                   Prob_loss < -thresh ~ 'Loss',
                                   Prob_neutral > thresh ~ 'Neutral',
                              .default = 'Ambiguous')) %>%
    pivot_wider(id_cols = c('Type', 'cast', 'pressure_def'), names_from = 'var', values_from = 'Change', names_prefix = 'Change_') %>% 
    inner_join(select(spatial_dat, Type, pressure_def, land_net_change, sea_net_change), by = c('Type', 'pressure_def')) %>% 
    mutate(ambig_threshold = thresh)
}
class <- do.call(rbind, tmp)
write.csv(class, 'outputs/validation/validation-results.csv', row.names = F)

# calculate overall prediction/classification accuracy
# and commission (users accuracy) and omission (producers accuracy) for each class

# function
calc_accuracy <- function(x, x2){ # x is vector of predictions, x2 is reference vector
cont.table <- confusionMatrix(factor(x), factor(x2))$table # contingency table
commission <- diag(cont.table)/rowSums(cont.table)*100
omission <- diag(cont.table)/colSums(cont.table)*100
overall.accuracy <- sum(diag(cont.table))/sum(cont.table)*100
class.df <- data.frame(class = levels(factor(x2)), Overall_accuracy = overall.accuracy, Producers_accuracy = omission, Users_accuracy = commission)
accuracy_list <- list(class.df, cont.table)
names(accuracy_list) <- c('accuracy.results', 'contingency.table')
return(accuracy_list)
}

# remove ambiguous responses and units where commodities and erosion are dominant drivers of loss, can't validate with our model
land_validate <- class %>% 
  left_join(drivers, by = 'Type') %>% 
  filter(Change_LandwardMang != 'Ambiguous' & Commodities < 0.1 & Erosion < 0.1) 
sea_validate <- class %>%
  left_join(drivers, by = 'Type') %>% 
  filter(Change_SeawardMang != 'Ambiguous' & Commodities < 0.1 & Erosion < 0.1) 

# net change validation
tmp2 <- list()
for(j in seq_along(unique(dat$pressure_def))){
  land <- land_validate %>% filter(pressure_def == j)
  sea <- sea_validate %>% filter(pressure_def == j)
  tmp <- list()
for(i in seq_along(threshold)){
  thresh <- threshold[i]
  landval <- land %>% filter(ambig_threshold == thresh)
  seaval <- sea %>% filter(ambig_threshold == thresh)
  results <- data.frame(validation = 'net',
                        mangrove = 'Landward',
                        pressure_def = j,
                        ambig_threshold = thresh, 
                        calc_accuracy(landval$Change_LandwardMang, landval$land_net_change)$accuracy.results)
  results2 <- data.frame(validation = 'net',
                        mangrove = 'Seaward',
                        pressure_def = j,
                        ambig_threshold = thresh, 
                        calc_accuracy(seaval$Change_SeawardMang, seaval$sea_net_change)$accuracy.results)
  tmp[[i]] <- rbind(results, results2)
}
  tmp2[[j]] <- do.call(rbind, tmp)
}
accuracy <- do.call(rbind, tmp2) %>%  mutate(mangrove = factor(mangrove, levels = c('Seaward', 'Landward')))
write.csv(accuracy, 'outputs/validation/accuracy-threshold-vary.csv', row.names = F)

# heatmap of overall accuracy for combinations of pressure definition and ambiguity thresholds

accuracy %>% 
  filter(validation == 'net') %>% 
  pivot_longer(cols = Overall_accuracy:Users_accuracy, names_to = 'metric', values_to = 'accuracy') %>% 
  mutate(class = ifelse(metric == 'Overall_accuracy', 'Loss & Gain', class)) %>% 
  distinct() %>% 
  ggplot() +
  aes(x = ambig_threshold, y = pressure_def, fill = accuracy) +
  geom_tile() +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1, name = 'Accuracy') +
  facet_nested_wrap(~factor(mangrove) + factor(class) + factor(metric), nrow = 2) +
  ylab('Pressure definition') +
  xlab('Ambiguity probability threshold') +
  theme_classic()

ggsave('outputs/validation/accuracy-heatmap.png', width = 10, height = 4.5)

# plot accuracy results

# End here

# net
a <- accuracy %>% 
  filter(validation == 'net' & pressure_def == 1) %>% 
  pivot_longer(cols = c(overall_accuracy:commission_accuracy), names_to = 'accuracy_cat', values_to = 'accuracy') %>% 
  mutate(remove = ifelse(accuracy_cat == 'overall_accuracy' & class == 'Loss', 1, 0)) %>% 
  filter(remove == 0) %>% 
  mutate(cat = ifelse(accuracy_cat != 'overall_accuracy', paste0(class, '_', accuracy_cat), accuracy_cat)) %>% 
  mutate(cat = factor(cat, levels = c('overall_accuracy', 'Gain_commission_accuracy', 'Loss_commission_accuracy',
                                      'Gain_omission_accuracy', 'Loss_omission_accuracy'))) %>% 
  filter(!is.na(accuracy)) %>% 
  mutate(accuracy =  ifelse(cat == 'overall_accuracy', accuracy + 2, accuracy)) %>%  # here adding small number to avoid over plotting
  ggplot() +
  geom_line(aes(x = ambig_threshold, y = accuracy, col = cat),
            linewidth = 1, alpha = 0.7) +
  #coord_flip() +
  facet_wrap(~mangrove) +
  #ylim(c(0,100)) +
  ylab('Accuracy') +
  xlab('Ambiguity probability threshold (%)') +
  scale_color_manual(values = c("black",'darkslategray3', 'darkolivegreen3',
                                'goldenrod3', 'lightcoral'),
                     labels = c('Overall accuracy', 'Gain commission accuracy',
                                'Loss commission accuracy', 'Gain omission accuracy',
                                'Loss omission accuracy')) +
  geom_vline(xintercept = 60, linetype = 'dashed', alpha = 0.4) +
  ggtitle('C) Net change') +
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
  filter(!is.na(accuracy)) %>% 
  mutate(accuracy =  ifelse(cat == 'overall_accuracy', accuracy - 2, accuracy)) %>%  # here adding small number to avoid over plotting
  ggplot() +
  geom_line(aes(x = ambig_threshold, y = accuracy, col = cat), linewidth = 1, alpha = 0.7) +
  facet_wrap(~mangrove) +
  ylim(c(0,100)) +
  ylab('') +
  xlab('Ambiguity probability threshold (%)') +
  scale_color_manual(values = c("black",'darkslategray3', 'darkolivegreen3',
                                'goldenrod3', 'lightcoral'),
                     labels = c('Overall accuracy', 'Loss commission accuracy',
                                'No Loss commission accuracy', 'Loss omission accuracy',
                                'No Loss omission accuracy')) +
  geom_vline(xintercept = 60, linetype = 'dashed', alpha = 0.4) +
  ggtitle('D) Gross losses') +
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
  filter(!is.na(accuracy)) %>% 
  mutate(accuracy =  ifelse(cat == 'overall_accuracy', accuracy + 2, accuracy),
         accuracy =  ifelse(cat == 'Gain_commission_accuracy', accuracy - 2, accuracy)) %>%  # here adding small number to avoid over plotting
  ggplot() +
  geom_line(aes(x = ambig_threshold, y = accuracy, col = cat), linewidth = 1, alpha = 0.7) +
  facet_wrap(~mangrove) +
  ylim(c(0,100)) +
  ylab('') +
  xlab('Ambiguity probability threshold (%)') +
  scale_color_manual(values = c("black",'darkslategray3', 'darkolivegreen3',
                                'goldenrod3', 'lightcoral'),
                     labels = c('Overall accuracy', 'Gain commission accuracy',
                                'No Gain commission accuracy', 'Gain omission accuracy',
                                'No Gain omission accuracy')) +
  geom_vline(xintercept = 60, linetype = 'dashed', alpha = 0.4) +
  ggtitle('E) Gross gains') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(size = 11)) +
  guides(color=guide_legend(ncol=1))
c

a+b+c
ggsave('outputs/validation/accuracy-plot.png', width = 11, height = 4)

# gross gains/losses validation
tmp2 <- list()
for(j in seq_along(unique(dat$pressure_def))){
  land <- land_validate %>% filter(pressure_def == j)
  sea <- sea_validate %>% filter(pressure_def == j)
  tmp <- list()
  for(i in seq_along(threshold)){
    thresh <- threshold[i]
    landval <- land %>% filter(ambig_threshold == thresh)
    seaval <- sea %>% filter(ambig_threshold == thresh)
    results <- data.frame(validation = 'gross_gain',
                          mangrove = 'Landward',
                          pressure_def = j,
                          ambig_threshold = thresh, 
                          calc_accuracy(landval$Land_Gain, landval$land_gain_obs)$accuracy.results)
    results2 <- data.frame(validation = 'gross_loss',
                           mangrove = 'Landward',
                           pressure_def = j,
                           ambig_threshold = thresh, 
                           calc_accuracy(landval$Land_Loss, landval$land_loss_obs)$accuracy.results)
    results3 <- data.frame(validation = 'gross_gain',
                           mangrove = 'Seaward',
                           pressure_def = j,
                           ambig_threshold = thresh, 
                           calc_accuracy(seaval$Sea_Gain, seaval$sea_gain_obs)$accuracy.results)
    results4 <- data.frame(validation = 'gross_loss',
                           mangrove = 'Seaward',
                           pressure_def = j,
                           ambig_threshold = thresh, 
                           calc_accuracy(seaval$Sea_Loss, seaval$sea_loss_obs)$accuracy.results)
    tmp[[i]] <- rbind(results, results2, results3, results4)
  }
  tmp2[[j]] <- do.call(rbind, tmp)
}
gross_val <- do.call(rbind, tmp2)

