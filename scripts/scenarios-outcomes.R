# determine probability of landward and seaward mangrove increase under different action, geomorphic, and pressure settings

library(igraph)
library(QPress)
library(tidyverse)
library(patchwork)
theme_set(theme_classic())
source('scripts/models.R')
source('scripts/helpers.R')

#TODO: figure out why the 'none' scenario is messing things up

# set up scenario simulations

numsims <- 10000
model <- modelA_driv

#actions = c('None', 'Restoration', 'Protection', 'PermeableDam')
settings =  c('TidalAmp', 'TidalFreq', 'HydroEnergy', 'SedSupply')
drivers = c('SeaLevelRise', 'SeaLevelFall', 'Cyclones', 'Dams', 'CoastalDev', 'BasementUp', 'Subsidence')

all <- data.frame(#action = rep(actions, each = length(drivers)*length(settings)),
                  setting = rep(settings, times = length(drivers)),
                  driver = rep(drivers, times = length(settings)))

# loop through scenarios with system.sim.press and store outcomes

out <- list()

for(i in 1:nrow(all)){
  p <- rep(1, ncol(all))
  names(p) <-  as.character(all[i,])
  if(length(which(names(p) == 'None')) != 0){
    p <- p[-which(names(p) == 'None')]
  }
  sim <- system.sim_press(numsims, model, perturb = p)
  out[[i]] <- sim$stableoutcome
}

# wrangle outcomes

outcomes <- do.call(rbind, out)
outcomes$scnr <- rep(apply(all, 1, paste, collapse = ' & '), each = c(numsims*length(node.labels(model))))
outcomes$setting <- rep(all$setting, each = c(numsims*length(node.labels(model))))
outcomes$driver <- rep(all$driver, each = c(numsims*length(node.labels(model))))
outcomes$setting_label <- recode(outcomes$setting, 
                                 'HydroEnergy' = 'Hydrodynamic Energy',
                                 'SedSupply' = 'Sediment Supply',
                                 'TidalAmp' = 'Tidal Amplitude',
                                 'TidalFreq' = 'Tidal Frequency')
outcomes$driver_label <- recode(outcomes$driver,
                                'SeaLevelRise' = 'Sea-level Rise',
                                'SeaLevelFall' = 'Sea-level Fall',
                                'CoastalDev' = 'Coastal Development',
                                'BasementUp' = 'Basement Uplift')
# calculate proportion of stable models that have positive, negative, or neutral outcome in landward/seaward mangrove response

seaward <- outcomes %>% 
  filter(var == 'SeawardMang') %>% 
  group_by(scnr, setting_label, driver_label) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

landward <- outcomes %>% 
  filter(var == 'LandwardMang') %>% 
  group_by(scnr, setting_label, driver_label) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

landsea <- rbind(data.frame(seaward, mangrove = 'seaward'), data.frame(landward, mangrove = 'landward'))

# Note, as potential TODO could boostrap resample here to get estimate of uncertainty around that probability

a <- ggplot(seaward) +
  geom_bar(aes(y = driver_label, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  xlab('Proportion of outcomes') +
  ylab('') +
  theme(legend.position = 'none') +
  facet_wrap(~setting_label) +
  ggtitle('Seaward mangroves')
a

b <- ggplot(landward) +
  geom_bar(aes(y = driver_label, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  xlab('Proportion of outcomes') +
  ylab('') +
  facet_wrap(~setting_label) +
  ggtitle('Landward mangroves') +
  theme(legend.title = element_blank(),
        axis.text.y =  element_blank())

a+b

ggsave('outputs/outcomes.png', width = 10, height = 5)

