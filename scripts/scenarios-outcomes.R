# determine probability of landward and seaward mangrove increase under different action, geomorphic, and pressure settings

library(igraph)
library(QPress)
library(tidyverse)
library(patchwork)
theme_set(theme_classic())
source('scripts/models.R')
source('scripts/helpers.R')

# run scenario simulations

numsims <- 1000

sim_SLR <- system.sim_press(numsims, modelA_driv, perturb=c(SeaLevelRise=1))
sim_SedSupp <- system.sim_press(numsims, modelA_driv, perturb=c(SedSupply=1))
sim_SLR_SedSupp <- system.sim_press(numsims, modelA_driv, perturb=c(SeaLevelRise=1,
                                                                 SedSupply=1))
sim_SLF_SedSupp <- system.sim_press(numsims, modelA_driv, perturb=c(SeaLevelFall=1,
                                                                 SedSupply=1))

#TODO: write a function to make all combinations of scenarios want to run, 
# based on actions, geomorphology, and drivers

actions <- c('Restoration', 'Protection', 'PermeableDam')
settings <- c('TidalAmp', 'TidalFreq', 'HydroEnergy', 'SedSupply')
drivers <- c('')

# then loop through scenarios with system.sim.press and store outcomes

# wrangle outcomes

outcomes <- data.frame(rbind(sim_SLR$stableoutcome, sim_SedSupp$stableoutcome, 
                             sim_SLR_SedSupp$stableoutcome, sim_SLF_SedSupp$stableoutcome))

outcomes$scnr <- rep(c("SLR", "Sediment supply", "Sediment supply & SLR", "Sediment supply & SLF"), 
                     each = c(numsims*length(node.labels(modelA_driv))))

# calculate proportion of stable models that have positive, negative, or neutral outcome in landward/seaward mangrove response

seaward <- outcomes %>% 
  filter(var == 'SeawardMang') %>% 
  group_by(scnr) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

landward <- outcomes %>% 
  filter(var == 'LandwardMang') %>% 
  group_by(scnr) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

# Note, as potential TODO could boostrap resample here to get estimate of uncertainty around that probability

a <- ggplot(seaward) +
  geom_bar(aes(y = scnr, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  xlab('Proportion of outcomes') +
  ylab('') +
  theme(legend.position = 'none') +
  ggtitle('Seaward mangroves')

b <- ggplot(landward) +
  geom_bar(aes(y = scnr, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  xlab('Proportion of outcomes') +
  ylab('') +
  ggtitle('Landward mangroves') +
  theme(legend.title = element_blank(),
        axis.text.y =  element_blank())

a+b

# now get parameter weights for increasing vs. decreasing
