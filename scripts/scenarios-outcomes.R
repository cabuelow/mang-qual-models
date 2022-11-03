# determine probability of landward and seaward mangrove increase under different action, geomorphic, and pressure settings

library(igraph)
library(QPress)
library(tidyverse)
library(patchwork)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
theme_set(theme_classic())
source('scripts/models.R')
source('scripts/helpers.R')

# set up scenario simulations

numsims <- 10000
model <- modelA_driv

# visual check and save

grViz(grviz.digraph(model))
grViz(grviz.digraph(model)) %>%
  save_png("outputs/model-signed-digraph.png")

# check for stability where all weights are equal (i.e., = 1)

stable.community(adjacency.matrix(model))

#actions = c('None', 'Restoration', 'Protection', 'PermeableDam')
settings =  c('TidalAmp', 'TidalFreq', 'HydroEnergy', 'SedSupply')
drivers = c('SeaLevelRise', 'SeaLevelFall', 'Cyclones', 'CoastalDev', 'BasementUp', 'Subsidence')

all <- data.frame(#action = rep(actions, each = length(drivers)*length(settings)),
                  setting = rep(settings, each = length(drivers)),
                  driver = rep(drivers, times = length(settings)))

# loop through scenarios with system.sim.press and store outcomes

out <- list()
stability <- list()
weights <- list()

for(i in 1:nrow(all)){
  p <- rep(1, ncol(all))
  names(p) <-  as.character(all[i,])
  if(length(which(names(p) == 'None')) != 0){
    p <- p[-which(names(p) == 'None')]
  }
  sim <- system.sim_press(numsims, model, perturb = p)
  out[[i]] <- sim$stableoutcome
  stability[[i]] <- sim$stability.df
  weights[[i]] <- sim$stableweights
}

# calculate potential stability in each scenario

scenario_stability <- data.frame(all, potential_stability = do.call(rbind, stability))

# wrangle outcomes

outcomes <- do.call(rbind, out) %>% 
  mutate(scnr = rep(apply(all, 1, paste, collapse = ' & '), each = c(numsims*length(node.labels(model)))),
         setting = rep(all$setting, each = c(numsims*length(node.labels(model)))),
         driver = rep(all$driver, each = c(numsims*length(node.labels(model)))))

outcomes <- outcomes %>% 
  mutate(setting_label = recode(outcomes$setting, 
                                 'HydroEnergy' = 'Hydrodynamic Energy',
                                 'SedSupply' = 'Sediment Supply',
                                 'TidalAmp' = 'Tidal Amplitude',
                                 'TidalFreq' = 'Tidal Frequency'),
         driver_label = recode(outcomes$driver,
                                'SeaLevelRise' = 'Sea-level Rise',
                                'SeaLevelFall' = 'Sea-level Fall',
                                'CoastalDev' = 'Coastal Development',
                                'BasementUp' = 'Basement Uplift'))

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
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  theme(legend.position = 'none') +
  facet_wrap(~setting_label) +
  ggtitle('A) Seaward mangroves')
a

b <- ggplot(landward) +
  geom_bar(aes(y = driver_label, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  facet_wrap(~setting_label) +
  ggtitle('B) Landward mangroves') +
  theme(legend.title = element_blank(),
        axis.text.y =  element_blank())

a+b

ggsave('outputs/outcomes.png', width = 10, height = 5)

# wrangle matrix weights

scnr_weights <- do.call(rbind, weights) %>% 
  mutate(scnr = rep(apply(all, 1, paste, collapse = ' & '), each = c(numsims*length(edge.labels(model)))),
         setting = rep(all$setting, each = c(numsims*length(edge.labels(model)))),
         driver = rep(all$driver, each = c(numsims*length(edge.labels(model))))) %>% 
  left_join(filter(outcomes, var == 'SeawardMang'))

scnr_weights$from_var <- sapply(strsplit(as.character(scnr_weights$param), "\\ "), `[`, 1)
scnr_weights$to_var <- sapply(strsplit(as.character(scnr_weights$param), "\\ "), `[`, 3)

# extract for ambiguous scenarios

ambig <- scnr_weights %>% 
  filter(setting == 'SedSupply', driver == 'Cyclones') %>% 
  mutate(outcome = ifelse(outcome > 0, 1, outcome)) %>% 
  mutate(outcome = ifelse(outcome < 0, 0, outcome))

paramsvis <- subset(ambig, from_var != to_var)

ggplot() +
  geom_density(data = filter(paramsvis, outcome == 1), 
               aes(x = weight), colour = 'blue', fill = 'blue', alpha = 0.2) +
  geom_density(data = filter(paramsvis, outcome == 0), 
               aes(x = weight), colour = 'yellow', fill = 'yellow', alpha = 0.2) +
  facet_wrap(~param, scales = 'free') +
  theme(strip.text.x = element_text(size = 5))

nrow(filter(paramsvis, outcome == 1))
nrow(filter(paramsvis, outcome == 0))
