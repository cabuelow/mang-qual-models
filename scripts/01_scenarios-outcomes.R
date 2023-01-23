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

set.seed(123)
numsims <- 10000
model <- modelA
model

# check for stability where all weights are equal (i.e., = 1)

stable.community(adjacency.matrix(model))

# set-up perturbation scenarios

pressures <- c('SeaLevelRise', 'Cyclones', 'GroundSubsid', 'CoastalDev', 'Erosion', 'Dams')
scenarios <- lapply(pressures, function (x) {
 p <- rep(1, length(x))
 names(p) <- as.character(x)
 return(p)
})

# loop through scenarios with system.sim.press and store outcomes

out <- list()
stability <- list()
weights <- list()

for(i in seq_along(pressures)){
  sim <- system.sim_press(numsims, model, perturb = scenarios[[i]])
  out[[i]] <- sim$stableoutcome
  stability[[i]] <- sim$stability.df
  weights[[i]] <- sim$stableweights
}

# calculate potential stability in each scenario (i.e. proportion of stable matrices out of total simulated)

scenario_stability <- data.frame(pressures, do.call(rbind, stability))

# wrangle outcomes

outcomes <- do.call(rbind, out) %>% 
  mutate(pressure = rep(pressures, each = c(numsims*length(node.labels(model)))))
  #mutate(scnr = rep(apply(all, 1, paste, collapse = ' & '), each = c(numsims*length(node.labels(model)))),
   #      setting = rep(all$setting, each = c(numsims*length(node.labels(model)))),
    #     driver = rep(all$driver, each = c(numsims*length(node.labels(model)))))

outcomes <- outcomes %>% 
  mutate(pressure_label = recode(outcomes$pressure,
                                'SeaLevelRise' = 'Sea-level Rise',
                                'CoastalDev' = 'Coastal Development',
                                'GroundSubsid' = 'Groundwater Subsidence'))

# calculate proportion of stable models that have positive, negative, or neutral outcome in landward/seaward mangrove response

seaward <- outcomes %>% 
  filter(var == 'SeawardMang') %>% 
  group_by(pressure, pressure_label) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')),
         pressure_label = factor(pressure_label, levels = c('Erosion', 'Sea-level Rise',
                                                          'Cyclones', 'Coastal Development',
                                                        'Groundwater Subsidence', 'Dams')))

landward <- outcomes %>% 
  filter(var == 'LandwardMang') %>% 
  group_by(pressure, pressure_label) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>%  
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')),
         pressure_label = factor(pressure_label, levels = c('Erosion', 'Sea-level Rise',
                                                            'Cyclones', 'Coastal Development',
                                                            'Groundwater Subsidence', 'Dams')))

landsea <- rbind(data.frame(seaward, mangrove = 'seaward'), data.frame(landward, mangrove = 'landward'))

# Note, as potential TODO could boostrap resample here to get estimate of uncertainty around that probability

a <- ggplot(seaward) +
  geom_bar(aes(y = pressure_label, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values=c("darkslategray3", 'darkolivegreen2', "brown3")) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  theme(legend.position = 'none') +
  #facet_wrap(~setting_label) +
  ggtitle('A) Seaward mangroves')

b <- ggplot(landward) +
  geom_bar(aes(y = pressure_label, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values=c("darkslategray3", 'darkolivegreen2', "brown3")) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  #facet_wrap(~setting_label) +
  ggtitle('B) Landward mangroves') +
  theme(legend.title = element_blank(),
        axis.text.y =  element_blank())

a+b

ggsave('outputs/outcomes.png', width = 10, height = 5.3)

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
