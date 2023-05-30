# determine probability of landward and seaward mangrove increase under different scenarios
# devtools::install_github("SWotherspoon/QPress",ref="Constrain")
# TODO: make scenario simulation more efficient and easier to implement
# Get rid of dams scenarios from function code

library(igraph)
library(QPress)
library(tidyverse)
library(patchwork)
#library(DiagrammeR)
#library(DiagrammeRsvg)
#library(rsvg)
theme_set(theme_classic())
source('scripts/helpers/models.R')
source('scripts/helpers/helpers.R')

# set up scenario simulations

set.seed(123)
numsims <- 1000

# perturbation scenarios

pressures <- c('Sea-level rise', 'Cyclones', 'Groundwater extraction', 'Coastal development', 'Erosion', #'Drought or Dams',
               'Drought', 'Dams', #'Extreme rainfall', 'Dams & Extreme Rainfall',
               'Sea-level rise & Cyclones', 'Sea-level rise & Groundwater extraction',
               'Sea-level rise  & Coastal development', 'Sea-level rise & Erosion', #'Sea-level rise & Drought or Dams')
               'Sea-level rise & Drought', 'Sea-level rise & Sediment', 'Sea-level rise & Dams') #, 'Sea-level rise & Extreme rainfall')
press.scenarios <- list(c(SeaLevelRise=1), c(Cyclones=1), c(GroundSubsid=1), c(CoastalDev=1), 
                        c(Erosion=1), #c(Sediment=-1),
                        c(Drought=1), c(Dams=1),
                        #c(Precipitation=1),
                        #c(Precipitation=1, Dams=1),
                        c(SeaLevelRise=1, Cyclones=1),
                        c(SeaLevelRise=1, GroundSubsid=1),
                        c(SeaLevelRise=1, CoastalDev=1), 
                        c(SeaLevelRise=1, Erosion=1),
                        c(SeaLevelRise=1, Drought=1),
                        c(SeaLevelRise=1, Sediment=1),
                        c(SeaLevelRise=1, Dams=1))
                        #c(SeaLevelRise=1, Precipitation=1))
                        # c(SeaLevelRise=1, Sediment=-1))

# edge constraint scenarios
#TODO: make this easier by changing how the function takes these constraints....

con.scenarios <- list(c('H', 'H', 'H', 'L'),
                    c('M', 'H', 'H', 'L'),
                    c('L', 'H', 'H', 'L'),
                    c('H', 'L', 'L', 'L'),
                    c('M',  'L', 'L', 'L'),
                    c('L', 'L', 'L', 'L'),
                    c('H', 'H', 'H', 'M'),
                    c('M', 'H', 'H', 'M'),
                    c('L', 'H', 'H', 'M'),
                    c('H', 'L', 'L', 'M'),
                    c('M',  'L', 'L', 'M'),
                    c('L', 'L', 'L', 'M'),
                    c('H', 'H', 'H', 'H'),
                    c('M', 'H', 'H', 'H'),
                    c('L', 'H', 'H', 'H'),
                    c('H', 'L', 'L', 'H'),
                    c('M',  'L', 'L', 'H'),
                    c('L', 'L', 'L', 'H'))
names(con.scenarios) <- c('Microtidal, High Hydro-connectivity, High Coastal squeeze', 
                          'Mesotidal, High Hydro-connectivity, High Coastal squeeze',
                          'Macrotidal, High Hydro-connectivity, High Coastal squeeze',
                          'Microtidal, Low Hydro-connectivity, High Coastal squeeze',
                          'Mesotidal, Low Hydro-connectivity, High Coastal squeeze',
                          'Macrotidal, Low Hydro-connectivity, High Coastal squeeze',
                          'Microtidal, High Hydro-connectivity, Medium Coastal squeeze',
                          'Mesotidal, High Hydro-connectivity, Medium Coastal squeeze',
                          'Macrotidal, High Hydro-connectivity, Medium Coastal squeeze',
                          'Microtidal, Low Hydro-connectivity, Medium Coastal squeeze',
                          'Mesotidal, Low Hydro-connectivity, Medium Coastal squeeze',
                          'Macrotidal, Low Hydro-connectivity, Medium Coastal squeeze',
                          'Microtidal, High Hydro-connectivity, Low Coastal squeeze',
                          'Mesotidal, High Hydro-connectivity, Low Coastal squeeze',
                          'Macrotidal, High Hydro-connectivity, Low Coastal squeeze',
                          'Microtidal, Low Hydro-connectivity, Low Coastal squeeze',
                          'Mesotidal, Low Hydro-connectivity, Low Coastal squeeze',
                          'Macrotidal, Low Hydro-connectivity, Low Coastal squeeze')

# set up relative edge constraint scenarios
# high sed supply model, sediment -> subVol will be greater than SLR neg interactions
# vice versa for low sed supply model

model.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), modelB),
                        parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), modelB))
names(model.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply')

# loop through scenarios with system.sim.press and store outcomes
# do a nested for loop, iterate over model, then class scenarios, then perturbation scenarios
# BUT do this better later! look into purrr

stability.ls1 <- list()
outcomes.ls1 <- list()

for(j in seq_along(model.scenarios)){
stability.ls <- list()
outcomes.ls <- list()
model <- model.scenarios[[j]]
for(k in seq_along(con.scenarios)){
  out <- list()
  stability <- list()
  weights <- list()
  class.con <- con.scenarios[[k]]
for(i in seq_along(pressures)){
  sim <- system.sim_press(numsims, constrainedigraph = model, 
                          from = c('SeaLevelRise', #'SeaLevelRise', #'SeaLevelRise', 
                                   'LandwardAvailableProp', 'SeawardAvailableProp',
                                   'SeaLevelRise'),
                          to = c('SeawardMang', #'SeawardPropag', #'SeawardEstabSpace', 
                                 'LandwardMang', 'SeawardMang',
                                 'LandwardMang'),
                          class = class.con,
                          perturb = press.scenarios[[i]])
  out[[i]] <- sim$stableoutcome
  stability[[i]] <- sim$stability.df
  weights[[i]] <- sim$stableweights
}
  # calculate potential stability in each scenario (i.e. proportion of stable matrices out of total simulated)
  stability.ls[[k]] <- data.frame(pressures, do.call(rbind, stability))
  # wrangle outcomes
  outcomes.ls[[k]] <- do.call(rbind, out) %>% 
    mutate(pressure = rep(pressures, each = c(numsims*length(node.labels(model$edges)))))
}
# calculate potential stability in each scenario (i.e. proportion of stable matrices out of total simulated)
stability.ls1[[j]] <- data.frame(constraint_scenario = names(con.scenarios), do.call(rbind, stability.ls))
# wrangle outcomes
outcomes.ls1[[j]] <- do.call(rbind, outcomes.ls) %>% 
  mutate(constraint_scenario = rep(names(con.scenarios), each = c(numsims*length(node.labels(model$edges))*length(pressures)))) 
}

# calculate potential stability in each scenario (i.e. proportion of stable matrices out of total simulated)
stability <- data.frame(constraint_scenario = names(model.scenarios), do.call(rbind, stability.ls1))
write.csv(stability, 'outputs/stability.csv', row.names = F)

# wrangle outcomes
outcomes <- do.call(rbind, outcomes.ls1) %>% 
  mutate(model_scenario = rep(names(model.scenarios), each = c(numsims*length(node.labels(modelB))*length(pressures)*length(con.scenarios)))) %>% 
  mutate(outcome = ifelse(pressure == 'Dams', round(outcome, 2), outcome))
saveRDS(outcomes, 'outputs/outcomes.rds')

# calculate proportion of stable models that have positive, negative, or neutral outcome in landward/seaward mangrove response

seaward_sedH <- outcomes %>% 
  filter(var == 'SeawardMang' & model_scenario == 'High Sediment Supply') %>% 
  group_by(constraint_scenario, pressure) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

landward_sedH <- outcomes %>% 
  filter(var == 'LandwardMang' & model_scenario == 'High Sediment Supply') %>% 
  group_by(constraint_scenario, pressure) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>%  
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

seaward_sedL <- outcomes %>% 
  filter(var == 'SeawardMang' & model_scenario == 'Low Sediment Supply') %>% 
  group_by(constraint_scenario, pressure) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>% 
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))

landward_sedL <- outcomes %>% 
  filter(var == 'LandwardMang' & model_scenario == 'Low Sediment Supply') %>% 
  group_by(constraint_scenario, pressure) %>% 
  summarise(Increase = sum(outcome>0)/n(),
            Neutral = sum(outcome==0)/n(),
            Decrease = sum(outcome<0)/n()) %>% 
  pivot_longer(Increase:Decrease ,names_to = 'outcome', values_to = 'prop') %>%  
  mutate(outcome = factor(outcome, levels = c('Increase', 'Neutral', 'Decrease')))


#landsea <- rbind(data.frame(seaward, mangrove = 'seaward'), data.frame(landward, mangrove = 'landward'))

# Note, as potential TODO could boostrap resample here to get estimate of uncertainty around that probability

a <- ggplot(seaward_sedH) +
  geom_bar(aes(y = pressure, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values=c("darkslategray3", 'darkolivegreen2', "brown3")) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  facet_wrap(~constraint_scenario) +
  ggtitle('A) Seaward mangroves, High sediment supply') +
  theme(legend.position = 'none', strip.text.x = element_text(size = 8))
a

b <- ggplot(landward_sedH) +
  geom_bar(aes(y = pressure, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values=c("darkslategray3", 'darkolivegreen2', "brown3")) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  facet_wrap(~constraint_scenario) +
  ggtitle('B) Landward mangroves, High sediment supply') +
  theme(legend.title = element_blank(),
        #axis.text.y =  element_blank(),
        strip.text.x = element_text(size = 8))
b

c <- ggplot(seaward_sedL) +
  geom_bar(aes(y = pressure, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values=c("darkslategray3", 'darkolivegreen2', "brown3")) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  facet_wrap(~constraint_scenario) +
  ggtitle('A) Seaward mangroves, Low sediment supply') +
  theme(legend.position = 'none', strip.text.x = element_text(size = 8))
c

d <- ggplot(landward_sedL) +
  geom_bar(aes(y = pressure, x = prop, fill = outcome),
           position = 'stack', stat = 'identity') +
  scale_fill_manual(values=c("darkslategray3", 'darkolivegreen2', "brown3")) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  xlab('Proportion of outcomes') +
  ylab('') +
  facet_wrap(~constraint_scenario) +
  ggtitle('B) Landward mangroves, Low sediment supply') +
  theme(legend.title = element_blank(),
        #axis.text.y =  element_blank(),
        strip.text.x = element_text(size = 8))
d

(a+b)/(c+d)

ggsave('outputs/outcomes.png', width = 17, height = 13)

# wrangle matrix weights

scnr_weights <- do.call(rbind, weights) %>% 
  mutate(pressure = rep(pressures, each = c(numsims*length(edge.labels(modelA))))) #%>% 
  #left_join(filter(outcomes, var == 'SeawardMang'))

scnr_weights$from_var <- sapply(strsplit(as.character(scnr_weights$param), "\\ "), `[`, 1)
scnr_weights$to_var <- sapply(strsplit(as.character(scnr_weights$param), "\\ "), `[`, 3)

# visualise edge weight distributions

paramsvis <- subset(scnr_weights, from_var != to_var)

ggplot() +
  geom_density(data = paramsvis, 
               aes(x = weight), colour = 'blue', fill = 'blue', alpha = 0.2) +
  #geom_density(data = filter(paramsvis, outcome == 0), 
   #            aes(x = weight), colour = 'yellow', fill = 'yellow', alpha = 0.2) +
  facet_wrap(~param, scales = 'free_y') +
  theme(strip.text.x = element_text(size = 5))

