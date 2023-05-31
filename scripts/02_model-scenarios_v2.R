# determine probability of landward and seaward mangrove increase under different scenarios
# devtools::install_github("SWotherspoon/QPress",ref="Constrain")

library(QPress)
library(tidyverse)
library(patchwork)
source('scripts/helpers/models_v2.R')
source('scripts/helpers/helpers.R')

# set up scenario simulations

set.seed(123) # set random number generator to make results reproducible
numsims <- 1000 # number of model simulations to run

# define perturbations scenarios
node.labels(model) # get model nodes for reference
press.scenarios <- list(c(SeaLevelRise=1), c(Cyclones=1), c(GroundSubsid=1), c(CoastalDev=1), 
                        c(Erosion=1), c(Drought=1), c(ExtremeRainfall=1), c(SeaLevelRise=1, Cyclones=1),
                        c(SeaLevelRise=1, GroundSubsid=1), c(SeaLevelRise=1, CoastalDev=1), c(SeaLevelRise=1, Erosion=1),
                        c(SeaLevelRise=1, Drought=1), c(SeaLevelRise=1, ExtremeRainfall=1))
# define names of above scenarios for plot labelling later
press.labels <- c('Sea-level rise', 'Cyclones', 'Groundwater extraction', 'Coastal development', 'Erosion','Drought',
               'Extreme rainfall', 'Sea-level rise & Cyclones', 'Sea-level rise & Groundwater extraction',
               'Sea-level rise  & Coastal development', 'Sea-level rise & Erosion', 'Sea-level rise & Drought', 
               'Sea-level rise & Extreme rainfall')

# define edge constraints - whether edge interaction strengths should be 'high', 'medium' or 'low'
edge.labels(model) # get model edges for reference
# create grid all combinations of high, medium and low for tidal range, propagule establishment capacity and coastal squeeze
edge.cons.grid <- expand.grid(TidalRange = c('H', 'M', 'L'),
                              PropEstab = c('H', 'M', 'L'),
                              CoastalSqueeze = c('H', 'M', 'L'))
# turn into a list
edge.cons.scenarios <- as.list(as.data.frame(t(edge.cons.grid)))
# label the list
labels.grid <- expand.grid(TidalRange = c('Microtidal', 'Mesotidal', 'Macrotidal'),
                           PropEstab = c('High propagule establishment capacity', 'Medium propagule establishment capacity', 'Low propagule establishment capacity'),
                           CoastalSqueeze = c('High coastal squeeze', 'Medium coastal squeeze', 'Low coastal squeeze'))
names(edge.cons.scenarios) <- apply(d, 1,function(row) paste(row, collapse = ", ")) # label the list of edge constraint scenarios


# define relative edge constraints - which edge interaction strengths are greater than other
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

