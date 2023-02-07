# determine probability of landward and seaward mangrove increase under different scenarios
#devtools::install_github("SWotherspoon/QPress",ref="Constrain")
# TODO: make scenario simulation more efficient and easier to implement

library(igraph)
library(QPress)
library(tidyverse)
library(patchwork)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(scales)
theme_set(theme_classic())
source('scripts/models.R')
source('scripts/helpers.R')

# set up scenario simulations

set.seed(123)
numsims <- 1000

# perturbation scenarios

pressures <- c('Sea-level rise', 'Cyclones', 'Groundwater extraction', 'Coastal development', 'Erosion', 'Drought or Dams',
               'Sea-level rise & Cyclones', 'Seal-level rise & Groundwater extraction',
               'Sea-level rise  & Coastal development', 'Sea-level rise & Erosion', 'Sea-level rise & Drought or Dams')
press.scenarios <- list(c(SeaLevelRise=1), c(Cyclones=1), c(GroundSubsid=1), c(CoastalDev=1), 
                        c(Erosion=1), c(Sediment=-1),
                        c(SeaLevelRise=1, Cyclones=1),
                        c(SeaLevelRise=1, GroundSubsid=1),
                        c(SeaLevelRise=1, CoastalDev=1), 
                        c(SeaLevelRise=1, Erosion=1),
                        c(SeaLevelRise=1, Sediment=-1))

# edge constraint scenarios
# **TODO: make this easier by changing how the function takes these constraints....

con.scenarios <- list(c('H', 'H', 'H'),
                      c('M', 'H', 'H'),
                      c('L', 'H', 'H'),
                      c('H', 'L', 'L'),
                      c('M',  'L', 'L'),
                      c('L', 'L', 'L'))
names(con.scenarios) <- c('Microtidal, High Hydro-connectivity',
                          'Mesotidal, High Hydro-connectivity',
                          'Macrotidal, High Hydro-connectivity',
                          'Microtidal, Low Hydro-connectivity',
                          'Mesotidal, Low Hydro-connectivity',
                          'Macrotidal, Low Hydro-connectivity')

# set up relative edge constraint scenarios
# high sed supply model, sediment -> subVol will be greater than SLR neg interactions
# vice versa for low sed supply model

model.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), modelB.1),
                        parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), modelB.1))
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
                                       'LandwardAvailableProp', 'SeawardAvailableProp'),
                              to = c('SeawardMang', #'SeawardPropag', #'SeawardEstabSpace', 
                                     'LandwardMang', 'SeawardMang'),
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

# wrangle outcomes
outcomes <- do.call(rbind, outcomes.ls1) %>% 
  mutate(model_scenario = rep(names(model.scenarios), each = c(numsims*length(node.labels(modelB))*length(pressures)*length(con.scenarios)))) 
saveRDS(outcomes, 'outputs/outcomes_alt-model.rds')

outcomes2 <- readRDS('outputs/outcomes.rds')

cyclone_pos <- outcomes %>% 
  filter(var %in% c('SeawardMang', 'LandwardMang') & 
           constraint_scenario %in% c('Macrotidal, High Hydro-connectivity',
                                      'Mesotidal, High Hydro-connectivity',
                                      'Microtidal, High Hydro-connectivity')) %>% 
  # pressure %in% c('Sea-level rise', 'Sea-level rise & Cyclones', 
  #               'Seal-level rise & Groundwater extraction', 'Sea-level rise  & Coastal development', 
  #               'Sea-level rise & Erosion', 'Sea-level rise & Drought or Dams')) %>% 
  group_by(model_scenario, constraint_scenario, pressure, var) %>% 
  summarise(Prob_gain_neutral = (sum(outcome>0) + sum(outcome==0))/n()) %>%
  mutate(constraint_scenario = recode(constraint_scenario, 'Macrotidal, High Hydro-connectivity' = 'Macrotidal',
                                      'Mesotidal, High Hydro-connectivity' = 'Mesotidal',
                                      'Microtidal, High Hydro-connectivity' = 'Microtidal'),
         var = recode(var, 'LandwardMang' = 'Landward mangrove', 'SeawardMang' = 'Seaward mangrove')) %>% 
  mutate(constraint_scenario = factor(constraint_scenario, levels = c('Microtidal', 'Mesotidal', 'Macrotidal')))

cyclone_pos$Prob_change <- rescale(cyclone_pos$Prob_gain_neutral, to = c(-100, 100))

cyclone_neg <- outcomes2 %>% 
  filter(var %in% c('SeawardMang', 'LandwardMang') & 
           constraint_scenario %in% c('Macrotidal, High Hydro-connectivity',
                                      'Mesotidal, High Hydro-connectivity',
                                      'Microtidal, High Hydro-connectivity')) %>% 
  # pressure %in% c('Sea-level rise', 'Sea-level rise & Cyclones', 
  #               'Seal-level rise & Groundwater extraction', 'Sea-level rise  & Coastal development', 
  #               'Sea-level rise & Erosion', 'Sea-level rise & Drought or Dams')) %>% 
  group_by(model_scenario, constraint_scenario, pressure, var) %>% 
  summarise(Prob_gain_neutral = (sum(outcome>0) + sum(outcome==0))/n()) %>%
  mutate(constraint_scenario = recode(constraint_scenario, 'Macrotidal, High Hydro-connectivity' = 'Macrotidal',
                                      'Mesotidal, High Hydro-connectivity' = 'Mesotidal',
                                      'Microtidal, High Hydro-connectivity' = 'Microtidal'),
         var = recode(var, 'LandwardMang' = 'Landward mangrove', 'SeawardMang' = 'Seaward mangrove')) %>% 
  mutate(constraint_scenario = factor(constraint_scenario, levels = c('Microtidal', 'Mesotidal', 'Macrotidal')))

cyclone_neg$Prob_change <- rescale(cyclone_neg$Prob_gain_neutral, to = c(-100, 100))

# compare predictions with cyclones pos vs. negative

cyclone_neg$Prob_change_pos <- cyclone_pos$Prob_change

ggplot(cyclone_neg) +
  geom_point(aes(x = Prob_change, y = Prob_change_pos, col = var))
