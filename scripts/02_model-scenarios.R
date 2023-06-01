# determine probability of landward and seaward mangrove increase under different scenarios
# devtools::install_github("SWotherspoon/QPress",ref="Constrain")

library(QPress)
library(tidyverse)
library(patchwork)
source('scripts/helpers/models_v2.R')
source('scripts/helpers/helpers_v2.R')

# which model do you want to run?

names(models) # names of available models
chosen_model <- models$mangrove_model
chosen_model_name <- 'mangrove_model'
  
# set up scenario simulations

set.seed(123) # set random number generator to make results reproducible
numsims <- 1000 # number of model simulations to run

# define perturbations scenarios
node.labels(chosen_model) # get model nodes for reference
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
edge.labels(chosen_model) # get model edges for reference
# create a grid of all combinations of high, medium and low for tidal range, propagule establishment capacity and coastal squeeze
edge.cons.grid <- expand.grid(TidalRange = c('H', 'M', 'L'),
                              PropEstab = c('H', 'M', 'L'),
                              CoastalSqueeze = c('H', 'M', 'L'))
# turn into a list
edge.cons.scenarios <- as.list(as.data.frame(t(edge.cons.grid)))
# label the list
labels.grid <- expand.grid(TidalRange = c('Microtidal', 'Mesotidal', 'Macrotidal'),
                           PropEstab = c('High propagule establishment capacity', 'Medium propagule establishment capacity', 'Low propagule establishment capacity'),
                           CoastalSqueeze = c('Low coastal squeeze', 'Medium coastal squeeze', 'High coastal squeeze')) # note that a low value for the interaction strength = High coastal squeeze and vice-versa
names(edge.cons.scenarios) <- apply(labels.grid, 1,function(row) paste(row, collapse = ", ")) # label the list of edge constraint scenarios
# for each scenario, set up 'from' 'to' vectors specifying edges they apply to
from_vec <-  c('SeaLevelRise','LandwardAvailableProp', 'SeaLevelRise')
to_vec <- c('SeawardMang', 'LandwardMang', 'LandwardMang')

# define relative edge constraints - which edge interaction strengths are greater than other
# in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
# under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
# the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                        parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 

# loop through scenarios defined above and store model outcomes
# do a nested for loop, iterate over relative edge constraint scenarios, absolute edge constraint scenarios, then perturbation scenarios

stability.ls1 <- list() # list for storing model stability results
outcomes.ls1 <- list() # list for storing model outcome results

for(k in seq_along(rel.edge.cons.scenarios)){
  stability.ls <- list()
  outcomes.ls <- list()
  model <- rel.edge.cons.scenarios[[k]]
  for(j in seq_along(edge.cons.scenarios)){
    out <- list()
    stability <- list()
    weights <- list()
    class.con <- edge.cons.scenarios[[j]]
    for(i in seq_along(press.scenarios)){
      sim <- system.sim_press(numsims, constrainedigraph = model, 
                              from = from_vec,
                              to = to_vec,
                              class = class.con,
                              perturb = press.scenarios[[i]])
      out[[i]] <- sim$stableoutcome
      stability[[i]] <- sim$stability.df
      weights[[i]] <- sim$stableweights
    } # end i loop
    # potential stability (i.e. proportion of stable matrices out of total simulated)
    stability.ls[[j]] <- data.frame(press.labels, do.call(rbind, stability))
    # outcomes
    outcomes.ls[[j]] <- do.call(rbind, out) %>% 
      mutate(pressure = rep(press.labels, each = c(numsims*length(node.labels(model$edges)))))
  } # end j loop
  # potential stability (i.e. proportion of stable matrices out of total simulated)
  stability.ls1[[k]] <- data.frame(constraint_scenario = names(edge.cons.scenarios), do.call(rbind, stability.ls))
  # outcomes
  outcomes.ls1[[k]] <- do.call(rbind, outcomes.ls) %>% 
    mutate(constraint_scenario = rep(names(edge.cons.scenarios), each = c(numsims*length(node.labels(model$edges))*length(press.scenarios)))) 
} # end k loop

# save potential stability in each scenario (i.e. proportion of stable matrices out of total simulated)
stability <- data.frame(constraint_scenario = names(rel.edge.cons.scenarios), do.call(rbind, stability.ls1))
write.csv(stability, paste0('outputs/simulation-outcomes/stability_', chosen_model_name, '.csv'), row.names = F)

# save outcomes in each scenario
outcomes <- do.call(rbind, outcomes.ls1) %>% 
  mutate(model_scenario = rep(names(rel.edge.cons.scenarios), each = c(numsims*length(node.labels(chosen_model))*length(press.scenarios)*length(edge.cons.scenarios))))
saveRDS(outcomes, paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '.rds'))

