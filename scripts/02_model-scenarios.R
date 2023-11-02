# determine probability of landward and seaward mangrove increase under different scenarios
# devtools::install_github("SWotherspoon/QPress",ref="Constrain")

library(QPress)
library(tidyverse)
library(patchwork)
source('scripts/helpers/models.R')
source('scripts/helpers/helpers_v2.R')

# loop through available model structures and store results in tmp files
tmp <- list()
tmp2 <- list()

for(l in seq_along(models)){
  
chosen_model <- models[[l]]
chosen_model_name <- names(models)[l]
  
# set up scenario simulations

set.seed(123) # set random number generator to make results reproducible
numsims <- 1000 # number of model simulations to run

# define perturbations scenarios
node.labels(chosen_model) # get model nodes for reference
press.scenarios <- list(c(SeaLevelRise=1), c(Cyclones=1), c(GroundSubsid=1), c(CoastalDev=1), c(Erosion=1), c(Drought=1), 
                        c(ExtremeRainfall=1), c(SeaLevelRise=1, Cyclones=1), c(SeaLevelRise=1, GroundSubsid=1), 
                        c(SeaLevelRise=1, CoastalDev=1), c(SeaLevelRise=1, Erosion=1),
                        c(SeaLevelRise=1, Drought=1), c(SeaLevelRise=1, ExtremeRainfall=1), c(SeaLevelRise=1, Cyclones=1, CoastalDev=1))
# define names of above scenarios for plot labelling later
press.labels <- c('Sea-level rise', 'Intense storms', 'Groundwater extraction', 'Coastal development', 'Erosion',
                  'Drought', 'Extreme rainfall', 'Sea-level rise & Intense storms', 'Sea-level rise & Groundwater extraction',
               'Sea-level rise & Coastal development', 'Sea-level rise & Erosion', 
               'Sea-level rise & Drought', 'Sea-level rise & Extreme rainfall', 'Sea-level rise & Intense storms & Coastal development')

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
names(edge.cons.scenarios) <- apply(labels.grid, 1, function(row) paste(row, collapse = ", ")) # label the list of edge constraint scenarios
# in each list element, repeat the middle character
edge.cons.scenarios <- lapply(edge.cons.scenarios, function(x) c(x[c(1,2)], x[2], x[3]))
# for each scenario, set up 'from' 'to' vectors specifying edges they apply to
from_vec <-  c('SeaLevelRise','LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')

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

system.time(
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
                              perturb = press.scenarios[[i]],
                              spatial = 'N',
                              arid = 'N',
                              prob = c(0.1, 0.5), # her assume probability of propagule establishment is 50% in all environments, and only 10% for landward propagules in arid environments
                              cdev = c(1, 1)) # here assume that coastal development to landward link is not uncertain (has a 100% probability)
      out[[i]] <- sim$stableoutcome
      stability[[i]] <- sim$stability.df
      #weights[[i]] <- sim$stableweights
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
)

# save potential stability in each scenario (i.e. proportion of stable matrices out of total simulated)
tmp[[l]] <- data.frame(model = chosen_model_name, constraint_scenario = names(rel.edge.cons.scenarios), do.call(rbind, stability.ls1))

# save outcomes in each scenario
tmp2[[l]] <- do.call(rbind, outcomes.ls1) %>% 
  mutate(model_scenario = rep(names(rel.edge.cons.scenarios), each = c(numsims*length(node.labels(chosen_model))*length(press.scenarios)*length(edge.cons.scenarios)))) %>% 
  mutate(model = chosen_model_name)
}

stability <- do.call(rbind, tmp)
outcomes <- do.call(rbind, tmp2)

write.csv(stability, 'outputs/simulation-outcomes/stability.csv', row.names = F)
saveRDS(outcomes, 'outputs/simulation-outcomes/outcomes.rds')

