# functions for simulating spatial models and running validation
source('scripts/helpers/helpers_v2.R')

# calculate accuracy
calc_accuracy <- function(x, x2){ # x is vector of predictions, x2 is reference vector
  cont.table <- confusionMatrix(factor(x), factor(x2))$table # contingency table
  users <- diag(cont.table)/rowSums(cont.table)*100
  producers <- diag(cont.table)/colSums(cont.table)*100
  overall.accuracy <- sum(diag(cont.table))/sum(cont.table)*100
  class.df <- data.frame(class = levels(factor(x2)), Overall_accuracy = overall.accuracy, Producers_accuracy = producers, Users_accuracy = users)
  accuracy_list <- list(class.df, cont.table)
  names(accuracy_list) <- c('accuracy.results', 'contingency.table')
  return(accuracy_list)
}

# simulate model matrices and store outcomes

sim_mod <- function(x, # spatial unit
                    numsims # number of simulations
){ # define model relative edge constraints - which edge interaction strengths are greater than other
  # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
  # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
  # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
  rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                  parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
  names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 

    # edge constraint scenarios
    
    if(x$csqueeze == 'None'){
      datselect2 <- dplyr::select(x, Tidal_Class, prop_estab) %>% 
        pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
      from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
      to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
      con.scenario <- c(datselect2$vals, datselect2$vals[2])
    }else{
      datselect2 <- dplyr::select(x, Tidal_Class, prop_estab, csqueeze) %>% 
        pivot_longer(Tidal_Class:csqueeze, names_to = 'press', values_to = 'vals') 
      from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
      to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
      con.scenario <- c(datselect2$vals, datselect2$vals[2], datselect2$vals[3])
      con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
    }
    
    # select model for sediment supply
    
    if(x$sed_supp == 'H'){
      model <- rel.edge.cons.scenarios[[1]]
    }else{
      model <- rel.edge.cons.scenarios[[2]]
    }
    
    # simulate outcomes
    
    sim <- system.sim(numsims, constrainedigraph = model, 
                            from = from_vec,
                            to = to_vec,
                            class = con.scenario, spatial = 'Y')
    return(sim)
  }