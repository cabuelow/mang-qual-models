# functions for simulating spatial models and running validation

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

# make spatial hindcast or forecast

cast <- function(x, # spatial unit
                 numsims, # number of simulations
                 type # hindcast or forecast?
                 ){
  if(type == 'forecast'){
    
    # define model relative edge constraints - which edge interaction strengths are greater than other
    # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
    # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
    # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
    rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                    parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
    names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
    
    # perturbations for forecasting
    datselect <- dplyr::select(x, fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
      pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
      filter(vals == 1) %>% 
      mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                            'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones'))
    
    if(nrow(datselect) != 0){ # if there are no perturbations, go to next typology
      
    press.scenario <- rep(1, nrow(datselect))
    names(press.scenario) <- datselect$press
    
    # edge constraint scenarios
    
    if(x$fut_csqueeze == 'None'){
      datselect2 <- dplyr::select(x, Tidal_Class, prop_estab) %>% 
        pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
      from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
      to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
      con.scenario <- c(datselect2$vals, datselect2$vals[2])
    }else{
      datselect2 <- dplyr::select(x, Tidal_Class, prop_estab, fut_csqueeze) %>% 
        pivot_longer(Tidal_Class:fut_csqueeze, names_to = 'press', values_to = 'vals') 
      from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
      to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
      con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
    }
    
    # select model for sediment supply
    
    if(x$fut_dams == 'L'){
      model <- rel.edge.cons.scenarios[[2]]
    }else if(x$sed_supp == 'H'){
      model <- rel.edge.cons.scenarios[[1]]
    }else{
      model <- rel.edge.cons.scenarios[[2]]
    }
    
    # simulate outcomes
    
    sim <- system.sim_press(numsims, constrainedigraph = model, 
                            from = from_vec,
                            to = to_vec,
                            class = con.scenario,
                            perturb = press.scenario,
                            spatial = 'Y')
    
    out <- sim$stableoutcome %>% 
      filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
      group_by(var) %>% 
      summarise(Prob_gain = (sum(outcome>0)/n())*100,
                Prob_neutral = (sum(outcome==0)/n())*100,
                Prob_loss = (sum(outcome<0)/n())*-100) %>% 
      mutate(Type = x$Type, cast = 'forecast')
    
    return(out)
    }
  }else{
    
    # define model relative edge constraints - which edge interaction strengths are greater than other
    # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
    # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
    # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
    rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                    parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
    names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
    
    # perturbations for hindcasting
    datselect <- dplyr::select(x, csqueeze_1, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
      pivot_longer(csqueeze_1:storms, names_to = 'press', values_to = 'vals') %>% 
      filter(vals == 1) %>% 
      mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                            'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
    
    if(nrow(datselect) != 0){ # if there are no perturbations, go to next typology

    press.scenario <- rep(1, nrow(datselect))
    names(press.scenario) <- datselect$press
    
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
    
    sim <- system.sim_press(numsims, constrainedigraph = model, 
                            from = from_vec,
                            to = to_vec,
                            class = con.scenario,
                            perturb = press.scenario,
                            spatial = 'Y')
    
    out <- sim$stableoutcome %>% 
      filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
      group_by(var) %>% 
      summarise(Prob_gain = (sum(outcome>0)/n())*100,
                Prob_neutral = (sum(outcome==0)/n())*100,
                Prob_loss = (sum(outcome<0)/n())*-100) %>% 
      mutate(Type = x$Type, cast = 'hindcast')
    
    return(out)
    }
  }
}

# train model and obtain calibrated (i.e. valid/invalid) interaction coefficients/strengths

train <- function(x, numsims){ # x is spatial training units
  
  # define model relative edge constraints - which edge interaction strengths are greater than other
  # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
  # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
  # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
  rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                  parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
  names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
  
  datselect <- dplyr::select(x, csqueeze_1, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
    pivot_longer(csqueeze_1:storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                          'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
  datmonitor <- dplyr::select(x, land_net_change_obs, sea_net_change_obs) %>% 
    pivot_longer(land_net_change_obs:sea_net_change_obs, names_to = 'monitor', values_to = 'vals') %>% 
    mutate(monitor = recode(monitor, 'land_net_change_obs' = 'LandwardMang', 'sea_net_change_obs' = "SeawardMang"))
  if(nrow(datselect) != 0){ # if there are no perturbations, go to next typology
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  monitor.scenario <- datmonitor$vals
  names(monitor.scenario) <- datmonitor$monitor
  
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
  
  sim <- system.sim_press_valid(numsims, constrainedigraph = model, 
                                from = from_vec,
                                to = to_vec,
                                class = con.scenario,
                                perturb = press.scenario,
                                monitor = monitor.scenario,
                                spatial = 'Y')
  
  outcomes <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(valid, pressures, var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    pivot_wider(id_cols = c('pressures'), names_from = 'var', values_from = c('Prob_gain', 'Prob_neutral', 'Prob_loss'))
  
  return(data.frame(sim$stableweights %>% 
                                  mutate(Type = x$Type,
                                         cast = 'training',
                                        # kfold = kfold, 
                                  ), outcomes[,-1], 
                                pivot_wider(datmonitor, names_from = 'monitor', values_from = 'vals', names_prefix = 'Observed_')))
}
  }

# test model using calibrated (i.e. valid/invalid) interaction coefficients/strengths

test <- function(x, numsims, params){ # params is calibrated parameters
  
  # define model relative edge constraints - which edge interaction strengths are greater than other
  # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
  # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
  # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
  rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                  parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
  names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
  
  datselect <- dplyr::select(x, csqueeze_1, ant_slr, gwsub, hist_drought, hist_ext_rain, storms) %>% 
    pivot_longer(csqueeze_1:storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'csqueeze_1' = 'CoastalDev', 'ant_slr' = "SeaLevelRise", 'gwsub' = "GroundSubsid", 
                          'hist_drought' = 'Drought', 'hist_ext_rain' = 'ExtremeRainfall', 'storms' = 'Cyclones'))
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario), collapse = '_') & Geomorphology == x$Geomorphology)
  
  if(nrow(datselect) != 0 & nrow(validweights) != 0){ # if there are no valid weights for this combination of pressures, go to next

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
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                           from = from_vec,
                           to = to_vec,
                           class = con.scenario,
                           perturb = press.scenario,
                           weights = validweights,
                           spatial = 'Y')
  
  out <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = x$Type)
  
  return(out)
  }
}

# make a calibrated forecast

forecast_calibrated <- function(x, numsims, params){
  
  # define model relative edge constraints - which edge interaction strengths are greater than other
  # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
  # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
  # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
  rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                  parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
  names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
  
  datselect <- select(x, fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
    pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                          'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones'))
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario), collapse = '_') & Geomorphology == x$Geomorphology)
  if(nrow(datselect) != 0 & nrow(validweights) != 0){ # if there are no valid weights for this combination of pressures, go to next

  # edge constraint scenarios
  
  if(x$fut_csqueeze == 'None'){
    datselect2 <- select(x, Tidal_Class, prop_estab) %>% 
      pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
    from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
    to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
    con.scenario <- c(datselect2$vals, datselect2$vals[2])
  }else{
    datselect2 <- select(x, Tidal_Class, prop_estab, fut_csqueeze) %>% 
      pivot_longer(Tidal_Class:fut_csqueeze, names_to = 'press', values_to = 'vals') 
    from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
    to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
    con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
  }
  
  # select model for sediment supply
  
  if(x$fut_dams == 'L'){
    model <- rel.edge.cons.scenarios[[2]]
  }else if(x$sed_supp == 'H'){
    model <- rel.edge.cons.scenarios[[1]]
  }else{
    model <- rel.edge.cons.scenarios[[2]]
  }
  
  # simulate outcomes
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                           from = from_vec,
                           to = to_vec,
                           class = con.scenario,
                           perturb = press.scenario,
                           weights = validweights,
                           spatial = 'Y')
  
  out <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = x$Type)
  
  return(out)
  }
}

# make a calibrated forecast with increased propagules (i.e., management/conservation/restoration)

forecast_calibrated_manage <- function(x, numsims, params){
  
  # define model relative edge constraints - which edge interaction strengths are greater than other
  # in all models the seaward mangrove -> substrate vol interaction strength is greater than the landward mangrove -> substrate vol interaction strength
  # under a high sediment supply scenario, the sediment -> subVol interaction strengths will be greater than 
  # the negative interaction between sea level rise -* and seaward mangroves; vice versa for the low sediment supply model
  rel.edge.cons.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model),
                                  parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), chosen_model))
  names(rel.edge.cons.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply') # label the list of relative edge constraint scenarios 
  
  datselect <- select(x, fut_csqueeze_1, fut_slr, fut_gwsub, fut_drought, fut_ext_rain, fut_storms) %>% 
    pivot_longer(fut_csqueeze_1:fut_storms, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'fut_csqueeze_1' = 'CoastalDev', 'fut_slr' = "SeaLevelRise", 'fut_gwsub' = "GroundSubsid", 
                          'fut_drought' = 'Drought', 'fut_ext_rain' = 'ExtremeRainfall', 'fut_storms' = 'Cyclones')) %>% 
    rbind(data.frame(press = c('LandwardAvailableProp', 'SeawardAvailableProp'), vals = 1))
  press.scenario <- rep(1, nrow(datselect))
  names(press.scenario) <- datselect$press
  validweights <- params %>% 
    filter(pressures == paste0(names(press.scenario[-c(length(press.scenario)-1, length(press.scenario))]), collapse = '_') & Geomorphology == x$Geomorphology)
  if(nrow(datselect) != 0 & nrow(validweights) != 0){ # if there are no valid weights for this combination of pressures, go to next
    
  # edge constraint scenarios
  
  if(x$fut_csqueeze == 'None'){
    datselect2 <- select(x, Tidal_Class, prop_estab) %>% 
      pivot_longer(Tidal_Class:prop_estab, names_to = 'press', values_to = 'vals') 
    from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp')
    to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang')
    con.scenario <- c(datselect2$vals, datselect2$vals[2])
  }else{
    datselect2 <- select(x, Tidal_Class, prop_estab, fut_csqueeze) %>% 
      pivot_longer(Tidal_Class:fut_csqueeze, names_to = 'press', values_to = 'vals') 
    from_vec <- c('SeaLevelRise', 'LandwardAvailableProp', 'SeawardAvailableProp', 'SeaLevelRise')
    to_vec <- c('SeawardMang', 'LandwardMang', 'SeawardMang', 'LandwardMang')
    con.scenario <- c(datselect2$vals[c(1,2)], datselect2$vals[2], datselect2$vals[3])
  }
  
  # select model for sediment supply
  
  if(x$fut_dams == 'L'){
    model <- rel.edge.cons.scenarios[[2]]
  }else if(x$sed_supp == 'H'){
    model <- rel.edge.cons.scenarios[[1]]
  }else{
    model <- rel.edge.cons.scenarios[[2]]
  }
  
  # simulate outcomes
  
  sim <- system.sim_press2(numsims, constrainedigraph = model, 
                           from = from_vec,
                           to = to_vec,
                           class = con.scenario,
                           perturb = press.scenario,
                           weights = validweights,
                           spatial = 'Y')
  
  out <- sim$stableoutcome %>% 
    filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
    group_by(var) %>% 
    summarise(Prob_gain = (sum(outcome>0)/n())*100,
              Prob_neutral = (sum(outcome==0)/n())*100,
              Prob_loss = (sum(outcome<0)/n())*-100) %>% 
    mutate(Type = x$Type)
  
  return(out)
  }
}