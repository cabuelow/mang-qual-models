SedSupply_SLF <- press.impact(modelA, perturb=c(LandwardLatAccom=-1,
                                 SeawardLatAccom=1,
                                 LandwardPropag=-1,
                                 LandwardMang=-1,
                                 SubVol=1))
SLR <- press.impact(modelA, perturb=c(LandwardLatAccom=1,
                                 SeawardLatAccom=-1,
                                 SeawardMang=-1))
SedSupply <- press.impact(modelA, perturb=c(SubVol=1))
SedSupply_SLR <- press.impact(modelA, perturb=c(LandwardLatAccom=1,
                                 SeawardLatAccom=-1,
                                 SeawardPropag=-1,
                                 SeawardMang=-1,
                                 SubVol=1))

simulate_scenarios <- function(i){
  
  simB <- community.sampler(modelA)
  varnames <- colnames(adjacency.matrix(modelA, labels = T))
  
  stable_comm <- FALSE
  
  while(!stable_comm){
    w <- simB$community()
    stable_comm <- stable.community(w)
  }
  
  modout <- data.frame(vars = rep(varnames, 4),
                       isim = i,
                       scnr = rep(c("Sediment supply", 
                                "SLR", 
                                "Sediment supply & SLF",
                                "Sediment supply & SLR"),
                                each = length(varnames)),
                       outcomes = c(SedSupply(w),
                                    SLR(w),
                                    SedSupply_SLF(w),
                                    SedSupply_SLR(w))
  )
  
  return(modout)
}

sim_valid_scenarios <- function(i){
  
  simB <- community.sampler(modelA)
  varnames <- colnames(adjacency.matrix(modelA, labels = T))
  
  stable_comm <- FALSE
  
  while(!stable_comm){
    w <- simB$community()
    stable_comm <- stable.community(w)
  }

ws <- simB$weights(w)
names(ws) <- simB$weight.labels

params <- data.frame(isim = i,
                     links = simB$weight.labels, 
                     weight = simB$weights(w))

}