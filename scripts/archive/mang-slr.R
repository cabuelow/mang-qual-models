# mangrove response to SLR

library(igraph)
library(QPress)
library(dplyr)
library(ggplot2)

# create conceptual model - all scenario models together
# make 4 different models from Kerrylee, see if they predict what we expect

# conceptual model with uncertain edges
# only negative links, because an increase in one variable cannot have 
# a positive effect on another - unless have mangrove dieback in the front
# front and back should just be one variable?

model <- parse.digraph(c("Landward--*Back", # landward can encroach on back with mangrove dieback, but uncertain
                'Back--*Landward', # back can encroach on landward with sea level rise, but uncertain
                'Back--*Front', # back can encroach on front
                'Front--*Back', # front can encroach on back
                'Front--*Seaward', # front can prograde
                'Seaward--*Front' # seaward can move to front
                # self regulation, back and front can't grow forever
                #'Back-*Back',
                #'Front-*Front'
                ))

model <- enforce.limitation(model)
wA <- adjacency.matrix(model, labels=TRUE)
wA
stable.community(wA) #stable model?

# plot model

plotmodel <- function(model){
  w <- adjacency.matrix(model, labels=TRUE)
  gdat <- graph_from_adjacency_matrix(t(w),
                                      mode = "directed",
                                      weighted = "b")
  E(gdat)$color <- ifelse(E(gdat)$b < 0, "red", "black")
  V(gdat)$color <- 'white'
  plot(gdat, 
       vertex.size = 50,
       vertex.label.cex = 1,
       vertex.label.color = "black",
       edge.arrow.size = 0.5,
       #edge.label = letters[1:length(E(gdat))],
       edge.label.cex = 0.9,
       edge.label.color = "black",
       edge.label.pos = 1)
}
plotmodel(model)

# does the model behave as expected under different scenarios?
# maybe we do need an 'accommodation space variable" 
# that represents the lateral and vertical?

# check scenarios work as expected

# scenario 1: high sediment supply
# expect increase in front mangroves, decrease in seaward

impactSedSupply <- press.impact(model, perturb = c(Front=1))
round(impactSedSupply(wA),2)

# scenario 2: reduction in hydroperiod (drought or sea level fall) 
# expect dieback in front mangroves only, landward increase

impactSeaLevelFall <- press.impact(model, perturb = c(Landward=1))
round(impactSeaLevelFall(wA),2)

# scenario 3: combined high sediment supply and reduction in hydroperiod
# expect increase decrease in back mangroves, increase in front

impactSeaLevelFallSedSupply <- press.impact(model, perturb = c(Landward = 1,
                                                               Front = 1))
round(impactSeaLevelFallSedSupply(wA),2)

# scenario 4: sea level rise
# expect back mangroves increase, landward decreases

impactSeaLevelRise <- press.impact(model, perturb = c(Back=1))
round(impactSeaLevelRise(wA),2)

# scenario 5: SLR and low sediment supply

impactSeaLevelRise <- press.impact(model, perturb = c(Back=1))
round(impactSeaLevelRise(wA),2)

# the problem is that, depending on the driver, there will be some links or there won't be,
# or some links will be heavier than others
# should really just get rid of front and back of mangroves?
# maybe can't make one model to rule them all? Just to locally and scenario-specific models?
# so need different conceptual models for different locations?
# but whats the use in that?

# or need to think about the 3-D more?
# look at Kerrylee's figures, 
# add in accomm. space variable, SubVol variable?

# do/can we do conditional weighting on the variables, so if sea level rise is the driver,
# then weighting on negative arrow to landward is higher than arrow from back to front?

simulate_scenarios <- function(i){
  
  simB <- community.sampler(model)
  
  wA <- adjacency.matrix(model, labels=TRUE)
  varnames <- names(impactSedSupply(wA))
  
  stable_comm <- FALSE
  
  while(!stable_comm){
    w <- simB$community()
    stable_comm <- stable.community(w)
  }
  
  modout <- data.frame(vars = rep(varnames, 4), 
                       isim = i,
                       scnr = rep(c("SedSupply"), 
                                  each = 4),
                       outcomes = c(impactSedSupply(w))
  )
  return(modout)
}
modout <- lapply(1:1000, simulate_scenarios)
modout <- do.call("rbind", modout)

#summarize results

modsum <- modout %>% group_by(vars, scnr) %>%
  summarize(net_change = sum(outcomes>0)/n())
ggplot(modsum) +
  aes(x = scnr, y = net_change) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~vars) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  ylab("Proportion positive change") +
  xlab("Scenario") +
  geom_hline(yintercept = 0.5)

mang_dens <- modout %>% filter(vars == "Seaward")
ggplot(mang_dens) +
  aes(x = outcomes) + 
  geom_density(color = "red", fill = "pink") +
  facet_wrap(~scnr, scales= "free") + 
  xlim(-10, 10) + 
  xlab("Net outcome") +
  geom_vline(xintercept = 0)

mang_dens <- modout %>% filter(vars == "Front")
ggplot(mang_dens) +
  aes(x = outcomes) + 
  geom_density(color = "red", fill = "pink") +
  facet_wrap(~scnr, scales= "free") + 
  xlim(-10, 10) + 
  xlab("Net outcome") +
  geom_vline(xintercept = 0)

mang_dens <- modout %>% filter(vars == "Back")
ggplot(mang_dens) +
  aes(x = outcomes) + 
  geom_density(color = "red", fill = "pink") +
  facet_wrap(~scnr, scales= "free") + 
  xlim(-10, 10) + 
  xlab("Net outcome") +
  geom_vline(xintercept = 0)

mang_dens <- modout %>% filter(vars == "Landward")
ggplot(mang_dens) +
  aes(x = outcomes) + 
  geom_density(color = "red", fill = "pink") +
  facet_wrap(~scnr, scales= "free") + 
  xlim(-10, 10) + 
  xlab("Net outcome") +
  geom_vline(xintercept = 0)

