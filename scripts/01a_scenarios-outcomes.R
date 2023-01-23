#install.packages(c("tcltk2","XML","devtools"))
#devtools::install_github("SWotherspoon/QPress",ref="Constrain")
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

# constrain edge weights for micro, meso, macro tides
# high, med, low connectivity
# high, med, low sediment delivery
# don't need to parse constraints - that just says one edge WITHIN model is greater than another
# instead paramaterise the models differently - instead of uniform dist from 0-1, make is 0-0.3, etc, etc

modconst <- parse.constraints(c(
  ""
))
# then make sure dealing with uncertain edges properly

# parse.constraints
# community.ordering.sampler
# system.simulate0

set.seed(17)
model2 <- parse.digraph(c(
  "DP*->GR", #trophic
  "DP->DMS", #phytoplankton produce DMS
  "GR->DMS", #grazing enhances release of DMS
  "GR*->PR", #trophic
  "DMS->PR" #predators are attracted to DMS
))
model2<- enforce.limitation(model2) #all nodes are self-limited

#constrain DMS attraction to contribute less to predators than consumption of grazers
constrained <- parse.constraints(c(
  "DMS -> PR < GR -> PR"),
  model2)
s <- community.ordering.sampler(constrained)
