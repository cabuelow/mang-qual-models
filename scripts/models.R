library(QPress)
library(dagitty)
library(ggdag)
library(igraph)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# build model signed digraphs

# mangrove system + local factors + anthropogenic pressures + climate pressures
modelA <- parse.digraph(c('LandwardEstabSpace->LandwardPropag',
                          'LandwardPropag->LandwardMang',
                          'LandwardMang->OM',
                          'SeawardEstabSpace->SeawardPropag',
                          'SeawardPropag->SeawardMang',
                          'SeawardMang->OM',
                          'OM->SubVol',
                          'SubVol->SeawardEstabSpace', 
                          'LandwardAvailableProp->LandwardPropag',
                          'SeawardAvailableProp->SeawardPropag',
                          'AccommSpace->LandwardEstabSpace',
                          'AccommSpace->SeawardEstabSpace',
                          #'HydroEnergy-*SeawardEstabSpace',
                          #'HydroEnergy->Erosion',
                          'Erosion-*SubVol',
                          'Erosion-*SeawardPropag',
                          'Erosion-*SeawardMang',
                          #'Precipitation->Sediment',
                          'Dams-*Sediment',
                          'Sediment->SubVol',
                          #'Autocompaction-*SubVol',
                          #'TidalFreq-*SeawardPropag',
                          #'TidalFreq->SubVol',
                          'CoastalDev-*LandwardEstabSpace',
                          'CoastalDev-*LandwardPropag',
                          'CoastalDev-*LandwardMang',
                          'CoastalDev--*SeawardMang',
                          'CoastalDev--*SeawardPropag',
                          'GroundSubsid->LandwardEstabSpace',
                          'GroundSubsid-*SeawardEstabSpace',
                          'GroundSubsid-*SeawardPropag',
                          'GroundSubsid-*SeawardMang',
                          'SeaLevelRise->LandwardEstabSpace',
                          'SeaLevelRise-*SeawardEstabSpace',
                          'SeaLevelRise-*SeawardPropag',
                          'SeaLevelRise-*SeawardMang',
                          'Cyclones--*LandwardEstabPropag',
                          'Cyclones--*LandwardMang',
                          'Cyclones-*SeawardMang',
                          'Cyclones-*SeawardEstabPropag',
                          'Cyclones-*SubVol',
                          'Cyclones-->LandwardAvailableProp'
)) %>% 
  enforce.limitation()

grViz(grviz.digraph(modelA)) %>%
export_svg %>% charToRaw %>% rsvg_png("outputs/modelA_signed-digraph.png")
