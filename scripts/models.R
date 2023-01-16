library(QPress)
library(dagitty)
library(ggdag)
library(igraph)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# build model signed digraphs

# mangrove system model
modelA <- parse.digraph(c('LandwardEstabSpace->LandwardPropag',
                'LandwardPropag->LandwardMang',
                'LandwardMang->OM',
                'SeawardEstabSpace->SeawardPropag',
                'SeawardPropag->SeawardMang',
                'SeawardMang->OM',
                'OM->SubVol',
                'SubVol->SeawardEstabSpace')) %>% 
  enforce.limitation()
grViz(grviz.digraph(modelA)) #%>%
  #export_svg %>% charToRaw %>% rsvg_png("outputs/graph.png")

# mangrove system + local factors                
modelB <- parse.digraph(c('LandwardEstabSpace->LandwardPropag',
                          'LandwardPropag->LandwardMang',
                          'LandwardMang->OM',
                          'SeawardEstabSpace->SeawardPropag',
                          'SeawardPropag->SeawardMang',
                          'SeawardMang->OM',
                          'OM->SubVol',
                          'SubVol->SeawardEstabSpace', 
                          'EcologicalConn->LandwardPropag',
                          'EcologicalConn->SeawardPropag',
                          'TidalRange->LandwardEstabSpace',
                          'TidalRange->SeawardEstabSpace',
                          'HydroEnergy-*SeawardEstabSpace',
                          'HydroEnergy->Erosion',
                          'Erosion-*SubVol',
                          'Erosion-*SeawardPropag',
                          'Erosion-*SeawardMang',
                          'Precipitation->Sediment',
                          'Sediment->SubVol',
                          'Autocompaction-*SubVol',
                          'TidalFreq-*SeawardPropag',
                          'TidalFreq->SubVol')) %>% 
  enforce.limitation()

grViz(grviz.digraph(modelB)) #%>%
#export_svg %>% charToRaw %>% rsvg_png("outputs/graph.png")

# mangrove system + local factors + anthropogenic pressures
modelC <- parse.digraph(c('LandwardEstabSpace->LandwardPropag',
                          'LandwardPropag->LandwardMang',
                          'LandwardMang->OM',
                          'SeawardEstabSpace->SeawardPropag',
                          'SeawardPropag->SeawardMang',
                          'SeawardMang->OM',
                          'OM->SubVol',
                          'SubVol->SeawardEstabSpace', 
                          'EcologicalConn->LandwardPropag',
                          'EcologicalConn->SeawardPropag',
                          'TidalRange->LandwardEstabSpace',
                          'TidalRange->SeawardEstabSpace',
                          'HydroEnergy-*SeawardEstabSpace',
                          'HydroEnergy->Erosion',
                          'Erosion-*SubVol',
                          'Erosion-*SeawardPropag',
                          'Erosion-*SeawardMang',
                          'Precipitation->Sediment',
                          'Dams-*Sediment',
                          'Sediment->SubVol',
                          'Autocompaction-*SubVol',
                          'TidalFreq-*SeawardPropag',
                          'TidalFreq->SubVol',
                          'Dams->Sediment',
                          'CoastalDev-*LandwardEstabSpace',
                          'CoastalDev-*LandwardPropag',
                          'CoastalDev-*LandwardMang',
                          'CoastalDev--*SeawardMang',
                          'CoastalDev--*SeawardPropag',
                          'GroundSubsid->LandwardEstabSpace',
                          'GroundSubsid-*SeawardEstabSpace',
                          'GroundSubsid-*SeawardPropag',
                          'GroundSubsid-*SeawardMang'
                          )) %>% 
  enforce.limitation()

grViz(grviz.digraph(modelC)) #%>%
#export_svg %>% charToRaw %>% rsvg_png("outputs/graph.png")

# mangrove system + local factors + anthropogenic pressures + climate pressures
# alternate model stucture 1: cyclones have neg. effect on landward prop, and neg. eff on substrate vol
modelD.1 <- parse.digraph(c('LandwardEstabSpace->LandwardPropag',
                          'LandwardPropag->LandwardMang',
                          'LandwardMang->OM',
                          'SeawardEstabSpace->SeawardPropag',
                          'SeawardPropag->SeawardMang',
                          'SeawardMang->OM',
                          'OM->SubVol',
                          'SubVol->SeawardEstabSpace', 
                          'EcologicalConn->LandwardPropag',
                          'EcologicalConn->SeawardPropag',
                          'TidalRange->LandwardEstabSpace',
                          'TidalRange->SeawardEstabSpace',
                          'HydroEnergy-*SeawardEstabSpace',
                          'HydroEnergy->Erosion',
                          'Erosion-*SubVol',
                          'Erosion-*SeawardPropag',
                          'Erosion-*SeawardMang',
                          'Precipitation->Sediment',
                          'Dams-*Sediment',
                          'Sediment->SubVol',
                          'Autocompaction-*SubVol',
                          'TidalFreq-*SeawardPropag',
                          'TidalFreq->SubVol',
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
                          'Cyclones-*SubVol'
)) %>% 
  enforce.limitation()

grViz(grviz.digraph(modelD.1)) #%>%
#export_svg %>% charToRaw %>% rsvg_png("outputs/graph.png")

# mangrove system + local factors + anthropogenic pressures + climate pressures
# alternate model srtucture 2: cyclones have pos. effect on landward prop, and pos. eff on substrate vol
modelD.2 <- parse.digraph(c('LandwardEstabSpace->LandwardPropag',
                            'LandwardPropag->LandwardMang',
                            'LandwardMang->OM',
                            'SeawardEstabSpace->SeawardPropag',
                            'SeawardPropag->SeawardMang',
                            'SeawardMang->OM',
                            'OM->SubVol',
                            'SubVol->SeawardEstabSpace', 
                            'EcologicalConn->LandwardPropag',
                            'EcologicalConn->SeawardPropag',
                            'TidalRange->LandwardEstabSpace',
                            'TidalRange->SeawardEstabSpace',
                            'HydroEnergy-*SeawardEstabSpace',
                            'HydroEnergy->Erosion',
                            'Erosion-*SubVol',
                            'Erosion-*SeawardPropag',
                            'Erosion-*SeawardMang',
                            'Precipitation->Sediment',
                            'Dams-*Sediment',
                            'Sediment->SubVol',
                            'Autocompaction-*SubVol',
                            'TidalFreq-*SeawardPropag',
                            'TidalFreq->SubVol',
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
                            'Cyclones-->LandwardEstabPropag',
                            'Cyclones--*LandwardMang',
                            'Cyclones-*SeawardMang',
                            'Cyclones-*SeawardEstabPropag',
                            'Cyclones->SubVol'
)) %>% 
  enforce.limitation()

grViz(grviz.digraph(modelD.2)) #%>%
#export_svg %>% charToRaw %>% rsvg_png("outputs/graph.png")