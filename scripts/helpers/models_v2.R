# function for saving models as images

save_model_image <- function(x, directory, name){
  require(DiagrammeR)
  require(DiagrammeRsvg)
  require(rsvg)
grViz(grviz.digraph(x)) %>%
  export_svg %>% charToRaw %>% rsvg_png(paste0(directory, '/', name, '.png'))
}

# build models as signed digraphs
# mangrove model
model <- parse.digraph(c( 'LandwardMang->SubVol',
                          'SeawardMang<->SubVol',
                          'LandwardAvailableProp->LandwardMang',
                          'SeawardAvailableProp->SeawardMang',
                          'Erosion-*SubVol',
                          'Erosion-*SeawardMang',
                          'ExtremeRainfall-->Sediment',
                          'Sediment->SubVol',
                          'CoastalDev-*LandwardMang',
                          'CoastalDev--*SeawardMang',
                          'GroundSubsid->LandwardMang',
                          'GroundSubsid-*SeawardMang',
                          'SeaLevelRise->LandwardMang',
                          'SeaLevelRise-*SeawardMang',
                          'Cyclones--*LandwardMang',
                          'Cyclones-*SeawardMang',
                          'Cyclones-*SubVol',
                          'Cyclones-->LandwardAvailableProp',
                          'Drought--*Sediment'
)) %>% enforce.limitation() %>% save_model_image('outputs/model-images', 'model')

# alternative model structure where cyclones increase sediment volume
model_cyc_pos <- parse.digraph(c( 'LandwardMang->SubVol',
                          'SeawardMang<->SubVol',
                          'LandwardAvailableProp->LandwardMang',
                          'SeawardAvailableProp->SeawardMang',
                          'Erosion-*SubVol',
                          'Erosion-*SeawardMang',
                          'ExtremeRainfall-->Sediment',
                          'Sediment->SubVol',
                          'CoastalDev-*LandwardMang',
                          'CoastalDev--*SeawardMang',
                          'GroundSubsid->LandwardMang',
                          'GroundSubsid-*SeawardMang',
                          'SeaLevelRise->LandwardMang',
                          'SeaLevelRise-*SeawardMang',
                          'Cyclones--*LandwardMang',
                          'Cyclones-*SeawardMang',
                          'Cyclones->SubVol',
                          'Cyclones-->LandwardAvailableProp',
                          'Drought--*Sediment'
)) %>% enforce.limitation() %>% save_model_image('outputs/model-images', 'model_cyc_pos')

# alternative model structure where rainfall is not uncertain
model_rain <- parse.digraph(c(  'LandwardMang->SubVol',
                                'SeawardMang<->SubVol',
                                'LandwardAvailableProp->LandwardMang',
                                'SeawardAvailableProp->SeawardMang',
                                'Erosion-*SubVol',
                                'Erosion-*SeawardMang',
                                'ExtremeRainfall->Sediment',
                                'Sediment->SubVol',
                                'CoastalDev-*LandwardMang',
                                'CoastalDev--*SeawardMang',
                                'GroundSubsid->LandwardMang',
                                'GroundSubsid-*SeawardMang',
                                'SeaLevelRise->LandwardMang',
                                'SeaLevelRise-*SeawardMang',
                                'Cyclones--*LandwardMang',
                                'Cyclones-*SeawardMang',
                                'Cyclones-*SubVol',
                                'Cyclones-->LandwardAvailableProp',
                                'Drought--*Sediment'
)) %>% enforce.limitation() %>% save_model_image('outputs/model-images', 'model_rain')

# alternative model structure where drought is not uncertain
model_drought <- parse.digraph(c(  'LandwardMang->SubVol',
                                   'SeawardMang<->SubVol',
                                   'LandwardAvailableProp->LandwardMang',
                                   'SeawardAvailableProp->SeawardMang',
                                   'Erosion-*SubVol',
                                   'Erosion-*SeawardMang',
                                   'ExtremeRainfall-->Sediment',
                                   'Sediment->SubVol',
                                   'CoastalDev-*LandwardMang',
                                   'CoastalDev--*SeawardMang',
                                   'GroundSubsid->LandwardMang',
                                   'GroundSubsid-*SeawardMang',
                                   'SeaLevelRise->LandwardMang',
                                   'SeaLevelRise-*SeawardMang',
                                   'Cyclones--*LandwardMang',
                                   'Cyclones-*SeawardMang',
                                   'Cyclones-*SubVol',
                                   'Cyclones-->LandwardAvailableProp',
                                   'Drought-*Sediment'
)) %>% enforce.limitation() %>% save_model_image('outputs/model-images', 'model_drought')

# alternative model structure where cyclone seaward is uncertain
model_cyc_seaward <- parse.digraph(c(  'LandwardMang->SubVol',
                                   'SeawardMang<->SubVol',
                                   'LandwardAvailableProp->LandwardMang',
                                   'SeawardAvailableProp->SeawardMang',
                                   'Erosion-*SubVol',
                                   'Erosion-*SeawardMang',
                                   'ExtremeRainfall-->Sediment',
                                   'Sediment->SubVol',
                                   'CoastalDev-*LandwardMang',
                                   'CoastalDev--*SeawardMang',
                                   'GroundSubsid->LandwardMang',
                                   'GroundSubsid-*SeawardMang',
                                   'SeaLevelRise->LandwardMang',
                                   'SeaLevelRise-*SeawardMang',
                                   'Cyclones--*LandwardMang',
                                   'Cyclones--*SeawardMang',
                                   'Cyclones-*SubVol',
                                   'Cyclones-->LandwardAvailableProp',
                                   'Drought--*Sediment'
)) %>% enforce.limitation() %>% save_model_image('outputs/model-images', 'model_cyc_seaward')

