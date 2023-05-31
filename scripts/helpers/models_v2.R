library(QPress)

# function for saving models as images

save_model_image <- function(x, directory, name){
  require(DiagrammeR)
  require(DiagrammeRsvg)
  require(rsvg)
grViz(grviz.digraph(x)) %>%
  export_svg %>% charToRaw %>% rsvg_png(paste0(directory, '/', name, '.png'))
}

models <- vector(mode = "list", length = 5)
names(models) <- c('mangrove_model', 'model_cyc_pos', 'model_rain', 'model_drought', 'model_cyc_seaward')

# build models as signed digraphs
# mangrove model
models[[1]] <- parse.digraph(c( 'LandwardMang->SubVol',
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
)) %>% enforce.limitation() 
models[[1]] %>% save_model_image('outputs/model-images', 'model')

# alternative model structure where cyclones increase sediment volume
models[[2]] <- parse.digraph(c( 'LandwardMang->SubVol',
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
)) %>% enforce.limitation() 
models[[2]] %>% save_model_image('outputs/model-images', 'model_cyc_pos')

# alternative model structure where rainfall is not uncertain
models[[3]] <- parse.digraph(c(  'LandwardMang->SubVol',
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
)) %>% enforce.limitation() 
models[[3]] %>% save_model_image('outputs/model-images', 'model_rain')

# alternative model structure where drought is not uncertain
models[[4]] <- parse.digraph(c(  'LandwardMang->SubVol',
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
)) %>% enforce.limitation() 
models[[4]] %>% save_model_image('outputs/model-images', 'model_drought')

# alternative model structure where cyclone seaward is uncertain
models[[5]]<- parse.digraph(c(  'LandwardMang->SubVol',
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
)) %>% enforce.limitation() 
models[[5]] %>% save_model_image('outputs/model-images', 'model_cyc_seaward')

