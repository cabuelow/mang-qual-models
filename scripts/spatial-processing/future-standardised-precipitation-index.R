# calculate projected future spi (standardised precipitation index) for each mangrove forest patch
# SSP5-8.5

library(sf)
library(terra)
library(exactextractr)
library(tidyverse)
library(tmap)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg') %>% st_buffer(0)
fils <- list.files('data/future-spi/', pattern = '.tif', full.name = T)
spi <- lapply(fils, rast)
spi <- c(spi[[1]], spi[[2]], spi[[3]])
names(spi) <- c('spi_change_percent_2081_2100', 'spi_change_percent_2041_2060', 'spi_change_percent_2021_2040')

results <- exact_extract(unlist(spi), typ, 'mean')
results <- data.frame(Type = typ$Type, results)

write.csv(results, 'outputs/processed-data/future-stand-precip-index.csv', row.names = F)


