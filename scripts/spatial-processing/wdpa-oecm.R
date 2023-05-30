## typology categarized by protected or unprotected area

library(tidyverse)
library(sf)
library(tmap)

typo <- st_read("data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg")

pa <- st_read("data/WDPA_May2023_Public_shp/WDPA_May2023_Public_shp_0/WDPA_May2023_Public_shp-polygons.shp")
