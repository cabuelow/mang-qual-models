library(sf)
library(tmap)
library(tidyverse)
library(rmapshaper)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg')
hydro <- read.csv('data/Hydro_Dat.csv')
slr <- read.csv('data/SLR_Data.csv')

# subset for Australia

# join attributes to typologies

typ2 <- typ %>% 
  left_join(hydro) %>% 
  left_join(slr)

# explore

tmap_mode('view')
qtm(typ2[1:20,], fill = 'ID')
