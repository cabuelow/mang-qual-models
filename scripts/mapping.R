library(sf)
library(tmap)
library(tidyverse)

#typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite.shp') %>% 
 # st_make_valid() %>% st_write('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg')
typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg')
hydro <- read.csv('data/Hydro_Dat.csv')
slr <- read.csv('data/SLR_Data.csv')

# join attributes to typologies

typ2 <- typ %>% 
  left_join(hydro) %>% 
  left_join(slr)

# explore

tmap_mode('view')

qtm(typ2)
