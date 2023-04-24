# Future Hydropower Reservoirs and Dams (FHReD) upstream of mangrove forest patches
# For each mangrove patch, just want the rivers that drain into each forest patch, i.e., are upstream
# then buffer those, and determine if a dam is within a certain radius
# Any mangroves without a river get 0

library(sf)
library(tidyverse)
library(tmap)

riv <- st_read('data/HydroRIVERS_v10_shp/mangrove_typ_rivers_dissolve.shp') %>% st_buffer(10000) # 10km buffer
typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')

fut_dams <- read.csv('data/FHReD_2015_future_dams_Zarfl_et_al_beta_version/future_dams.csv') %>% 
 st_as_sf(coords = c('Lon_Cleaned', 'LAT_cleaned'), crs = 4326)

riv_fut_dams <- riv %>% 
  st_join(fut_dams) %>% 
  st_drop_geometry() 

fut_dams.df <- riv_fut_dams %>% 
  mutate(dam = ifelse(is.na(DAM_ID), 0, 1)) %>% 
  group_by(Type) %>% 
  summarise(number_future_dams = sum(dam))

typ2 <- typ %>% 
  left_join(fut_dams.df) %>% 
  mutate(number_future_dams = ifelse(is.na(number_future_dams), 0, number_future_dams)) # if NA set to 0

# check
tmap_mode('view')
qtm(filter(riv, Type == 'Delta_70010')) +  qtm(fut_dams) + qtm(typ2, dots.col = 'number_future_dams')

write.csv(st_drop_geometry(typ2), 'outputs/processed-data/future-dams.csv', row.names = F)

