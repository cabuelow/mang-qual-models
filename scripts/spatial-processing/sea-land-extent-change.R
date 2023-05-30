# seaward and landward gains in mangrove extent from 1996 -2020 were obtained with following steps in ArcGIS:
# 1. intersect and erase tools used to isolate mangrove losses and gains in each time interval between 96-2020
# 2. all time intervals of losses or gains merged into one gross loss and gross gain layer 96-2020
# 3. intersect with eez polygons = seaward gains/losses
# 4. erase eez polygons = landward gains/losses

# here we calculate area of seaward/landward loss/gain in each mangrove typology

library(sf)
library(tidyverse)
# s2 is on so assuming a sherical approximation of earth when doing area calcs

typ_cent <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')

sea_gain <- st_read('data/typologies/extent-change-gmwv3/seaward_gross_gains_96_20_eez_typ.shp') %>% 
  select(Type, geometry) %>% 
  st_make_valid() %>% 
  mutate(sea_gain_ha = as.numeric(st_area(.)/10000)) %>% 
  st_drop_geometry()

sea_loss <- st_read('data/typologies/extent-change-gmwv3/seaward_gross_loss_96_20_eez_typ.shp') %>% 
  select(Type, geometry) %>% 
  st_make_valid() %>% 
  mutate(sea_loss_ha = as.numeric(st_area(.)/10000)) %>% 
  st_drop_geometry()

land_gain <- st_read('data/typologies/extent-change-gmwv3/landward_gross_gains_96_20_eez_typ.shp') %>% 
  select(Type, geometry) %>% 
  st_make_valid() %>% 
  mutate(land_gain_ha = as.numeric(st_area(.)/10000)) %>% 
  st_drop_geometry()

land_loss <- st_read('data/typologies/extent-change-gmwv3/landward_gross_loss_96_20_eez_typ.shp') %>% 
  select(Type, geometry) %>% 
  st_make_valid() %>% 
  mutate(land_loss_ha = as.numeric(st_area(.)/10000)) %>% 
  st_drop_geometry()

# bind all together in dataframe and save

df <- select(st_drop_geometry(typ_cent), Type) %>% 
  left_join(sea_gain) %>% 
  left_join(sea_loss) %>% 
  left_join(land_gain) %>% 
  left_join(land_loss) %>% 
  mutate_at(vars(sea_gain_ha:land_loss_ha), as.numeric)

write.csv(df, 'outputs/processed-data/sea-land-extent-change.csv', row.names = F)
