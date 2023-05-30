library(sf)
sf_use_s2(FALSE)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite.shp')
eez <- st_read('data/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp') 

# spatial join whereby if multiple countries intersect a typology, the one with largest area gets assigned

typ_eez <- typ %>% 
  st_join(eez, largest = T)

typ_eez_df <- typ_eez %>% 
  st_drop_geometry() %>% 
  select(Type:UNION, TERRITORY1, ISO_TER1, SOVEREIGN1, ISO_SOV1)

write.csv(typ_eez_df, 'outputs/processed-data/typology-land-eez.csv', row.names = F)
