# hydrobasins and rivers that intersect with mangrove typological units are identified and dissolved
# to create one catchment for each typology
# intersections and dissolves are completed in arcGIS for speed where necessary

library(sf)
library(tidyverse)
library(tmap)
sf_use_s2(TRUE)

# get river networks that intersect with mangrove forest patches

mang_riv <- st_read('data/HydroRivers_v10_shp/Mangrove_Typolgy_HRiver_Intersect.shp')
st_read('data/HydroRivers_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp') %>% 
filter(MAIN_RIV %in% mang_riv$MAIN_RIV) %>% 
left_join(select(st_drop_geometry(mang_riv), Type, MAIN_RIV), by = 'MAIN_RIV') %>% 
st_make_valid() %>% 
st_write('data/HydroRivers_v10_shp/mangrove_typ_rivers.gpkg', append = F, overwrite = T)

mang_riv2 <- st_read('data/HydroRivers_v10_shp/mangrove_typ_rivers.gpkg')
mang_riv2_diss <- st_read('data/HydroRivers_v10_shp/mangrove_typ_rivers_dissolve.shp') # dissolved in arcgis

# filter for basins that belong to each river network, export, then dissolve by Type in arcgis (faster)
# include basins that directly intersect with typologies (and no river network) and export

mang_bas <- st_read('data/BasinATLAS_Data_v10_shp/Mangrove_Typology_HBasin_L12.shp')
bas <- st_read('data/BasinATLAS_Data_v10_shp/BasinATLAS_v10_shp/BasinATLAS_v10_lev12.shp') 

bas_mang <- bas %>% 
filter(HYBAS_ID %in% c(mang_riv2$HYBAS_L12, mang_bas$HYBAS_ID))
st_write(bas_mang, 'data/BasinATLAS_Data_v10_shp/Mangrove_catchment_basins.gpkg')
bas_mang_df <- st_drop_geometry(bas_mang)
write.csv(bas_mang_df, 'data/BasinATLAS_Data_v10_shp/Mangrove_catchment_basins.csv', row.names = F)

mang_catchments <- bas_mang_df %>%   # river basins
left_join(select(st_drop_geometry(mang_riv2), Type, HYBAS_L12), by = c('HYBAS_ID' = 'HYBAS_L12')) 

mang_catchments2 <- bas_mang_df %>% # non-river basins that intersect mangroves
left_join(select(st_drop_geometry(mang_bas), Type, HYBAS_ID), by = 'HYBAS_ID')

# now bind and get distinct basins, then dissolve in arcgis

mang_catchments3 <- rbind(mang_catchments, mang_catchments2) %>% distinct() 
length(unique(mang_catchments3$Type)) # should be 3918
mang_catchments4 <- mang_catchments3 %>% filter(!is.na(Type))
length(unique(mang_catchments4$Type)) 
write.csv(mang_catchments4, 'data/BasinATLAS_Data_v10_shp/Mangrove_catchment_basins_final.csv', row.names = F)

mang_catch_Final <- bas %>% # select the unique basins belonging to mangroves and their river networks
  inner_join(select(mang_catchments4, HYBAS_ID, Type), by = 'HYBAS_ID') %>% st_make_valid() %>% st_wrap_dateline()
st_write(mang_catch_Final, 'data/BasinATLAS_Data_v10_shp/Mangrove_catchment_basins_final.shp', overwrite = T, append = F)

# dissolve in argis, read in to check

dat <- st_read('data/BasinATLAS_Data_V10_shp/Mangrove_catchment_basins_final_dissolve.shp') 
