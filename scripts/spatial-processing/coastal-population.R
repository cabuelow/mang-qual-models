# calculate current population density (2020) and future pop density (SSP5, 2100) as proxy for 'coastal squeeze'
# in the lower elevation coastal zone (lecz) of typology catchments. Or where don't have catchments, use a 50 km buffer
# rasters represent population count per cell - so will need to divide by lecz land area to get density

library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(tmap)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg')
catch_lecz <- st_read('data/BasinATLAS_Data_v10_shp/mangrove_catchment_lecz.shp')
catch <- st_read('data/BasinATLAS_Data_v10_shp/Mangrove_catchment_basins_final_dissolve.shp') # note there are 8 catchments that don't have a lecz
lecz <- st_read('data/lecz-v3/lecz-mask.gpkg')
pop_2020 <- rast('data/gridded-pop-world/gpw_v4_population_count_rev11_2020_30_sec.tif')
pop_2040 <- rast('data/pop_SSP/DIVA_SSP5_2040.tif')
pop_2060 <- rast('data/pop_SSP/DIVA_SSP5_2060.tif')
pop_2100 <- rast('data/pop_SSP/DIVA_SSP5_2100.tif')
pop <- c(pop_2020, pop_2040, pop_2060, pop_2100)
names(pop) <- c('pop_size_lecz_2020', 'pop_size_lecz_2040', 'pop_size_lecz_2060', 'pop_size_lecz_2100')

# buffer mangrove patches without a catchment by 50 km and intersect with lower elevation coastal zone
#typ_miss <- typ %>% filter(!Type %in% catch_lecz$Type) %>% st_buffer(50000) %>% st_wrap_dateline() #st_intersection(lecz) 
#st_write(typ_miss, 'miss2.shp') # do intersection with lecz in arcgis
typ_buff <- st_read('data/BasinATLAS_Data_v10_shp/typ-buff-lecz.shp') # missing 4 where lecz doesn't intersect with buffer
# add to the catchment lecz
catch_lecz <- rbind(select(catch_lecz, Type, geometry), select(typ_buff, Type, geometry)) %>% st_wrap_dateline() %>% st_make_valid() %>%  mutate(area_ha = as.numeric(st_area(.))/10000)

# calculate total population size in each typological units lecz, then divide by total area of the lecz

results <- exact_extract(pop, catch_lecz, 'sum')
results2 <- data.frame(Type = catch_lecz$Type, results)

# join

results_final <- select(st_drop_geometry(typ), Type) %>% 
  left_join(select(st_drop_geometry(catch_lecz), Type, area_ha)) %>% 
  left_join(results2)
results_final[is.na(results_final)] <- 0

write.csv(results_final, 'outputs/processed-data/coastal-population.csv', row.names = F)

# below code makes the lecz mask used above
#lecz <- rast('data/lecz-v3/lecz_v3_spatial_data/data/merit_leczs.tif')
# reclassify lecz so that any pixel not classified as within 0 to 10m of the coast is NA,
# and any pixel within 0 to 10 m is 1
#m <- rbind(c(5,1), c(10,1))
#lecz.rc <- classify(lecz, m, others=NA)
#lecz.poly <- as.polygons(lecz.rc, dissolve = T)
#writeVector(lecz.poly, 'data/lecz-v3/lecz-mask.gpkg')

