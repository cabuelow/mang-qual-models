# calculate  projected future SLR (metres increase from baseline scenario) for each forest patch
# under SSP5-8.5 pathway

library(sf)
library(terra)
library(exactextractr)
library(tidyverse)
library(tmap)

typ_cent <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
fils <- list.files('data/future-slr/', pattern = '.tif', full.name = T)
slr <- lapply(fils, rast)
names(slr) <- c('slr_m_2081_2100', 'slr_m_2041_2060', 'slr_m_2021_2040')

# identify slr raster pixels that are closest to centroid of each mangrove forest patch
# loop through the different slr rasters to get estimates for each different time periods

results <- data.frame(Type = typ_cent$Type)

system.time( # takes 74 mins to run
for(i in seq_along(slr)){
  
slr.p <- st_as_sf(as.polygons(slr[[i]], dissolve = F, values = T, na.rm = T))
system.time(dist <- data.frame(st_distance(typ_cent, slr.p))) # s2 is on so calculates great circle distances on sphere (instead of ellipsoid)
# above takes 16 minutes

# filter values from raster that corresponds to column name and add to results dataframe
results[,i+1] <- na.omit(values(slr[[i]]))[as.numeric(substr(colnames(dist)[apply(dist,1,which.min)], 2, nchar(colnames(dist)[apply(dist,1,which.min)])))]
colnames(results)[i+1] <- paste0(names(slr)[i])

}) # end loop

write.csv(results, 'outputs/processed-data/future-slr.csv', row.names = F)

tmap_mode('view')

# do some checks the above works
sub <- crop(slr[[1]], vect(st_buffer(filter(typ_cent, Type == 'Estuary_7177'), dist = 200000)))
qtm(slr[[1]]) + qtm(sub) + qtm(filter(typ_cent, Type == 'Estuary_7177'))
filter(results, Type == 'Estuary_7177')

# plot all to check

typ_slr <- typ_cent %>% 
  left_join(results)

qtm(slr[[1]]) + qtm(typ_slr, dots.col = 'slr_m_2081_2100')
