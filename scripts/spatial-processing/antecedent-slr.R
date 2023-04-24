# regional mean sea level trends
# truncate extreme values to +or- 5mm/yr

library(sf)
library(terra)
library(exactextractr)
library(tidyverse)
library(tmap)

slr <- rast('data/regional-mean-sea-level-trends/ESACCI-SEALEVEL-IND-MSLTR-MERGED-20161202000000-fv02.nc')
typ_cent <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg') %>% st_transform(crs(slr[[1]]))

# identify slr raster pixels that are closest to centroid of each mangrove forest patch
# loop through the different slr rasters to get estimates for each different time periods

results <- data.frame(Type = typ_cent$Type)

slr.p <- st_as_sf(as.polygons(slr[[1]], dissolve = F, values = T, na.rm = T))
    system.time(dist <- data.frame(st_distance(typ_cent, slr.p))) # s2 is on so calculates great circle distances on sphere (instead of ellipsoid)
    # above takes 16 minutes
    
    # filter values from raster that corresponds to column name and add to results dataframe
    results[,2] <- na.omit(values(slr[[1]]))[as.numeric(substr(colnames(dist)[apply(dist,1,which.min)], 2, nchar(colnames(dist)[apply(dist,1,which.min)])))]
    colnames(results) <- c('Type', 'local_msl_trend')

# truncate to plus or minus 5mm/yr

    results2 <- results %>% 
      mutate(local_msl_trend = ifelse(local_msl_trend >  5, 5, local_msl_trend)) %>% 
      mutate(local_msl_trend = ifelse(local_msl_trend <  -5, -5, local_msl_trend))

write.csv(results2, 'outputs/processed-data/antecedent-slr.csv', row.names = F)
