# Calculate mangrove clumpiness base on typologies units

library(tidyverse)
library(terra)
library(sf)
library(landscapemetrics)
library(tmap)

terraOptions(tempdir = "e:/KC/temp")

## read and merge global mangrove watch raster data
# fils <- list.files("data/GMW_v3/gmw_v3_2020_gtiff/gmw_v3_2020", full.names = TRUE)
# gmw_rast <- lapply(fils, rast)
# gmw_sprc <- sprc(gmw_rast)
# gmw <- merge(gmw_sprc) 

# Warning messages:                         
# 1: In doTryCatch(return(expr), name, parentenv, handler) :
#   restarting interrupted promise evaluation
# 2: [merge] Estimated disk space needed without compression: 488GB. Available: 67 GB. 

# merged and compressed, output tif retrieved from temp folder and saved in data folder

gmw <- rast("data/GMW_v3/gmw_v3_2020_gtiff/gmw_v3_2020_merged.tif")

# read typologies data as SpatVector for masking
typo <- vect("data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg")

## Method 1: using GMW 2020 raster layer mask over typologies polygons
results <- data.frame(Type = NA, clumpy = NA)

for (i in 1:nrow(typo)) {
  gmw_mask <- crop(gmw, typo[i,], mask = TRUE) # crop and mask by typology unit
  gmw_mask <- project(gmw_mask, "+proj=cea", method = "near")
  # replace NA to 0, to make classes of 0(non-mangrove) and 1(mangrove)
  gmw_mask <- subst(gmw_mask, NA, 0) 

  clumpy <- lsm_c_clumpy(gmw_mask)
  
  # There are 2 units fully covered by mangrove(only contain class 1), clumpiness 
  # returns NA. We decided to replace the NA to 0 for now to be conservative,
  # meaning it is neither aggregated nor disaggregated.
  if (nrow(clumpy) == 1 & 1 %in% clumpy$class) {
    clumpy$value <- 0
  } 
  
  results[i,] <- c(values(typo[i,1]), clumpy$value[clumpy$class == 1])
}

write.csv(results, 'outputs/processed-data/clumpiness.csv', row.names = F)

## Method 2: using just rasterized typologies vector layer
results2 <- data.frame(Type = NA, clumpy = NA)

for (i in 1:nrow(typo)) {
  typo1 <- typo[i,]
  
  # make raster layer that has the same resolution as the GMW raster 
  r <- rast(ext(typo1), res = 0.0002222222, crs = crs(typo1)) 
  r <- init(r, 1)
  
  typo1_r <- rasterize(typo1, r)
  typo1_r <- project(typo1_r, "+proj=cea", method = "near")
  typo1_r <- subst(typo1_r, NA, 0)
  # replace NA to 0, to make classes of 0(non-mangrove) and 1(mangrove)
  
  clumpy <- lsm_c_clumpy(typo1_r)
  
  # There are a few units fully covered by mangrove(only contain class 1), clumpiness 
  # returns NA. We decided to replace the NA to 0 for now to be conservative,
  # meaning it is neither aggregated nor disaggregated.
  if (nrow(clumpy) == 1 & 1 %in% clumpy$class) {
    clumpy$value <- 0
  }  
  
  results2[i,] <- c(values(typo[i,1]), clumpy$value[clumpy$class == 1])
}

write.csv(results2, 'outputs/processed-data/clumpiness2.csv', row.names = F)


## test: to see if landscapemetrics works, compare two projection cea vs moll
# gmw2 <- rast("data/GMW_v3/gmw_v3_2020_gtiff/gmw_v3_2020/GMW_N00E008_2020_v3.tif")
# 
# # cea
# gmw2_cea <- project(gmw2, "+proj=cea", method = "near")
# gmw2_cea <- subst(gmw2_cea, NA, 0)
# 
# check_landscape(gmw2_cea)
# lsm_c_clumpy(gmw2_cea)
#
# # cea by location
# library(crsuggest)
# crs <- suggest_crs(gmw2, gcs = 4326)[[grep("cea", suggest_crs(gmw2, gcs = 4326)[[6]]),6]]
# 
# gmw2_cea <- project(gmw2, crs, method = "near")
# gmw2_cea <- subst(gmw2_cea, NA, 0)
# 
# check_landscape(gmw2_cea)
# lsm_c_clumpy(gmw2_cea)
# 
# # mollwied
# gmw2_moll <- project(gmw2, "+proj=moll", method = "near")
# gmw2_moll <- subst(gmw2_moll, NA, 0)
# 
# check_landscape(gmw2_moll)
# lsm_c_clumpy(gmw2_moll)

