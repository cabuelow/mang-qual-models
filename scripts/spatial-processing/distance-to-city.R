# Extract raster data for typology of mangroves data

library(tidyverse)
library(terra)
library(exactextractr)
library(sf)
library(tmap)

dist <- rast("data/Global_Travel_Time_to_Cities_2015/201501_Global_Travel_Time_to_Cities_2015.tif")

typo <- st_read("data/typologies/Mangrove_Typology_v3_Composite_valid_EEZ_qgis.gpkg")

value <- exact_extract(dist, typo, c('mean', 'min', 'max', 'stdev'))
results <- data.frame(Type_MRGID = typo$Type_MRGID, value)

## Method 1: replace na with max distance in the entire dist raster
# results1 <- results %>% 
#   mutate(mean = ifelse(is.na(mean), max(max, na.rm = TRUE), mean))
# 
# write.csv(results1, 'outputs/processed-data/distance-to-city-na-max-all.csv', row.names = F)


## Method 2: replace na with max distance of the nearest polygon

#subset typology sf
typo_non_na <- typo[!is.na(results$mean),]
typo_na <- typo[is.na(results$mean),]

# st_write(typo_na, "outputs/typo_na.gpkg")

# find max value from nearest features 
nearest <- data.frame(Type_MRGID = typo_na$Type_MRGID, 
                      Type_nearest = typo_non_na$Type_MRGID[st_nearest_feature(typo_na, typo_non_na)])

typo_na_results <- left_join(nearest, results, by = join_by(`Type_nearest` == `Type_MRGID`)) %>% 
  select(Type_MRGID, max)
colnames(typo_na_results)[2] <- "max.y"

results2 <- left_join(results, typo_na_results, by = join_by(Type_MRGID)) %>% 
  mutate(mean = ifelse(is.na(mean), max.y, mean))
results2 <- select(results2, -c(max.y))

write.csv(results2, 'outputs/processed-data/distance-to-city.csv', row.names = F)



# NA polygons distance to nearest pixel
# try with smaller area
# typo_na <- read_sf("outputs/typo_na.gpkg")
# 
# dist_sub <- crop(dist, c(107.46260551557143,108.93916801557143,-3.3845122872966567,-2.3855289825021364))
# dist_sub_v <- as.polygons(dist_sub, values = FALSE) # too heavy to run
# 
# grid <- rast(ext(dist_sub), nrow = nrow(dist_sub), ncol = ncol(dist_sub))
# 
# grid_dist <- distance(grid, dist_sub_v)
# 
# result <- exact_extract(grid_dist, typo_na, 'mean')
# result2 <- data.frame(Type = typo_na$Type, result)
# typo_na_dist <- result2[!is.na(result),]
# 
# typo_result <- right_join(typo_na, typo_na_dist, by = "Type")
# 
# plot(typo_result[5])

     