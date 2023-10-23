# process median and average aridity index values to mangrove typology units

library(terra)
library(sf)
library(exactextractr)
library(tidyverse)
library(tmap)

typ <- st_read('data/typologies/v3.14/Mangrove_Typology_v3.14_Composite_valid.gpkg') %>% st_buffer(0)
typ_points <- st_read('data/typologies/v3.14/Mangrove_Typology_v3.14_Composite_valid_centroids.gpkg')
ai <- rast('data/Global-AridityIndex_ET0_v3_annual/ai_v3_yr.tif')*0.0001 # aridity index annual average (1970-2000), to obtain correct index values, need to multiply by 0.0001

# reclassify raster so that 0 is NA

ai <- classify(ai, cbind(0,NA))

# calculate summary stats (mean, median, mode, min, max) in each forest patch
system.time(summary <- exact_extract(ai, typ, c('mean', 'median', 'mode', 'min', 'max')))
results <- data.frame(Type = typ$Type, summary)

# where results are 'NA'. buffer forest patch by 50 km and calculate summary stats
typ_sub <- typ %>% filter(Type %in% filter(results, is.na(mean))$Type) %>% st_buffer(50000)
summary_buff <- exact_extract(ai, typ_sub, c('mean', 'median', 'mode', 'min', 'max'))
results_buff <- data.frame(Type = typ_sub$Type, summary_buff)

# one unit still with missing value, provide nearest value
typ_sub2 <- typ %>% filter(Type %in% filter(results_buff, is.na(mean))$Type) %>% st_buffer(200000)
summary_buff2 <- exact_extract(ai, typ_sub2, c('mean', 'median', 'mode', 'min', 'max'))
results_buff2 <- data.frame(Type = typ_sub2$Type, summary_buff2)
results_buff2

# combine into single dataframe
results_final <- results %>% 
  filter(!Type %in% results_buff$Type) %>% 
  rbind(filter(results_buff, !Type %in% results_buff2$Type)) %>%
  rbind(results_buff2)

# check
typ_points <- typ_points %>% left_join(results_final, by = 'Type')

# save
write.csv(results_final, paste0('outputs/processed-data/aridity.csv'), row.names = F)
