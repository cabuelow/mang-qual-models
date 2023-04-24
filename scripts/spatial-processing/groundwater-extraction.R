# Potential global subsidence from groundwater extraction in 2010 and by 2040
# extract values in each mangrove forest patch
# where forest patches don't intersect with subsidence raster, use buffer of 10km to calculate nearest values
# otherwise provide a value of 1 for 'very low probability  of subsidence

library(terra)
library(sf)
library(exactextractr)
library(tidyverse)
library(tmap)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg') %>% st_buffer(0)
gw_2040 <- rast('data/subsidence/GSH_2040/GSH_2040.tif')
gw <- rast('data/subsidence/GSH/GSH.tif')
gw3 <- list(gw, gw_2040)
names(gw3) <- c('gw_subsid_2010', 'gw_subsid_2040')
nam <- c('gw_subsid_2010', 'gw_subsid_2040')

for(i in seq_along(gw3)){
  
gw <- gw3[[i]]
# calculate summary stats (mean, median, mode, min, max) in each forest patch
system.time(gw_summary <- exact_extract(gw, typ, c('mean', 'median', 'mode', 'min', 'max')))
results <- data.frame(Type = typ$Type, gw_summary)

# where results are 'NA'. buffer forst patch by 10 km and calculate summary stats
typ_sub <- typ %>% filter(Type %in% filter(results, is.na(mean))$Type) %>% st_buffer(10000)
gw_summary_buff <- exact_extract(gw, typ_sub, c('mean', 'median', 'mode', 'min', 'max'))
results_buff <- data.frame(Type = typ_sub$Type, gw_summary_buff)

results_final <- results %>% 
  left_join(results_buff, by = 'Type') %>% 
  mutate(across(c(mean.x:max.x), 
                ~ ifelse(is.na(.), get(paste0(substr(cur_column(), 1, nchar(cur_column())-2), '.y')), .))) %>% 
  select(-c(mean.y:max.y)) %>% 
  mutate(mode.x = ifelse(is.na(mode.x), 1, mode.x))

colnames(results_final)[-1] <- c(substr(colnames(results_final[,-1]), 1, nchar(colnames(results_final[,-1]))-2))

write.csv(results_final, paste0('outputs/processed-data/', nam[i] ,'.csv'), row.names = F)
}

tmap_mode('view')

# do some checks the above works
sub <- crop(gw, vect(st_buffer(filter(typ_sub, Type == 'Estuary_30019'), 100000)))
qtm(sub) + qtm(filter(typ_sub, Type == 'Estuary_30019')) + qtm(filter(typ, Type == 'Lagoon_20036'))
filter(results, Type == 'Estuary_30019')

# plot all to check

typ_gw <- typ_cent %>% 
  left_join(results)

qtm(gw3[[1]]) + qtm(typ_gw, dots.col = 'gw_subsid_2010')
