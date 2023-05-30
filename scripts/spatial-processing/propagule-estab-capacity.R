# find the average nearest neighbour distance between extant and lost forest within each typological unit
# if there is no loss patch it will get a 'no loss' character string

library(sf)
library(tidyverse)

extant <- st_read('data/typologies/mangrove_typology_2020_dissolve.shp')
loss <- st_read('data/typologies/gmw_v3_f96_t19_typ.shp') 
loss2 <- loss %>% filter(chng_type == 'loss')
typ_cent <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')

df <- data.frame(Type = NA, extant_loss_dist_mean_m = NA, extant_loss_dist_min_m = NA, extant_loss_dist_max_m = NA)

system.time(
for(i in 1:nrow(extant)){
  extant_sub <- extant[i,]
  loss_sub <- filter(loss2, Type == extant_sub$Type)
  if(nrow(loss_sub)>0){
  nn <- st_nearest_points(st_make_valid(extant_sub), st_make_valid(loss_sub))
  nn2 <- st_distance(nn)
  df[i,1] <- extant_sub$Type
  df[i,2] <- as.numeric(mean(nn2))
  df[i,3] <- as.numeric(min(nn2))
  df[i,4] <- as.numeric(max(nn2))
  }else{
    df[i,1] <- extant_sub$Type
    df[i,2] <- 'no_loss'
    df[i,3] <- 'no_loss'
    df[i,4] <- 'no_loss'
  }
  write.csv(df, 'outputs/processed-data/extant-loss-distances.csv', row.names = F)
  print(i)
})

datfinal <- typ_cent %>% 
  st_drop_geometry() %>% 
  select(Type) %>% 
  left_join(df) %>% 
  mutate(extant_loss_dist_mean_m = ifelse(is.na(extant_loss_dist_mean_m), 'no_extant', extant_loss_dist_mean_m),
         extant_loss_dist_min_m = ifelse(is.na(extant_loss_dist_min_m), 'no_extant', extant_loss_dist_min_m),
         extant_loss_dist_max_m = ifelse(is.na(extant_loss_dist_max_m), 'no_extant', extant_loss_dist_max_m))


write.csv(datfinal, 'outputs/processed-data/propagule-establishment-distances.csv', row.names = F)

