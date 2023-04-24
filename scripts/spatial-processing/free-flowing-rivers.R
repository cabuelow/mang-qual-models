# in arcgis, intersect ffr dataset with mangrove catchments (made in 'make-mangrove-catchments.R')
# read in here, find most downstream, coastal outlet for each typology catchment
# to get sediment trapping index for that typology
# do weighted average of sed. trapping index if there are multiple, weighted by the FFR's average longterm river discharge (DIS_AV_CMS)
# with discharge set to the minimum value for segments with no flow

library(sf)
library(tidyverse)
library(tmap)
sf_use_s2(TRUE)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
ffr <- st_read('data/free-flowing-rivers/mangrove_catchment_basins_final_dissolve_ffr.shp')

# identify coastal segment of each ffr intersecting with each forest catchment

ffr_outlet <- ffr %>% filter(NDOID == 0)

# check these are all coastal and not going to a terminal inland sink

ffr_outlet.pt <- st_centroid(ffr_outlet)
tmap_mode('view')
qtm(ffr_outlet.pt) # looks good!

ffr_outlet.df <- st_drop_geometry(ffr_outlet)
write.csv(ffr_outlet.df, 'data/free-flowing-rivers/mangrove-typology-ffr-outlet.csv', row.names = F)

# now calculate weighted average for each forest patch

ffr_sum <- ffr_outlet.df %>% 
  group_by(Type) %>% 
  summarise(SED_weighted_average = sum(SED*DIS_AV_CMS)/sum(DIS_AV_CMS),
            DOF_weighted_average = sum(DOF*DIS_AV_CMS)/sum(DIS_AV_CMS),
            DOR_weighted_average = sum(DOR*DIS_AV_CMS)/sum(DIS_AV_CMS),
            USE_weighted_average = sum(USE*DIS_AV_CMS)/sum(DIS_AV_CMS),
            URB_weighted_average = sum(URB*DIS_AV_CMS)/sum(DIS_AV_CMS),
            RDD_weighted_average = sum(RDD*DIS_AV_CMS)/sum(DIS_AV_CMS),
            FLD_weighted_average = sum(FLD*DIS_AV_CMS)/sum(DIS_AV_CMS),
            CSI_weighted_average = sum(CSI*DIS_AV_CMS)/sum(DIS_AV_CMS))

# check missing basins

miss <- filter(typ, !Type %in% ffr_sum$Type) %>% st_buffer(20000) %>% st_write('miss.shp') # all looks okay

# if don't have a ffr, set to 0

ffr_final <- st_drop_geometry(typ) %>% 
  select(Type) %>% 
  left_join(ffr_sum) %>% 
  mutate_at(vars(SED_weighted_average:CSI_weighted_average), ~ifelse(is.na(.), 0, .))

#save

write.csv(ffr_final, 'outputs/processed-data/free-flowing-rivers.csv', row.names = F)

