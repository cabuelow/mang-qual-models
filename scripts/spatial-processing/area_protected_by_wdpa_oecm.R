library(sf)
library(tidyverse)
library(mapview)

sf_use_s2(FALSE)

pa <- st_read("data/WDPA/wdpa_woecm_jun2023_clean.gpkg") 
st_agr(pa) <- "constant"

typo <- st_read("data/typologies/v3.14/Mangrove_Typology_v3.14_Composite_valid_EEZ_qgis.gpkg") %>% 
  st_transform(crs = st_crs(pa)) 
st_agr(typo) <- "constant"

# calculate protected area for each typology unit
    
df <- data.frame(Type_MRGID = NA, wdpa_ha = NA, woecm_ha = NA, pa_oecm_ha =  NA)
    
for(i in 1:nrow(typo)){
    
  typo_sub <- typo[i,]
  typo_wdpa <- st_intersection(typo_sub, pa[pa$PA_DEF == "PA",]) %>% 
    st_union() # merge in case intersects with multiple pa
  typo_oecm <- st_intersection(typo_sub, pa[pa$PA_DEF == "OECM",]) %>% 
    st_union() # merge in case intersects with multiple oecm 
  pa_oecm <- st_union(typo_oecm, typo_wdpa) # merge pa and oecm to get total area
      
  wdpa_area <- st_area(typo_wdpa)/10000 # pa area in typo
  oecm_area <- st_area(typo_oecm)/10000 # oecm area in typo
  pa_oecm_area <- st_area(pa_oecm)/10000 # pa + oecm area in typo
      
  df[i,1] <- typo_sub$Type_MRGID
  df[i,2] <- ifelse(length(wdpa_area) == 0, NA, wdpa_area) # replace zeros with na
  df[i,3] <- ifelse(length(oecm_area) == 0, NA, oecm_area)
  df[i,4] <- ifelse(length(pa_oecm_area) == 0, NA, pa_oecm_area)
        
  write.csv(df, 'outputs/processed-data/wdpa_oecm_area.csv', row.names = F)
  }
    

