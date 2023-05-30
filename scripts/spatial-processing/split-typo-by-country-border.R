# split typology unit by country EEZ border

library(tidyverse)
library(sf)
library(mapview)

sf_use_s2(FALSE)

typo <- st_read("data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg") %>%  
  st_set_agr("constant")

eez <- st_read("data/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp")%>% 
  select(c(ISO_TER1, ISO_SOV1, geometry)) %>% 
  st_set_agr("constant")

datalist <- list()

for (i in 1:nrow(typo)) {
  typo_sub <- typo[i, ]
  
  j <- st_intersects(typo_sub, eez)[[1]] # eez features that intersect with typo_sub
  
  if (length(j) >= 2) {
    
    eez_sub <- eez[j,]
    
    typo_eez_sub <- st_intersection(typo_sub, eez_sub) %>%
      mutate(Type_TER1 = paste0(Type, "_", ISO_TER1))
    
    datalist[[i]] <- typo_eez_sub
    
  } else{
    datalist[[i]] <- typo_sub %>% 
      mutate(ISO_TER1 = NA,
             ISO_SOV1 = NA,
             Type_TER1 = NA)
  }
}

typo_eez <- do.call(rbind, datalist)

st_write(typo_eez, "data/typologies/Mangrove_Typology_v3_Composite_valid_EEZ.gpkg")


# typo_sub <- typo[typo$Type %in% "Estuary_6103", ]

# typo_eez <- st_intersection(typo, eez) %>%
#   mutate(Type_TER1 = paste0(Type, "_", ISO_TER1))
