# split typology unit by country EEZ border
# edited KC 20230927

library(tidyverse)
library(sf)
library(mapview)
library(tmap)

# issue when st_make_valid typo causing fiji (typo[typo$ID == 60299, ]) to 
# expand across the world, similar issue was reported for sf at https://github.com/r-spatial/sf/issues/1985
# In our case, make valid with s2 turned off solved the issue,
# see below comparison between s2 on and off

# comparing reading typologies and make valid using s2 on and off
# sf - GEOS
sf_use_s2(FALSE)

typo_sf <- st_read("data/typologies/v3.14/Mangrove_Typology_v3_Composite.shp")

# fiji_sf <- typo_sf[typo_sf$Type == "OpenCoast_60299",] %>% 
#   st_make_valid() # make valid without did not create weird rectangle
# 
# mapview(fiji_sf)

# # S2
# sf_use_s2(TRUE)
# 
# typo_s2 <- st_read("data/typologies/v3.14/Mangrove_Typology_v3_Composite.shp") #%>%  
# 
# fiji_s2 <- typo_s2[typo_s2$Type == "OpenCoast_60299",] %>% # before make valid, fiji mapview is 'normal'
#   st_make_valid() # make valid with S2 turned on will create long rectangle spanning across the globe
# 
# mapview(fiji_s2)

# make valid typo shp and save as gpkg
typo <- typo_sf %>% 
  st_make_valid() #%>%
  # st_set_agr("constant") %>% # to turn off warning about planar 
  # st_set_precision(1e8) # issue with intersection taking too long https://github.com/r-spatial/sf/issues/1510
st_write(typo, "data/typologies/v3.14/Mangrove_Typology_v3.14_Composite_valid_sf.gpkg")

# read eez_land
eez <- st_read("data/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp")%>% 
  #select(c(ISO_TER1, ISO_SOV1, geometry)) %>% 
  st_set_agr("constant") %>% 
  st_set_precision(1e8)

# intersect typo with eez
typo_eez <- st_intersection(typo, eez) %>% 
  mutate(Type_MRGID = paste0(Type, "_", MRGID))

# check Fiji
fiji <- typo_eez[typo_eez$Type == "OpenCoast_60299",]
mapview(fiji)

st_write(typo_eez, "data/typologies/v3.14/Mangrove_Typology_v3.14_Composite_valid_eez.gpkg")
