# calculate patch connectivity within each typological unit
# for the most recent year (i.e. 2020)

library(sf)
library(tidyverse)
library(landscapemetrics)

patches <- st_read('data/GMW_v3/gmw_v3_1996-2020_Merge_Typology_join.shp') %>% st_make_valid()
st_write(patches)