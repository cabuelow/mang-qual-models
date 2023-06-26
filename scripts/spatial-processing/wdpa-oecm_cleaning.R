## typology categarized by protected or unprotected area

library(wdpar)
library(tidyverse)
library(sf)
library(tmap)

# check what is the latest version of data available 
wdpa_latest_version() 

pa_clean <- st_read("data/WDPA/wdpa_woecm_jun2023_clean.gpkg")

# #if you need the latest version, run below step 1-3 to download and clean data
# #1. download latest wdpa and woecm data, for the latest data
# pa <- wdpa_fetch("global", wait = TRUE, download_dir = "data/WDPA/")
# 
# #2. clean data without erasing overlaps
# pa_clean <- wdpa_clean(pa, erase_overlaps = FALSE)
# st_write(pa_clean, "data/WDPA/wdpa_woecm_jun2023_clean.gpkg")

# 
# #3. separate wdpa and wocem, then dissolve overlaps
woecm_dis <- wdpa_dissolve(pa_clean[pa_clean$PA_DEF == "OECM",]) %>%   
  mutate(PA_DEF = "OECM")
write_sf(woecm_dis, "data/WDPA/woecm_jun2023_clean_dissolved.gpkg", )

wdpa_dis <- wdpa_dissolve(pa_clean[pa_clean$PA_DEF == "PA",]) %>%   
  mutate(PA_DEF = "PA")
write_sf(wdpa_dis, "data/WDPA/wdpa_jun2023_clean_dissolved.gpkg")

pa_dis <- rbind(wdpa_dis, woecm_dis)

write_sf(pa_dis, "data/WDPA/wdpa_woecm_jun2023_clean_dissolved.gpkg")

# # 4. Recommendation for erase overlaps, use st_erase_overlaps by IUCN_CAT order
iucn_order <- c("Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Assigned",
                "Not Reported", "Not Applicable")

# subset overlapping polygons
overlap <- st_intersects(pa_clean)
pa_sub <- pa_clean[lengths(overlap) >=2,]

write_sf(pa_sub, "data/WDPA/wdpa_woecm_jun2023_clean_overlap_only.gpkg")

# erase overlaps, keep higher IUCN categories
pa_sub <- pa_sub[order(factor(pa_sub$IUCN_CAT, level = iucn_order)),]
pa_sub2 <- st_erase_overlaps(pa_sub)

# combine overlaps with non-overlaps
pa_clean2 <- rbind(pa_clean[!pa_sub], pa_sub2)

write_sf(pa_clean2, "data/WDPA/wdpa_woecm_jun2023_clean_no_overlap.gpkg")
