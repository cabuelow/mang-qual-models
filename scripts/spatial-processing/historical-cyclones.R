# estimate tropical storm frequency at each mangrove unit between 1996 and 2020
# number of tropical cyclone occurrences within a 200km buffer of mangrove typology centroids

library(sf)
library(tidyverse)
library(tmap)

typ_cent <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
typ_cent_200buff <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg') %>% st_buffer(200000) %>% st_wrap_dateline()
cyc <- st_read('data/IBTrACs-database/IBTrACS.since1980.list.v04r00.lines/IBTrACS.since1980.list.v04r00.lines.shp') %>% 
  filter(year %in% c(1996:2020))

df <- data.frame(Type = NA, cyclone_tracks = NA, cyclone_wind_knots = NA, cyclone_speed_knots = NA) # df for storing results

for(i in 1:nrow(typ_cent_200buff)){
  
  unit <- typ_cent_200buff[i,]
  
  tracks <- cyc %>% 
    mutate(unit = lengths(st_intersects(cyc, unit, sparse = T)) > 0) %>% 
    filter(unit == 'TRUE')
  
  df[[i,1]] <- unit$Type
  df[[i,2]] <- length(unique(tracks$SID))
  df[[i,3]] <- mean(tracks$USA_WIND, na.rm = T)
  df[[i,4]] <- mean(tracks$STORM_SPD, na.rm = T)
  
}

df2 <- df %>% 
  mutate(cyclone_wind_knots = ifelse(cyclone_tracks == 0, 0, cyclone_wind_knots)) %>% 
  mutate(cyclone_speed_knots = ifelse(cyclone_tracks == 0, 0, cyclone_speed_knots))

# note becuase some buffers are split over the antimeridian, will need to group by typology ID and summarise
# will use mean because it could be the same cyclone in each unit

df.final <- df2 %>% 
  group_by(Type) %>% 
  summarise(cyclone_tracks_1996_2020 = mean(cyclone_tracks),
            cyclone_wind_knots = mean(cyclone_wind_knots),
            cyclone_speed_knots = mean(cyclone_wind_knots))

write.csv(df.final, 'outputs/processed-data/cyclone-tracks-wind_1996_2020.csv', row.names = F)

# check above has worked

sub <- filter(typ_cent_200buff, Type == 'OpenCoast_40178')
cyc_sub <- cyc %>% st_crop(st_bbox(sub))
qtm(sub) + qtm(cyc_sub)

typ_cyc <- typ_cent %>% 
  left_join(df.final)

qtm(typ_cyc, dots.col = 'cyclone_tracks_1996_2020')
