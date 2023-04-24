# future cyclone risk under climate change
# SSP5-8.5, 2015-2050
# estimate tropical storm frequency at each mangrove unit between in 10000 simulation years
# number of tropical cyclone occurrences within a 200km buffer of mangrove typology centroids

library(sf)
library(tidyverse)
library(tmap)

typ_cent <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
typ_cent_200buff <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg') %>% st_buffer(200000) %>% st_wrap_dateline()
GCM <- 'CMCC'
fils <- list.files(paste0('data/STORM-tracks/', GCM, '/'), full.names = T)
dat <- lapply(fils, read.table, fill = T, header = F, sep = ',')
dat <- do.call(rbind, dat)
colnames(dat) <- c('Year', 'Month', 'TC_number', 'Time_step', 'Basin_ID',
                   'Latitude', 'Longitude', 'Min_pressure_hpa', 'Max_wind_speed_m_s', 
                   'Max_wind_radius_km', 'Category', 'Landfall', 'Dist_to_land_km')

# turn into spatialdataframe

dat.sf <- dat %>% 
  mutate(Year_TC_number_Basin = paste0(Year, '_', TC_number, '_', Basin_ID)) %>% 
  #filter(Year == 1) %>% 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) #%>% 
  #group_by(Year_TC_number_Basin) %>% 
  #summarise() %>% 
  #st_cast("LINESTRING")

#st_write(dat.sf, 'data/STORM-tracks/CMCC.gpkg', append = F, overwrite = T)

#dat.sf <- st_read('data/STORM-tracks/CMCC.gpkg')

# get number of unique cyclone occurrences per year, in 10000 year period, forest patches 200 km buffer

df <- data.frame(Type = NA, cyclone_occurrences_10000yrs = NA, cyclone_max_wind_speed_m_s_mean = NA) # df for storing results

for(i in 1:nrow(typ_cent_200buff)){
  
  unit <- typ_cent_200buff[i,]
  
  #system.time(occurrences <- unit %>% st_intersection(dat.sf))
  
  occurrences <- dat.sf %>% 
    mutate(unit = lengths(st_intersects(dat.sf, unit, sparse = T)) > 0) %>% 
    filter(unit == 'TRUE')
  
  df[[i,1]] <- unit$Type
  df[[i,2]] <- length(unique(occurrences$Year_TC_number_Basin))
  df[[i,3]] <- mean(occurrences$Max_wind_speed_m_s, na.rm = T)
  
}

df2 <- df %>% 
  mutate(cyclone_max_wind_speed_m_s_mean = ifelse(cyclone_occurrences_10000yrs == 0, 0, cyclone_max_wind_speed_m_s_mean))

# note becuase some buffers are split over the antimeridian, will need to group by typology ID and summarise
# will use mean because it could be the same cyclone in each unit

df.final <- df2 %>% 
  group_by(Type) %>% 
  summarise(cyclone_occurrences_10000yrs = mean(cyclone_occurrences_10000yrs),
            cyclone_max_wind_speed_m_s_mean = mean(cyclone_max_wind_speed_m_s_mean))

write.csv(df.final, paste0('future-cyclone-occurrences_', GCM, '_10000yrs.csv'), row.names = F)
