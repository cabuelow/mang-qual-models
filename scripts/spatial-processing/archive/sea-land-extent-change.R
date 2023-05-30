# determine whether shoreline is in 100m search radius of each patch of mangrove loss or gain
# if it isn't, classify as landward loss/gain. if it is, classify as seaward loss/gain
# to speed up could dissolve patches of gain or loss by typology - but then higher risk of misclassification, so won't do
# NOTE** this approach ended up taking way too long, could not use
library(nngeo)
library(sf)
library(dplyr)

# set HPC working directory
setwd('/export/home/s2988833')

file <- 'gains_96_07_typology'

shore <- st_read('data/shoreline/USGS_WCMC_GlobalShoreline.shp')
shore2 <- st_read('data/shoreline/USGS_WCMC_GlobalIslandsv1_BigIslands_lines.shp')
shore3 <- st_read('data/shoreline/USGS_WCMC_GlobalIslandsv1_SmallIslands_lines.shp')
shore4 <- st_read('data/shoreline/USGS_WCMC_GlobalIslandsv1_VerySmallIslands_lines.shp')
shoreall <- rbind(shore, shore2, shore3, shore4)
change <- st_read(paste0('extent-change-gmwv3/', file, '.shp')) %>% st_transform(st_crs(shore)) %>% st_make_valid()
change$id <- 1:nrow(change)

# find the 'nearest neighbour' shoreline to each mangrove patch within a 100 metre search radius
# if there is no shoreline within 100m it will get a 0

system.time(
nn <- st_nn(change, shoreall, k = 1, maxdist = 100))
nn <- lapply(nn, function(x) ifelse(length(x) == 0, 0, x))

# turn into dataframe with seaward vs. landward classification

df <- data.frame(id = 1:length(nn), class = do.call(rbind, nn)) %>% 
  mutate(class = ifelse(class == 0, 'landward', 'seaward'))

# join to spatial dataframe by id and save

change <- change %>% 
  left_join(df, by = 'id')

st_write(change, paste0(file, '_class.gpkg'))