# estimate historical drought conditions in each forest patch
# extract average of all pixel values in each year (1996-2020)
# using a 12-month time scale

library(terra)
library(exactextractr)
library(lubridate)
library(tidyverse)
library(sf)
library(tmap)

typ <- st_read('data/typologies/v3.14/Mangrove_Typology_v3.14_Composite_valid.gpkg') %>% st_buffer(0)
dat <- rast('data/SPEI/spei12.nc')
dat1 <- dat[[50]]
#writeRaster(dat1, 'data/SPEI/exampleraste.tiff')
plot(dat1)

# subset data for every month from 1996-2020

dates <- data.frame(date = seq(as.Date("1901-01-16"), as.Date("2020-12-16"), by = "months")) %>% 
  mutate(month = month(date),
         year = year(date), 
         index = row.names(.)) %>%
  filter(year %in% c(1996:2020)) #& month == 12)
  #filter(year %in% c(1996:2020) & month == 12)

datsub <- dat[[as.numeric(dates$index)]]

# extract

drought <- exact_extract(datsub, typ, 'mean')
colnames(drought) <- paste0(rep('mean_monthly_spei_', 300), paste0(dates$month, '_', dates$year))

final.df <- data.frame(st_drop_geometry(typ), drought)

# 10km buffer to gapfill those near land, others get 0 for drought - no land

# which are missing values - there are 319

miss <- filter(typ, Type %in% filter(final.df, is.na(mean_monthly_spei_11_2020))$Type) %>% st_buffer(20000)
#st_write(miss, 'miss.shp', overwrite = T, append = F)

drought2 <- exact_extract(datsub, miss, 'mean')
colnames(drought2) <- paste0(rep('mean_monthly_spei_', 300), paste0(dates$month, '_', dates$year))
drought2 <- data.frame(st_drop_geometry(miss), drought2)
drought2[is.na(drought2)] <- 0

final.df2 <- final.df %>%
  filter(!Type %in% drought2$Type) %>% 
  rbind(drought2)

final.df2$mean_spei_1996_2020 = rowMeans(select(final.df2, mean_monthly_spei_1_1996:mean_monthly_spei_12_2020))
final.df2$max_spei_1996_2020 = apply(select(final.df2, mean_monthly_spei_1_1996:mean_monthly_spei_12_2020), 1, max)
final.df2$min_spei_1996_2020 = apply(select(final.df2, mean_monthly_spei_1_1996:mean_monthly_spei_12_2020), 1, min)

write.csv(final.df2, 'outputs/processed-data/historical-drought.csv', row.names = F)

# check Estuary_30003

tmap_mode('view')
sub <- typ %>% filter(Type == 'Estuary_30003') %>% st_buffer(20000)
qtm(sub)
subr <- crop(dat1, vect(sub), mask = F)
qtm(subr) + qtm(sub) + qtm(filter(typ, Type == 'Estuary_30003'))
exact_extract(subr, filter(typ, Type == 'Estuary_30003'), 'mean')
exact_extract(subr, sub, 'mean') # works for the buffer

sub <- typ %>% filter(Type == 'Estuary_30003') %>% st_buffer(0)
qtm(sub)
exact_extract(subr, sub, 'mean') # works for the buffer
