# join exposure factors for madagascar/belize

library(tidyverse)
library(sf)
library(tmap)
tmap_mode('view')

typ <- st_read('data/typologies/Mangrove_Typology_v3.14_Composite_valid_centroids.gpkg')
dat <- list()

#### tidal range

dat[[1]] <- read.csv('data/typologies/SLR_Data.csv') %>% 
  filter(Type %in% typ$Type) %>% 
  select(Type, Tidal_Class)

#### future sea level rise

dat[[2]] <- read.csv('outputs/processed-data/future-slr.csv') %>% 
  filter(Type %in% typ$Type) %>% 
  select(Type, slr_m_2041_2060) %>% 
  rename(future_slr = slr_m_2041_2060)

#### antecedent sea level rise

dat[[3]] <- read.csv('outputs/processed-data/antecedent-slr.csv') %>% 
  filter(Type %in% typ$Type) %>% 
  rename(antecedent_slr = local_msl_trend)

#### future drought/extreme rainfall

dat[[4]] <- read.csv('outputs/processed-data/future-stand-precip-index.csv') %>% 
  filter(Type %in% typ$Type) %>% 
  mutate(future_drought = abs(ifelse(mean.spi_change_percent_2041_2060 > 0, 0, mean.spi_change_percent_2041_2060)),
         future_extreme_rainfall = ifelse(mean.spi_change_percent_2041_2060 < 0, 0, mean.spi_change_percent_2041_2060)) %>% 
  select(Type, future_drought, future_extreme_rainfall)

#### historical drought/extreme rainfall

dat[[5]] <- read.csv('outputs/processed-data/historical-drought.csv') %>% 
  filter(Type %in% typ$Type) %>% 
  mutate(historical_drought = abs(ifelse(min_spei_1996_2020 > 0, 0, min_spei_1996_2020)),
         historical_extreme_rainfall = ifelse(max_spei_1996_2020 < 0, 0, max_spei_1996_2020)) %>% 
  select(Type, historical_drought, historical_extreme_rainfall)
  
#### historical cyclones

dat[[6]] <- read.csv('outputs/processed-data/cyclone-tracks-wind_1996_2020.csv') %>% 
  filter(Type %in% typ$Type) %>% 
  select(Type, cyclone_tracks_1996_2020) %>% 
  rename(historical_storms = cyclone_tracks_1996_2020)

#### future cyclones

fils <- list.files('outputs/processed-data/', pattern = '10000yrs.csv', full.names = T)
dat[[7]] <- lapply(fils, read.csv) %>% 
  lapply(select, Type, cyclone_occurrences_10000yrs) %>% 
  reduce(left_join, by = 'Type') %>% 
  pivot_longer(cols = cyclone_occurrences_10000yrs.x:cyclone_occurrences_10000yrs.y.y, names_to = 'model', values_to = 'cyclone_occurrences_10000yrs') %>% 
  group_by(Type) %>% 
  summarise(cyclone_occurrences_10000yrs = median(cyclone_occurrences_10000yrs)) %>% 
  mutate(future_storms = 1 - (1 - (cyclone_occurrences_10000yrs/10000))^(2050-2023)) %>% 
  filter(Type %in% typ$Type) %>% 
  select(Type, future_storms)

# merge into final master database

mast_dat <- data.frame(Reduce(full_join, dat))

write.csv(mast_dat, 'outputs/climate-exposure_global.csv', row.names = F)

# classify as low, medium or high according to percentile thresholds

# percentile thresholds for classifying as low, medium or high (anything inbetween 0.3 and 0.7 is medium)
low <- 0.3
high <- 0.7

mast_dat2 <- mast_dat %>% 
  mutate(future_slr_class = case_when(future_slr <= quantile(.$future_slr, low) ~ 'Low',
                                future_slr >= quantile(.$future_slr, high) ~ 'High',
                                .default = 'Medium'),
         antecedent_slr_class = case_when(antecedent_slr <= quantile(.$antecedent_slr, low) ~ 'Low',
                                    antecedent_slr >= quantile(.$antecedent_slr, high) ~ 'High',
                                   .default = 'Medium'),
         future_drought_class = case_when(future_drought <= quantile(.$future_drought, low) ~ 'Low',
                                          future_drought >= quantile(.$future_drought, high) ~ 'High',
                                          .default = 'Medium'),
         future_extreme_rainfall_class = case_when(future_extreme_rainfall == 0 ~ 'Low', # set 0 as low and exclude from percentile classification
                                          future_extreme_rainfall <= quantile(.$future_extreme_rainfall, low) ~ 'Low',
                                          future_extreme_rainfall >= quantile(.$future_extreme_rainfall, high) ~ 'High',
                                          .default = 'Medium'),
         historical_drought_class = case_when(historical_drought <= quantile(.$historical_drought, low) ~ 'Low',
                                          historical_drought >= quantile(.$historical_drought, high) ~ 'High',
                                          .default = 'Medium'),
         historical_extreme_rainfall_class = case_when(historical_extreme_rainfall == 0 ~ 'Low', # set 0 as low and exclude from percentile classification
                                                   historical_extreme_rainfall <= quantile(.$historical_extreme_rainfall, low) ~ 'Low',
                                                   historical_extreme_rainfall >= quantile(.$historical_extreme_rainfall, high) ~ 'High',
                                                   .default = 'Medium'),
         historical_storms_class = case_when(historical_storms <= quantile(.$historical_storms, low) ~ 'Low',
                                            historical_storms >= quantile(.$historical_storms, high) ~ 'High',
                                           .default = 'Medium'),
         future_storms_class = case_when(future_storms <= quantile(.$future_storms, low) ~ 'Low',
                                             future_storms >= quantile(.$future_storms, high) ~ 'High',
                                             .default = 'Medium'))

typdat <- typ %>% 
  left_join(mast_dat2)

qtm(typdat, dots.col = 'future_slr_class')
qtm(typdat, dots.col = 'antecedent_slr_class')
qtm(typdat, dots.col = 'future_drought_class')
qtm(typdat, dots.col = 'future_extreme_rainfall_class')
qtm(typdat, dots.col = 'historical_drought_class')
qtm(typdat, dots.col = 'historical_extreme_rainfall_class')
qtm(typdat, dots.col = 'historical_storms_class')
qtm(typdat, dots.col = 'future_storms_class')

