# take spatial data processed to mangrove typological units
# convert to meaningful input for climate risk model
# and combine into master dataframe

library(tidyverse)
library(sf)
library(tmap)
tmap_mode('view')

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
typ_area <- read.csv('outputs/processed-data/typology-area.csv')
dat <- list() # list to store wrangled dat

# do pressure classifications within a range for sensitivity analysis
# mid-range in sensitvity ranges below is used for main results
sens <- c(30, 40, 50, 60, 70)
sens_fdrought <- c(2, 3, 4, 5, 6)
sens_hdrought <- c(0.5, 1, 1.5, 2, 2.5)

# loop through pressure classifications and save in master dataframe
tmp <- list()
for(i in seq_along(sens)){
  
#### coastal squeeze
# TODO: replace below ifelses with casewhens where possible
coast <- read.csv('outputs/processed-data/coastal-population.csv') %>% 
  select(-sum.pop_size_lecz_2020) %>% 
  mutate_at(vars(sum.pop_size_lecz_2000:sum.pop_size_lecz_2100), ~./area_ha) %>% # divide by lecz area to get population density
  mutate(sum.pop_size_lecz_2000 = ifelse(area_ha == 0, max(.$sum.pop_size_lecz_2000, na.rm = T), sum.pop_size_lecz_2000),# assign max pop. density to units with 0 lecz area (indicates no space for mangroves to grow)
         sum.pop_size_lecz_2040 = ifelse(area_ha == 0, max(.$sum.pop_size_lecz_2040, na.rm = T), sum.pop_size_lecz_2040),
         sum.pop_size_lecz_2060 = ifelse(area_ha == 0, max(.$sum.pop_size_lecz_2060, na.rm = T), sum.pop_size_lecz_2060),
         sum.pop_size_lecz_2100 = ifelse(area_ha == 0, max(.$sum.pop_size_lecz_2100, na.rm = T), sum.pop_size_lecz_2100)) %>% 
  mutate(csqueeze = ifelse(log(sum.pop_size_lecz_2000) > quantile(log(.$sum.pop_size_lecz_2000[.$sum.pop_size_lecz_2000>0]), 0.66), 'High', NA)) %>% 
  mutate(csqueeze = ifelse(log(sum.pop_size_lecz_2000) < quantile(log(.$sum.pop_size_lecz_2000[.$sum.pop_size_lecz_2000>0]), 0.66) & log(sum.pop_size_lecz_2000) > quantile(log(.$sum.pop_size_lecz_2000[.$sum.pop_size_lecz_2000>0]), 0.33), 'Medium', csqueeze)) %>% 
  mutate(csqueeze = ifelse(log(sum.pop_size_lecz_2000) < quantile(log(.$sum.pop_size_lecz_2000[.$sum.pop_size_lecz_2000>0]), 0.33), 'Low', csqueeze)) %>% 
  mutate(csqueeze = ifelse(sum.pop_size_lecz_2000 == 0, 'None', csqueeze)) %>% 
  mutate(fut_csqueeze = ifelse(log(sum.pop_size_lecz_2100) > quantile(log(.$sum.pop_size_lecz_2100[.$sum.pop_size_lecz_2100>0]), 0.66), 'High', NA)) %>% 
  mutate(fut_csqueeze = ifelse(log(sum.pop_size_lecz_2100) < quantile(log(.$sum.pop_size_lecz_2100[.$sum.pop_size_lecz_2100>0]), 0.66) & log(sum.pop_size_lecz_2100) > quantile(log(.$sum.pop_size_lecz_2100[.$sum.pop_size_lecz_2100>0]), 0.33), 'Medium', fut_csqueeze)) %>% 
  mutate(fut_csqueeze = ifelse(log(sum.pop_size_lecz_2100) < quantile(log(.$sum.pop_size_lecz_2100[.$sum.pop_size_lecz_2100>0]), 0.33), 'Low', fut_csqueeze)) %>% 
  mutate(fut_csqueeze = ifelse(sum.pop_size_lecz_2100 == 0, 'None', fut_csqueeze)) %>% 
  select(Type, csqueeze, fut_csqueeze) %>% 
  mutate(fut_csqueeze = ifelse(csqueeze == 'Low' & fut_csqueeze == 'None', 'Low', fut_csqueeze),
         fut_csqueeze = ifelse(csqueeze == 'Medium' & fut_csqueeze == 'None', 'Medium', fut_csqueeze),
         fut_csqueeze = ifelse(csqueeze == 'High' & fut_csqueeze == 'None', 'High', fut_csqueeze),
         fut_csqueeze = ifelse(csqueeze == 'Medium' & fut_csqueeze == 'Low', 'Medium', fut_csqueeze),
         fut_csqueeze = ifelse(csqueeze == 'High' & fut_csqueeze == 'Low', 'High', fut_csqueeze),
         fut_csqueeze = ifelse(csqueeze == 'High' & fut_csqueeze == 'Medium', 'High', fut_csqueeze))

# map to check

typ2 <- typ %>% 
  left_join(coast)
qtm(typ2, dots.col = 'csqueeze') 
qtm(typ2, dots.col = 'fut_csqueeze')
dat[[1]] <- coast # if happy add to dat list

#### sediment supply

sed <- read.csv('outputs/processed-data/free-flowing-rivers.csv') %>%
  left_join(select(st_drop_geometry(typ), Type, Class), by = 'Type') %>% 
  mutate(SED_weighted_average = ifelse(Class != 'Delta', 100, SED_weighted_average)) %>%
  mutate(sed_supp = ifelse(SED_weighted_average >= sens[i], 'Low', NA)) %>% 
  mutate(sed_supp = ifelse(SED_weighted_average < sens[i], 'High', sed_supp)) %>% 
  select(Type, sed_supp)
  
# map to check

typ2 <- typ %>% 
  left_join(sed)
qtm(typ2, dots.col = 'sed_supp')
dat[[2]] <- sed # if happy add to dat list

#### future dams

dams <- read.csv('outputs/processed-data/future-dams.csv') %>% 
  mutate(fut_dams = factor(ifelse(number_future_dams >= 1, 1, 0))) %>%
  select(Type, fut_dams)

# map to check

typ2 <- typ %>% 
  left_join(dams)
qtm(typ2, dots.col = 'fut_dams')
dat[[3]] <- dams # if happy add to dat list

#### future sea level rise

fslr <- read.csv('outputs/processed-data/future-slr.csv') %>% 
  mutate(fut_slr = ifelse(slr_m_2081_2100 > quantile(.$slr_m_2081_2100, sens[i]/100), 1, 0)) %>% 
  select(Type, fut_slr)
  
# map to check

typ2 <- typ %>% 
  left_join(fslr)
qtm(typ2, dots.col = 'fut_slr')
dat[[4]] <- fslr # if happy add to dat list

#### antecedent sea level rise

aslr <- read.csv('outputs/processed-data/antecedent-slr.csv') %>% 
  mutate(ant_slr = ifelse(local_msl_trend > quantile(.$local_msl_trend, sens[i]/100), 1, 0)) %>% 
  select(Type, ant_slr)

# map to check

typ2 <- typ %>% 
  left_join(aslr)
qtm(typ2, dots.col = 'ant_slr')
dat[[5]] <- aslr # if happy add to dat list

#### future subsidence from groundwater extraction

fgw <- read.csv('outputs/processed-data/gw_subsid_2040.csv') %>% 
  mutate(fut_gwsub = ifelse(mode >= sens_fdrought[i], 1, 0)) %>% 
  select(Type, fut_gwsub)

# map to check

typ2 <- typ %>% 
  left_join(fgw)
qtm(typ2, dots.col = 'fut_gwsub')
dat[[6]] <- fgw # if happy add to dat list

#### current subsidence from groundwater extraction

gw <- read.csv('outputs/processed-data/gw_subsid_2010.csv') %>% 
  mutate(gwsub = ifelse(mode >= sens_fdrought[i], 1, 0)) %>% 
  select(Type, gwsub)

# map to check

typ2 <- typ %>% 
  left_join(gw)
qtm(typ2, dots.col = 'gwsub') 
dat[[7]] <- gw # if happy add to dat list

#### future drought

fdro <- read.csv('outputs/processed-data/future-stand-precip-index.csv') %>% 
  mutate(fut_drought = ifelse(mean.spi_change_percent_2081_2100 < -sens[i], 1, 0)) %>% 
  select(Type, fut_drought)

# map to check

typ2 <- typ %>% 
  left_join(fdro)
qtm(typ2, dots.col = 'fut_drought') 
dat[[8]] <- fdro # if happy add to dat list

#### future extreme rainfall

frain <- read.csv('outputs/processed-data/future-stand-precip-index.csv') %>% 
  mutate(fut_ext_rain = ifelse(mean.spi_change_percent_2081_2100 > sens[i], 1, 0)) %>% 
  select(Type, fut_ext_rain)

# map to check

typ2 <- typ %>% 
  left_join(frain)
qtm(typ2, dots.col = 'fut_ext_rain') 
dat[[9]] <- frain # if happy add to dat list

#### historical cyclones

cyc <- read.csv('outputs/processed-data/cyclone-tracks-wind_1996_2020.csv') %>% 
  mutate(storms = ifelse(cyclone_tracks_1996_2020 > 1, 1, 0)) %>% 
  select(Type, storms)

# map to check

typ2 <- typ %>% 
  left_join(cyc)
qtm(typ2, dots.col = 'storms') 
dat[[10]] <- cyc # if happy add to dat list

#### historical drought

hdro <- read.csv('outputs/processed-data/historical-drought.csv') %>% 
  mutate(hist_drought = ifelse(min_spei_1996_2020 < -sens_hdrought[i], 1, 0)) %>% 
  select(Type, hist_drought)

# map to check

typ2 <- typ %>% 
  left_join(hdro)
qtm(typ2, dots.col = 'hist_drought') 
dat[[11]] <- hdro # if happy add to dat list

#### historical extreme rainfall

hrain <- read.csv('outputs/processed-data/historical-drought.csv') %>% 
  mutate(hist_ext_rain = ifelse(max_spei_1996_2020 > sens_hdrought[i], 1, 0)) %>% 
  select(Type, hist_ext_rain)

# map to check

typ2 <- typ %>% 
  left_join(hrain)
qtm(typ2, dots.col = 'hist_ext_rain') 
dat[[12]] <- hrain # if happy add to dat list

#### tidal range

tide <- read.csv('data/typologies/SLR_Data.csv') %>% 
  select(Type, Tidal_Class)

# map to check

typ2 <- typ %>% 
  left_join(tide)
qtm(typ2, dots.col = 'Tidal_Class') 
dat[[13]] <- tide # if happy add to dat list

#### future cyclones

fils <- list.files('outputs/processed-data/', pattern = '10000yrs.csv', full.names = T)
fstorms <- lapply(fils, read.csv) %>% 
  lapply(select, Type, cyclone_occurrences_10000yrs) %>% 
  reduce(left_join, by = 'Type') %>% 
  pivot_longer(cols = cyclone_occurrences_10000yrs.x:cyclone_occurrences_10000yrs.y.y, names_to = 'model', values_to = 'cyclone_occurrences_10000yrs') %>% 
  group_by(Type) %>% 
  summarise(cyclone_occurrences_10000yrs = median(cyclone_occurrences_10000yrs)) %>% 
  mutate(fut_storms = 1 - (1 - (cyclone_occurrences_10000yrs/10000))^(2050-2023)) %>% 
  mutate(fut_storms = ifelse(fut_storms > sens[i]/100, 1, 0)) %>% 
  select(Type, fut_storms)

# map to check

typ2 <- typ %>% 
  left_join(fstorms)
qtm(typ2, dots.col = 'fut_storms') 
dat[[14]] <- fstorms # if happy add to dat list

#### propagule establishment distances

propest <- read.csv('outputs/processed-data/propagule-establishment-distances.csv') %>% 
  left_join(typ_area) %>% 
  mutate(prop_establishment = ifelse(extant_loss_dist_mean_m == 'no_loss', NA, extant_loss_dist_mean_m)) %>% 
  mutate(prop_establishment = as.numeric(ifelse(prop_establishment == 'no_extant', NA, prop_establishment))) %>% 
  mutate(prop_establishment = prop_establishment/area_ha) %>% 
  mutate(prop_establishment = ifelse(is.na(prop_establishment), max(.$prop_establishment, na.rm = T), prop_establishment)) %>% 
  mutate(prop_estab = ifelse(prop_establishment > quantile(.$prop_establishment, 0.66), 'Low', NA),
         prop_estab = ifelse(prop_establishment < quantile(.$prop_establishment, 0.33), 'High', prop_estab),
         prop_estab = ifelse(is.na(prop_estab), 'Medium', prop_estab)) %>% 
  select(Type, prop_estab)
  

# map to check

typ2 <- typ %>% 
  left_join(propest)
qtm(typ2, dots.col = 'prop_estab') 
dat[[15]] <- propest # if happy add to dat list

#### landward vs. seaward loss

land <- read.csv('data/typologies/gmw_v3_change_96_20_land_typologies.csv') %>% 
  group_by(chng_type, Type) %>% 
  summarise(area_ha = sum(area_ha)) %>% 
  pivot_wider(names_from = 'chng_type', values_from = 'area_ha') %>% 
  mutate(gain = ifelse(is.na(gain), 0, gain),
         loss = ifelse(is.na(loss), 0, loss)) %>% 
  mutate(land_net_change_ha = gain - loss) %>% 
  select(Type, land_net_change_ha)
  
sea <- read.csv('data/typologies/gmw_v3_change_96_20_sea_typologies.csv') %>% 
  group_by(chng_type, Type) %>% 
  summarise(area_ha = sum(area_ha)) %>% 
  pivot_wider(names_from = 'chng_type', values_from = 'area_ha') %>% 
  mutate(gain = ifelse(is.na(gain), 0, gain),
         loss = ifelse(is.na(loss), 0, loss)) %>% 
  mutate(sea_net_change_ha = gain - loss) %>% 
  select(Type, sea_net_change_ha)

# map to check

typ2 <- typ %>% 
  left_join(land) %>% 
  left_join(sea) %>% 
  mutate(sea_net_change = case_when(sea_net_change_ha > 0 ~ 'Gain',
                                    sea_net_change_ha < 0 ~ 'Loss',
                                    is.na(sea_net_change_ha) ~ 'Neutral'),
         land_net_change = case_when(land_net_change_ha > 0 ~ 'Gain',
                                    land_net_change_ha < 0 ~ 'Loss',
                                    is.na(land_net_change_ha) ~ 'Neutral')) %>% 
  select(Type, sea_net_change:land_net_change)

qtm(typ2, dots.col = 'sea_net_change') 
qtm(typ2, dots.col = 'land_net_change') 
dat[[16]] <- st_drop_geometry(typ2) # if happy add to dat list

# merge into final master database

tmp[[i]] <- data.frame(pressure_def = i, Reduce(full_join, dat))

}

mast.dat <- do.call(rbind, tmp)

# save

write.csv(mast.dat, 'outputs/master-dat.csv', row.names = F)

# End here
