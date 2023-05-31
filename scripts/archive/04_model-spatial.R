library(sf)
library(tmap)
library(tidyverse)
library(igraph)
library(QPress)
library(patchwork)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(scales)
source('scripts/models_v2.R')
source('scripts/helpers_v2.R')
sf_use_s2(FALSE)

# convert typologies to points for faster mapping
#typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg')
#typ_points <- st_centroid(typ)
#st_write(typ_points, 'data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
hydro <- read.csv('data/Hydro_Dat.csv') # only have for areas of restorable loss
slr <- read.csv('data/SLR_Data.csv')
world <- data("World")

# join attributes to typologies

typ2 <- typ_points %>% 
  st_drop_geometry() %>% 
  left_join(slr, 'Type') %>% 
  inner_join(hydro) # note missing some hydro data

# classify pressure and biophysical scenarios for each typology
# deciles for each pressure, and anything in top 10th decile the pressure is present
# also tide, hydro connectivity, and sediment supply

dat <- typ2 %>% 
  mutate(Erosion = ntile(Erosion, 10),
         Extreme_Weather = ntile(Extreme_Weather, 10),
         Settlement = ntile(Settlement, 10),
         Connectivity = ntile(DOF_Mean, 3), 
         SedSupply = ntile((100-SED_Mean),2)) %>% 
  mutate(Erosion = ifelse(Erosion == 10, 1, 0),
         Extreme_Weather = ifelse(Extreme_Weather == 10, 1, 0),
         Settlement = ifelse(Settlement == 10, 1, 0),
         SLR = ifelse(SLR_Class == 'High', 1, 0),
         Connectivity = ifelse(Connectivity == 3, 'H', Connectivity),
         Connectivity = ifelse(Connectivity == 2, 'M', Connectivity),
         Connectivity = ifelse(Connectivity == 1, 'L', Connectivity),
         SedSupply = ifelse(SedSupply == 2, 'H', SedSupply),
         SedSupply = ifelse(SedSupply== 1, 'L', SedSupply),
         Tidal_Class = ifelse(Tidal_Class == 'Micro', 'H', Tidal_Class),
         Tidal_Class = ifelse(Tidal_Class == 'Meso', 'M', Tidal_Class),
         Tidal_Class = ifelse(Tidal_Class == 'Macro', 'L', Tidal_Class))

# run a model for each typlogy, based on pressures and biophysical settings

# set up relative edge constraint scenarios
# high sed supply model, sediment -> subVol will be greater than SLR neg interactions
# vice versa for low sed supply model

model.scenarios <- list(parse.constraints(c('SeaLevelRise -* SeawardMang < Sediment -> SubVol', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), modelB),
                        parse.constraints(c('Sediment -> SubVol < SeaLevelRise -* SeawardMang', 'LandwardMang -> SubVol < SeawardMang -> SubVol'), modelB))
names(model.scenarios) <- c('High Sediment Supply', 'Low Sediment Supply')

set.seed(123)
numsims <- 1000

tmp <- list()

for(i in 1:nrow(dat)){

  # perturbation scenarios
  
datselect <- select(dat[i,], SLR, Extreme_Weather, Settlement, Erosion) %>% 
  pivot_longer(SLR:Erosion, names_to = 'press', values_to = 'vals') %>% 
  filter(vals == 1) %>% 
  mutate(press = recode(press, 'SLR' = "SeaLevelRise", 
                        'Extreme_Weather' = "Cyclones", 
                        'Settlement' = 'CoastalDev'))

if(nrow(datselect) == 0){ # if there are no perturbations, go to next typology
  next
}

press.scenario <- rep(1, nrow(datselect))
names(press.scenario) <- datselect$press

# edge constraint scenarios
# **TODO: make this easier by changing how the function takes these constraints....

datselect2 <- select(dat[i,], Tidal_Class, Connectivity) %>% 
  pivot_longer(Tidal_Class:Connectivity, names_to = 'press', values_to = 'vals') 

con.scenario <- c(datselect2$vals, datselect2$vals[2])

# select model for sediment supply

if(dat[i,]$SedSupply == 'H'){
model <- model.scenarios[[1]]
}else{
  model <- model.scenarios[[2]]
}

# simulate outcomes

      sim <- system.sim_press(numsims, constrainedigraph = model, 
                              from = c('SeaLevelRise', #'SeaLevelRise', #'SeaLevelRise', 
                                       'LandwardAvailableProp', 'SeawardAvailableProp'),
                              to = c('SeawardMang', #'SeawardPropag', #'SeawardEstabSpace', 
                                     'LandwardMang', 'SeawardMang'),
                              class = con.scenario,
                              perturb = press.scenario)

      out <- sim$stableoutcome %>% 
        filter(var %in% c('SeawardMang', 'LandwardMang')) %>% 
        group_by(var) %>% 
        summarise(Prob_gain_neutral = (sum(outcome>=0))/n())

out$Type <- rep(dat[i, 'Type'], nrow(out))
      
tmp[[i]] <- out
}

allout <- do.call(rbind, tmp)
allout$Prob_change <- rescale(allout$Prob_gain_neutral, to = c(-100, 100))

write.csv(allout, 'outputs/typology_outcomes.csv', row.names = F)

allout <- read.csv('outputs/typology_outcomes.csv')

# join outcomes to spatial typologies

landward <- typ_points %>% 
  left_join(filter(allout, var == 'LandwardMang'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

seaward <- typ_points %>% 
  left_join(filter(allout, var == 'SeawardMang'), by = 'Type') %>% 
  st_crop(xmin = -150, ymin = -40, xmax = 180, ymax = 33)

world_mang <- st_crop(World, xmin = -150, ymin = -40, xmax = 180, ymax = 33)

# map

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(landward, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(landward, !is.na(Prob_change))) +
  tm_dots('Prob_change', 
          palette = 'Spectral',
          breaks = c(-100,-50,0,50,100),
          midpoint = 0,
          title = '',
          legend.is.portrait = T) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'A) Landward mangrove',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
lmap

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray95') +
  tm_shape(filter(seaward, is.na(Prob_change))) +
  tm_dots('darkgrey') +
  tm_shape(filter(seaward, !is.na(Prob_change))) +
  tm_dots('Prob_change', 
          palette = 'Spectral',
          breaks = c(-100,-50,0,50,100),
          midpoint = 0,
          title = '',
          legend.is.portrait = T) +
  tm_layout(legend.outside = F,
            #legend.outside.position = 'bottom',
            #legend.position = c(0.35, 0.6),
            title.size = 0.8,
            title.position = c(0.01,0.45),
            legend.title.size = 0.9,
            main.title = 'B) Seaward mangrove',
            title = 'Probability of Loss (red) \nor Neutrality/Gain (blue)',
            main.title.size = 1,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0.8)
smap

maps <- tmap_arrange(lmap, smap, ncol = 1)

tmap_save(maps, 'outputs/map_change.png', width = 10, height = 5)





