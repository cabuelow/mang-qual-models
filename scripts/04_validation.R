# validation
# for now use goldberg drivers to say whether a pressure is present or not (might need a threshold?)
# simulate with those combination of pressures
# compare probability of predicted loss and gain with gross loss and gain
# eventually use Lagomasino et al. drivers 2.0 - gain and loss

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
library(patchwork)
library(caret)
source('scripts/models.R')
source('scripts/helpers.R')
sf_use_s2(FALSE)

typ_points <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid_centroids.gpkg')
hydro <- read.csv('data/Hydro_Dat.csv')
slr <- read.csv('data/SLR_Data.csv')
world <- data("World")

# join attributes to typologies

typ2 <- typ_points %>% 
  st_drop_geometry() %>% 
  left_join(slr, 'Type') %>% 
  inner_join(hydro) # note missing some hydro data

# classify pressure and biophysical scenarios for each typology for validation
# if proportion of loss due to a pressure is >0, make it a 1 (i.e., it is present)
# also tide, hydro connectivity, and sediment supply

dat <- typ2 %>% 
  mutate(Connectivity = ntile(DOF_Mean, 3), 
         SedSupply = ntile((100-SED_Mean),2)) %>% 
  mutate(Erosion = ifelse(Erosion > 0.1, 1, 0),
         Extreme_Weather = ifelse(Extreme_Weather > 0.1, 1, 0),
         Settlement = ifelse(Settlement > 0.1, 1, 0),
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
  
  datselect <- select(dat[i,], Extreme_Weather, Settlement, Erosion) %>% 
    pivot_longer(Extreme_Weather:Erosion, names_to = 'press', values_to = 'vals') %>% 
    filter(vals == 1) %>% 
    mutate(press = recode(press, 'Extreme_Weather' = "Cyclones", 
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

write.csv(allout, 'outputs/typology_outcomes_validation.csv', row.names = F)

allout <- read.csv('outputs/typology_outcomes_validation.csv')

# join outcomes to typologies, and compare probability of loss/gain with gross loss/gain
# for landward and seaward

landsea <- allout %>% 
  pivot_wider(id_cols = -Prob_gain_neutral, names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Land_Gain = ifelse(LandwardMang >= 50, 1, 0),
         Land_Ambig = ifelse(LandwardMang <50 & LandwardMang > -50, 1, 0),
         Land_Loss = ifelse(LandwardMang <= -50, 1, 0),
         Sea_Gain = ifelse(SeawardMang >= 50, 1, 0),
         Sea_Ambig = ifelse(SeawardMang <50 & SeawardMang > -50, 1, 0),
         Sea_Loss = ifelse(SeawardMang <= -50, 1, 0)) %>% 
  filter(Land_Ambig != 1 & Sea_Ambig != 1) %>% # filter out ambiguous predictions
  mutate(SeaLand_Gain = ifelse(Land_Gain == 1 & Sea_Loss == 0 | Land_Loss == 0 & Sea_Gain == 1, 1, 0),
         SeaLand_Loss = ifelse(Land_Gain == 0 & Sea_Loss == 1 | Land_Loss == 1 & Sea_Gain == 0, 1, 0)) %>% 
  filter(SeaLand_Gain == 1 | SeaLand_Loss == 1) %>% # filter out predictions where land cancels sea
  inner_join(typ2, by = 'Type') %>%
  mutate(Net_Gain = ifelse(Net_Change >= 0, 1, 0),
         Net_Loss = ifelse(Net_Change < 0, 1, 0)) %>% 
  as.data.frame()
head(landsea)

# now do confusion matrix of non-ambiguous predictions of gain or loss, compared to net gain or loss
# also remove ambiguous predictions where land-cancels sea

conMatrix <- confusionMatrix(factor(c(landsea$SeaLand_Loss,landsea$SeaLand_Gain), levels = c(0,1)), 
                             factor(c(landsea$Net_Loss, landsea$Land_Gain)))
conMatrix

# plot probability of loss against gross loss

a <- ggplot(land) +
  geom_point(aes(x = Prob_change, y = Net_Change), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ylab('Net change in mangrove area (km2)') +
  xlab('') +
  ggtitle('A) Landward mangrove') +
  theme_classic() +
  theme(plot.title = element_text(size = 11))
a

b <- ggplot(sea) +
  geom_point(aes(x = Prob_change, y = Net_Change), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ylab('Net change in mangrove area (km2)') +
  xlab('Probability of loss (-) and gain (+)') +
  ggtitle('B) Seaward mangrove') +
  theme_classic() +
  theme(plot.title = element_text(size = 11))
b

gg_axis <- cowplot::get_plot_component(ggplot() +
                                         labs(y = "Net change in mangrove area (km2)"), "ylab-b")

c <- (a/b & labs(y = NULL))
c

blanklabelplot<-ggplot()+labs(y="Net change in mangrove area (km2)")+theme_classic()+ 
  guides(x = "none", y = "none")

blanklabelplot+c+plot_layout(widths=c(1,1000))

ggsave('outputs/validation_loss.png', width = 4, height = 4)

# calculate proportion of points that are valid vs. invalid

land <- land %>% 
  mutate(valid_loss = ifelse(Net_Change < 0 & Prob_change < 0, 1, 0),
         valid_gain = ifelse(Net_Change > 0 & Prob_change > 0, 1, 0))
(sum(land$valid_gain) + sum(land$valid_loss))/nrow(land) * 100

sea <- sea %>% 
  mutate(valid_loss = ifelse(Net_Change < 0 & Prob_change < 0, 1, 0),
         valid_gain = ifelse(Net_Change > 0 & Prob_change > 0, 1, 0))
(sum(sea$valid_gain) + sum(sea$valid_loss))/nrow(sea) * 100

# join outcomes to spatial typologies for mapping

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
            main.title = 'B) Landward mangrove',
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

tmap_save(maps, 'outputs/map_validation.png', width = 10, height = 5)
