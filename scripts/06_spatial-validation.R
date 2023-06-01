# compare hindcast probability of mangrove loss and gain with historical loss and gain to validate network model

library(tidyverse)
source('scripts/helpers/models_v2.R')

# wrangle historical SRS observations of mangrove loss and gain into categories of change (no change, loss, gain, loss and gain)

spatial_dat <- read.csv('outputs/master-dat.csv') %>% 
  mutate(sea_change = sea_gain + sea_loss,
        land_change = land_gain + land_loss) %>% 
  mutate(sea_change_c = ifelse(sea_change == 2, 'Loss & Gain', NA),
         sea_change_c = ifelse(sea_change != 2 & sea_gain ==1, 'Gain', sea_change_c),
         sea_change_c = ifelse(sea_change != 2 & sea_loss ==1, 'Loss', sea_change_c),
         sea_change_c = ifelse(sea_change ==0, 'No change', sea_change_c),
         land_change_c = ifelse(land_change == 2, 'Loss & Gain', NA),
         land_change_c = ifelse(land_change != 2 & land_gain ==1, 'Gain', land_change_c),
         land_change_c = ifelse(land_change != 2 & land_loss ==1, 'Loss', land_change_c),
         land_change_c = ifelse(land_change ==0, 'No change', land_change_c))

# which model outcomes to validate? Get outcomes for that model

names(models) # names of available models
chosen_model_name <- 'mangrove_model'
dat <- read.csv(paste0('outputs/simulation-outcomes/outcomes_', chosen_model_name, '_spatial.csv')) %>% 
  mutate(Prob_change = ifelse(Prob_gain_neutral > 50, Prob_gain_neutral, Prob_loss))

# join outcomes to typologies and compare probability of loss/gain with historical gross loss/gain (1996-2020)

landsea <- dat %>% 
  pivot_wider(id_cols = -c(Prob_loss, Prob_gain_neutral), names_from = 'var', values_from = 'Prob_change') %>% 
  mutate(Land_Gain = ifelse(LandwardMang > 75, 1, 0),
         Land_Ambig = ifelse(LandwardMang <= 50 & LandwardMang >= -75, 1, 0),
         Land_Loss = ifelse(LandwardMang < -75, 1, 0),
         Sea_Gain = ifelse(SeawardMang > 75, 1, 0),
         Sea_Ambig = ifelse(SeawardMang <= 50 & SeawardMang >= -75, 1, 0),
         Sea_Loss = ifelse(SeawardMang < -75, 1, 0)) %>% 
  filter(Land_Ambig != 1 & Sea_Ambig != 1) %>% # filter out ambiguous predictions
  inner_join(select(spatial_dat, Type, sea_gain:land_loss), by = 'Type') %>%
  as.data.frame()
head(landsea)

# now do confusion matrix of non-ambiguous predictions of gain or loss, compared to net gain or loss
# also remove ambiguous predictions where land-cancels sea

conMatrix <- confusionMatrix(factor(landsea$SeaLand_Loss, levels = c(0,1)), 
                             factor(landsea$Net_Loss, levels = c(0,1)))
conMatrix

conMatrix2 <- confusionMatrix(factor(landsea$SeaLand_Gain, levels = c(0,1)), 
                             factor(landsea$Net_Gain, levels = c(0,1)))
conMatrix2

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
