library(tidyverse)
library(scales)
library(patchwork)

# plot scenario outcomes

names(models) # names of available models
model_name <- 'mangrove_model'
pthfil <- paste0('outputs/simulation-outcomes/outcomes_', model_name, '.rds') # which outcomes to plot
dat <- readRDS(pthfil) 

# calculate proportion of stable model outcomes that have a positive, negative, or neutral landward/seaward mangrove response

dat2 <- dat %>% 
  filter(var %in% c('SeawardMang', 'LandwardMang') & 
           constraint_scenario %in% c('Macrotidal, High propagule establishment capacity, High coastal squeeze',
                                      'Microtidal, High propagule establishment capacity, High coastal squeeze',
                                      'Macrotidal, High propagule establishment capacity, Low coastal squeeze',
                                      'Microtidal, High propagule establishment capacity, Low coastal squeeze')) %>% #& 
          # !pressure %in% c('Cyclones', 'Dams', 'Drought', 'Erosion', 'Groundwater extraction',
           #                 'Sea-level rise & Cyclones', 'Sea-level rise & Dams', 'Sea-level rise & Drought',
            #                'Sea-level rise & Erosion', 'Sea-level rise & Groundwater extraction')) %>% 
  group_by(model_scenario, constraint_scenario, pressure, var) %>% 
  summarise(Prob_gain_neutral = ((sum(outcome>0) + sum(outcome==0))/n())*100,
            Prob_loss = (sum(outcome<0)/n())*-100) %>%
  mutate(tide = strsplit(constraint_scenario, ', ')[[1]][1],
         coastalsqueeze = strsplit(constraint_scenario, ', ')[[1]][3],
         var = recode(var, 'LandwardMang' = 'Landward mangrove', 'SeawardMang' = 'Seaward mangrove')) %>% 
  mutate(tide = factor(tide, levels = c('Microtidal', 'Mesotidal', 'Macrotidal')))

dat2$Prob_change <- rescale(dat2$Prob_gain_neutral, to = c(-100, 100))

# plot

a <- ggplot(filter(dat2, var == 'Landward mangrove'), 
            aes(tide, pressure, fill = Prob_change)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', 
                       name = 'Probability of Loss (red) or Neutrality/Gain (blue)', 
                       direction = 1,
                       limits = c(-100, 100),
                       breaks = c(-100, -50, 0, 50, 100),
                       labels = c("-100", "-75", "50", '75', '100')) +
  #facet_wrap(~factor(coastalsqueeze)) +
  facet_wrap(vars(factor(model_scenario), factor(coastalsqueeze)), ncol = 2) +
  theme_classic() +
  theme(legend.position = 'none',
        legend.justification = 'left',
        #axis.text.y =  element_blank(),
        strip.text.x = element_text(size = 9),
        axis.title = element_blank()) +
  ggtitle('A) Landward mangrove') +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,))
a

b <- ggplot(filter(dat2, var == 'Seaward mangrove'), 
            aes(tide, pressure, fill = Prob_change)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', 
                       name = 'Probability of Loss (red) or Neutrality/Gain (blue)',
                       direction = 1,
                       limits = c(-100, 100),
                       breaks = c(-100, -50, 0, 50, 100),
                       labels = c("-100", "-75", "50", '75', '100')) +
  facet_wrap(vars(factor(model_scenario), factor(coastalsqueeze)), ncol = 2) +
  #facet_wrap(~factor(model_scenario)) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.justification = 'left',
        #axis.text.y =  element_blank(),
        strip.text.x = element_text(size = 9),
        axis.title = element_blank()) +
  ggtitle('B) Seaward mangrove') +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,))
b

c <- a/b #+ plot_layout(heights = c(0.8, 2))
c 
ggsave(paste0('outputs/heatmap_outputs/', model_name ,'_heatmap.png'), width = 7, height = 10)

