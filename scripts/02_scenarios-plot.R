library(tidyverse)
library(scales)

# plot scenario outcomes

dat <- readRDS('outputs/outcomes.rds')

# calculate proportion of stable models that have positive, negative, or neutral outcome in landward/seaward mangrove response

dat2 <- dat %>% 
  filter(var %in% c('SeawardMang', 'LandwardMang') & 
           constraint_scenario %in% c('Macrotidal, High Hydro-connectivity',
                                      'Mesotidal, High Hydro-connectivity',
                                      'Microtidal, High Hydro-connectivity') &
           pressure %in% c('Sea-level rise', 'Cyclones', 'Sea-level rise & Cyclones', 
                           'Seal-level rise & Groundwater extraction', 'Sea-level rise  & Coastal development', 
                           'Sea-level rise & Erosion', 'Sea-level rise & Drought or Dams')) %>% 
  group_by(model_scenario, constraint_scenario, pressure, var) %>% 
  summarise(Prob_gain = (sum(outcome>0) + (sum(outcome==0))/2)/n())  # divide neutral outcomes by two and add to increase, will just show as ambiguous
dat2$Prob_change <- rescale(dat2$Prob_gain, to = c(-1, 1))

a <- ggplot(filter(dat2, model_scenario == 'High Sediment Supply'), aes(var, pressure, fill = Prob_change)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', 
                       name = 'Probability of Gain (blue) or Loss (red)', 
                       direction = 1) +
  facet_wrap(~constraint_scenario, nrow = 3, ncol = 1) +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        strip.text.x = element_text(size = 8),
        axis.title = element_blank()) +
  ggtitle('High sediment supply') +
  guides(fill = guide_colorbar(title.position = "top"))
a  

b <- ggplot(filter(dat2, model_scenario == 'Low Sediment Supply'), aes(var, pressure, fill = Prob_change)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', name = 'Probability of Gain (red) or Loss (blue)', direction = 1) +
  facet_wrap(~constraint_scenario, nrow = 3, ncol = 1) +
  ggtitle('Low sediment supply') +
  theme(legend.position = 'none',
    axis.text.y =  element_blank(),
    strip.text.x = element_text(size = 8),
    axis.title = element_blank())

a + b

ggsave('outputs/outcomes-heatmap.png', width = 7, height = 10)
