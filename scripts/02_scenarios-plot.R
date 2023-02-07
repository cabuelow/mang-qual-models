library(tidyverse)
library(scales)

# plot scenario outcomes

dat <- readRDS('outputs/outcomes_alt-model.rds')

# calculate proportion of stable models that have positive, negative, or neutral outcome in landward/seaward mangrove response

dat2 <- dat %>% 
  filter(var %in% c('SeawardMang', 'LandwardMang') & 
           constraint_scenario %in% c('Macrotidal, High Hydro-connectivity',
                                      'Mesotidal, High Hydro-connectivity',
                                      'Microtidal, High Hydro-connectivity')) %>% 
          # pressure %in% c('Sea-level rise', 'Sea-level rise & Cyclones', 
           #               'Seal-level rise & Groundwater extraction', 'Sea-level rise  & Coastal development', 
            #               'Sea-level rise & Erosion', 'Sea-level rise & Drought or Dams')) %>% 
  group_by(model_scenario, constraint_scenario, pressure, var) %>% 
  summarise(Prob_gain_neutral = (sum(outcome>0) + sum(outcome==0))/n()) %>%
  mutate(constraint_scenario = recode(constraint_scenario, 'Macrotidal, High Hydro-connectivity' = 'Macrotidal',
         'Mesotidal, High Hydro-connectivity' = 'Mesotidal',
         'Microtidal, High Hydro-connectivity' = 'Microtidal'),
         var = recode(var, 'LandwardMang' = 'Landward mangrove', 'SeawardMang' = 'Seaward mangrove')) %>% 
  mutate(constraint_scenario = factor(constraint_scenario, levels = c('Microtidal', 'Mesotidal', 'Macrotidal')))

dat2$Prob_change <- rescale(dat2$Prob_gain_neutral, to = c(-100, 100))

# plot

a <- ggplot(dat2, aes(constraint_scenario, pressure, fill = Prob_change)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', 
                       name = 'Probability of Loss (red) or Neutrality/Gain (blue)', 
                       direction = 1) +
  facet_wrap(vars(factor(var), factor(model_scenario))) +
  #facet_wrap(~factor(model_scenario)) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.justification = 'left',
        #axis.text.y =  element_blank(),
        strip.text.x = element_text(size = 9),
        axis.title = element_blank()) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,))
a

ggsave('outputs/outcomes-heatmap_alt-model.png', width = 6.5, height = 8)
