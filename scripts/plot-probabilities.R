library(tidyverse)
library(scales)
library(patchwork)
library(ggh4x)

naive_outcomes <- read.csv(paste0('outputs/validation/naive_outcomes.csv'))

out <- naive_outcomes %>% 
  mutate_at(vars(LandwardMang, SeawardMang), ~ifelse(. == -1, 0, 1)) %>% 
  group_by(scenario,press) %>% 
  summarise_at(vars(LandwardMang, SeawardMang), ~(sum(.)/1000)*100) %>% 
  separate_wider_delim(scenario, delim = '.', names = c('CoastalDev', 'Sediment', 'Tide', 'Propagule', 'Climate', 'CoastalDev2')) %>% 
  filter(Propagule == 'Propestab_H' & Sediment == 'Sedsupp_H', Tide != 'TidalClass_M')

# turn into a table of probabilities for each scenario

a <- ggplot(out, 
            aes(Tide, press, fill = SeawardMang)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', 
                       name = 'Probability of Loss (red) or Gain/Neutrality (blue)',
                       direction = 1,
                       breaks = c(0, 25, 50, 75, 100),
                       labels = c("-100", "-75", "50", '75', '100')) +
  #facet_nested(~factor(model_scenario)) +
  facet_nested(~factor(CoastalDev) + factor(Climate) + factor(CoastalDev2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.justification = 'left',
        strip.text.x = element_text(size = 9),
        title = element_text(size = 8),
        axis.title = element_blank()) +
  ggtitle('A) Seaward mangrove') +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
a

b <- ggplot(out, 
            aes(Tide, press, fill = LandwardMang)) +
  geom_tile(color = 'black') +
  scale_fill_distiller(palette = 'Spectral', 
                       name = 'Probability of Loss (red) or Gain/Neutrality (blue)',
                       direction = 1,
                       breaks = c(0, 25, 50, 75, 100),
                       labels = c("-100", "-75", "50", '75', '100')) +
  #facet_nested(~factor(model_scenario)) +
  facet_nested(~factor(CoastalDev) + factor(Climate)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.justification = 'left',
        strip.text.x = element_text(size = 9),
        title = element_text(size = 8),
        axis.title = element_blank()) +
  ggtitle('B) Landward mangrove') +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
b
