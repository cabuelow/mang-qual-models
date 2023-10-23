# plot results of sensitivity analysis

library(tidyverse)
library(ggh4x)
#library(patchwork)

dat <- readRDS('outputs/simulation-outcomes/outcomes.rds') # sensitivity to structural model assumptions
dat_press <- read.csv('outputs/simulation-outcomes/outcomes_mangrove_model_spatial.csv') # sensitivity to pressure definition

##### sensitivity to structural model assumptions ######
# calculate proportion of stable model outcomes that have a positive, negative, or neutral landward/seaward mangrove response

dat2 <- dat %>% 
  filter(var %in% c('SeawardMang', 'LandwardMang') & 
           constraint_scenario %in% c('Macrotidal, High propagule establishment capacity, High coastal squeeze',
                                      'Microtidal, High propagule establishment capacity, High coastal squeeze',
                                      'Macrotidal, High propagule establishment capacity, Low coastal squeeze',
                                      'Microtidal, High propagule establishment capacity, Low coastal squeeze')) %>% #& 
  group_by(model, model_scenario, constraint_scenario, pressure, var) %>% 
  summarise(Prob_gain_neutral = ((sum(outcome>0) + sum(outcome==0))/n())*100,
            Prob_loss = (sum(outcome<0)/n())*-100) %>%  
  mutate(tide = strsplit(constraint_scenario, ', ')[[1]][1],
         coastalsqueeze = strsplit(constraint_scenario, ', ')[[1]][3],
         var = recode(var, 'LandwardMang' = 'Landward mangrove', 'SeawardMang' = 'Seaward mangrove')) %>% 
  pivot_wider(id_cols = c(model_scenario, constraint_scenario, pressure, var, tide, coastalsqueeze, var), names_from = 'model', values_from = 'Prob_gain_neutral') %>% 
  mutate(model_cyc_pos = abs(model_cyc_pos - mangrove_model),
         model_cyc_seaward = abs(model_cyc_seaward - mangrove_model),
         model_drought = abs(model_drought - mangrove_model),
         model_rain = abs(model_rain - mangrove_model)) %>% 
  pivot_longer(cols = c(model_cyc_pos:model_rain), names_to = 'model', values_to = 'absolute_deviation') %>% 
  mutate(model = recode(model, 'model_cyc_pos' = 'A) Storms -> substrate', 'model_cyc_seaward' = 'B) Storms --* seaward', 
                        'model_drought' = 'C) Drought -* Sediment', 'model_rain' = 'D) Rain -> Sediment')) %>% 
  mutate(#pressure = ifelse(pressure == 'Sea-level rise & Coastal development', paste0(pressure, ' (', coastalsqueeze, ')'), pressure),
         #pressure = ifelse(pressure == 'Sea-level rise & Intense storms & Coastal development', paste0(pressure, ' (', coastalsqueeze, ')'), pressure),
         pressure = ifelse(pressure == 'Sea-level rise & Groundwater extraction', 'Sea-level rise & Subsidence (Groundwater extraction)', pressure),
         pressure = ifelse(pressure == 'Groundwater extraction', 'Subsidence (Groundwater extraction)', pressure)) %>% 
  filter(!pressure %in% c('Sea-level rise & Intense storms & Coastal development'))

# plot

ggplot(filter(dat2, tide == 'Macrotidal' & model_scenario == 'High Sediment Supply' & coastalsqueeze == 'High coastal squeeze')) +
  geom_point(aes(x = absolute_deviation, y = pressure, col = var, alpha = 0.5)) +
  facet_nested(~model) +
  ylab('') +
  xlab('Percent change in probability of mangrove gain/neutrality from baseline model') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom') +
  guides(alpha = 'none')

ggsave('outputs/sensitivity/sensitivity_structural-model-assumptions.png', height = 3, width = 9)

##### sensitivity to pressure definition ######
# pressure definition = 4 is the baseline

dat_press2 <- dat_press %>% 
  filter(cast == 'forecast') %>% 
  mutate(Prob_gain_neutral = Prob_gain + Prob_neutral) %>% 
  pivot_wider(id_cols = c(var, Type), names_from = 'pressure_def', values_from = 'Prob_gain_neutral', names_prefix = 'press_') %>% 
  mutate(`1_4` = abs(press_1-press_4),
         `2_4` = abs(press_2-press_4),
         `3_4` = abs(press_3-press_4),
         `5_4` = abs(press_5-press_4)) %>% 
  pivot_longer(cols = c(`1_4`, `2_4`, `3_4`, `5_4`), names_to = 'sensitivity', values_to = 'change') %>% 
  mutate(var = recode(var, 'LandwardMang' = 'A) Landward mangrove', 'SeawardMang' = 'B) Seaward mangrove'))

ggplot(dat_press2) +
  geom_jitter(aes(x = sensitivity, y = change), width = 0.2, alpha = 0.1) +
  #geom_violin(aes(x = sensitivity, y = change)) +
  facet_wrap(~var) +
  ylab('% change from baseline') +
  xlab('') +
  theme_classic()

ggsave('outputs/sensitivity/sensitivity_pressure-definition.png', width = 3.5, height = 2.2)

