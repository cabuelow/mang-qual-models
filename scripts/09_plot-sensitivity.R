# plot results of sensitivity analysis

library(tidyverse)
library(ggh4x)
library(sf)
library(tmap)
library(RColorBrewer)
library(QPress)
source('scripts/helpers/models.R')
model_names <- names(models) # which model do you want to run?
sf_use_s2(FALSE)
press <- 4
thresh <- 75

typ_points <- st_read('data/typologies/Mangrove_Typology_v3.14_Composite_valid_centroids.gpkg')
world <- data("World")
# sensitivity to model assumptions - non-spatial
dat <- readRDS('outputs/simulation-outcomes/outcomes.rds') 
# sensitivity to model assumptions - spatial
spatial_pred_fit <- lapply(model_names, function(i){read.csv(paste0('outputs/predictions/forecast-predictions', press, '_', thresh, '_', 'SeaLevelRise_', i, '_fit.csv'))})
names(spatial_pred_fit) <- model_names

##### non-spatial ######
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
  xlab('Percent change in probability of mangrove gain/stability from baseline model') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom') +
  guides(alpha = 'none')

ggsave('outputs/sensitivity/sensitivity_structural-model-assumptions.png', height = 3, width = 9)

##### spatial ######

spatial_dat <- spatial_pred_fit %>% 
  map(~ select(.x, Type, Landward, Seaward)) %>% 
  imap(~ rename(.x, "{.y}_Landward" := Landward)) %>% 
  imap(~ rename(.x, "{.y}_Seaward" := Seaward)) %>% 
  reduce(left_join) %>% 
  mutate(model_cyc_pos_Landward = ifelse(mangrove_model_Landward != model_cyc_pos_Landward, 'Change', 'No change'),
         model_cyc_pos_Seaward = ifelse(mangrove_model_Seaward != model_cyc_pos_Seaward, 'Change', 'No change'),
         model_rain_Landward = ifelse(mangrove_model_Landward != model_rain_Landward, 'Change', 'No change'),
         model_rain_Seaward = ifelse(mangrove_model_Seaward != model_rain_Seaward, 'Change', 'No change'),
         model_drought_Landward = ifelse(mangrove_model_Landward != model_drought_Landward, 'Change', 'No change'),
         model_drought_Seaward = ifelse(mangrove_model_Seaward != model_drought_Seaward, 'Change', 'No change'),
         model_cyc_seaward_Landward = ifelse(mangrove_model_Landward != model_cyc_seaward_Landward, 'Change', 'No change'),
         model_cyc_seaward_Seaward = ifelse(mangrove_model_Seaward != model_cyc_seaward_Seaward, 'Change', 'No change'))

preds_fit <- typ_points %>% 
  left_join(spatial_dat) %>%
  st_crop(xmin = -180, ymin = -40, xmax = 180, ymax = 33)
world_mang <- st_crop(World, xmin = -180, ymin = -40, xmax = 180, ymax = 33)  

lmap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(preds_fit, model_cyc_pos_Landward == 'Change')) +
  tm_symbols('model_cyc_pos_Landward', size = 0.06, border.lwd = 0, col = 'grey10', legend.shape.show = F) +
  tm_shape(filter(preds_fit, model_cyc_seaward_Landward == 'Change')) +
  tm_symbols('model_cyc_seaward_Landward', size = 0.04, border.lwd = 0, col = 'cyan4', legend.shape.show = F) +
  tm_shape(filter(preds_fit, model_drought_Landward == 'Change')) +
  tm_symbols('model_drought_Landward', size = 0.02, border.lwd = 0, col = 'orange', legend.shape.show = F) +
  tm_shape(filter(preds_fit, model_rain_Landward == 'Change')) +
  tm_symbols('model_rain_Landward', size = 0.005, border.lwd = 0, col = 'darkblue', legend.shape.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.07),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.35,
            main.title = 'B) Change in landward forecast from baseline model',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = c('grey10', 'cyan4', 'orange', 'darkblue'),
                labels =  c('A) Storms -> substrate', 'B) Storms --* seaward', 'C) Drought -* Sediment', 'D) Rain -> Sediment'), border.alpha = 0, size = 0.25)
lmap
tmap_save(lmap, paste0('outputs/maps/landward-forecast_map_', press, '_', thresh, '_all-data', '_', chosen_model_name, '_sensitivity.png'), width = 5, height = 1, dpi = 5000)

smap <- tm_shape(world_mang) +
  tm_fill(col = 'gray88') +
  tm_shape(filter(preds_fit, model_cyc_pos_Seaward == 'Change')) +
  tm_symbols('model_cyc_pos_Seaward', size = 0.06, border.lwd = 0, col = 'grey10', legend.shape.show = F) +
  tm_shape(filter(preds_fit, model_cyc_seaward_Seaward == 'Change')) +
  tm_symbols('model_cyc_seaward_Seaward', size = 0.04, border.lwd = 0, col = 'cyan4', legend.shape.show = F) +
  tm_shape(filter(preds_fit, model_drought_Seaward == 'Change')) +
  tm_symbols('model_drought_Seaward', size = 0.02, border.lwd = 0, col = 'orange', legend.shape.show = F) +
  tm_shape(filter(preds_fit, model_rain_Seaward == 'Change')) +
  tm_symbols('model_rain_Seaward', size = 0.005, border.lwd = 0, col = 'darkblue', legend.shape.show = F) +
  tm_layout(legend.outside = F,
            legend.position = c(0.13, 0.07),
            title.position = c(0.01,0.45),
            legend.title.size = 0.35,
            legend.text.size = 0.35,
            main.title = 'A) Change in seaward forecast from baseline model',
            main.title.size = 0.4,
            frame = T,
            legend.bg.color = 'white',
            legend.bg.alpha = 0) +
  tm_add_legend('symbol', col = c('grey10', 'cyan4', 'orange', 'darkblue'),
                labels =  c('A) Storms -> substrate', 'B) Storms --* seaward', 'C) Drought -* Sediment', 'D) Rain -> Sediment'), border.alpha = 0, size = 0.25)
smap
tmap_save(smap, paste0('outputs/maps/seaward-forecast_map_', press, '_', thresh, '_all-data', '_', chosen_model_name, '_sensitivity.png'), width = 5, height = 1, dpi = 5000)

# summarise percent of units with change in forecasted outcome

summary <- spatial_dat %>% 
  select(model_cyc_pos_Landward:model_cyc_seaward_Seaward) %>% 
  mutate(across(starts_with('model'), ~ifelse(. == 'Change', 1, 0))) %>% 
  mutate(model_cyc_pos = model_cyc_pos_Landward + model_cyc_pos_Seaward,
         model_rain = model_rain_Landward + model_rain_Seaward,
         model_drought = model_drought_Landward + model_drought_Seaward,
         model_cyc_seaward = model_cyc_seaward_Landward + model_cyc_seaward_Seaward) %>% 
  select(model_cyc_pos:model_cyc_seaward) %>% 
  mutate(across(starts_with('model'), ~ifelse(. > 0, 1, 0))) %>% 
  summarise(across(starts_with('model'), ~(sum(., na.rm=T)/n()))*100)
summary

write.csv(summary, 'outputs/summary-stats/spatial-sensitivity.csv', row.names = F)
