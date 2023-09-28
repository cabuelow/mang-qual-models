# classify typologies as low, medium or high according to percentile thresholds

library(tidyverse)

dat <- read.csv('outputs/climate-exposure_global.csv')

#TODO: here filter for Types that intersect with Madagascar land-eez before doing the percentile classification below

# set the percentile thresholds for classifying as low, medium or high
# anything between low and high will be medium - check with Jaramar this sounds okay, can easily change
low <- 0.3
high <- 0.7

# classify according to thresholds

dat <- dat %>% 
  mutate(future_slr_class = case_when(future_slr <= quantile(.$future_slr, low) ~ 'Low',
                                      future_slr >= quantile(.$future_slr, high) ~ 'High',
                                      .default = 'Medium'),
         antecedent_slr_class = case_when(antecedent_slr <= quantile(.$antecedent_slr, low) ~ 'Low',
                                          antecedent_slr >= quantile(.$antecedent_slr, high) ~ 'High',
                                          .default = 'Medium'),
         future_drought_class = case_when(future_drought <= quantile(.$future_drought, low) ~ 'Low',
                                          future_drought >= quantile(.$future_drought, high) ~ 'High',
                                          .default = 'Medium'),
         future_extreme_rainfall_class = case_when(future_extreme_rainfall == 0 ~ 'Low', # set 0 as low and exclude from percentile classification
                                                   future_extreme_rainfall <= quantile(.$future_extreme_rainfall, low) ~ 'Low',
                                                   future_extreme_rainfall >= quantile(.$future_extreme_rainfall, high) ~ 'High',
                                                   .default = 'Medium'),
         historical_drought_class = case_when(historical_drought <= quantile(.$historical_drought, low) ~ 'Low',
                                              historical_drought >= quantile(.$historical_drought, high) ~ 'High',
                                              .default = 'Medium'),
         historical_extreme_rainfall_class = case_when(historical_extreme_rainfall == 0 ~ 'Low', # set 0 as low and exclude from percentile classification
                                                       historical_extreme_rainfall <= quantile(.$historical_extreme_rainfall, low) ~ 'Low',
                                                       historical_extreme_rainfall >= quantile(.$historical_extreme_rainfall, high) ~ 'High',
                                                       .default = 'Medium'),
         historical_storms_class = case_when(historical_storms <= quantile(.$historical_storms, low) ~ 'Low',
                                             historical_storms >= quantile(.$historical_storms, high) ~ 'High',
                                             .default = 'Medium'),
         future_storms_class = case_when(future_storms <= quantile(.$future_storms, low) ~ 'Low',
                                         future_storms >= quantile(.$future_storms, high) ~ 'High',
                                         .default = 'Medium'))