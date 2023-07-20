### Qualitative network analysis of mangrove extent change

This repository provides code for reproducing results and figures in Buelow et al. (*in prep*), 'Forecasting mangrove futures under climate change'.

Link to documents describing all spatial data processing:

-   [Part 1 here](https://mangrove-climate-risk-mapping.netlify.app/)
-   [Part 2 here](https://mangrove-climate-risk-mapping-2.netlify.app/)

#### Scripts

1.  01_wrangle-dat.R: wrangles processed data into a master dataframe

2.  02_model-scenarios.R: a script to simulate scenarios and plot probability of predicted outcomes

3.  03_plot-scenarios.R: plot the scenario results

4.  04_model-spatial.R: map forecasts and hindcasts to mangrove typologies

5.  05_spatial-validation.R: compare hindcast to historical mangrove loss and gain to validate

6.  06_map-spatial.R: map spatial outcomes

7.  07_plot-sensitivity.R: plot sensitivity analyses

8.  helpers/models.R: a script that builds different models

9.  helpers/helpers.R: a script with helper functions for simulating responses

##### TODO

 - Validation
  1. Remove locations where commodities & erosion are major drivers of loss (agriculture/aquaculture) - we can't validate that
  2. Grid search - ambiguity prob. threshold vs. pressure threshold - develop some criteria for picking the optimal on this grid, e.g., maximizing producer's accuracy - then try forecasts
  3. Model calibration/fitting - use the subset of the 10000 models that match observed responses. Then have the range of parameter values for each interaction strength in the model to make forecasts
    - Challenging part is need to match the range of interaction strength parameter ranges associated with a given set of stressors combinations to what we see in forecasts
  4. Global k-fold x-validation with random removal - e.g., 5-fold
  
 
-   [X] Decide whether to map all typological units, or just those with extant forest (as of 2020)
        - Will use all units (extant or only historical loss) as areas with historical loss could be restored