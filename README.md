### Qualitative forecasts of mangrove extent change

This repository provides code for reproducing results and figures in Buelow et al. (*in prep*), 'Forecasting mangrove futures under climate change'.

The codebase draws heavily from the excellent R package [{QPress}](https://github.com/SWotherspoon/QPress).

Link to documents describing all spatial data processing:

-   [Part 1 here](https://mangrove-climate-risk-mapping.netlify.app/)
-   [Part 2 here](https://mangrove-climate-risk-mapping-2.netlify.app/)

##### TODO

-   [ ] Consider removing slr and subsidence and other longterm threats from validation

#### Scripts

1.  01_wrangle-dat.R: wrangles processed data into a master dataframe

2.  02_model-scenarios.R: a script to simulate scenarios and plot probability of predicted outcomes

3.  03_plot-scenarios.R: plot the scenario results

4.  04_model-spatial.R: map forecasts and hindcasts to mangrove typologies under different pressure definition thresholds

5.  05_spatial-accuracy.R: compare hindcasts to historical observations of mangrove loss and gain across pressure and ambiguity thresholds and quantify accuracy

6.  06_spatial-crossvalidation.R: 5-fold crossvalidation training models to historical observations of mangrove change and testing on unseen data

7.  07_map-spatial.R: map spatial hindcasts and forecasts

8.  07_plot-sensitivity.R: plot sensitivity analyses

9.  helpers/models.R: a script that builds different models

10. helpers/helpers.R: a script with helper functions for simulating responses

