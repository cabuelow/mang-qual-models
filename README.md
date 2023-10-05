### Qualitative forecasts of mangrove extent change with network models

This repository provides code for reproducing results and figures in Buelow et al. (*in prep*), 'Forecasting mangrove futures under climate change'.

The codebase draws heavily from the R package [{QPress}](https://github.com/SWotherspoon/QPress).

Link to documents describing all spatial data processing:

-   [Part 1 here](https://mangrove-climate-risk-mapping.netlify.app/)
-   [Part 2 here](https://mangrove-climate-risk-mapping-2.netlify.app/)

## TODO

- [ ] Clean up all scripts, especially helpers

#### Scripts

1.  01_wrangle-dat.R: wrangles processed data into a master dataframe

2.  02_model-scenarios.R: a script to simulate scenarios and plot probability of predicted outcomes

3.  03_plot-scenarios.R: plot the scenario results

4.  04_spatial-model-hindcast-validation.R: fit and cross-validate spatial hindcasts

5.  05_model-spatial-hindcast-uncertainty.R: quantify uncertainty in cross-validation accuracy estimators

6.  06_map-hindcast-validation.R: map the hindcasts

7.  07_spatial-model-forecast.R: make calibrated forecasts and map

8.  08_plot-sensitivity.R: plot sensitivity analyses

9.  helpers/models.R: a script that builds different models

10.  helpers/helpers.R: a script with helper functions for simulating models

11. helpers/spatial-helpers.R: a script with helper functions for simulating models spatially

12. misc-plotting.R: miscellaneous plotting and mapping

