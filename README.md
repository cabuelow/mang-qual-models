### Probabilistic forecasts of the direction of future change in mangrove extent using network models

This repository provides code for reproducing results and figures in Buelow et al. (*in prep*), 'Mangrove persistence under climate change doubles with management and restoration'.

The codebase draws heavily from the R package [{QPress}](https://github.com/SWotherspoon/QPress).

Link to webpages where you can interactively view all spatial data underpinning the analysis, and read the steps taken to process the spatial data:

-   [Part 1 here](https://mangrove-climate-risk-mapping.netlify.app/)
-   [Part 2 here](https://mangrove-climate-risk-mapping-2.netlify.app/)

## TODO

- [ ] Clean up all scripts, especially helpers

#### Scripts

1.  01_wrangle-dat.R: wrangles processed data into a master dataframe

2.  02_model-scenarios.R: a script to simulate scenarios and plot probability of predicted outcomes

3.  03_plot-scenarios.R: plot the scenario results

4.  04_spatial-model-hindcast-validation.R: fit and cross-validate spatial hindcasts

5.  05_spatial-model-hindcast-uncertainty.R: quantify uncertainty in cross-validation accuracy estimators

6.  06_map-hindcast-validation.R: map the hindcasts

7.  07_spatial-model-forecast.R: make forecasts and map

8.  08_spatial-model-forecast-uncertainty.R: quantify 95% confidence intervals around forecast classes

9.  09_plot-sensitivity.R: plot sensitivity analyses

10.  helpers/models.R: a script that builds different models

11.  helpers/helpers.R: a script with helper functions for simulating models

12. helpers/spatial-helpers.R: a script with helper functions for simulating models spatially

13. misc-plotting.R: miscellaneous plotting and mapping

