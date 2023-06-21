### Qualitative hindcasts and forecasts of mangrove extent change

This repository provides code for reproducing results and figures in Buelow et al. (*in prep*), 'Resolving ambiguous mangrove futures'.

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

7.  helpers/models.R: a script that builds different conceptual models

8.  helpers/helpers.R: a script with helper functions for simulating responses

##### TODO

-   [X] Decide whether to map all typological units, or just those with extant forest (as of 2020)
        - Will use all units (extant or only historical loss) as areas with historical loss could be restored
-   [ ] What happens to maps when we include local factors, e.g., erosion, hydrodynamic energy, etc
-   [ ] Maple/loop analyst, which loops drive ambiguity? o Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes
