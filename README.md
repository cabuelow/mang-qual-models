### Qualitative analysis of mangrove response to drivers of change

Link to documents describing all spatial data processing:

-   [Part 1 here](https://mangrove-climate-risk-mapping.netlify.app/)
-   [Part 2 here](https://mangrove-climate-risk-mapping-2.netlify.app/)

#### Scripts

1.  01_wrangle-dat.R: wrangles processed data into a master dataframe

2.  02_model-scenarios.R: a script to simulate scenarios and plot probability of predicted outcomes

3.  03_plot-scenarios.R: plot the scenario results

4.  04_model-spatial.R: map models to mangrove typologies

5.  05_model-spatial-validation.R: hindcast and validate spatial model

6.  models.R: a script that builds different conceptual models

7.  helpers.R: a script with helper functions for simulating responses

##### TODO

-   Maple/loop analyst, which loops drive ambiguity? o Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes

##### Development notes

-   removed density dependence from sea level rise - matrix isn't sign stable. so will keep
