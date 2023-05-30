### Qualitative analysis of mangrove response to drivers of change

Link to documents describing all spatial data processing:

-   [Part 1 here](https://mangrove-climate-risk-mapping.netlify.app/)
-   [Part 2 here](https://mangrove-climate-risk-mapping-2.netlify.app/)

#### Scripts

1.  01_wrangle-dat.R: wrangles processed data into a master dataframe

2.  02_model-scenarios.R: a script to simulate scenarios and plot probability of predicted outcomes

3.  03_plot-scenarios.R: plot the scenario results

4.  04_model-sptial.R: map models to mangrove typologies

5.  05_model-spatial-validation.R: hindcast and validate spatial model

6.  helpers.R: a script with helper functions

7.  models.R: a script that builds different conceptual models for testing

##### TODO

-   Clean up code and make scenario simulation more efficient! purrrr....
-   Maple/loop analyst, which loops drive ambiguity? o Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes

##### Development notes

-   high vs. low sed supply (clarify that interaction strenght between sediment variable and sub volume = supply) - supply rate can be modified by hydrological connectivity
-   Sediment variable = amount of sediment available in lower catchment, can be modified by drought, dams, precipitation
-   model is saying that under slr and high sed supply mangroves have 50/50 chance of increasing (prograding) or decreasing
-   if we have an increase in sed to the system (i.e., via extreme rainfall event), model says mangroves have 100% chance of increasing
-   removed density dependence from sea level rise - matrix isn't sign stable. so will keep
