### Qualitative analysis of mangrove response to drivers of change

TODO:: 
- Clean up code and make scenario simulation more efficient! purrrr....
-	Process data for missing for typologies
-	Maple/loop analyst, which loops drive ambiguity?
o	Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes

- 1. Finish spatial processing
- 1.a Decide whether to get mangrove loss/gain to validate

- 2. Ask co-authors what they think about new results - final logic check
- Do we need a 'coastal-squeeze' paramaterisation - if coastal development = 1 (is perturbed),
then SLR to Landward Mang interaction is 0, or close to 0 - or let strength of SLR->LandwardMang interaction vary by pop.size?
- also futher explain what is happening with high vs. low sed supply (clarify that interaction strenght between sediment variable and sub volume = supply) - supply rate can be modified by hydrological connectivity
- Sediment variable = amount of sediment available in lower catchment, can be modified by drought, dams, precipitation
- model is saying that under slr and high sed supply mangroves have 50/50 chance of increasing (prograding) or decreasing
- if we have an increase in sed to the system (i.e., via extreme rainfall event), model says mangroves have 100% chance of increasing

Notes:
- coastal squeeze isn't adequately represented in the model - coastal development just directly reduces mangrove forest, doesn't limit its ability to migrate landward with SLR
- if we are assuming the mangroves are directly reduced by coastal development, should be assuming that mangroves can't move landward? 
- so we could do a scenario fo coastal develoment like roads, that doesn't stop landward migration,
- and one scenario for coastal squeeze, where SLR can't allow landward mangroves to migrate (thats just the coastal development scenario)
- So this is done - just explain difference between the scenarios
- except this won't work for spatial predictions. in that way the model isn't adequately capturing coastal squeeze. Unless we do a model structure where if there is coastal development, interaction between SLR and mangroves is 0.
- so for spatial predictions, we should have both scenarios, and identify sites where this changes predictions substantially. another source of ambiguity.
- removed density dependence from sea level rise - matrix isn't sign stable. so will keep

In the scripts folder you will find:

1. 01_scenarios-outcomes.R: a script to simulate scenarios and plot probability of predicted outcomes

2. 02_scenarios-plots.R: plot the scenario results

3. 03_mapping.R: map models to mangrove typologies

4. helpers.R: a script with helper functions

5. models.R: a script that builds different conceptual models for testing
  
