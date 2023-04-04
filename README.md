### Qualitative analysis of mangrove response to drivers of change

TODO:: 
- Clean up code and make scenario simulation more efficient! purrrr....
-	Process data for missing for typologies
-	Maple/loop analyst, which loops drive ambiguity?
o	Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes

Step 1
- add in an alternative model structure for coastal squeeze - mangroves can't migrate landward with SLR - edge between SLR and landward mangroves removed
- remove density dependence from sea level rise
- see what co-authors think about SLR high sediment supply as seaward mangroves can either increase or decrease, but if we increase sediment under SLR, still could increase or decrease - but higher probability of gain
- alternative is to constrain paramaterisation of high sed supply even more, and decrease infludence of SLR - put to coauthors

Step 2: finish global predictions
- process spatial data
- Decide whether to get mangrove loss/gain to validate

Notes:
- coastal squeeze isn't adequately represented in the model - coastal development just directly reduces mangrove forest, doesn't limit its ability to migrate landward with SLR
- if we are assuming the mangroves are directly reduced by coastal development, should be assuming that mangroves can't move landward? 
- so we could do a scenario fo coastal develoment like roads, that doesn't stop landward migration,
- and one scenario for coastal squeeze, where SLR can't allow landward mangroves to migrate (thats just the coastal development scenario)
- So this is done - just explain difference between the scenarios
- except this won't work for spatial predictions. in that way the model isn't adequately capturing coastal squeeze. Unless we do a model structure where if there is coastal development, interaction between SLR and mangroves is 0.
- so for spatial predictions, we should have both scenarios, and identify sites where this changes predictions substantially. another source of ambiguity.

In the scripts folder you will find:

1. 01_scenarios-outcomes.R: a script to simulate scenarios and plot probability of predicted outcomes

2. 02_scenarios-plots.R: plot the scenario results

3. 03_mapping.R: map models to mangrove typologies

4. helpers.R: a script with helper functions

5. models.R: a script that builds different conceptual models for testing
  
