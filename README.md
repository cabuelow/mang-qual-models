### Qualitative analysis of mangrove response to drivers of change

TODO:: 
- Clean up code and make scenario simulation more efficient! purrrr....
-	Process data for missing for typologies
-	Maple/loop analyst, which loops drive ambiguity?
o	Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes

Step1: Fix the base model:
- coastal squeeze isn't adequately represented in the model - coastal development just directly reduces mangrove forest, doesn't limit its ability to migrate landward with SLR
- if we are assuming the mangroves are directly reduced by coastal development, should be assuming that mangroves can't move landward? 
- so we could do a scenario fo coastal develoment like roads, that doesn't stop landward migration,
- and one scenario for coastal squeeze, where SLR can't allow landward mangroves to migrate (thats just the coastal development scenario)
- So this is done - just explain difference between the scenarios
- except this won't work for spatial predictions. in that way the model isn't adequately capturing coastal squeeze. Unless we do a model structure where if there is coastal development, interaction between SLR and mangroves is 0.
- so for spatial predictions, we should have both scenarios, and identify sites where this changes predictions substantially. another source of ambiguity.

- add in an alternative model structure for coastal squeeze - mangroves can't migrate landward with SLR - edge between SLR and landward mangroves removed
- add in extreme rainfall - what do do - include as variables in model or directly perturb?
- remove density dependence from sea level rise

- alternative model paramaterisation where added sediment means mangroves prograde seaward
- but actually, I think model results make sense. Seaward mangrove may increase or decrease with SLR and high sediment supply
- add in scenario of increased sediment - see what happens

- If we are certain that there are some scenarios where seaward mangroves will prograde, we can parameterise the model to capture this. This is common practice with network models, where interaction strengths are parameterised so that they meet some condition of validity, e.g., seaward mangroves increase
- CL: I would do it - high sediment supply locations that get colonised by mangroves - identify locations for validation

- CL: Management scenarios are nice idea. Compare directly with Schuerch et al. maps? Those are +/- coastal development? That might be interesting to see where they don’t agree. But I'm not sure we would know why.

Step 2: finish global predictions
- process spatial data
- Decide whether to get mangrove loss/gain to validate

In the scripts folder you will find:

1. 01_scenarios-outcomes.R: a script to simulate scenarios and plot probability of predicted outcomes

2. 02_scenarios-plots.R: plot the scenario results

3. 03_mapping.R: map models to mangrove typologies

4. helpers.R: a script with helper functions

5. models.R: a script that builds different conceptual models for testing
  
