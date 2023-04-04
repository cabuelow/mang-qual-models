### Qualitative analysis of mangrove response to drivers of change

TODO:: 
- Clean up code and make scenario simulation more efficient! purrrr....
-	Process data for missing for typologies
-	Maple/loop analyst, which loops drive ambiguity?
o	Although Ward (2021) says symbolic analysis is powerful for identifying feedbacks, and is difficult to carry out for models with more than five nodes

- add in extreme rainfall - what do do - include as variables in model or directly perturb?
- remove density dependence from sea level rise
- make SLR + coastal development more realistic - mangroves will not migrate landward, so no chance of increase/retreat: 'coastal squeeze' phenomenon in the model, i.e., retreat is not possible with coastal development.  
- If we are certain that there are some scenarios where seaward mangroves will prograde, we can parameterise the model to capture this. This is common practice with network models, where interaction strengths are parameterised so that they meet some condition of validity, e.g., seaward mangroves increase
- CL: I would do it - high sediment supply locations that get colonised by mangroves - identify locations for validation
- CL: Management scenarios are nice idea. Compare directly with Schuerch et al. maps? Those are +/- coastal development? That might be interesting to see where they don’t agree. But I'm not sure we would know why.

In the scripts folder you will find:

1. 01_scenarios-outcomes.R: a script to simulate scenarios and plot probability of predicted outcomes

2. 02_scenarios-plots.R: plot the scenario results

3. 03_mapping.R: map models to mangrove typologies

4. helpers.R: a script with helper functions

5. models.R: a script that builds different conceptual models for testing
  
