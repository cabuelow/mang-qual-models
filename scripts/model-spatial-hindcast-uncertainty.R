# quantify uncertainty in model hindcast accuracy
# bootstrap resample data (ie random sampling with replacement) into 5 folds 40 times (which is 5*40 = 200 test-training pairs, recommended as stable by Lyons et al. 2018)
# and obtain sampling distribution of each accuracy estimate (overall, producers, users)
# can use median and +/- 0.025 and 0.975 percentiles to estimate 95% confidence intervals around accuracy estimate
# propagate uncertainty to estimates of number of units of loss, gain/neutrality, or ambiguity
# by taking the 95th percentiles of the commission and omission errors for each class, e.g.,
# num_gain_units(upper) = num_gain_units + (num_gain_units * omission_error(95 percentile)) # omission error is the rate of false negatives, so predictions of a class we are missing, so is the upper
# num_gain_units(lower) = num_gain_units - (num_gain_units * commission_error(95 percentile)) # commission error is the rate of false positives, so predictions of a class we shouldn't have, so is the lower