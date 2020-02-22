# PLTMLE

This repository contains sample R code to execute the pooled LTMLE algorithm, developped by Petersen *et al.* (2014). This code is similar to the one implemented in the *ltmleMSM* function in the *ltmle* R package version 1-1.0 developped by Schwab *et al.* (2013). However, this code is specific to the discretization problem presented in the manuscript "Impact of discretization of the timeline for longitudinal causal inference methods". The provided code is for a two time-point example and is presented as raw steps to illustrate the corresponding pooled LTMLE algorithm.

The function PLTMLE2TP_CV serves to compute the cross-validated variance of the pooled LTMLE estimate. Note that one can obtain the estimate by simply running the first part of the code, i.e. on the training data only.

### References

Petersen M, Schwab J, Gruber S, BlaserN, Schomaker M, Laan v. dM. Targeted maximum likelihood estimation for dynamic
and static longitudinal marginal structural working models. Journal of causal inference 2014; 2(2): 147–185.

Schwab J, Lendle S, Petersen M, Laan v. dM. LTMLE: longitudinal targeted maximum likelihood estimation. 2013. R
package.
