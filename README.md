# apc_IE
Intrinsic Estimator for Age, Period, Cohort Analysis

This R routine estimates age, period, and cohort effects using the Intrinsic Estimator as described in Fu (2000) [Fu, W.J., 2000. "The ridge estimator in singular design with application to age, period, cohort analysis of disease rates." Communication in Statistics – Theory and Methods 29, 263–27], and Y. Yang and K.C. Land (2013), [Age-Period-Cohort Analysis
New Models, Methods, and Empirical Applications, Chapman & Hall.]

The input data frame must contain variables named: Y, N, age, period, cohort, based on an age by period table of counts Y and exposures N. Variable names are checked and the data are passed to a function that constructs a centered-effects design matrix for the age, period, and cohort factors. Although the "last" factor level is commonly used as a referent in existing software (i.e., the Stata module apc_ie)--and is the default normalization here, this routine also allows the "first" level to be used as a reference (see, e.g., L. Luo, et al.'s commentary in the American Journal of Sociology). 

The resulting design matrix is passed to routines for the linear model (lnear log rate), 
or nonlinear models (Poisson or Binomial) regression. Results are then passed to a function that normalizes the
apc parameters so that the full set of results is output to a data frame.

### Example:.
Read data and source the function script.
```
dat <- read.csv(file='APCimrAP.csv')
source("apc.functions.R")
```
Fit a linear (log rate) model to get starting values.
```
linmod <- IE_linear(dat, ref="last", out="raw") 
b <- linmod$estimate   # use for start values
```
Here we have requested the last category normalization of the apc factors (default), 
the option ```out="raw"```
returns a list of estimates that can be passed as start values as shown below.
```
poismod <- IE_rate(dat, family="Poisson", bstart=b) 
print(poismod, digits=3)
```
The ```family``` argument defaults to "Poisson"
Similarly, we can request a logit model using ``family="Binomial"
```
logitmod <- IE_rate(dat, family="Binomial", bstart=b)
print(logitmod, digits=3)
```
