# apc_IE
Intrinsic Estimator for Age, Period, Cohort Analysis

This R routine estimates age, period, and cohort (apc) effects using the Intrinsic Estimator (IE) as described in Fu (2000) 

Fu, W.J., 2000. "The ridge estimator in singular design with application to age, period, cohort analysis of disease rates." _Communication in Statistics – Theory and Methods_ 29, 263–27, 

and in 

Y. Yang and K.C. Land (2013), _Age-Period-Cohort Analysis
New Models, Methods, and Empirical Applications_, Chapman & Hall.

The input data frame must contain variables named: Y, N, age, period, cohort, based on an age by period table of counts Y and exposures N. Variable names are checked and the data are passed to a function that constructs a centered-effects design matrix for the age, period, and cohort factors. Although the "last" factor level is commonly used as a referent in existing software (i.e., the Stata module apc_ie)--and is the default normalization here, this routine also allows the "first" level to be used as a reference (see, e.g., L. Luo, et al.'s (2016) commentary in the _American Journal of Sociology_). 

The resulting design matrix is passed to routines for the linear model (log rate), 
or nonlinear models (Poisson or Binomial). Results are then passed to a function that normalizes the
apc parameters so that the full set of results is output to a data frame. 

This approach is described in more detail in the attached paper, 

D.A Powers (2013) "Black–white differences in maternal age, maternal birth cohort, and period effects on infant mortality in the US (1983–2002)," _Social Science Research_, 1033–1045.

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
### Manipulating Results
The returned results are in a dataframe (unless the ```out="raw"``` option is used). This is handy for processing in ```ggplot```.
As an illustration we fit both last- and first-category normalized Poisson regresssions using data from the attached article.
```
#
# generate two sets of apc results (first and last category reference) 
#
# last category level as reference
#
p.modL <- IE_rate(dat, ref="last", family="Poisson")
p.modL[grep('age', rownames(p.modL)),] %>% 
    mutate(x = seq(15,45,5), apcfact=factor(rep("Age"))) -> agedfL
p.modL[grep('period', rownames(p.modL)),] %>%
    mutate(x= seq(1985,2000,5), apcfact=factor(rep("Period"))) -> perdfL
p.modL[grep('cohort', rownames(p.modL)),] %>%
    mutate(x = seq(40,85,5), apcfact=factor(rep("Cohort")))  -> cohdfL

apcdfL <- rbind(agedfL, perdfL, cohdfL)

#
# first category level as reference
#
p.modF <- IE_rate(dat, ref="first", family="Poisson")
p.modF[grep('age', rownames(p.modF)),] %>% 
  mutate(x = seq(15,45,5), apcfact=factor(rep("Age"))) -> agedfF
p.modF[grep('period', rownames(p.modF)),] %>%
  mutate(x= seq(1985,2000,5), apcfact=factor(rep("Period"))) -> perdfF
p.modF[grep('cohort', rownames(p.modF)),] %>%
  mutate(x = seq(40,85,5), apcfact=factor(rep("Cohort")))  -> cohdfF

apcdfF <- rbind(agedfF, perdfF, cohdfF)
```
Next we average these sets of results to get the ingredients for a ```ggplot```.
```
# average results
average.apc(apcdfF, apcdfL) -> apcdfAve
```
Then compare the plots (for sensitivity to reference category)
```
# make plots
apcdfL %>%
  ggplot(aes(y=estimate, x=x), color=apcfact) +
  facet_wrap(apcfact ~., scales="free_x") +
  geom_line(col="darkslategrey", alpha=6.0) +
  geom_ribbon(aes(ymin=lower, ymax=upper), 
              fill="darkorange", alpha=.6) +
  xlab("") + ylim(-.5, 1) +
  theme(strip.text = element_text(size=15)) -> pL

apcdfF %>%
  ggplot(aes(y=estimate, x=x), color=apcfact) +
  facet_wrap(apcfact ~., scales="free_x") +
  geom_line(col="darkslategrey", alpha=6.0) +
  geom_ribbon(aes(ymin=lower, ymax=upper), 
              fill="darkorange", alpha=.6) +
  xlab("") +  ylim(-.5, 1) +
  theme(strip.text = element_text(size=15)) -> pF

apcdfAve %>%
  ggplot(aes(y=y, x=x), color=apcfact) +
  facet_wrap(apcfact ~., scales="free_x") +
  geom_line(col="darkslategrey", alpha=6.0) +
  geom_ribbon(aes(ymin=lower, ymax=upper), 
              fill="darkorange", alpha=.6) +
  xlab("") + ylim(-.5, 1) +
  theme(strip.text = element_text(size=15)) -> pA

gridExtra::grid.arrange(pL, pF, pA, nrow=1, ncol=3 )

```
