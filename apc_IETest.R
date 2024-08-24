
rm(list=ls())
library(MASS)
library(tidyverse)
library(tidyverse)
dat <- read.csv(file='APCimrAP.csv')
source("apc.functions.R")

linmod <- IE_linear(dat, ref="last", out="raw") 
b <- linmod$estimate
# use for start values

poismod <- IE_rate(dat, family="Poisson", bstart=b) 
print(poismod, digits=3)

logitmod <- IE_rate(dat, family="Binomial", bstart=b)
print(logitmod, digits=3)

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

# average results
average.apc(apcdfF, apcdfL) -> apcdfAve

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


