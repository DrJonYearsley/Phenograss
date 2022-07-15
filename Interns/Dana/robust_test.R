# Test out robust lmm and quantile regression
#
# Jon Yearsley
# Sept 2020
# **********************************************

rm(list=ls())

library(lme4)
library(robustlmm)
library(performance)
library(lqmm)


nGroups = 10
sdGroups = 3        # Variation between groups (variaion within groups=1)
nDataGroup = 100      # Number of data points per group
slope = 10

# Create expected values
yIntercept = 2+rep((rnorm(nGroups, sd=sdGroups)), each=nDataGroup)
x = rep(runif(nDataGroup), each=nGroups, min=0, max=1)
groups = rep(letters[c(1:nGroups)], each=nDataGroup)

# Draw y from different distributions
# y has a normal error distribution
y1 = abs(yIntercept + slope*x + rnorm(length(yIntercept)))

# y has a gamma distribution 
# (expectation = scale * shape, variance = shape*scale^2, skew = 2/sqrt(shape))
shape = 1
y2 = array(NA, dim=nGroups*nDataGroup)
for (i in 1:length(y2)) {
  y2[i] =  rgamma(1,shape=shape, scale=abs(yIntercept[i] + slope*x[i]))
}


d = data.frame(y1, y2, x, as.factor(groups))

# Fit some mixed models
m1 = lmer(y2~1+x+(1|groups), data=d)
m2 = glmer(y2~1+x+(1|groups), data=d, 
           family=Gamma(link="identity"))
m3 = rlmer(y2~1+x+(1|groups), data=d)
m4 = lqmm(y2~1+x, random=~1, group=groups, data=d, tau=0.5)

check_model(m1, check=c('qq','homogeneity','reqq'))
check_model(m2, check=c('qq','homogeneity','reqq'))

summary(m1)
summary(m2)
summary(m3)
summary(m4)

ranef(m4)
VarCorr(m4)
