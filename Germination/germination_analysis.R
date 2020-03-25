# Script to perform survival analysis on the 
# PHENOGRASS germination data
# 
# Data are preprocessed by the read_germination.R script
# and these are saved to file germination_data.RData
#
# 
# Jon Yearsley (24th March 2020)
# Jon.Yearsley@ucd.ie
# **************************************************

rm(list=ls())
setwd('./Germination')

library(survival)
library(survminer)
library(coxme)
library(broom)

# Load processed data ----
load('germination_data.RData')

# Calculate Kaplan-Meier survival curve ----
# This is a non-parametric approach, and doesn't 
# have the same assumptions as Cox PH model

# The hazard rate normally corresponds to the
# probability of dying. In our case it is the
# probability of germinating 


# Effect of Treatment for each Variety
germ_fit = survfit(Surv(Day, germinated) ~ Variety + Treatment, 
                             data=germination)

# Print estimates of time to germination
germ_fit



# Effect of treatment averaging over variety
germ_fit_treatment = survfit(Surv(Day, germinated) ~ Treatment, 
                             data=germination)

# Fit using Nelson-Aalen estimator (rather than Kaplan-Meier)
germ_fit_treatment_fh = survfit(Surv(Day, germinated) ~ Treatment, 
                             data=germination, type='fh')


# Display results of average time to germination across varieties
germ_fit_treatment

ggsurvplot(fit=germ_fit_treatment,
           xlab='Days',
           ylab='Prob not germinated',
           conf.int=TRUE,
           pval=TRUE)


# Log-ratio test of difference between survival 
# curves for control - treatment
m_logratio = survdiff(Surv(Day, germinated) ~ Treatment, 
         data=germination)
m_logratio


# Fit Cox proportional hazards models -----
# This is a parametric model that must be validated


# Fit the full model with main effects and interaction
coxph.fit <- coxph(Surv(Day, germinated) ~ Treatment+Variety, 
                   data=germination, 
                   method="breslow")  # Could use efron

# # Fit the full model with main effects and random term for Variety
coxph.rand <- coxph(Surv(Day, germinated) ~ Treatment + frailty(Variety),
                   data=germination)


# Validate hazard function assumption
test.ph <- cox.zph(coxph.fit)
test.ph

par(mfrow=c(1,2))
plot(test.ph)
par(mfrow=c(1,1))
# Looks like proportional hazard decreases over time for 
# the treatment effect

# Look at model with Variety as a random effect
test2.ph <- cox.zph(coxph.rand)
test2.ph

plot(test2.ph)
# Same message. Looks like proportional hazard decreases over time for 
# the treatment effect


# Need to find extension to Cox PH model




# Hypothesis tests --------

# Fit full model testing for Variety on control data
coxph.fit_full <- coxph(Surv(Day, germinated) ~ Variety, 
                         data=subset(germination, Treatment=='CON'), 
                         method="breslow")  # Could use efron

# Fit null model testing for Variety on control data
coxph.fit_null <- coxph(Surv(Day, germinated) ~ 1, 
                         data=subset(germination, Treatment=='CON'), 
                         method="breslow")  # Could use efron

anova(coxph.fit_full, coxph.fit_null, test='ChiSq')  

# Strong effect of Variety on germination

# Plot the hazard ratio for the model with main effect variety
ggforest(coxph.fit_full, data=germination)

survfit(coxph.fit_full)


# Time varying coefficients -----

# Modify Cox PH model -----
# Fit a Cox PH model but with several time bins

# Bin data into 3 bins t<14, 14<t<18, t>18
germin2 <- survSplit(Surv(Day, germinated) ~ ., 
                  data= germination, 
                  cut=c(14,18),
                  episode= "tgroup", id="id")

# Fit model
coxph.fit_full2 <- coxph(Surv(Day, germinated) ~ Treatment:strata(tgroup), 
                        data=germin2, 
                        method="breslow")  # Could use efron

# Look at mean fitted values
coxph.fit_full2$means

# # This isn't working
# coxph.fit_full3 <- coxph(Surv(Day, germinated) ~ Treatment + tt(Treatment), 
#                          data=germination, 
#                          method="breslow",
#                          tt = function(x, t, ...) x * (15*t/(15+t)))  # Could use efron





# Fit model with time varying coefficients for Treatment -----
library(timereg)

# Treat Variety as constant and allow Treatment to time vary
m = timecox(Surv(Day, germinated) ~ const(Variety)+Treatment, 
            data=germination, 
            max.time=25, 
            residuals=TRUE)


par(mfrow=c(1,2))
plot(m)
plot(m,score=TRUE)
par(mfrow=c(1,1))

summary(m)


