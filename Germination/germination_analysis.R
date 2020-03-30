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
#setwd('./Germination')

library(survival)
library(survminer)
library(coxme)
library(broom)
library(ggplot2)
library(emmeans)
library(DHARMa)

# Load processed data ----
load('./Germination/germination_data.RData')


# Look at distribution of germination times
ggplot(data=subset(germination, Day<30), 
       aes(x=Day, fill=Treatment)) +
  geom_histogram(bins=20) +
  theme_bw() +
  scale_fill_brewer(palette ='Dark2') +
  labs(x='Day post sowing',
         y='Number Germinated')


# Test whether there's a diference in overall germination chance
germin_tab = table(germination$Treatment, 
                   germination$germinated)

chisq.test(germin_tab)

# Test for difference between species and treatments
m=glm(germinated~Treatment*Variety, data=germination, family='binomial')

summary(m)
# Residual deviance is about equal to degrees of freedom


# Validate the GLM using DHARMa
val = simulateResiduals(m, n=500, refit=FALSE)
plot(val)
testResiduals(val)
testOverdispersion(val)

# Logistic model look valid

m0 = update(m, .~.-Treatment:Variety)

anova(m0,m, test='Chisq')

# No evidence of an interaction for total germination proportion

# Test for effect of variety
m0b = update(m, .~.-Treatment:Variety - Variety)
anova(m0b,m0, test='Chisq')

# A strong effect of Variety

# Test for effect of Treatement
m0c = update(m, .~.-Treatment:Variety - Treatment)
anova(m0c,m0, test='Chisq')

# No effect of Treatment


# Posthoc test for the effect of Variety
m_eff = emmeans(m0, specs = 'Variety')

# Compare each Variety to the average
posthoc = contrast(m_eff, method='eff', type='contrast')

posthoc

plot(posthoc)

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

survplot = ggsurvplot(fit=germ_fit_treatment,
           xlab='Days',
           ylab='Prob not germinated',
           conf.int=TRUE,
           pval=TRUE,
           risk.table=FALSE)
survplot$plot = survplot$plot + 
  scale_color_brewer(palette = 'Dark2') + 
  scale_fill_brewer(palette = 'Dark2') + 
  # theme(axis.title = element_text(size=20),
  #       axis.text = element_text(size=20)) +
  theme_bw()
survplot

ggsave('KaplanMeier_curve.png')

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

tmp <- coxph(Surv(Day, germinated) ~ 1+Treatment, 
                        data=germination, 
                        method="breslow")  # Could use efron
r = residuals(tmp,'schoenfeld')

plot(germination$Day, r)

# Plot expected and predicted
plot(germ_fit_treatment)
newdat  <- data.frame(Treatment = levels(germination$Treatment))         # newdata to get expected curves
lines(survfit(tmp, newdata = newdat),
      col = "red", lty = 1:2)



# Look at a time varying model
coxph.tt <- coxph(Surv(Day, germinated) ~ Treatment + tt(Treatment),
                         data=germination,
                         method="breslow",
                         tt = function(x, t, ...) (x=='eCO2') * (1-t))  # Could use efron



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



est = CsmoothB(m$cum, c(7:20), b=2)

plot(est[,1],exp(est[,3]),pch=20,xlab='Hazard Rate', type='b')
plot(est[,1],exp(est[,2]),pch=20,xlab='Baseline Hazard Rate', type='b')


library(MRH)
m2 = estimateMRH(Surv(Day, germinated) ~ Variety + nph(Treatment), 
            data=germination, 
            M=2,
            maxStudyTime=20,
            maxIter=1000)
