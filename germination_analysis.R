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

# Load processed data ----
load('germination_data.RData')

# Calculate Kaplan-Meier survival curve ----

# Effect of treatment averaging over variety
germ_fit_treatment = survfit(Surv(Day, germinated) ~ Treatment, 
                             data=germination)

# Effect of variety for control treatment
germ_fit_CON = survfit(Surv(Day, germinated) ~ Variety, 
                       se.fit=TRUE,
                       data=subset(germination, Treatment=='CON'))

# Effect of variety for control treatment
germ_fit_TRT = survfit(Surv(Day, germinated) ~ Variety, 
                       se.fit=TRUE,
                       data=subset(germination, Treatment=='eCO2'))


# Plot survival curves -----
ggsurvplot(fit=germ_fit_treatment,
           xlab='Days',
           ylab='Prob not germinated',
           conf.int=TRUE)

ggsurvplot(fit=germ_fit_CON,
           xlab='Days',
           ylab='Prob not germinated',
           title='Effect of variety for Treatment=CON',
           conf.int=TRUE)

ggsurvplot(fit=germ_fit_TRT,
           xlab='Days',
           ylab='Prob not germinated',
           title='Effect of variety for Treatment=eCO2',
           conf.int=TRUE)





# Fit Cox proportional hazards models -----
# Note these models are not yet validated

# Fit the full model with main effects and interaction
coxph.fit <- coxph(Surv(Day, germinated) ~ Variety*Treatment, 
                   data=germination, 
                   method="breslow")  # Could use efron

# Fit null model testing interaction between Variety and Treatment
coxph.fit_null1 <- coxph(Surv(Day, germinated) ~ Variety+Treatment, 
                   data=germination, 
                   method="breslow")  # Could use efron

anova(coxph.fit, coxph.fit_null1, test='ChiSq')  # No evidence of an interaction

# Fit null model testing effect of Treatment
coxph.fit_null2 <- coxph(Surv(Day, germinated) ~ Variety, 
                         data=germination, 
                         method="breslow")  # Could use efron

anova(coxph.fit, coxph.fit_null2, test='ChiSq')  # Effect of treatment


# Fit null model testing effect of Variety
coxph.fit_null3 <- coxph(Surv(Day, germinated) ~ Treatment, 
                         data=germination, 
                         method="breslow")  # Could use efron

anova(coxph.fit, coxph.fit_null3, test='ChiSq')  # Effect of variety



