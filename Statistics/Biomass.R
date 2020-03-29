#analysis of biomass 
#stats to account for difference between CON and WAT Treatment
#Script:
  #1. Biomass second cut
    #data exploration
    #assumptions test
    #kruskal wallis test
    #posthoc tests
    #GLM
    #boxcox
  #2. Biomass first cut

setwd("C:/00 Dana/Uni/Internship/Work") #set working directory
library(agricolae)
library(PMCMRplus)
library(MASS)
library(PMCMR)
library(tidyverse)
#for biomass
biomass=read.table("Dried Biomass Data Collection.csv", sep=";", dec=".", header=T)
#delete empty cloumns off data.frame
biomass=Filter(function(x)!all(is.na(x)), biomass)
str(biomass)#check 

#data exploration
boxplot(biomass$Biomass.Cut.2~biomass$Variety)
boxplot(biomass$Biomass.Cut.2~biomass$Treatment)
#Anova Assumptions
#normality
# test for whole data
biomass_all=biomass$Biomass.Cut.2
qqnorm(biomass_all)
qqline(biomass_all)
shapiro.test(biomass_all) #p value = <2.2*10^-16
#data is not normally distributed
#try to transform the dependent variable
#log transformation
biomass_log=log(biomass_all)
qqnorm(biomass_log)
qqline(biomass_log)
shapiro.test(biomass_log) #1.210^-9

#square root transformation
biomass_sqrt=sqrt(biomass_all)
qqnorm(biomass_sqrt)
qqline(biomass_sqrt)
shapiro.test(biomass_sqrt) #p value 6.59*10^-7
#independence
#equality of variance
#as assumption of normality is not met use kruskal wallis test instead
biomass_non_param=cbind.data.frame(biomass$Variety, biomass$Treatment, biomass$Biomass.Cut.2)
colnames(biomass_non_param)=c("Variety", "Treatment", "Weight")
kruskal.test(Weight~Treatment, data=biomass_non_param) #p value 1.08*10^-12
kruskal.test(Weight~Variety, data=biomass_non_param) #p value 2.06*10^-12
interTV=interaction(biomass_non_param$Treatment, biomass_non_param$Variety)
kruskal.test(biomass_non_param$Weight~interTV) #p value 1.48*10^-12
#install.packages("PMCMRplus")
?posthoc.kruskal.nemenyi.test
#Calculate pairwise multiple comparisons between group levels. 
#These tests are sometimes referred to as Nemenyi-tests for multiple 
#comparisons of (mean) rank sums of independent samples
#--> not appropriate for unequal samples
attach(biomass_non_param)
posthoc.kruskal.nemenyi.test(x=Weight, g=Treatment, dist="Chisquare")
posthoc.kruskal.nemenyi.test(x=Weight, g=Variety, dist="Chisquare")
options(max.print=99999)
posthoc.kruskal.nemenyi.test(x=Weight, g=interTV, dist="Chisquare")

#Calculate pairwise multiple comparisons between group levels according to Conover
?posthoc.kruskal.dunn.test
posthoc.kruskal.dunn.test(x=Weight, g=Treatment, dist="Chisquare")
posthoc.kruskal.dunn.test(x=Weight, g=Variety, dist="Chisquare")
#same as bonferroni?
#Calculate pairwise multiple comparisons between group levels according to Dunn

#linear model with interaction
sum.lm=summary(lm(formula=biomass$Biomass.Cut.2~biomass$Treatment*biomass$Variety))
lm=lm(formula=biomass$Biomass.Cut.2~biomass$Treatment*biomass$Variety)
plot(lm)

#lm without interaction
summary(lm(formula=biomass$Biomass.Cut.2~biomass$Treatment+biomass$Variety))
lm.add=lm(formula=biomass$Biomass.Cut.2~biomass$Treatment+biomass$Variety)
plot(lm)

glm=glm(biomass$Biomass.Cut.2~biomass$Treatment*biomass$Variety)
plot(glm)

#try boxcox
lm.biom=lm(biomass$Biomass.Cut.2~biomass$Variety*biomass$Treatment)
plot(lm.biom)
bx.bio=boxcox(lm.biom)
lm.biom.bx=lm((biomass$Biomass.Cut.2*0.4)~biomass$Variety*biomass$Treatment)
biomass2_log=log(biomass$Biomass.Cut.2)
shapiro.test(biomass2_log)
qqnorm(biomass2_log)
qqline(biomass2_log)
#boxcox return 0.3 -> try log transformation

biomass2_test=biomass$Biomass.Cut.2^0.5
shapiro.test(biomass2_test)
qqnorm(biomass2_test)
qqline(biomass2_test)

biomass2_rezi=biomass$Biomass.Cut.2*-1
shapiro.test(biomass2_rezi)
#*****************************************************************************
#2. first biomass cut
#data exploration
boxplot(biomass$Biomass_Cut_1.~biomass$Variety)
boxplot(biomass$Biomass_Cut_1.~biomass$Treatment)
boxplot(biomass$Biomass_Cut_1.~biomass$type)
boxplot(biomass$Biomass_Cut_1.~biomass$sum_treatments)

#Anova Assumptions
#normality
# test for whole data
qqnorm(biomass$Biomass_Cut_1.)
qqline(biomass$Biomass_Cut_1.)
shapiro.test(biomass$Biomass_Cut_1.) #p value = <2.2*10^-16
#data is not normally distributed
#try to transform the dependent variable

#square root transformation
biomass_sqrt=sqrt(biomass$Biomass_Cut_1.)
qqnorm(biomass_sqrt)
qqline(biomass_sqrt)
shapiro.test(biomass_sqrt) #p value 4.1*10^-7
#log doesn't work because of zeros in data
#reciproce transformation
biomass_rezi=biomass$Biomass_Cut_1.*(-1)
qqnorm(biomass_rezi)
qqline(biomass_rezi)
shapiro.test(biomass_rezi) #p value: 2.2*10^-16

#boxcox transformation
lm.biom1=lm(biomass$Biomass_Cut_1.~biomass$Treatment)
plot(lm.biom1)
boxcox(lm.biom1)

#independence
#equality of variance
#as assumption of normality is not met use kruskal wallis test instead
kruskal.test(biomass$Biomass_Cut_1.~biomass$Variety) #p value 8.2*10^-15
kruskal.test(biomass$Biomass_Cut_1.~biomass$Treatment) #p value 1.81*10^-13
kruskal.test(biomass$Biomass_Cut_1.~biomass$type) #2.03*10^-6

wilcox.test(biomass$Biomass_Cut_1.~biomass$sum_treatments) #p value 3.07*10^-15

posthoc.kruskal.nemenyi.test(biomass$Biomass_Cut_1.~biomass$Variety, dist="Chisquare")
posthoc.kruskal.nemenyi.test(biomass$Biomass_Cut_1.~biomass$Treatment, dist="Chisquare")
posthoc.kruskal.nemenyi.test(biomass$Biomass_Cut_1.~biomass$type, dist="Chisquare")

#GLM
summary(glm(biomass$Biomass_Cut_1.~biomass$Variety*biomass$Treatment))