#analysis of biomass 
#stats to account for difference between treatments
#Script:
  #1. Biomass second cut
    #data exploration
    #assumptions test
    #mann whitney test
    #posthoc tests
    #GLM
    #boxcox
  #2. Biomass first cut

setwd("C:/00_Dana/Uni/Internship/Work/Data Rosemount") #set working directory
library(agricolae)
library(PMCMRplus)
library(MASS)
library(PMCMR)
library(tidyverse)
#for biomass (read in excel data on biomass from google drive)
biomass=read.table("Biomass Data.csv", sep=";", dec=",", header=T)
str(biomass)
#delete empty cloumns off data.frame
biomass=Filter(function(x)!all(is.na(x)), biomass)
str(biomass)#check 
names(biomass)=c("Plant_ID", "Variety", "Treatment","Chamber", "Cut_1", "Cut_2")

#data exploration
boxplot(biomass$Cut_2~biomass$Variety)
boxplot(biomass$Cut_2~biomass$Treatment)
hist(biomass$Cut_2)
#t.test
#normality
# test for whole data
biomass_all=biomass$Cut_2
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
#as assumption of normality is not met use mann whitney test
wilcox.test(biomass$Cut_2 ~ biomass$Treatment) #p-value = 2.179e-14

#do kurskal wallis test for Varieties
kruskal.test(biomass$Cut_2~biomass$Variety) #p-value = 2.056e-12

#Calculate pairwise multiple comparisons between group levels. 
#These tests are sometimes referred to as Nemenyi-tests for multiple 
#comparisons of (mean) rank sums of independent samples
#--> not appropriate for unequal samples

posthoc.kruskal.nemenyi.test(x=biomass$Cut_2, g=biomass$Variety, dist="Chisquare")

#linear model with interaction
sum.lm=summary(lm(formula=biomass$Cut_2~biomass$Treatment*biomass$Variety))
lm=lm(formula=biomass$Cut_2~biomass$Treatment*biomass$Variety)
plot(lm)

#lm without interaction
summary(lm(formula=biomass$Cut_2~biomass$Treatment+biomass$Variety))
lm.add=lm(formula=biomass$Cut_2~biomass$Treatment+biomass$Variety)
plot(lm)

glm=glm(biomass$Cut_2~biomass$Treatment*biomass$Variety)
plot(glm)

#try boxcox
lm.biom=lm(biomass$Cut_2~biomass$Variety*biomass$Treatment)
plot(lm.biom)
bx.bio=boxcox(lm.biom)
lm.biom.bx=lm((biomass$Cut_2*0.4)~biomass$Variety*biomass$Treatment)
biomass2_log=log(biomass$Cut_2)
shapiro.test(biomass2_log)
qqnorm(biomass2_log)
qqline(biomass2_log)
#boxcox return 0.3 -> try log transformation

biomass2_test=biomass$Cut_2^0.5
shapiro.test(biomass2_test)
qqnorm(biomass2_test)
qqline(biomass2_test)

biomass2_rezi=biomass$Cut_2*-1
shapiro.test(biomass2_rezi)

#*****************************************************************************
#2. first biomass cut
#data exploration
boxplot(biomass$Cut_1~biomass$Variety)
boxplot(biomass$Cut_1~biomass$Treatment)
hist(biomass_Cut_1)
#Assumptions test
#normality
# test for whole data
qqnorm(biomass$Cut_1)
qqline(biomass$Cut_1)
shapiro.test(biomass$Cut_1) #p value = <2.2*10^-16
#data is not normally distributed
#try to transform the dependent variable

#square root transformation
biomass_sqrt=sqrt(biomass$Cut_1)
qqnorm(biomass_sqrt)
qqline(biomass_sqrt)
shapiro.test(biomass_sqrt) #p value 4.1*10^-7
#log doesn't work because of zeros in data
#reciproce transformation
biomass_rezi=biomass$Cut_1*(-1)
qqnorm(biomass_rezi)
qqline(biomass_rezi)
shapiro.test(biomass_rezi) #p value: 2.2*10^-16

#boxcox transformation
lm.biom1=lm(biomass$Cut_1~biomass$Treatment)
plot(lm.biom1)
boxcox(lm.biom1)

#independence
#equality of variance
#as assumption of normality is not met use wilcox test instead

wilcox.test(biomass$Cut_1~biomass$Treatment) #p-value = 9.145e-16

#GLM
summary(glm(biomass$Cut_1~biomass$Variety*biomass$Treatment))