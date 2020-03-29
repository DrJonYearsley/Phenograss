#significance test for germination data
#by Dana Looschelders
#the script performs the following tasks
  #data exploration (boxplot, histogram)
  #assumption test
  #kruskal test 
  #posthoc test (nemenyi test)
  #fit glm with quasibinomial distribution

#for cumulative germination 
#to do: analysis of D50
library(agricolae)
library(PMCMR)
library(PMCMRplus)
library(tidyverse)

setwd("C:/00 Dana/Uni/Internship/Work")
germ=read.table("germination_4_stats.csv", sep=";", dec=".", header=T)
str(germ) #use fraction of germination
summary(germ)

#cumulative germination
boxplot(germ$Fraction~germ$Variety)
boxplot(germ$Fraction~germ$Treatment)
boxplot(germ$Fraction~germ$type)
boxplot(germ$Fraction~germ$sum_treatments)
#test for normal distribution
hist(germ$Fraction)
qqnorm(germ$Fraction)
qqline(germ$Fraction)
shapiro.test(germ$Fraction) 

#as values only can be between 0 and 1 it is a binomial distribution
#not normally distributed, therefore non-parametric test

kruskal.test(germ$Fraction~germ$Variety) #significant: p-value: 2.2*10^-16
kruskal.test(germ$Fraction~germ$Treatment) #not significant p-value: 0.82
kruskal.test(germ$Fraction~germ$type) #significant: p-value: 2.3*10^-9
kruskal.test(germ$Fraction~germ$sum_treatments) #not significant: p-value: 0.78

#posthoc tests
posthoc.kruskal.nemenyi.test(germ$Fraction~germ$Variety, dist = "Chisquare")
posthoc.kruskal.nemenyi.test(germ$Fraction~germ$type, dist="Chisquare")

#fir glm with binomial distribution
#for 
model_var.treat=glm(germ$Fraction~germ$Variety*germ$Treatment, family="quasibinomial")
summary(model_variety)

model_typ.comtreat=glm(germ$Fraction~germ$type*germ$sum_treatments, family="quasibinomial")
summary(model_typ.comtreat)

model_var.typ=glm(germ$Fraction~germ$Variety*germ$sum_treatments, family="quasibinomial")
summary(model_var.typ)

#try boxcox 
lm.germ=lm(germ$Fraction~germ$Variety*germ$Treatment)
plot(lm.germ)
bx=boxcox(lm.germ)
#stats analysis for D50 Germination

