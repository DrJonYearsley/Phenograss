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

setwd("C:/00_Dana/Uni/Internship/Work/Data Rosemount/")
germ=read.table("Raw Data Germination Trial.csv", sep=";", dec=".", header=T, skip=1)
str(germ) #use fraction of germination
#delete empty cloumns off data.frame
germ=Filter(function(x)!all(is.na(x)), germ)
str(germ)#check 

summary(germ)

germ$sum=rowSums(germ[,5:34])
germ$total=rep(NA)
germ$Fraction=rep(NA)
#no of seeds per pot
#5 seeds per pot: Moy, Semi-natural, wild
germ$total[germ$Variety=="Moy"|
              germ$Variety=="Semi-natural6"|
              germ$Variety=="Semi-natural7"|
              germ$Variety=="Semi-natural11"|
              germ$Variety=="Wild4"|
              germ$Variety=="Wild6"|
              germ$Variety=="Wild7"]=5

#10 seeds per pot: cultivars
germ$total[germ$Variety=="Lilora"|
                germ$Variety=="Solomon"|
                germ$Variety=="Carraig"|
                germ$Variety=="Dunluce"|
                germ$Variety=="Aberchoice"|
                germ$Variety=="Abergain"|
                germ$Variety=="Aspect"]=10
germ$Fraction=germ$sum/germ$total

#QAQC
any(germ$Fraction[germ$fraction>1|germ$Fraction<0]) #none smaller than 0 or bigger than 1

#data exploration 
boxplot(germ$Fraction~germ$Variety)
boxplot(germ$Fraction~germ$Treatment)



#test for normal distribution
hist(germ$Fraction) 
qqnorm(germ$Fraction)
qqline(germ$Fraction)
shapiro.test(germ$Fraction) 

#as values only can be between 0 and 1 it is a binomial distribution
#not normally distributed, therefore non-parametric test

kruskal.test(germ$Fraction~germ$Variety) #significant: p-value: < 2.2e-16
kruskal.test(germ$Fraction~germ$Treatment) #not significant p-value: 0.7829

#posthoc tests
posthoc.kruskal.nemenyi.test(germ$Fraction~germ$Variety, dist = "Chisquare")

#fir glm with binomial distribution
#for 
model_var.treat=glm(germ$Fraction~germ$Variety*germ$Treatment, family="quasibinomial")
summary(model_variety)


model_var.typ=glm(germ$Fraction~germ$Variety*germ$Treatment, family="quasibinomial")
summary(model_var.typ)

#try boxcox 
lm.germ=lm(germ$Fraction~germ$Variety*germ$Treatment)
plot(lm.germ)
bx=boxcox(lm.germ)

#stats analysis for D50 Germination

plot(glm(germ$Fraction~germ$Treatment*germ$Variety))

