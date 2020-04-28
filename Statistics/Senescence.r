#Senescence significance test
#by Dana Looschelders
#the script performs the following analysis
  #data exploration (boxplot, histogram)
  #assumptions test for ttest
  #ttest 
  #wilcoxon test and kruskal test (with posthoc test)
#use significance level of 0.05
#stats to account for difference between Treatments
setwd("C:/00_Dana/Uni/Internship/Work/Data Rosemount/") #set working directory

library(agricolae)
library(PMCMR)
#for Sene.data
Sene=read.table("Time to Leaf Senescence from Germination.csv", sep=";", dec=".", header=T)
#delete empty cloumns off data.frame
Sene=Filter(function(x)!all(is.na(x)), Sene)
str(Sene)#check 

#data exploration
summary(Sene)

#boxplot
hist(Sene$Senescence)
boxplot(Sene$Senescence~Sene$Treatment)
boxplot(Sene$Senescence~Sene$Variety)

#ttest Assumptions
#normality
qqnorm(Sene$Senescence)
qqline(Sene$Senescence)
shapiro.test(Sene$Senescence) #p value = 1.505e-09 -> normalty cannot be assumed as p>0.05
#independence -> can be assumed
#equality of variance 
  #use wilcox test as normality cannot be assumed
wilcox.test(Sene$Senescence~Sene$Treatment) #p-value = 0.07853 -> not significant

#for Variety
kruskal.test(Sene$Senescence~Sene$Variety) #p-value = 0.0003811
#post hoc test
posthoc.kruskal.nemenyi.test(Sene$Senescence~Sene$Variety)

