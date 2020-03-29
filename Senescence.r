#Senescence significance test
#by Dana Looschelders
#the script performs the following analysis
  #data exploration (boxplot, histogram)
  #assumptions test for ANOVA
  #ANOVA (additive and interaction model)
  #LSD test as post hoc test (and scheffe test)
#use significance level of 0.05
#stats to account for difference between CON and WAT Treatment
setwd("C:/00 Dana/Uni/Internship/Work") #set working directory
#install.packages("agricolae")
library(agricolae)
#for Sene.data
Sene.data=read.table("leaf_phenology.csv", sep=";", dec=".", header=T)
#delete empty cloumns off data.frame
Sene.data=Filter(function(x)!all(is.na(x)), Sene.data)
str(Sene.data)#check 
#set all negative values to NA (NA means that the first leaf did not yet senesce)
Sene.data$Phyllochron.1.Senesced[Sene.data$Phyllochron.1.Senesced<0]=NA
str(Sene.data$Phyllochron.1.Senesced)
#create new data.frame
Sene.all=cbind.data.frame(Sene.data$Treatment, Sene.data$Variety, Sene.data$Phyllochron.1.Senesced)
str(Sene.all) #check
colnames(Sene.all)=c("Treatment", "Variety","Days") #set column names
#data exploration
stats=summary(Sene.data)
summary(Sene.all)
str(Sene.all)
#boxplot
boxplot(Days~Treatment, Sene.all)
#Anova Assumptions
#normality
Sene.data$Treatment=as.character(Sene.data$Treatment)
qqnorm(Sene.all$Days)
qqline(Sene.all$Days)
shapiro.test(Sene.all$Days) #p value = 0.069 -> normalty can be assumed as p>0.05
#independence -> can be assumed
#equality of variance 
  #use bartlett test as normality can be assumed
bartlett.test(Days~Treatment, Sene.all)
#p value =0.35 -> homogenity of variance can be assumed
#ANOVA 
aov.test=aov(Days~Treatment, Sene.all)
plot(aov.test)
summary(aov.test) #p value = 0.791 
#-> no significant difference bewteen groups
#-> if significant run least significant difference

#two way anova
#Assumptions
  #dependant variable continious
  #independent variables categorical
  #sample independence
  #normality
    #assumptions were tested for one way anova already
  #variance equality
  bartlett.test(Days ~ interaction(Treatment,Variety), data=Sene.all)
  #p value = 0.34 so assumption is met
two.aov.test=aov(Days~Treatment*Variety, Sene.all)
summary(two.aov.test)
  #Treatment: P value = 0.76
  #Variety: p value = 2.23*10^-5 --> highly significant
  #Interaction: p value = 0.37
 
#interaction is not significant -> use additive model
add.two.aov.test=aov(Days~Treatment+Variety, Sene.all)
summary(add.two.aov.test)
  #Treatment: p value = 0.76
  #Variety: p value = 2.2*10^-5 -> highly significant

#perform least significance difference
df.error=df.residual(add.two.aov.test)
ms.error=deviance(add.two.aov.test)/df.error
lsd.test=LSD.test(add.two.aov.test, "Variety", console=T, group=T, DFerror = df.error, MSerror = ms.error, p.adj = "holm")

#maybe use turkey test instead?
#behavioural ecology suggests scheffe post hoc as best practice
?scheffeTest
scheffeTest(Days~Treatment, Sene.all)
scheffeTest(Days~Variety, Sene.all)

