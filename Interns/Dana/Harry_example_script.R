#DW modes
library(car)
library(lattice)
library(lme4)
library(influence.ME)
library(robustlmm)
library(lmerTest)

setwd("C:/Users/chgio/OneDrive/Work folder/R/IRC")
raw<-read.csv("data.csv", header = TRUE)

#Select columns and omit rows with empty cells
sub_raw_DW<-subset(raw, select=c(Treatment, Chamber,Variety_abbr,Weight_first_harvest))
sub_raw_DW<-na.omit(sub_raw_DW)  

#Assign Chambers and Varieties abbreviations as factors, also rename to "Variety"
sub_raw_DW$Chamber <- as.factor(sub_raw_DW$Chamber)
sub_raw_DW$Variety <- as.factor(sub_raw_DW$Variety_abbr)

#compare treatments
boxplot(Weight_first_harvest~Treatment, data=sub_raw_DW)

#compare chambers
boxplot(Weight_first_harvest~Chamber, data=sub_raw_DW)

#Mixed effects model Model 1 (you will get the error: fixed-effect model matrix is rank deficient so dropping 1 column / coefficient. This happens because chambers are nested in treatments and cannot be avoided, I think)
modelDW1 <- lmer(Weight_first_harvest ~ Treatment+Chamber+(1|Variety) ,data = sub_raw_DW, REML=T)   #REML=F when I want to do model selection between models differing in their fixed effects
summary(modelDW1)

#Check normality and heteroscedacity of residuals
plot(modelDW1)  #if there is not heteroscedacity points shoud be randomly distributed - looks bad
qqmath(modelDW1, id=0.05)    # looks bad
shapiro.test(resid(modelDW1))  #test residuals for normal distribution - not normal
leveneTest(residuals(modelDW1) ~ sub_raw_DW$Treatment)   #test for equal variances - not homogenus

#Check outliers
inf=influence(modelDW1,obs=T)
plot(inf,which="cook")  #no outliers however is depends on the limit that you will set

#New robust mixed effects Model
modelDW2 <- rlmer(Weight_first_harvest ~ Treatment+Chamber+(1|Variety),data = sub_raw_DW)
summary(modelDW2)

#get p values by using t-values from the robust model and the Satterthwaite approximated degrees of freedom 
#from the equivalent regular mixed effect model. You can site this for this method: 
#Geniole SN, Proietti V, Bird BM, Ortiz TL, Bonin PL, Goldfarb B, Watson N V., Carré JM. 2019. Testosterone reduces the threat premium in competitive resource division. Proceedings of the Royal Society B: Biological Sciences 286.
coefs <- data.frame(coef(summary(modelDW1)))
coefs.robust <- coef(summary(modelDW2))
p.values <- 2*pt(abs(coefs.robust[,3]), coefs$df, lower=FALSE)
p.values

#Get random effects intercepts
ranef(modelDW2) 
