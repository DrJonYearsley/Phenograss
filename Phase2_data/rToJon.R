#setwd("C:/Users/FRIC2/Desktop/DublinStuff/Research/pData-New/Modelling")

library(emmeans)
library(car)
library(nlme)
library(tidyr)
library(ggplot2)

setwd('./Expt')
datC <- read.delim("FULLchambergrowth.csv",sep=",",header=TRUE)

datMCM1 <- subset(datC, DataSet=="MCM")
datMCM <- subset(datMCM1, Month!="May")
datCF1 <- subset(datC, DataSet=="CF")
datCF <- subset(datCF1, Month!="May")
str(datCF)


# nam <- names(datCF)
# names(datCF) <- nam
datCF$BiomassGrowthRate <- datCF$BiomassGrowthRate*1000
# datCF$BiomassSc <- scale(datCF$BiomassGrowthRate)
# datCF$HeightSc <- scale(datCF$HeightGrowthRate)
datCF$WaterStatus <- as.factor(datCF$WaterStatus)
datCF$Chamber <- as.factor(datCF$Chamber)
datCF$ClimateTreatment <- as.factor(datCF$ClimateTreatment)
datCF$Month <- as.factor(datCF$Month)
datCF$Ploidy <- as.factor(datCF$Ploidy)
datCF$Variety <- as.factor(datCF$Variety)


datCF_long <- pivot_longer(datCF, cols= c('BiomassGrowthRate','HeightGrowthRate'),
			names_to = "names", values_to = "response")

datCF_long$whichresponse=1
datCF_long$whichresponse[datCF_long$names=="BiomassGrowthRate"] = 2
datCF_long$whichresponse = as.factor(datCF_long$whichresponse)
datCF_long$ID = as.factor(datCF_long$ID)	# ID only within month
datCF_long$X = as.factor(datCF_long$X)	# Could be ID between months
datCF_long$Harvest = as.factor(datCF_long$Month)
head(datCF_long)

datCF_harvestJune = subset(datCF_long, Harvest%in%c("June"))

m <- lme(	response ~ -1 + whichresponse, 
          random = ~ -1 + whichresponse|Chamber,
		weights=varIdent(form=~1|whichresponse), 
		corr=corSymm(form=~as.numeric(whichresponse)|Chamber/ID),
		data=datCF_harvestJune, 
		method="ML",
		control = list(maxIter=500, msMaxIter=500, tolerance=1e-8, niterEM=250))

summary(m)
m_out = intervals(m)
VarCorr(m)
m_out$reStruct

m$modelStruct$varStruct
exp(coef(m$modelStruct$varStruct))

m_out$sigma

summary(m$modelStruct$corStruct)

####

# datCF_harvestJune = subset(datCF_long, Harvest=="June")
# datCF_harvestJuly = subset(datCF_long, Harvest=="July")
# datCF_harvestJuneJuly = subset(datCF_long, Harvest=="June" | Harvest=="July")

start.time <- Sys.time()
m2 <- lme(	response ~ -1 + whichresponse:((Month+WaterStatus+ClimateTreatment+Variety)^2), 
           random = ~ -1 + whichresponse|Chamber,
           weights=varIdent(form=~1|whichresponse), 
           corr=corSymm(form=~as.numeric(whichresponse)|Chamber/X),
           data=datCF_long, 
           method="ML",
           control = list(maxIter=500, msMaxIter=500, niterEM=250, sing.tol=1e-20, msMaxEval = 500))
end.time <- Sys.time()
time.taken <- end.time-start.time
time.taken # 12.92045 secs

# ++++++++++++++++++++++++++++++++++++++++++
# Look at different Anova approaches

Anova(m2, type=2)

drop1(m2, scope=~whichresponse:Month:WaterStatus, test="Chisq")


Anova(m2, type=3)
summary(m2)





e2 <- emmeans(m2, ~ WaterStatus * Month, nesting=NULL)
pwpp(e2, method="pairwise", sort=FALSE, adjust="sidak")
pwpm(e2, adjust="sidak")

e3 <- emmeans(m2, ~ WaterStatus * ClimateTreatment, nesting=NULL)
pwpp(e3, method="pairwise", sort=FALSE, adjust="sidak")
pwpm(e3, adjust="sidak")

e4 <- emmeans(m2, ~ ClimateTreatment * Month, nesting=NULL)
pwpp(e4, method="pairwise", sort=FALSE, adjust="sidak")
pwpm(e4, adjust="sidak")

e5 <- emmeans(m2, ~ Variety * Month, nesting=NULL)
pwpp(e5, method="pairwise", sort=FALSE, adjust="sidak")
pwpm(e5, adjust="sidak")

e6 <- emmeans(m2, ~ WaterStatus * Variety, nesting=NULL)
pwpp(e6, method="pairwise", sort=FALSE, adjust="sidak")
pwpm(e6, adjust="sidak")

e7 <- emmeans(m2, ~ ClimateTreatment * Variety, nesting=NULL)
pwpp(e7, method="pairwise", sort=FALSE, adjust="sidak")
pwpm(e7, adjust="sidak")

#### 


datCF_long2 <- pivot_longer(datCF, cols= c('Tiller','BiomassGrowthRate','HeightGrowthRate'),
			names_to = "names", values_to = "response")

datCF_long2$whichresponse=1
datCF_long2$whichresponse[datCF_long2$names=="BiomassGrowthRate"] = 2
datCF_long2$whichresponse[datCF_long2$names=="Tiller"] = 3
datCF_long2$whichresponse = as.factor(datCF_long2$whichresponse)
datCF_long2$ID = as.factor(datCF_long2$ID)	# ID only within month
datCF_long2$X = as.factor(datCF_long2$X)	# Could be ID between months
datCF_long2$Harvest = as.factor(datCF_long2$Month)
head(datCF_long2)


start.time <- Sys.time()
m3 <- lme(	response ~ -1 + whichresponse:((Month+WaterStatus+ClimateTreatment+Variety)^2), random = ~ -1 + whichresponse|Chamber,
		weights=varIdent(form=~1|whichresponse), 
		corr=corSymm(form=~as.numeric(whichresponse)|Chamber/X),
		data=datCF_long2, method="ML",
		control = list(maxIter=500, msMaxIter=500, tolerance=1e-8, niterEM=250))
end.time <- Sys.time()
time.taken <- end.time-start.time
time.taken # 4.525498 mins

# Error in lme.formula(response ~ -1 + whichresponse:((Month + WaterStatus +  : 
#  nlminb problem, convergence error code = 1
#  message = singular convergence (7)


start.time <- Sys.time()
m3 <- lme(	response ~ -1 + whichresponse:(Month+WaterStatus+ClimateTreatment+Variety), 
           random = ~ -1 + whichresponse|Chamber,
           weights=varIdent(form=~1|whichresponse), 
           corr=corSymm(form=~as.numeric(whichresponse)|Chamber/X),
           data=datCF_long2, method="ML",
           control = list(maxIter=500, msMaxIter=500, niterEM=250, sing.tol=1e-20, msMaxEval = 500))
end.time <- Sys.time()
time.taken <- end.time-start.time
time.taken # 2.217183 mins

## tolerance=1e-8 ## dunno, should be the same as sing.tol, but it isn't
## msVerbose=TRUE ## lets you see each iteration

Anova(m3, type=3)
summary(m3)

m3_out = intervals(m3)
m3_out$reStruct


e1 <- emmeans(m3, ~ WaterStatus)
pwpp(e1, method="pairwise", sort=FALSE, adjust="sidak", by="WaterStatus %in% whichresponse")
pwpp(e1, method="pairwise", sort=FALSE, adjust="sidak", by="whichresponse")
pwpm(e1, adjust="sidak")

str(e1)


