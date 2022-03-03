# Multivariate analysis of growth chamber results
#
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Feb 2022
# +++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd('~/git_repos/Phenograss/Expt/')

library(tidyr)
library(nlme)
library(ggplot2)


# Import data
d = read.csv('GrowthNewHarvest.csv')

# Rename Variables
nam = names(d)
nam[14] = "BiomassGR"
nam[15] = "HeightGR"
names(d) = nam

# Scale and centre response variables
d$BiomassGRsc = scale(d$BiomassGR, center=T, scale=T)
d$HeightGRsc = scale(d$HeightGR, center=T, scale=T)

# Visualise
ggplot(data=d,
       aes(x=BiomassGR,
           y=HeightGR)) +
  geom_point()

# +++++++++++++++++++++++++++++++++++++++++++++++++++
# Test data set ---------
path <- "https://content.sph.harvard.edu/fitzmaur/ala2e/cholesterol-data.txt"
dfW.data <- read.table(path, na.string = ".")
colnames(dfW.data) <- c("group","id","y0","y1","y2","y3","y4")
dfW.data$group <- factor(dfW.data$group,
                         levels = 1:2, labels = c("T","C"))
dfW.data$group <- relevel(dfW.data$group, ref = "C")
dfW.data$id <- as.factor(dfW.data$id)
str(dfW.data)
library(reshape2)
dfL.data <- melt(dfW.data, id.vars = c("group","id"),
                 value.name = "cholesterol", variable.name = "time")
dfL.data$time <- as.factor(gsub("y","visit",dfL.data$time))
dfL.data$time <- relevel(dfL.data$time, ref = "visit0")
dfL.data <- dfL.data[order(dfL.data$id,dfL.data$time),]

dfL.data$treatment <- as.character(dfL.data$group)
dfL.data[dfL.data$time == "visit0", "treatment"] <- "none"
dfL.data$treatment <- factor(dfL.data$treatment,
                             levels = c("none","C","T"),
                             labels = c("none","pl","tr"))


# Check... should give the same as melt
dfL.2 = pivot_longer(dfW.data, cols=c(3:7), names_to = "time", values_to="cholesterol")
dfL.2$time = as.factor(dfL.2$time)

summary(dfL.2)

#  Define correlation structures
e.glsCS <- gls(cholesterol~group*time,
               data = dfL.data,
               correlation = corCompSymm(form= ~1|id),
               na.action = na.omit)


summary(e.glsCS$modelStruct)
getVarCov(e.glsCS, individual=1)




# Unstructured correlation matrix
e.glsUN <- gls(model = cholesterol~group * time,
               data = dfL.data,
               correlation = corSymm(form =~as.numeric(time)|id),
               weights = varIdent(form =~1|group),
               na.action = na.omit)
summary(e.glsUN$modelStruct$corStruct)
summary(e.glsUN$modelStruct$varStruct)


Sigma <- getVarCov(e.glsUN, individual = 1)
Sigma
# +++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++


# Look at correlation between Biomass and Height growth rates
cov(d[,c(14,15)])

d_long =  pivot_longer(d, cols = c(14,15), 
                       names_to = "ResponseID", 
                       values_to = "Values" )

head(d_long)

d_long$ResponseID = as.factor(d_long$ResponseID)
d_long$ID = as.factor(d_long$ID)

d_long$heightID = d_long$ResponseID=="HeightGR"
d_long$biomassID = d_long$ResponseID=="BiomassGR"

table(d_long$ResponseID, d_long$ID)

m.gls <- gls(Values~1,
               data = subset(d_long, ID%in%c(1:10)),
               correlation = corCompSymm(form= ~1|ID),
               na.action = na.omit)

summary(m.gls)
getVarCov(m.gls)

# Broadly same model using lme
m.lme <- lme(Values~1,
             data = subset(d_long, ID%in%c(1:10)),
             random = ~1|ID,
             na.action = na.omit)

summary(m.lme)
getVarCov(m.lme, type="marginal")


m.gls <- gls(Values~1+ResponseID,
             data = subset(d_long, ID%in%c(1:10)),
             correlation = corSymm(form= ~1|ID),
             weights=varIdent(form=~1|ResponseID),
             na.action = na.omit)


# Fit a mixed model with a correlation structure
m = lme(Values~ResponseID-1, 
        random=list(ID=pdDiag(~ResponseID)), 
        data=subset(d_long, ID%in%c(1:100)))


m = lme(Values~ResponseID-1, 
        random=list(ID=pdBlocked(list(pdDiag(~heightID-1),pdDiag(~biomassID-1)))), 
        data=subset(d_long, ID%in%c(1:10)))


summary(m)
getVarCov(m, type="random.effects")
getVarCov(m, type="conditional")
getVarCov(m, type="marginal")
head(d)


# Look at latent variable model

library(lavaan)

myModel <- ' # regressions
             HeightGRsc + BiomassGRsc ~ Waterlogging + Conditions

             # # latent variable definitions 
             #   f1 =~ y1 + y2 + y3 
             #   f2 =~ y4 + y5 + y6 
             #   f3 =~ y7 + y8 + y9 + y10

             # variances and covariances 
               BiomassGRsc ~~ HeightGRsc 

             # intercepts 
               BiomassGRsc ~ 1 
               HeightGRsc ~ 1
           '

m = cfa(myModel, data=d)


summary(m, fit.measures=TRUE)
