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

# Set data types
d$Waterlogging = as.factor(d$Waterlogging)
d$Chamber = as.factor(d$Chamber)


# Visualise
ggplot(data=d,
       aes(x=BiomassGR,
           y=HeightGR)) +
  geom_point()



# +++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++

# Manova approach -------

library(MASS)
library(emmeans)


# +++++++++++++++++++++++++++++++++
# Base functions: user defined contrasts

# Create contrast matrix
con_inv = matrix(c(1/4,1/4,1/4,1/4,
               1/2,1/2,-1/2,-1/2,
               1,-1,0,0,
               0,0,1,-1), 
             ncol=4, 
             byrow=TRUE)

con = ginv(con)



# Specify contrast matrix in the lm command 
# (remove first column because it is the intercept)
m1 = lm(BiomassGRsc~Chamber+Waterlogging, data=d, contrasts=list(Chamber=con[,-1]))
summary(m1)


# +++++++++++++++++++++++++++++++++
# Use emmeans to produce contrasts

# Fit a linear model with the default contrasts (treatment contrasts)
m = lm(BiomassGRsc~Chamber+Waterlogging, data=d)
summary(m)

m_eff = emmeans(m, spec="Chamber")
m_contrast1 = contrast(m_eff, list('Chamber 1 & 2 - Chamber 3 & 4'=c(0.5,0.5,-0.5,-0.5),
                                   'Chamber 1 - 2'=c(1,-1,0,0),
                                   'Chamber 3 - 4'=c(0,0,1,-1)))

summary(m_contrast1)











# +++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++
# Create multivariate model using a mixed modelling approach ---------

library(nlme)

# Look at correlation between Biomass and Height growth rates
cov(d[,c('BiomassGR','HeightGR')])
cor(d[,c('BiomassGR','HeightGR')])


d_long =  pivot_longer(d, cols = c('BiomassGR','HeightGR'), 
                       names_to = "names", 
                       values_to = "response" )

d_long$whichresponse = 1
d_long$whichresponse[d_long$names=="BiomassGR"] = 2
d_long$whichresponse = as.factor(d_long$whichresponse)
d_long$ID = as.factor(d_long$ID)
d_long$Harvest = as.factor(d_long$Harvest)
head(d_long)



# Create a mixed model that accounts for correlations between the two response variables

# Model inspired by book https://books.google.ie/books?id=N1BQvcomDdQC&lpg=PP1&pg=PA282&redir_esc=y#v=onepage&q&f=false
# Code inspired by https://www.stats.ox.ac.uk/~snijders/

# This page was also helpful
# http://staff.pubhealth.ku.dk/~jufo/courses/rm2018/nlmePackage.pdf

# Having "- 1" as part of a formula drops the intercept.
# For weights=varIdent, see Pinheiro & Bates (2000), page 208-209.
# For corr=corSymm, see Pinheiro & Bates (2000), page 234-235.

# Look at one harvest (harvest 2)
d_harvest2 = subset(d_long, Harvest==2)

m <- lme(response ~ - 1 + whichresponse, 
         random = ~ -1 + whichresponse|Chamber, 
         weights=varIdent(form=~1|whichresponse),
         corr=corSymm(form=~as.numeric(whichresponse)|Chamber/ID),
         data=d_harvest2, method="ML",
         control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                        niterEM=250))
summary(m)

# Model output with 95% confidence intervals
m_out = intervals(m)

# The covariance matrix at the Chamber level and the residual variance for
VarCorr(m)
m_out$reStruct

# Variance at Chamber level (between chambers) for HeightGR = 0.00291
# Variance at Chamber level  (between chambers) for BiomassGR = 129.25
# Covariance between HeightGR and BiomassGR at Chamber level = 0.437   (correlation=0.712)


# Within chamber variances... 
# Residual variance weighting factors are given by 
m$modelStruct$varStruct
# These numbers give the proportionality constants between the standard deviations.
# The internal coefficients are the logarithms, and
exp(coef(m$modelStruct$varStruct))

# Variance within chambers for BiomassGR (residual variance) = 6821.68 
m_out$sigma^2

# Variance within chambers for HeightGR (Residual variance) =  6821.68 * 0.003829643^2 = 0.1000479

# Correlation between HeightGR and BiomassGR within chambers = 0.6553833
m_out$corStruct

# Covariance between HeightGR and BiomassGR within chambers = 0.6553833 * sqrt(6821.68 * 0.1000479) = 17.12



# Can also extract this from
getVarCov(m, type="random.effects")




# Observed correlation between HeightGR and BiomassGR is combination of within and between correlations
# Total covariance divided by square-root of total variances for HeightGR and BiomassGR
# (17.12 + 0.437) / sqrt((0.1000479+0.00291)*( 6821.68+129.25)) = 0.6562

(17.12 + 0.437) / sqrt((0.1000479+0.00291)*( 6821.68+129.25))





# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Add in some fixed effects 
m2 <- lme(response ~ - 1 + whichresponse + whichresponse:(Waterlogging+Conditions), 
          random = ~ -1 + whichresponse|Chamber, 
          weights=varIdent(form=~1|whichresponse),
          corr=corSymm(form=~as.numeric(whichresponse)|Chamber/ID),
          data=d_harvest2, method="ML",
          control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                         niterEM=250))
