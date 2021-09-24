# 
# Tests for improvement of the initial model with environmental covariates
#
# Claire-Marie Alla
# 09/08/2021
# ++++++++++++++++++++++++++++++++++++++++++++++


rm(list=ls())

library(sf)
library(stars)
library(segmented)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)
library(gstat)
library(dplyr)
library(rgdal)
library(emmeans)
library(corrplot)


dataDir = '~/Stage/Data/MODIS/Phenophase_estimates'


# Load data
filename1 = '/gls_improve_model_square1_to_21_in_2017.Rdata'
filename2 = '/gls_improve_model_square20_in_2013_to_2017.Rdata'

load(file.path(dataDir,filename1))


# Fit model for phenology dates of phase 1 -----

# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time

phenology_wide$dummy = interaction(phenology_wide$year, phenology_wide$square, sep='_')
# x_ITM_centre + y_ITM_centre + elevation + slope + class_aspect + CATEGORY + soilmoisture + cumul_temp + cumul_precip0 + cumul_precip1 + cumul_precip5

# If filename 1 add in gls : factor(year)
gls_phase1 = gls(t_phase1~1 + (x_ITM_centre + y_ITM_centre + elevation + slope + class_aspect + CATEGORY + soilmoisture + cumul_temp + cumul_precip0 + cumul_precip1 + cumul_precip5),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre| dummy,
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)


summary(gls_phase1)

# Correlation matrix
mcor = summary(gls_phase1)$corBeta

anova(gls_phase1, type='marginal')

acf(gls_phase1$residuals, lag.max=20, type='partial')  # Makes sense... spatial correlation roughly 1-2 km

# # If filename1
# m_eff = emmeans(gls_phase1, spec='year')
# sos = as.data.frame(summary(m_eff))
# contrast(m_eff, infer=c(T,T))
# 
# as.Date(paste(round(sos$emmean),sos$year),format="%j %Y")
# 
# # WOrk out date from day of year
# as.Date(paste('52',sos$year),format="%j %Y")


# Keep result and residuals
phenology_wide$phase1_residual = residuals(gls_phase1)
phenology_wide$phase1_fit = predict(gls_phase1)


# Calculate some R-squared between predictions and data -----
m = lm(t_phase1~factor(square) + phase1_fit,
       data=phenology_wide)

# # If filename2
#m = lm(t_phase1~phase1_fit, data=phenology_wide)

print(summary(m))

print('moy residuals :')
print(mean(abs(phenology_wide$phase1_residual),na.rm=T))





# For correlation matrix
# colnames(mcor)[colnames(mcor) == 'class_aspectN'] <- 'N'
# colnames(mcor)[colnames(mcor) == 'class_aspectNE'] <- 'NE'
# colnames(mcor)[colnames(mcor) == 'class_aspectW'] <- 'W'
# colnames(mcor)[colnames(mcor) == 'class_aspectSE'] <- 'SE'
# colnames(mcor)[colnames(mcor) == 'class_aspectSW'] <- 'SW'
# colnames(mcor)[colnames(mcor) == 'class_aspectNW'] <- 'NW'
# colnames(mcor)[colnames(mcor) == 'class_aspectW'] <- 'W'
# colnames(mcor)[colnames(mcor) == 'class_aspectE'] <- 'E'
# colnames(mcor)[colnames(mcor) == 'class_aspectS'] <- 'S'
# colnames(mcor)[colnames(mcor) == 'CATEGORYMade'] <- 'Made'
# colnames(mcor)[colnames(mcor) == 'CATEGORYPeat'] <- 'Peat'
# colnames(mcor)[colnames(mcor) == 'CATEGORYWater'] <- 'Water'
# colnames(mcor)[colnames(mcor) == 'CATEGORYPoorly Drained'] <- 'PoorlyDrained'
# colnames(mcor)[colnames(mcor) == 'CATEGORYWell Drained'] <- 'WellDrained'
# 
# rownames(mcor)[rownames(mcor) == 'class_aspectN'] <- 'N'
# rownames(mcor)[rownames(mcor) == 'class_aspectNE'] <- 'NE'
# rownames(mcor)[rownames(mcor) == 'class_aspectW'] <- 'W'
# rownames(mcor)[rownames(mcor) == 'class_aspectSE'] <- 'SE'
# rownames(mcor)[rownames(mcor) == 'class_aspectSW'] <- 'SW'
# rownames(mcor)[rownames(mcor) == 'class_aspectNW'] <- 'NW'
# rownames(mcor)[rownames(mcor) == 'class_aspectW'] <- 'W'
# rownames(mcor)[rownames(mcor) == 'class_aspectE'] <- 'E'
# rownames(mcor)[rownames(mcor) == 'class_aspectS'] <- 'S'
# rownames(mcor)[rownames(mcor) == 'CATEGORYMade'] <- 'Made'
# rownames(mcor)[rownames(mcor) == 'CATEGORYPeat'] <- 'Peat'
# rownames(mcor)[rownames(mcor) == 'CATEGORYWater'] <- 'Water'
# rownames(mcor)[rownames(mcor) == 'CATEGORYPoorly Drained'] <- 'PoorlyDrained'
# rownames(mcor)[rownames(mcor) == 'CATEGORYWell Drained'] <- 'WellDrained'
# 
# col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# 
# corrplot(mcor, method='color', col=col(11), type="upper", tl.col="black", tl.srt=45, diag=FALSE)



#save(phenology_wide, file=paste0(dataDir, '/gls_improve_model_with_res_square1_to_21_in_2017.Rdata'))


