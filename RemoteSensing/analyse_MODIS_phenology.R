# Analyse MODIS phenology
#
#
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Jan 2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(segmented)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)

dataDir = '~/Research/Phenograss/Data/PhenologyOutput'

input_file_preffix = 'phenology'

# Import data --------
squares = c(1:9,13:21)
years = c(2003:2004,2011:2014,2016:2017)


# Filename segmented data
for (i in 1:length(squares)) {
  for (y in 1:length(years)) {
    filename = paste0(input_file_preffix,'_square_',squares[i],'_',years[y],'.RData')
    load(file.path(dataDir,filename))
    
    output_smoothed$square = squares[i]
    
    if (y==1 & i==1) {
      phenology = output_smoothed
    } else {
      phenology = rbind(phenology,
                        output_smoothed)
    }
  }
}

# Create a wide version of phenology
phenology_wide = pivot_wider(data=subset(phenology, warning==FALSE),
                             id_cols = c('pixelID','year', 'x_MODIS','y_MODIS','x_ITM','y_ITM','square'),
                             names_from = 'phase',
                             names_prefix = 'phase',
                             values_from = c(t,slope),
                             values_fn = mean)

# Centre the x and y coordinates
phenology_wide$x_ITM_centre = phenology_wide$x_ITM - mean(phenology_wide$x_ITM, na.rm=TRUE)
phenology_wide$y_ITM_centre = phenology_wide$y_ITM - mean(phenology_wide$y_ITM, na.rm=TRUE)

# Coordinates
locations = subset(phenology_wide, 
                   subset=year==years[1],
                   select=c('pixelID','x_ITM','y_ITM'))

# Calculate the centre coordinates for each square
tmp = aggregate(cbind(x_ITM_centre, y_ITM_centre)~square+year, 
                data=phenology_wide, 
                FUN=mean, 
                na.rm=TRUE)


# # Try an INLA model
# library(inlabru)
# library(INLA)


# Fit model for phenology dates of phase 1 -----


# Mixed model for phase 1
m_phase1 = lme(t_phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)

summary(m_phase1)

# # Fit this same model using INLA
# m_inla = inla(t_phase1~1+as.factor(year)+f(pixelID, model='iid'),
#               data=phenology_wide,
#               family='gaussian',
#               control.compute = list(dic = TRUE))


# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time

phenology_wide$dummy = interaction(phenology_wide$year, phenology_wide$square, sep='_')
gls_phase1 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)

summary(gls_phase1)


gls_null1 = update(gls_phase1, .~.-x_ITM_centre)
gls_null2 = update(gls_phase1, .~.-y_ITM_centre)

acf(gls_phase1$residuals, lag.max=20, type='partial')  # Makes sense... spatial correlation roughly 1-2 km


library(emmeans)

m_eff = emmeans(gls_phase1, spec='year')
summary(m_eff)


gls_phase1_v2 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                 data=phenology_wide,
                 na.action=na.exclude)


# # INLA version of spatial model
# Mesh <- inla.mesh.2d(locations, 
#                      max.edge = c(10, 10))
