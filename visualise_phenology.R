# Imports data frame containing MODIS VI data
# and visualises phenology
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Nov 2019
# ************************************************


library(ggplot2)
library(mgcv)
library(gstat)

gstatrm(list=ls())
#setwd("F:/")
setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology/RScript/MODIS_Testing/')

load('modis_pasture_data_A2018.RData')


radius = 10  # Radius in m
focal = c(639598.2, 702988.7)
dsub = subset(d,  (x_ITM-focal[1])^2+(y_ITM-focal[2])^2<radius^2)

# Plot EVI time series
ggplot(dsub, aes(x=doy, y=ndvi)) + geom_point() + theme_bw()

# Look for saturation in NDVI
ggplot(d, aes(x=evi, y=ndvi)) + geom_point(alpha=0.5, shape='.', size=2) + theme_bw()


# Look at variogram of the EVI signal through space
variogram(evi~1, data=subset(d, is.finite(evi)), locations=d[,c('x_ITM','y_ITM')])
