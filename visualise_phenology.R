# Imports data frame containing MODIS VI data
# and visualises phenology
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Nov 2019
# ************************************************

# New comment
library(ggplot2)
library(mgcv)
library(gstat)
library(sp)

rm(list=ls())
#setwd("F:/")
setwd('~/WorkFiles/PeopleStuff/GrasslandPhenology/Data/MODIS')

load('modis_pasture_data_A20.RData')

# Only keep observations with a finite EVI
d = subset(d, is.finite(evi))

# Create a date for each observation
d$date = strptime(paste0(d$year,'-',d$doy), format="%Y-%j")

radius = 500  # Radius in m
focal = c(639598.2, 702988.7)
dsub = subset(d,  (x_ITM-focal[1])^2+(y_ITM-focal[2])^2<radius^2)

# Plot EVI time series (for subset)
ggplot(dsub, aes(x=doy, y=evi)) + geom_point() + theme_bw()

# Plot full time series (for subset)
ggplot(dsub, aes(x=as.Date(date), 
                 y=evi)) + 
  geom_smooth(method='loess', span=0.01) + 
  geom_point(shape=19, size=0.5, alpha=0.5) + 
  labs(x='Date', y='EVI') +
  theme_bw()

# Look for saturation in NDVI
ggplot(d, aes(x=ndvi, y=evi)) + geom_point(alpha=0.3, shape='.', size=2) + theme_bw()


# Create a spatial data frame
d_sp <- SpatialPointsDataFrame(coords = d[,c('x_ITM','y_ITM')], data=d, 
                          proj4string = CRS("+init=epsg:2157"))


# Look at variogram of the EVI signal through space
v = variogram(evi~1, data=subset(d_sp, year==2017))

# Now fit a variogram using a Spherical, Matern and an Exponential model and give best fit
# Warnings are normal for the Matern model
(vfit = fit.variogram(v, vgm(c("Exp", "Sph",  "Mat", "Ste")), fit.kappa = TRUE))

# Plot empirical variogram and the fitted model (vfit)
plot(v, model=vfit)
