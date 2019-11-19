# Imports data frame containing MODIS VI data
# and starts to model the phenology of the evi signal
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Nov 2019
# ************************************************


library(ggplot2)
library(mgcv)
library(gstat)

rm(list=ls())
#setwd("F:/")
setwd('~/WorkFiles/PeopleStuff/GrasslandPhenology/Data/MODIS')

load('modis_pasture_data_A20.RData')


# Only keep observations with a finite EVI
d = subset(d, is.finite(evi))

# Create a date for each observation
d$date = strptime(paste0(d$year,'-',d$doy), format="%Y-%j")
d$julian = as.numeric(julian(d$date, origin=as.POSIXct("2000-01-01", tz = "GMT")))

radius = 500  # Radius in m
focal = c(639598.2, 702988.7)
dsub = subset(d,  (x_ITM-focal[1])^2+(y_ITM-focal[2])^2<radius^2)

# Plot EVI time series
ggplot(dsub, aes(x=doy, y=evi)) + geom_point() + theme_bw()



# Model data using a GAM to model long-term changes 
# and seasonal trend
# Add in a tensor smooth to include spatial trend

dsub$period = NA
dsub$period[dsub$year%in%c(2001:2005)] = "A"
dsub$period[dsub$year%in%c(2013:2017)] = "B"
dsub$period = as.factor(dsub$period)

m = gam(evi~s(julian) + s(doy, by=period) + te(x_ITM, y_ITM), 
        data=dsub)

plot(m)
