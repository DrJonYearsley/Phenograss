# Imports data frame containing MODIS VI data
# and performs EDM analysis
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Nov 2019
# ************************************************

library(rEDM)
library(ggplot2)
library(mgcv)
library(gstat)

rm(list=ls())
#setwd("F:/")
setwd('~/WorkFiles/PeopleStuff/GrasslandPhenology/EDM')

load('modis_multi_landcover_data_A2017.RData')


# pasture = d$pasture
# urban = d$urban
# forest = d$forest


d_use = d$urban
# Create a date for each observation
d_use$date = strptime(paste0(d_use$year,'-',d_use$doy), format="%Y-%j")
d_use$julian = as.numeric(julian(d_use$date, origin=as.POSIXct("2000-01-01", tz = "GMT")))

# Find all unique locations
tmp = unique(d_use[,c('x_ITM','y_ITM')])
locations = data.frame(x=tmp[,1], y=tmp[,2],label=c(1:nrow(tmp)))
for (l in 1:nrow(locations)) {
  ind = d_use$x_ITM==locations$x[l] & d_use$y_ITM==locations$y[l]
  d_use$location[ind] = locations$label[l]
}

# Sort data frame in terms of julian day
d_use = d_use[order(d_use$julian),]


# Perform EDM analysis on each location
for (l in 1:nrow(locations)) {
  d_location = droplevels(subset(d_use, location==l))
  training = c(1:floor(2*nrow(d_location)/3))
  prediction = c(ceiling(2*nrow(d_location)/3):nrow(d_location))
  
  # Apply EDM method to estimate embedding dimension
  simplex_output = simplex(d_location$evi, 
                           lib = training, 
                           pred = prediction, 
                           E=c(1:10))
  
  if (l==1) {
    edm_out = data.frame(location=l, simplex_output)
  } else {
    edm_out = rbind(edm_out, data.frame(location=l, simplex_output))
  }
}

# Create a spatial subset
radius = 10000  # Radius in m around a point
focal = c(639598.2, 702988.7)  # The focal point

# Aggregate data across spatial region
#dsub = subset(d_use,  (x_ITM-focal[1])^2+(y_ITM-focal[2])^2<radius^2)
dsub = aggregate(evi ~ year+doy, data=d_use, FUN=median, na.rm=T)
  
# Plot EVI time series
ggplot(dsub, aes(x=doy, y=evi)) + geom_point() + theme_bw()



# Model data using a GAM to model long-term changes 
# and seasonal trend
# Add in a tensor smooth to include spatial trend


# Fit a GAM, not including any autocorrleation
m1 = gam(evi~te(x_ITM, y_ITM) +  s(doy) , 
        data=dsub)




plot(m1, rug=T, residuals=T)

gam.check(m1, pch='.')
summary(m1)

# Apply EDM method to estimate embedding dimension
simplex_output = simplex(dsub$evi, lib = c(1, 50), pred = c(51, 80), E=c(1:20))

plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")
