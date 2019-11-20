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

radius = 1000  # Radius in m
focal = c(639598.2, 702988.7)
dsub = subset(d,  (x_ITM-focal[1])^2+(y_ITM-focal[2])^2<radius^2)

# Plot EVI time series
ggplot(dsub, aes(x=doy, y=evi)) + geom_point() + theme_bw()



# Model data using a GAM to model long-term changes 
# and seasonal trend
# Add in a tensor smooth to include spatial trend

# Create a factor for two periods
dsub$period = NA
dsub$period[dsub$year%in%c(2001:2005)] = "A"
dsub$period[dsub$year%in%c(2015:2019)] = "B"
dsub$period = as.factor(dsub$period)

# Fit a GAM, not including any autocorrleation
m1 = gam(evi~te(x_ITM, y_ITM) + s(julian) + s(doy, by=period) , 
        data=dsub)


m2 = gam(evi~te(x_ITM, y_ITM) + year + s(doy, by=period) , 
         data=dsub)
#plot(m)

gam.check(m1, pch='.')
summary(m1)
summary(m2)


# Visualise model predictions, 
# correcting for spatial variation and a broad temporal trend
d_predict = data.frame(period=rep(c('A','B'), each=360), 
                       doy = rep(c(1:360), times=2),
                       julian = 3800,
                       year=2015,
                       x_ITM = focal[1],
                       y_ITM=focal[2],
                       evi1=NA,
                       evi2=NA)
tmp1 = predict(m1, newdata = d_predict, se.fit = TRUE)
d_predict$evi1 = tmp1$fit
d_predict$evi1_se = tmp1$se.fit
tmp2 = predict(m2, newdata = d_predict, se.fit = TRUE)
d_predict$evi2 = tmp2$fit
d_predict$evi2_se = tmp2$se.fit

# Plot predictions from EVI1
ggplot(data=d_predict, aes(x=doy, 
                           y=evi1, 
                           ymax=evi+2*evi_se,
                           ymin=evi-2*evi_se,
                           fill=period)) +
  geom_ribbon() +
  geom_line() + 
  scale_fill_brewer(palette = 'Set2') +
  theme_bw()
