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

# Create date going from July - July
d$doy_2 = d$doy
d$year_2 = as.numeric(as.character(d$year))
d$doy_2[d$doy_2>200] = d$doy_2[d$doy_2>200]-365  # NOTE: doesn't account for leap years
d$year_2[d$doy_2>200] = 1+d$year_2[d$doy_2>200]  # NOTE: doesn't account for leap years
d$year_2 = as.factor(d$year_2)



radius = 2000  # Radius in m
focal = c(639598.2, 702988.7)
dsub = subset(d,  (x_ITM-focal[1])^2+(y_ITM-focal[2])^2<radius^2)

# Plot EVI time series
#ggplot(dsub, aes(x=doy, y=evi)) + geom_point() + theme_bw()



# Model data using a GAM to model long-term changes 
# and seasonal trend
# Add in a tensor smooth to include spatial trend

# Create a factor for two periods
dsub$period = NA
dsub$period[dsub$year%in%c(2001:2004)] = "A"
dsub$period[dsub$year%in%c(2005:2008)] = "B"
dsub$period[dsub$year%in%c(2009:2012)] = "C"
dsub$period[dsub$year%in%c(2013:2016)] = "D"
dsub$period = as.factor(dsub$period)

# Create a factor for two periods
dsub$period_2 = NA
dsub$period_2[dsub$year_2%in%c(2001:2004)] = "A"
dsub$period_2[dsub$year_2%in%c(2005:2008)] = "B"
dsub$period_2[dsub$year_2%in%c(2009:2012)] = "C"
dsub$period_2[dsub$year_2%in%c(2013:2016)] = "D"
dsub$period_2 = as.factor(dsub$period_2)

# Fit a GAM, not including any autocorrleation
m1 = gam(evi~te(x_ITM, y_ITM) + s(julian) + s(doy, by=period) , 
        data=dsub)


m2 = gam(evi~te(x_ITM, y_ITM) + year + s(doy, by=period) , 
         data=dsub)


# Fit a model with year from July to July
m3 = gam(evi~te(x_ITM, y_ITM) + year_2 + s(doy_2, by=period_2) , 
         data=dsub)

#plot(m3, rug=T, residuals=T)

gam.check(m1, pch='.')
summary(m1)
summary(m2)


# Visualise model predictions, 
# correcting for spatial variation and a broad temporal trend
d_predict = data.frame(period=rep(c('A','B','C','D'), each=max(dsub$doy)), 
                       doy = rep(c(1:max(dsub$doy)), times=4),
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
                           y=evi2, 
                           ymax=evi2+2*evi2_se,
                           ymin=evi2-2*evi2_se,
                           fill=period)) +
  geom_ribbon(alpha=0.5) +
  geom_line(aes(colour=period)) + 
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  theme_bw()


# ***************************************************#
# Visualise model centred on winter 
# (i.e. year runs from July to July)

d_predict_2 = data.frame(period_2=rep(c('A','B','C','D'), each=1+diff(range(dsub$doy_2))), 
                       doy_2 = rep(c(min(dsub$doy_2):max(dsub$doy_2)), times=4),
                       julian = 3800,
                       year_2=2015,
                       x_ITM = focal[1],
                       y_ITM=focal[2],
                       evi=NA)
tmp = predict(m3, newdata = d_predict_2, se.fit = TRUE)
d_predict_2$evi = tmp$fit
d_predict_2$evi_se = tmp$se.fit

# Plot predictions from EVI1
ggplot(data=d_predict_2, aes(x=doy_2, 
                           y=evi, 
                           ymax=evi+2*evi_se,
                           ymin=evi-2*evi_se,
                           fill=period_2)) +
  geom_ribbon(alpha=0.5) +
  geom_line(aes(colour=period_2)) + 
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  theme_bw()
