# Import processed MERA data and look at time-series
#
#
#
# Jon Yearsley (Jon.Yearsley@ucd.ie)
# June 2021
# ++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd('/Volumes/MODIS_data/MERA/')

library(ggplot2)

file_processed = "DailyData_subsetted/SoilMoist_subset_2017.RData"

file2_processed = "DailyData_subsetted/TotalPrecip_subset_2017.RData"

# Import data ----
load(file_processed)
sm = mera_subset 

load(file2_processed)
precip = mera_subset 

# Data wrangling -----

names(precip) = c("Longitude","Latitude","validityDate","Value","square")

# Convert dates to a date object
sm$date = as.Date(as.character(sm$ValidityDate),
                           format="%Y%m%d")
precip$date = as.Date(as.character(precip$validityDate),
                  format="%Y%m%d")


# Look at one pixel and plot data -------
sm_sub = subset(sm, Latitude==sm$Latitude[1] & 
                  Longitude==sm$Longitude[1])

precip_sub = subset(precip, Latitude==sm$Latitude[1] & 
                  Longitude==sm$Longitude[1])



ggplot(data=sm_sub,
       aes(x=date,
           y=DailyMean)) + 
  geom_point()  + 
  geom_linerange(data=precip_sub,
                 aes(ymin=0.15,
                     ymax=0.15+Value/40*0.2,
                     y=0.15+Value/40)) +
  labs(y="Daily Soil Moisture",
    title="Average Daily Soil Moisture 2017 and precipitation (bars)") +
  lims(y=c(0.15,0.3)) +
  theme_bw()

ggplot(data=precip_sub,
       aes(x=date,
           ymax=Value)) + 
  geom_linerange(aes(ymin=0))  + 
  labs(title="Average Daily Soil Moisture 2017") +
  theme_bw()
