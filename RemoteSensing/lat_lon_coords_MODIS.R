#get the latitude and longitude of every MODIS pixel
#Dana Looschelders
#Jan 2020
#***********************************************************
setwd("C:/00 Dana/Uni/Internship/Work/remote sensing/MODIS_Data")

library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(gstat)
library(MODISTools)
library(MODIS)

#*************************************************************
#load workspace to get dataframe created in phenograss_process_MODIS
load("modis_pasture_data_A2018.RData")
#the x and y coordinates in the MODIS sinusoidal grid
x.coords=d$x_MODIS
y.coords=d$y_MODIS
#create new dataframe for lat and long cords
MODIS_coords=sin_to_ll(x = x.coords, y=y.coords)
#name columns and add column for MODIS coordinates
MODIS_coords=cbind("x_lat"=MODIS_coords$latitude_ll, "y_lon"= MODIS_coords$longitude_ll,"x_MODIS"=d$x_MODIS, "y_MODIS"=d$y_MODIS)
