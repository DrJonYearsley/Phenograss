# Read in MERA data and subset for the 20 10Km squares
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# June 2021
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd('~/git_repos/Phenograss/RemoteSensing/')

library(rgdal)
library(sp)


# Import MERA data
d = read.table('~/WorkFiles/Data/Climate/MERA/TPrecip/TotalPrecip_2012_01.txt',
               sep='',
               header=TRUE,
               colClasses = rep('numeric', 7),
               nrow=2*10^5)

d$Longitude.[d$Longitude.>180] =  d$Longitude.[d$Longitude.>180]-360


# Convert MERA data to be a spatial object
coordinates(d) = ~Longitude. + Latitude.
latlong = "+init=epsg:4326"
proj4string(d) = CRS(latlong)
str(d)


# Import quadrat data
squares = readOGR(dsn='~/WorkFiles/PeopleStuff/GrasslandPhenology/Data/Quadrats/agriclimate_quadrats_Ireland.shp',
                  layer='agriclimate_quadrats_Ireland')


squares_WGS84 = spTransform(squares, CRSobj = proj4string(d))

# Transform MERA data onto squares CRS (modis CRS)
d_MODIS = spTransform(d, CRSobj = proj4string(squares[1,]))



tmp = over(squares[2,], d_MODIS, returnList=TRUE)


plot(squares_WGS84)
plot(d, add=T)
