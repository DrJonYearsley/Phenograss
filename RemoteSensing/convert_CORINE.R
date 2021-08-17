# Convert CORINE to ITM and crop to Ireland
#
# Primarily used to calculate percentage of pasture within a MODIS pixel 
# using either CORINE or the semi-natural grassland survey polygons
#
# Jon Yearsley  Mar 2020
# jon.yearsley@ucd.ie
#
# **********************************************************************



library(sp)
library(rgdal)
library(raster)

rm(list=ls())

# Define top level data directory
setwd('~/Research/Phenograss/Data')

corineFile = './Corine_Europe/u2018_clc2018_v2020_20u1_fgdb/DATA/U2018_CLC2018_V2020_20u1.gdb/'

# Define extent of Irleand (roughly) in MODIS CRS 
ir = extent(-7.5E5, -3.3E5,5.7E6, 6.17E6)

# Define extent of Irleand (roughly) in CORINE CRS
#ir = extent(-7.5E5, -3.3E5,5.7E6, 6.17E6)


# Load CORINE data
ogrListLayers(corineFile)
ogrInfo(corineFile, "U2018_CLC2018_V2020_20u1")

corine = readOGR(dsn=corineFile,
                 layer='U2018_CLC2018_V2020_20u1',
                 require_geomType = "wkbPolygon",
                 use_iconv = TRUE)


# Load MODIS grid data
modis = raster('./MODIS/modis_grid_ireland.grd')
# Save MODIS CRS
modis_crs = crs(modis)


