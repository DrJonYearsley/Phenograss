# Calculate percentage of a MODIS pixel that has various land cover classes
#
# Primarily used to calculate percentage of pasture within a MODIS pixel 
# using either CORINE or the semi-natural grassland survey polygons
#
# Jon Yearsley  Jan 2020
# jon.yearsley@ucd.ie
#
# **********************************************************************

library(sp)
library(rgdal)
library(raster)

rm(list=ls())

# Define top level data directory
setwd('~/Research/Phenograss/Data')

# Load MODIS grid data
modis = raster('./MODIS/modis_grid_ireland.grd')
# Save MODIS CRS
modis_crs = crs(modis)

# Analyse semi-natural grassland  data ---------------

# Define location of semi-natural grassland data
grasslandPath = './IrishSeminaturalGrasslandSurvey/ISG13_GIS_Datasets/'

# Read in semi-natural grassland data. Make sure this is vector data (shapefile)
grassland = readOGR(dsn=file.path(grasslandPath,'ISGS13_Habitats_01a.shp'), layer='ISGS13_Habitats_01a')


# Transform grassland data into MODIS CRS
grassland_modis = spTransform(grassland, CRS=modis_crs)

# Subset sn grassland for certain Fossit habitats
grassland_subset = subset(grassland_modis, FOSS_HAB %in% c('GS','GS1','GS2','GS3','GS4'))

# Calculate fraction of MODIS pixel covered by grass
sn_grass_cover= rasterize(grassland_subset, modis, getCover=T)

writeRaster(sn_grass_cover,filename='sn_grasscover',format='raster')
writeRaster(sn_grass_cover,filename='sn_grasscover',format='GTiff')


# Analyse CORINE data ---------------

# Define location of CORINE data (2018)
corinePath = './CORINE_Ireland/clc18_All_Ireland'

# Read in CORINE data. make sure this is vector data (shapefile)
corine = readOGR(dsn=file.path(corinePath,'clc18_All_Ireland.shp'), layer='clc18_All_Ireland')

# Transform grassland data into MODIS CRS
corine_modis = spTransform(corine, CRS=modis_crs)

# Subset CORINE for agricultural pasture
corine_subset = subset(corine_modis, Code_18%in%231)

# Calculate fraction of MODIS pixel covered by grass
corine_pasture_cover = rasterize(corine_subset, modis, getCover=T)
  
writeRaster(corine_pasture_cover,filename='corine2018_pasturecover_All_Ireland',format='raster')
writeRaster(corine_pasture_cover,filename='corine2018_pasturecover_All_Ireland',format='GTiff')
