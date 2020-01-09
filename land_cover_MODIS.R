# Calculate percentage of a MODIS pixel that has various land cover classes
#
# Primarily used to calculate percentage of pasture within a MODIS pixel 
# using either CORINE or the semi-natural grassland survey polygons
#
# Jon Yearsley  Jan 2020
#
###########################################################################

library(sp)
library(rgdal)
library(raster)

rm(list=ls())

setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology/Data')

# Directories containing the input and output MODIS data 
outputDir = './MODIS'
outputSuffix = 'pasture'


# Define location of semi-natural grassland data
grasslandPath = './IrishSeminaturalGrasslandSurvey/ISG13_GIS_Datasets/'

# Read in semi-natural grassland data. Make sure this is vector data (shapefile)
grassland = readOGR(dsn=file.path(grasslandPath,'ISGS13_Habitats_01a.shp'), layer='ISGS13_Habitats_01a')

# Define location of CORINE data
corinePath = './CORINE_Ireland/CLC18_IE_ITM'
# Read in CORINE data. make sure this is vector data (shapefile)
corine = readOGR(dsn=file.path(corinePath,'CLC18_IE_ITM.shp'), layer='CLC18_IE_ITM')


# Load MODIS grid data
modis = raster('./MODIS/modis_grid_ireland.grd')
# Save MODIS CRS
modis_crs = crs(modis)

# Transform grassland data into MODIS CRS
grassland_modis = spTransform(grassland, CRS=modis_crs)
corine_modis = spTransform(corine, CRS=modis_crs)


# define a square within Ireland
# ####
# # Define a square quadrat and the size of the square
# square_centre_TM75 = data.frame(east=279084, north=211675) # Define in Irish Grid TM75 (espg:29903)
# square_size_m = 10000  # units = m. Each pixel is roughly 250 m
# 
# sq_df = data.frame(east = square_size_m * c(-0.5, -0.5, 0.5, 0.5) + square_centre_TM75$east, 
#                    north = square_size_m * c(-0.5, 0.5, 0.5, -0.5)+ square_centre_TM75$north)
# 
# sq = SpatialPoints(sq_df, proj4string=CRS("+init=epsg:29903"))
# 
# # Convert square to MODIS CRS and extract extent
# sq_modis = extent(spTransform(sq, CRS=modis_crs))
# grassland_sq = crop(grassland_modis, sq_modis)


# Subset sn grassland for certain Forrit habitats
grassland_subset = subset(grassland_modis, FOSS_HAB %in% c('GS','GS1','GS2','GS3','GS4'))
sn_grass_cover= rasterize(grassland_subset, modis, getCover=T)

# Subset CORINE for agricultural pasture
corine_subset = subset(corine_modis, CODE_18%in%231)
corine_pasture_cover = rasterize(corine_subset, modis, getCover=T)


writeRaster(sn_grass_cover,filename='sn_grasscover',format='raster')
writeRaster(sn_grass_cover,filename='sn_grasscover',format='GTiff')


writeRaster(corine_pasture_cover,filename='corine_pasturecover',format='raster')
writeRaster(corine_pasture_cover,filename='corine_pasturecover',format='GTiff')
