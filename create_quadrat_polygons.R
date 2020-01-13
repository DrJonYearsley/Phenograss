# create_quadrat_polygons.R
#
# Reads in the eatsings and northings (TM75) of centre of 
# each quadrat, converts coordinates to MODIS CRS and then 
# defines a 10 km (roughly) box around the centre. 
# The corners of this box are saved as a shapefile in the MODIS CRS.
# 
# Jon Yearsley  jon.yearsley@ucd.ie
# Jan 2020
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())


library(sp)
library(rgdal)
library(raster)

rm(list=ls())

setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology')
outputDir = './Data'

# Define a square quadrat and the size of the square
square_centre_TM75 = data.frame(east=116386, north=328904) # Define in Irish Grid TM75 (espg:29903)
square_size_m = 10000  # units = m. Each pixel is roughly 250 m


# Load MODIS grid data
modis = raster('./Data/MODIS/modis_grid_ireland.grd')

# Save MODIS CRS
modis_crs = crs(modis)

side = SpatialPoints(data.frame(east=square_centre_TM75$east+c(0.5,-0.5)*square_size_m,
                              north=square_centre_TM75$north)
                   , proj4string=CRS("+init=epsg:29903"))

# Convert square and side length into  MODIS
square_centre_modis = spTransform(square_centre_TM75, modis_crs)
side_modis = spTransform(side, modis_crs)

# define a square within Ireland
sq_df = data.frame(east = square_size_m * c(-0.5, -0.5, 0.5, 0.5) + square_centre_TM75$east, 
                   north = square_size_m * c(-0.5, 0.5, 0.5, -0.5)+ square_centre_TM75$north)

sq = SpatialPointsDataFrame(sq_df, data=sq_df, proj4string=CRS("+init=epsg:29903"))
sq_name = 'test'

# Convert square to MODIS CRS and extract extent
sq_modis = spTransform(sq, CRS=modis_crs)


writeOGR(sq_modis, dsn=file.path(outputDir,paste0(sq_name,'.shp')), 
         layer=sq_name, 
         driver='ESRI Shapefile') 



