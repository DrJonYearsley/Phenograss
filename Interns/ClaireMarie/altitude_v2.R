# Display of the altitude of the 21 squares
#
#
# Claire-Marie
# 31/05/2021
# 
# Edited Jon Yearsley (2/6/2021)
#
#1. Import squares
#2. Import GTOPO30 data (easier to start because it's 1 raster)
#3. Put data in the same reference MODIS
#4. Cut data by squares
#5. Divider GTOPO30 data in modis pixels
#6. Collect data
#
# ****************************************


rm(list=ls())
setwd('/Volumes/MODIS_data')

library(sf) # combine sp, rgeos et rgdal
library(stars)
library(ggplot2)



quadratPath = 'Quadrats/'
outputDir = 'Data_created/'

GTOPO30Path = 'Altitude/'
GTOPO30Filename = 'gt30w020n90.tif'


# List of squares to analyse
squareList = c(1:21)
# squareList = c(1:9,13:21)

# GTOPO30 data : WGS84
GTOPO30 = read_stars(file.path(GTOPO30Path,GTOPO30Filename), proxy = TRUE)
crs_GTOPO30 = st_crs(GTOPO30)
print(crs_GTOPO30)


# Import squares from shapefile
squares = st_read(file.path(quadratPath,"agriclimate_quadrats_Ireland.shp"))
crs_squares = st_crs(squares)
# contains easting and northing

# Read in MODIS grid
modis = read_stars('MODIS/modis_grid_ireland.tif')
st_crs(modis) =  crs_squares  # Make sure modis and squares have same CRS


# Elevation data by squares --- to review => MODIS pixel
for (s in squareList) {
  print(squares[s,])
  # Crop to one of the squares
  # https://r-spatial.github.io/stars/reference/st_as_sf.html
  
  # Crop modis grid to a square and then resample elevation onto this cropped grid
  crop_square = st_crop(modis, squares[s,])
  crop_modis = st_warp(GTOPO30, crop_square, method="near", use_gdal=T)
  
  
  # Save data to a GeoTiff file
  fname = paste0("GTOPO30_square_",s,".tif")
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  write_stars(adrop(crop_modis), file.path(outputDir,fname) ) 
}


# Data display
ggplot() + 
  geom_stars(data=crop_modis) +
  coord_equal() + 
  scale_fill_gradientn(colours = terrain.colors(20)) +
  theme_bw()




