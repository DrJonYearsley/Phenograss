# 
# Display of the elevation of the 21 squares
#
# Claire-Marie
# 31/05/2021
# 
# Edited Jon Yearsley (2/6/2021)
#
#1. Import squares
#2. Import GTOPO30 data (easier to start because it's 1 raster)
#3. Put squares in the same reference WGS84
#4. Cut modis data by squares
#5. Resample GTOPO30 data with modis squares data
#6. Collect data
#
# ****************************************


rm(list=ls())

library(sf) # combine sp, rgeos et rgdal
library(stars)
library(ggplot2)


quadratPath = '~/Stage/Data/Quadrats'
outputDir = '~/Stage/Data_created/elevation_GTOPO30'

GTOPO30Path = '~/Stage/Data/DigitalElevation/GTOPO30'
GTOPO30Filename = 'gt30w020n90.tif'


# List of squares to analyse
squareList = c(1:9,13:21)

# GTOPO30 data : WGS84
GTOPO30 = read_stars(file.path(GTOPO30Path,GTOPO30Filename), proxy = TRUE)
crs_GTOPO30 = st_crs(GTOPO30)

# Import squares from shapefile
squares = st_read(file.path(quadratPath,"agriclimate_quadrats_Ireland.shp"))
crs_squares = st_crs(squares)


# Convert squares to gtopo crs (because cropping is faster on a regular grid)
squares_gtopo = st_transform(squares, crs = crs_GTOPO30)

# Read in MODIS grid
modis = read_stars('~/Stage/Data/MODIS/modis_grid_ireland.tif')
st_crs(modis) =  crs_squares  # Make sure modis and squares have same CRS


# Elevation data by squares
for (s in squareList) {
  
  # Crop modis grid to a square and then resample elevation onto this cropped grid
  crop_square = st_crop(modis, squares[s,])
  crop_modis = st_warp(GTOPO30, crop_square, method="cubic", use_gdal=T)
  
  # Save data to a GeoTiff file
  fname = paste0("GTOPO30_square_",s,".tif")
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  write_stars(adrop(crop_modis), file.path(outputDir,fname))
  
  
  
  # Data display
  ggplot() + 
    geom_stars(data=crop_modis) +
    coord_equal() + 
    scale_fill_gradientn(colours = terrain.colors(20)) +
    labs(x='X Coord (MODIS CRS)',
         y='Y Coord (MODIS CRS)',
         title=paste('Elevation GTOPO30 square', s)) +
    theme_bw()
  
  # Save data for visualisation
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/GTOPO30_square_', s, '.png'))
}