# 
# Display of the soil moisture data from MERA for each square
#
# Claire-Marie
# 11/06/2021
#
# Edited to include two types of interpolation (Jon Yearsley 17th June)
#
#1. Import modis grid
#2. Import soil moisture MERA data
#3. Import phenophase square
#4. Put soil data in modis
#5. Transform dataframe soil data in stars class
#6. Extract soil data
#7. Collect data
#
# ****************************************


rm(list=ls())

library(sf) # combine sp, rgeos et rgdal
library(stars)
library(ggplot2)
library(gstat)
library(dplyr)


input_file_prefix = 'phenology'
soilPath = '~/Stage/Data/Climate/SoilMoisture_quadrats'
quadratPath = '~/Stage/Data/Quadrats'
modisPath = '~/Stage/Data/MODIS'

datadir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created/soil_moisture_data'


#soilPath = '/Volumes/MODIS_data/MERA/DailyData_subsetted/'
#quadratPath = '/Volumes/MODIS_data/Quadrats/'
#modisPath = '/Volumes/MODIS_data/MODIS/'

#datadir = '~/Research/Phenograss/Data/PhenologyOutput/'
#outputDir = '~/Research/Phenograss/Data/MERA_processed/'


# List of squares to analyse
squareList = c(1)
year = 2017


# Read in MODIS grid
modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
crs_modis = st_crs(modis)

# Import squares from shapefile
squares = st_read(file.path(quadratPath,"agriclimate_quadrats_Ireland.shp"))
squares_modis = st_transform(squares, crs_modis)  # Make sure CRS is identical (required for gstat interpolation)
crs_squares = st_crs(squares)

# Import soil data from shapefile
filename = paste0('SoilMoist_subset_2017.RData')
load(file.path(soilPath,filename))


# Soil moisture data by squares
for (s in squareList) {
  
  # Import phenophase estimate
  filename = paste0(input_file_prefix,'_square_',s,'_',year,'.RData')
  load(file.path(datadir,filename))
  
  # Select the square s in January / February / March
  mera_square = subset(mera_subset, 
                       mera_subset$SquareID == s & mera_subset$ValidityDate < paste0(year, '0401'))

  # Create a geometry column
  mera_square_sf = st_as_sf(mera_square, coords = c("Longitude", "Latitude"), 
                         crs = st_crs("EPSG:4326"))

  # Transform into the modis CRS
  mera_square_modis = st_transform(mera_square_sf, crs = crs_modis)

  # +++++++++++++++++++++++++++++++++++++++++++++
  # Interpolate MERA data

  crop_square = st_crop(modis, squares_modis[s,,])
  
  # Method 1: 
  # Inverse distance weighted interpolation of the MERA data
  g_idw = gstat::idw(formula = DailyMean ~ 1, 
                 locations=subset(mera_square_modis, ValidityDate==20170101), 
                 newdata=crop_square, 
                 idp=2)
  
  # Plot the result
  ggplot() +
    geom_sf(data=squares[s,]) +
    geom_stars(data=g_idw)
  
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/inverse_distance_', s, '.png'))
  
  # Method 2:
  # Interpolate using a model variogram (try a linear variogram)
  
  # Look at the empirical variogram
  v_emp = variogram(DailyMean ~ 1, data = subset(mera_square_modis, ValidityDate==20170101))
  plot(v_emp)
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/variogramme_empiric_', s, '.png'))
  
  # Fit variogram model (try linear)  use show.vgm() to display all possible models
  v_mod = fit.variogram(v_emp, model = vgm(NA,"Lin",0))

  # Look at fit of variogram model versus the empirical variiogram
  plot(v_emp,v_mod)
  
  # Now do some ordinary krigging to interpolate
  g_mod = gstat::krige(formula = DailyMean ~ 1, 
                 locations=subset(mera_square_modis, ValidityDate==20170101),
                 model=v_mod,
                 newdata=crop_square)

  # Plot the result 
  ggplot() +
    geom_sf(data=squares[s,]) +
    geom_stars(data=g_mod)
  
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/krigging_', s, '.png'))

  # +++++++++++++++++++++++++++++++++++++++++++++
}