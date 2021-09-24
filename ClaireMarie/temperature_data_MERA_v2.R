# 
# Display of the temperatures data from MERA for each square
#
# Claire-Marie
# 24/06/2021
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
tempPath = '~/Stage/Data/Climate'
quadratPath = '~/Stage/Data/Quadrats'
modisPath = '~/Stage/Data/MODIS'

datadir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created/temperature_data'

# List of squares to analyse
squareList = c(20)
year = 2013


# Read in MODIS grid
modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
crs_modis = st_crs(modis)

# Import squares from shapefile
squares = st_read(file.path(quadratPath,"agriclimate_quadrats_Ireland.shp"))
squares_modis = st_transform(squares, crs_modis)  # Make sure CRS is identical (required for gstat interpolation)
crs_squares = st_crs(squares)

# Import soil data from shapefile
filename = paste0('temperature_degrees_2012_2013_2015.RData')
load(file.path(tempPath,filename))


# Temperature data by squares
for (s in squareList) {
  
  # Import phenophase estimate
  filename = paste0(input_file_prefix,'_square_',s,'_',year,'.RData')
  load(file.path(datadir,filename))
  
  # Select the square s in January / February / March
  mera_square = subset(temperature_mera, 
                       temperature_mera$square == s 
                       & temperature_mera$validityDate > paste0(year, '0101')
                       & temperature_mera$validityDate < paste0(year, '0601'))
  
  # Create a geometry column
  mera_square_sf = st_as_sf(mera_square, coords = c("Longitude", "Latitude"), 
                            crs = st_crs("EPSG:4326"))
  
  # Transform into the modis CRS
  mera_square_modis = st_transform(mera_square_sf, crs = crs_modis)
  
  crop_square = st_crop(modis, squares_modis[s,,])
  
  # Create a geometry column
  output_smoothed = subset(output_smoothed,phase==1 & warning==FALSE & !is.na(t))
  output_smoothed = st_as_sf(output_smoothed, 
                             coords = c("x_MODIS", "y_MODIS"), crs = crs_modis)
  
  # Reduce day values to 0 decimal places
  output_smoothed$t = round(output_smoothed$t, digit=0)
  
  # Match the day of the year to the GPS day
  dates =  unique(mera_square_modis$validityDate)
  nDates = length(dates)
  days = c(1:nDates)
  print(dates)
  print(nDates)
  
  pixel_list = unique(output_smoothed$pixelID)
  nPixel = length(pixel_list)
  
  tmp = data.frame()
  
  for (i in pixel_list) {
    print(i)
    
    # Retrieve the index of one of the pixels
    ind = which(output_smoothed$pixelID == i)[1]

    day = output_smoothed$t[ind]
    
    if (day > 0 & day <= nDates) {

      test = subset(mera_square_modis, validityDate==dates[day])
      
      # +++++++++++++++++++++++++++++++++++++++++++++
      # Interpolate MERA data
      
      # Method 2:
      # Interpolate using a model variogram (try a linear variogram)
      
      # Look at the empirical variogram
      v_emp = variogram(Value ~ 1, data = test)
      
      if (length(which(v_emp$dist == 0.0)) != 0) {
        # Fit variogram model (try linear)  use show.vgm() to display all possible models
        v_mod = fit.variogram(v_emp, model = vgm(NA,"Lin",0))
        
        # Look at fit of variogram model versus the empirical variiogram
        plot(v_emp,v_mod)
        
        # Now do some ordinary krigging to interpolate
        g_mod = gstat::krige(formula = Value ~ 1, 
                             locations=test,
                             model=v_mod,
                             newdata=crop_square)
        
        # +++++++++++++++++++++++++++++++++++++++++++++
        
        # Create a temperature colomn thanks to geometry
        output_smoothed$temperature[ind] = as.numeric(st_extract(g_mod,
                                                                  output_smoothed$geometry[ind]))
        
        tmp = rbind(tmp, output_smoothed[ind,])
      }
    }
    
  }
  
  # For viewing
  tmp$x_MODIS = sapply(tmp$geometry,"[[",1)
  tmp$y_MODIS = sapply(tmp$geometry,"[[",2)
  
  # not the best solution for the moment
  tmp$temperature = round(tmp$temperature, digit=0)
  tmp$class[tmp$temperature <= 2] = '< 2'
  tmp$class[tmp$temperature > 2 & tmp$temperature <= 5] = '< 5'
  tmp$class[tmp$temperature > 5 & tmp$temperature <= 8] = '< 8'
  tmp$class[tmp$temperature > 8 & tmp$temperature <= 10] = '< 10'
  tmp$class[tmp$temperature > 10] = '> 10'
  
  # Display temperature per SOS
  ggplot(data=tmp,
         aes(x=t,
             y=temperature)) +
    geom_point() +
    labs(x='Day of year',
         y='Temperature',
         title=paste('Temperature data square for day of year 2013 square ', s)) +
    theme_bw()
  
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/temperature_date_square_sos_',s,'.png'))
  
  # Display
  ggplot(data=tmp,
         aes(x=x_MODIS,
             y=y_MODIS,
             fill= class)) +
    geom_tile(colour='darkblue') +
    coord_equal() +
    scale_fill_brewer(palette='Blues') +
    labs(x='X Coord (MODIS CRS)',
         y='Y Coord (MODIS CRS)',
         title=paste('Temperature data square', s)) +
    theme_bw()
  
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/temperature_square_sos_',s,'.png'))

}