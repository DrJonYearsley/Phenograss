# 
# Display of the precipitation data from MERA for each square
# (similarly same as temperature (temperature_data_MERA.R))
#
# Claire-Marie
# 25/06/2021
#
#
# ****************************************


rm(list=ls())

library(sf) # combine sp, rgeos et rgdal
library(stars)
library(ggplot2)
library(gstat)
library(dplyr)


input_file_prefix = 'phenology'
precipPath = '~/Stage/Data/Climate/Precip_quadrats'
quadratPath = '~/Stage/Data/Quadrats'
modisPath = '~/Stage/Data/MODIS'

datadir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created/precipitation_data'

# Save file
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

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
filename = paste0('TotalPrecip_subset_2013.RData')
load(file.path(precipPath,filename))


# Soil moisture data by squares
for (s in squareList) {
  
  # Import phenophase estimate
  filename = paste0(input_file_prefix,'_square_',s,'_',year,'.RData')
  load(file.path(datadir,filename))
  
  # Select the square s in January / February / March
  mera_square = subset(mera_subset, 
                       mera_subset$SquareID == s & mera_subset$ValidityDate < paste0(year, '0601'))

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
  
  # réduire les valeurs des jours à 0 chiffres après la virgule
  output_smoothed$t = round(output_smoothed$t, digit=0)
  
  # faire correspondre le jour de l'année au jour GPS
  dates =  unique(mera_square_modis$ValidityDate)
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
    print(day)
    
    if (day > 0 & day <= nDates) {
      
      test = subset(mera_square_modis, ValidityDate==dates[day])
      
      # +++++++++++++++++++++++++++++++++++++++++++++
      # Interpolate MERA data
      
      # Method 2:
      # Interpolate using a model variogram (try a linear variogram)
      
      # Look at the empirical variogram
      v_emp = variogram(DailyMean ~ 1, data = test)
      plot(v_emp)
      
      if (length(which(v_emp$dist == 0.0)) != 0) {

        # Fit variogram model (try linear)  use show.vgm() to display all possible models
        v_mod = fit.variogram(v_emp, model = vgm(NA,"Lin",0))
        
        # Look at fit of variogram model versus the empirical variiogram
        plot(v_emp,v_mod)
        
        # Now do some ordinary krigging to interpolate
        g_mod = gstat::krige(formula = DailyMean ~ 1, 
                             locations=test,
                             model=v_mod,
                             newdata=crop_square)
        
        # +++++++++++++++++++++++++++++++++++++++++++++
        
        # Create a precipitation colomn thanks to geometry
        output_smoothed$precipitation[ind] = as.numeric(st_extract(g_mod,
                                                                   output_smoothed$geometry[ind]))
        
        tmp = rbind(tmp, output_smoothed[ind,])
      }
    }
    
  }
  
  # For viewing
  tmp$x_MODIS = sapply(tmp$geometry,"[[",1)
  tmp$y_MODIS = sapply(tmp$geometry,"[[",2)
  
  # not the best solution for the moment
  tmp$precipitation = round(tmp$precipitation, digit=2)
  tmp$class[tmp$precipitation <= 0.04] = '< 0.04'
  tmp$class[tmp$precipitation > 0.04 & tmp$precipitation <= 0.25] = '< 0.25'
  tmp$class[tmp$precipitation > 0.25 & tmp$precipitation <= 1] = '< 1'
  tmp$class[tmp$precipitation > 1 & tmp$precipitation <= 4] = '< 4'
  tmp$class[tmp$precipitation > 4 & tmp$precipitation <= 7] = '< 7' 
  tmp$class[tmp$precipitation > 7 & tmp$precipitation <= 9] = '< 9' 
  tmp$class[tmp$precipitation > 9] = '> 9'
  
  # Display precipitation per SOS
  ggplot(data=tmp,
         aes(x=t,
             y=precipitation)) +
    geom_point() +
    labs(x='Day of year',
         y='Precipitation',
         title=paste('Precipitation data square for day of year 2013 square ', s)) +
    theme_bw()
  
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/precipitation_date_square_sos_',s,'.png'))
  
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
         title=paste('Precipitation data square', s)) +
    theme_bw()
  
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/precipitation_square_sos_',s,'.png'))

}