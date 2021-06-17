# 
# Display of the soil moisture data from MERA for each square
#
# Claire-Marie
# 11/06/2021
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
# soilPath = '~/Stage/Data/Climate/SoilMoisture_quadrats'
# quadratPath = '~/Stage/Data/Quadrats'
# 
# datadir = '~/Stage/Data/MODIS/Phenophase_estimates'
# outputDir = '~/Stage/Data_created/soil_moisture_data'


soilPath = '/Volumes/MODIS_data/MERA/DailyData_subsetted/'
quadratPath = '/Volumes/MODIS_data/Quadrats/'
modisPath = '/Volumes/MODIS_data/MODIS/'

datadir = '~/Research/Phenograss/Data/PhenologyOutput/'
outputDir = '~/Research/Phenograss/Data/MERA_processed/'


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

    # Transform in stars class => make an interpolation
  crop_square = st_crop(modis, squares_modis[s,,])
  
  # Inverse distance weighted interpolation of the MERA data
  g = gstat::idw(formula = DailyMean ~ 1, 
                 locations=subset(mera_square_modis, ValidityDate==20170101), 
                 newdata=crop_square, 
                 idp=2)
  # 
  # ggplot() +  
  #   geom_sf(data=squares[s,]) +
  #   geom_sf(data=d_modis) 
  
  # +++++++++++++++++++++++++++++++++++++++++++++
  
  # Group pixels and sum daily
  mera_square$grp = sapply(st_equals(mera_square$geometry), max)
  mera_square_grp = mera_square %>% 
    group_by(grp) %>% 
    summarize(sum_DailyMean = sum(DailyMean))

  
  
  
  mera_square_r = st_rasterize(mera_square_grp)
  mera_square_r = st_warp(mera_square_r, crop_square, method="cubic", use_gdal=T) # ---- doesn't work solution 1
  
  #crop_square[is.na(crop_square)] <- 0
  #print(crop_square[[1]])
  #plot(mera_square)
  
  
  ggplot() + 
    geom_stars(data=mera_square_r) +
    coord_equal() + 
    scale_fill_gradientn(colours = terrain.colors(20)) +
    labs(x='X Coord (MODIS CRS)',
         y='Y Coord (MODIS CRS)',
         title=paste('test', s)) +
    theme_bw()
  
  # Save data for visualisation
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/test_', s, '.png'))

  
  g = gstat(formula = DailyMean ~ 1, data = mera_square_grp)
  print('ok')
  print(st_crs(mera_square_grp))
  print(st_crs(crop_square))
  print(!identical(st_crs(crop_square), st_crs(mera_square_grp)))
  
  z = predict(g, crop_square) # ------------------------------------------------------- doesn t work solution 2
  print(z["var1.pred"])
  print('ok1')
  #names(z) = "DailyMean"
  
  b = seq(0, 1200, 100)
  plot(z, breaks = b, col = hcl.colors(length(b)-1, "Spectral"), reset = FALSE)
  #plot(st_geometry(mera_square), pch = 3, add = TRUE)
  #contour(z, breaks = b, add = TRUE)
  
  # solution 3 : cr?er un grille d'interpolation ?
  

  
  # Create a geometry column
  d_final = st_as_sf(d_final, coords = c("x_MODIS", "y_MODIS"), crs = crs_modis)
  
  pixel_list = unique(d_final$pixelID)
  nPixel = length(pixel_list)
  
  tmp = data.frame()

  for (i in pixel_list) {
    print(i)
    
    # Retrieve the index of one of the pixels
    ind = which(d_final$pixelID == i)[1]

    # Create a soilmoisture colomn thanks to geometry
    #print(mera_square[[1]])
    d_final$soilmoisture[ind] = as.numeric(st_extract(mera_square_r, 
                                                      d_final$geometry[ind]))
    
    tmp = rbind(tmp, d_final[ind,])
  }
  
  
  # average over every day by pixels ?
  # graph per day on the evolution ?
  
  # For viewing
  tmp$x_MODIS = mera_square$geometry[[1]]
  tmp$y_MODIS = mera_square$geometry[[1]][[2]]
  
  # Display
  ggplot(data=tmp,
         aes(x=x_MODIS,
             y=y_MODIS,
             fill= as.factor(soilmoisture))) +
    geom_tile(colour='darkblue') +
    coord_equal() +
   scale_fill_brewer(palette='Blues', na.value='darkgray') +
    labs(x='X Coord (MODIS CRS)',
         y='Y Coord (MODIS CRS)',
         title=paste('Soil data square', s)) +
    theme_bw()
  
  # Save file as png
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
   
  ggsave(width=11, height=6,
         filename = paste0(outputDir, '/soil_moisture_square_',s,'.png'))
}