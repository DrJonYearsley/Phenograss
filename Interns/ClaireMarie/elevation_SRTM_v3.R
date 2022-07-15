# 
# Display elevation and slope pasture for each square
# algorithm optimazed in improvement_model.R
#
# Claire-Marie
# 10/06/2021
#
#1. Import squares
#2. Import SRTM3 mosaic data
#3. Put squares in the same reference WGS84
#4. Cut modis data by squares
#5. Resample SRTM3 data with modis squares data
#6. Extract data pasture
#7. Display the number of pixels of same elevation
#8. Display slope
#8. Save data
#
# Adapted by Jon Yearsley
# 28th Sept 2021
#  Outputs elevation for all squares across all of Ireland
#  Only produces graphics on request (plot_on=TRUE)
# ****************************************


rm(list=ls())

library(sf) # combine sp, rgeos et rgdal
library(stars)
library(ggplot2)
library(raster)


# datadir = '~/Stage/Data/MODIS/Phenophase_estimates'
# quadratPath = '~/Stage/Data/Quadrats'
# outputDir = '~/Stage/Data_created/elevation_SRTM3_square'
# 
# SRTM3Path = '~/Stage/Data/DigitalElevation/SRTM_data/hgt'

#datadir = '~/Stage/Data/MODIS/Phenophase_estimates'
modis_dir = '/Volumes/MODIS_data/MODIS/'
quadratPath = '/Volumes/MODIS_data/Quadrats/'
outputDir = '/Volumes/MODIS_data/Elevation_SRTM_squares/'

SRTM3Path = '/Volumes/MODIS_data/Elevation_SRTM_squares/hgt/'

plot_on = FALSE


# input_file_prefix = 'phenology'

setwd(SRTM3Path)

# List of squares to analyse
squareList = c(10:12)
# year = 2019


# SRTM3 data : WGS84
files = list.files(pattern='.hgt')
SRTM = st_mosaic(files)

SRTM3 = st_mosaic(read_stars("N51W009.hgt", NA_value = 0), 
                  read_stars("N51W010.hgt", NA_value = 0),
                  read_stars("N52W007.hgt", NA_value = 0), 
                  read_stars("N52W008.hgt", NA_value = 0),
                  read_stars("N52W009.hgt", NA_value = 0), 
                  read_stars("N52W010.hgt", NA_value = 0),
                  read_stars("N53W007.hgt", NA_value = 0), 
                  read_stars("N53W008.hgt", NA_value = 0),
                  read_stars("N53W009.hgt", NA_value = 0), 
                  read_stars("N53W010.hgt", NA_value = 0),
                  read_stars("N54W007.hgt", NA_value = 0), 
                  read_stars("N54W008.hgt", NA_value = 0),
                  read_stars("N54W009.hgt", NA_value = 0), 
                  read_stars("N55W007.hgt", NA_value = 0),
                  read_stars("N55W008.hgt", NA_value = 0))

crs_SRTM3 = st_crs(SRTM3)

# Import squares from shapefile
squares = st_read(file.path(quadratPath,"agriclimate_quadrats_Ireland.shp"))
crs_squares = st_crs(squares)

# Convert squares to gtopo crs (because cropping is faster on a regular grid)
squares_gtopo = st_transform(squares, crs = crs_SRTM3)

# Read in MODIS grid
modis = read_stars(file.path(modis_dir,'modis_grid_ireland.tif'))
st_crs(modis) =  crs_squares  # Make sure modis and squares have same CRS

# nb_pixels = data.frame()

# Elevation data by squares
for (s in squareList) {
  
  # # Import phenophase estimate
  # filename = paste0(input_file_prefix,'_square_',s,'_',year,'.RData')
  # load(file.path(datadir,filename))
  # 
  # pixel_list = unique(d_final$pixelID)
  # nPixel = length(pixel_list)
  
  # Crop modis grid to a square and then resample elevation onto this cropped grid
  crop_square = st_crop(modis, squares[s,])
  crop_modis = st_warp(SRTM3, crop_square, method="cubic", use_gdal=T)
  
  
  
  # Save modis square - not find better to use a raster function (not available in stars)
  fname = paste0("square_",s,".tif")
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  write_stars(adrop(crop_modis), file.path(outputDir,fname))
  
  
  if (plot_on) {
    # Visualise some of the data ----------
    # Use the field function
    crop_slope = raster(file.path(outputDir,fname))
    crop_slope = terrain(crop_slope, opt='slope', unit='degrees')
    
    
    tmp = data.frame()
    
    for (i in pixel_list) {
      print(i)
      
      # Retrieve the index of one of the pixels
      ind = which(d_final$pixelID == i)[1]
      
      # Create a modis geometry column
      d_final$geometry[ind] = st_sfc(st_point(c(d_final$x_MODIS[ind], 
                                                d_final$y_MODIS[ind])))
      # Create a elevation colomn thanks to geometry
      d_final$elevation[ind] = as.numeric(st_extract(crop_modis, 
                                                     st_sfc(d_final$geometry[ind], 
                                                            crs = st_crs(modis))))
      # Create a slope colomn
      d_final$slope[ind] = as.numeric(st_extract(st_as_stars(crop_slope), 
                                                 st_sfc(d_final$geometry[ind], 
                                                        crs = st_crs(modis))))
      
      tmp = rbind(tmp, d_final[ind,])
      nb_pixels = rbind(nb_pixels, d_final[ind,])
      
    }
    
    # Min, max elevations values
    min = min(tmp$elevation, na.rm=T)
    max = max(tmp$elevation, na.rm=T)
    
    # Display elevation
    ggplot(data=tmp,
           aes(x=x_MODIS,
               y=y_MODIS,
               fill=elevation)) +
      geom_tile(colour='darkblue') +
      coord_equal() +
      scale_fill_gradientn(colours = terrain.colors(20),
                           na.value='darkgray') +
      labs(x='X Coord (MODIS CRS)',
           y='Y Coord (MODIS CRS)',
           title=paste('Elevation SRTM3 square', s),
           subtitle = paste('min = ', min, 'm max = ', max, 'm')) +
      theme_bw()
    
    # Save file as png
    #ggsave(width=11, height=6,
    #       filename = paste0(outputDir, '/elevation_SRTM3_square_',s,'.png'))
    
    
    # Min, max slopes values
    mins = min(tmp$slope, na.rm=T)
    maxs = max(tmp$slope, na.rm=T)
    
    # Display slope
    ggplot(data=tmp,
           aes(x=x_MODIS,
               y=y_MODIS,
               fill=slope)) +
      geom_tile(colour='darkblue') +
      coord_equal() +
      scale_fill_gradientn(colours = terrain.colors(20),
                           na.value='darkgray') +
      labs(x='X Coord (MODIS CRS)',
           y='Y Coord (MODIS CRS)',
           title=paste('Slope square', s),
           subtitle = paste('min = ', mins, '? max = ', maxs, '?')) +
      theme_bw()
    
    # Save file as png
    ggsave(width=11, height=6,
           filename = paste0(outputDir, '/slope_elevation_square_',s,'.png'))
    
    
    
    # Display the number of pixels of same elevation
    ggplot(data=tmp, aes(x = elevation)) +
      geom_histogram(bins=30) +
      labs(x='Elevation',
           y = 'Number of Pixels',
           title=paste('Elevation square ',s))
    
    #ggsave(width=11, height=6,
    #       filename = paste0(outputDir, '/nb_pixels_square_',s,'.png'))
  }
}

# Display the number total of pixels of same elevation
ggplot(data=nb_pixels, aes(x = elevation)) +
  geom_histogram(bins=30) +
  labs(x='Elevation',
       y = 'Number of Pixels',
       title=paste('Elevation square ',s))

ggsave(width=11, height=6,
       filename = paste0(outputDir, '/nb_pixels_squares.png'))