# Put all slope, aspect, soil data together into one data frame
#
# Claire-Marie Alla  & Jon Yearsley (Jon.Yearsley@ucd.ie)
# Started 29/06/2021
# Last edited 27/09/2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd("~/MEGAsync/Projects/GrasslandPhenology/Admin/SteeringCommittee/Sept2021/")

library(sf)
library(stars)
library(terra)
library(ggplot2)

elevation_dir = '/Volumes/MODIS_data/Elevation_SRTM_squares/'
soil_file =  '/Volumes/MODIS_data/Soils_IE_WetDry/Soils_IE_WetDry.shp'

outputDir = '/Volumes/MODIS_data/Data_created/'
squaresList = c(1:21) 
plot_on=FALSE


# Import data --------


# Read soil data as a SpatVector object (terra package)
soil = terra::vect(x=soil_file)

# Create a categorical raster layer for the soil types
soil_cat = unique(soil$CATEGORY)
cat_df = data.frame(ID = c(1:length(soil_cat)-1), 
                    CATEGORY = soil_cat)
soil$ID = cat_df$ID[match(soil$CATEGORY, cat_df$CATEGORY)]

# Reproject onto MODIS CRS
tmp = rast(file.path(elevation_dir, 'square_1.tif'))
soil_modis = project(soil, tmp)

# Method 1 (JY version)
for (e in 1:length(squaresList)) {
  filename = paste0('square_',squaresList[e],'.tif')
  elevation = rast(file.path(elevation_dir, filename))
  names(elevation) = "elevation"
  slope = terrain(elevation, v="slope", unit="degrees", neighbors=8)
  aspect = terrain(elevation, v="aspect", unit="degrees", neighbors=8)
  
  if (squaresList[e]%in%c(10:12)) { 
    # No soil data for Northern Ireland
    square_rast = c(elevation, slope, aspect)
    square_rast$SOIL_TYPE = -9
    
  } else {
    tmp = crop(soil_modis, elevation)
    soil_rast = rasterize(tmp, elevation, field="ID", fun="last")
    levels(soil_rast) = cat_df$CATEGORY
    names(soil_rast) = "SOIL_TYPE"

    square_rast = c(elevation, slope, aspect, soil_rast)
  }
  square_rast$SQUARE_ID = squaresList[e]
  
  if (e==1) {
    all_squares = square_rast
  } else {
    all_squares = merge(all_squares, square_rast)
  }
  
  writeRaster(square_rast, filename=file.path(outputDir,filename), overwrite=TRUE)
}


# Write all squares in one raster file
all_squares$SOIL_TYPE[all_squares$SOIL_TYPE==-9] = NA
writeRaster(all_squares, filename=file.path(outputDir,'combined_data.tif'), overwrite=TRUE)

# Write data as a data frame
all_squares_df = as.data.frame(all_squares, xy=TRUE)
names(all_squares_df) = c('x_MODIS','y_MODIS','ELEVATION','SLOPE','ASPECT','SOIL_TYPE','SQUARE_ID')
all_squares_df$SOIL_TYPE = as.factor(cat_df$CATEGORY[match(all_squares_df$SOIL_TYPE, cat_df$ID)])
all_squares_df$SQUARE_ID = as.factor(all_squares_df$SQUARE_ID)
all_squares_df$SLOPE_PERCENT = slope_percent = atan(all_squares_df$SLOPE*pi/180)*100

save(all_squares_df, file = file.path(outputDir,'all_envData.RData'))



# Plot some of the data ------------
if (plot_on) {
  # PLot elevation
  ggplot() + 
    geom_stars(data=st_as_stars(elevation)) +
    coord_equal() + 
    theme_void() +
    scale_fill_gradientn("Elevation (m)",colours = terrain.colors(10))
  ggsave("elevation_square18.png", width=6, height=5, units="cm")
  
  
  slope_percent = atan(slope*pi/180)*100
  # Plot slope
  ggplot() + 
    geom_stars(data=st_as_stars(slope_percent)) +
    coord_equal() + 
    theme_void() +
    scale_fill_viridis_c("Slope", 
                         option = "magma",
                         breaks=c(2,4,6,8), 
                         labels=c("2%","4%","6%","8%"))
  ggsave("slope_square18.png", width=6, height=5, units="cm")
  
  
  
  tmp = aspect
  tmp[slope_percent<1] = NA
  
  # Plot elevation
  ggplot() + 
    geom_stars(data=st_as_stars(tmp)) +
    coord_equal() + 
    scale_fill_gradient2("Aspect\n(slope>1%)", 
                         low = 'purple4', 
                         mid = "yellow", 
                         high = 'purple4', 
                         midpoint=180,
                         breaks=c(1,90,180,270,359),
                         labels = c("N","E","S","W","N")) +
    theme_void() 
  
  ggsave("aspect_square18.png", width=6, height=5, units="cm")
  
}
