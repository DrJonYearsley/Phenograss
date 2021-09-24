# Slope aspect of SRTM elevation
#
# Claire-Marie Alla
# 29/06/2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd("~/MEGAsync/Projects/GrasslandPhenology/Admin/SteeringCommittee/Sept2021/")

library(sf)
library(stars)
library(terra)
library(ggplot2)

elevation_dir = '/Volumes/MODIS_data/Elevation_SRTM_squares/'
# modisPath = '~/Stage/Data/MODIS'
# dataDir = '~/Stage/Data/MODIS/Phenophase_estimates'
# outputDir = '~/Stage/Data_created'

# input_file_preffix = 'phenology'

# Import data --------
squaresList = c(18)


# 
# # Method 1
# for (e in 1:length(squaresList)) {
#   filename = paste0('square_',squaresList[e],'.tif')
#   square = read_stars(file.path(elevation_dir, filename))
#   
#   filename2 = paste0(input_file_preffix,'_square_',squaresList[e],'_',year,'.RData')
#   load(file.path(dataDir,filename2))
#   
#   c = st_contour(square, breaks=c(100, 150))
#   print(c)
#   
#   # Create the gradient from the outline ? 
#   # But we don't have the expression of the fonction
# }

# Method 1 (JY version)
for (e in 1:length(squaresList)) {
  filename = paste0('square_',squaresList[e],'.tif')
  elevation = rast(file.path(elevation_dir, filename))
  slope = terrain(elevation, v="slope", unit="degrees", neighbors=8)
  aspect = terrain(elevation, v="aspect", unit="degrees", neighbors=8)
}


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
