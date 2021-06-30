# Slope aspect of SRTM elevation
#
# Claire-Marie Alla
# 29/06/2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(sf)
library(stars)
library(raster)
library(ggplot2)
library(dplyr)
library(biwavelet)
library(grid)

elevation_dir = '~/Stage/Data_created/elevation_SRTM3_square'
modisPath = '~/Stage/Data/MODIS'
dataDir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created'

input_file_preffix = 'phenology'

# Import data --------
squaresList = c(20)
year = 2019


# Method 1
for (e in 1:length(squaresList)) {
  filename = paste0('square_',squaresList[e],'.tif')
  square = read_stars(file.path(elevation_dir, filename))
  
  filename2 = paste0(input_file_preffix,'_square_',squaresList[e],'_',year,'.RData')
  load(file.path(dataDir,filename2))
  
  c = st_contour(square, breaks=c(100, 150))
  print(c)
  
  # Create the gradient from the outline ? 
  # But we don't have the expression of the fonction
}





# Method 2
# Sobel filter : calculate for each point the gradient, the direction and
# the norm of the gradient of the raster

square2 = raster(file.path(elevation_dir, filename))
square2 = as.matrix(square2)

# Create the 2 gradient filter
deriv_hor = matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), nrow = 3)
deriv_vert = matrix(c(1, 0, -1, 2, 0, -2, 1, 0, -1), nrow = 3)

# Replace NA by 0
#square2[is.na(square2)] = 0

print(class(square2))
print(length(square2)) 

Gx = convolve2D(square2, deriv_hor, type="open")
Gy = convolve2D(square2, deriv_vert, type="open")

# Norm of the gradient
norm_grad = sqrt(Gx**2 + Gy**2)

# Direction of the gradient
#arctan(Gy/Gx)
direction_grad = atan2(Gx, Gy)
direction_grad = direction_grad * (360/pi)

# doesn't work
direction_grad[direction_grad < 0] = - direction_grad

print(direction_grad)


plot(c, reset = FALSE)
#contour(square_1, add = TRUE) # contour plot


# Read in MODIS grid
# modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
# crs_modis = st_crs(modis)
# 
# #direction_grad = st_as_stars(direction_grad, crs = crs_modis)
# 
# #norm_grad = st_as_stars(norm_grad, crs = crs_modis)
# 
# output_smoothed = st_as_sf(output_smoothed, coords = c("x_MODIS", "y_MODIS"),
#                            crs = crs_modis)
# 
# output_smoothed$grp = sapply(st_equals(output_smoothed$geometry), max)
# output_smoothed = output_smoothed %>% group_by(grp, pixelID) %>% summarize(t = mean(t))


# Create a dataframe with coord / grad dir et grad norme
slope_aspect = data.frame()
c = 0
for (i in 1:44) {
  for (j in 1:44) {
    print(i)
    print(j)
    print(square2[[i,j]])
    slope_aspect$geometry[c] = square2[[i,j]]
    c = c +1
  }
}

# if between 45 and 135 degrees => North
# if between 135 and 225 degrees => West
# if between 225 and 315 degrées => South
# if between 315 and 360, 0 and 45 => East

# display with ggpolt2 + geom_segment and arrow (coord of origin point, norm and direction)
#ggplot(data=as.data.frame(direction_grad))


