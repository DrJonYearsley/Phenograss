# 
# Distribution of the temperature anomaly with the slope aspect
#
# 1. Retrieve slope aspect
# 2. Retrieve slope
# 3. Retrieve temperature data
# 4. Calcul temperature anomaly
# 5. separate high and slow slope
#
# Claire-Marie Alla
# 21/07/2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(sf)
library(stars)
library(segmented)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)
library(gstat)
library(dplyr)
library(raster)


elevation_dir = '~/Stage/Data_created/elevation_SRTM3_square'
aspect_dir = '~/Stage/Data_created/elevation_SRTM3_square/exposition_squares'
modisPath = '~/Stage/Data/MODIS'
soilPath = '~/Stage/Data/soil_data/Soils_IE_WetDry'
climatePath = '~/Stage/Data/Climate'

dataDir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created'

input_file_preffix = 'phenology'

# Import data --------
squaresList = c(1:9, 13:21)
years = c(2015)

# # Filename segmented data
# for (i in 1:length(squaresList)) {
#   
#   filename = paste0(input_file_preffix,'_square_',squaresList[i],'_',years[1],'.RData')
#   load(file.path(dataDir,filename))
#   
#   if (i==1) {
#     phenology = d_final
#   } else {
#     phenology = rbind(phenology,
#                       d_final)
#   }
# }
# 
# # Read in MODIS grid
# modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
# crs_modis = st_crs(modis)
# 
# IR = st_read('~/Stage/Data/Quadrats/country.shp')
# squares = st_read('~/Stage/Data/Quadrats/agriclimate_quadrats_Ireland.shp')
# squares_modis = st_transform(squares, crs_modis)  # Make sure CRS is identical (required for gstat interpolation)
# crs_squares = st_crs(squares)
# 
# 
# # Create a geometry column
# phenology = st_as_sf(phenology, coords = c("x_MODIS", "y_MODIS"), crs = crs_modis)
# 
# phenology$grp = sapply(st_equals(phenology$geometry), max)
# tmp = phenology %>% group_by(grp) %>% summarize(pixelID = first(pixelID), square = first(square))
# 
# # 1.bis Retrieve aspect slope data
# aspect_slope = st_mosaic(read_stars(file.path(aspect_dir, 'img_exposition_square1.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square2.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square3.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square4.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square5.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square6.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square7.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square8.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square9.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square13.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square14.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square15.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square16.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square17.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square18.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square19.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square20.tif')),
#                          read_stars(file.path(aspect_dir, 'img_exposition_square21.tif')))
# st_crs(aspect_slope) = crs_modis
# 
# # 1.bis 2 : Retrieve slope data
# 
# # Save modis square - not find better to use a raster function (not available in stars)
# fname = paste0("all_squares.tif")
# # if (!dir.exists(outputDir)) {
# #   dir.create(outputDir)
# # }
# # write_stars(adrop(elevation), file.path(outputDir,fname))
# 
# # Use the field function
# crop_slope = raster(file.path(outputDir,fname))
# crop_slope = terrain(crop_slope, opt='slope', unit='degrees')
# 
# 
# 
# filename2 = paste0('temperature_degrees_2012_2013_2015.RData')
# load(file.path(climatePath,filename2))
# 
# temperature_mera = subset(temperature_mera, temperature_mera$validityDate >= '20150101'
#                           & temperature_mera$validityDate < '20160101')
# 
# 
# # Function : Create variogram
# variogramme <- function(covariate_mera, date1, date2) {
# 
#   mera_square = subset(covariate_mera, covariate_mera$square == s
#                        & covariate_mera$validityDate > date1
#                        & covariate_mera$validityDate < date2)
# 
# 
#   # Create a geometry column
#   mera_square_sf = st_as_sf(mera_square, coords = c("Longitude", "Latitude"),
#                             crs = st_crs("EPSG:4326"))
# 
#   # Transform into the modis CRS
#   mera_square_modis = st_transform(mera_square_sf, crs = crs_modis)
# 
#   mera_square_modis$grp = sapply(st_equals(mera_square_modis$geometry), max)
#   test = mera_square_modis %>% group_by(grp) %>% summarize(Value = mean(Value,na.rm=T))
# 
#   crop_square = st_crop(modis, squares_modis[s,,])
# 
#   # Method 2:
#   # Interpolate using a model variogram (try a linear variogram)
# 
#   # Look at the empirical variogram
#   v_emp = variogram(Value ~ 1, data = test)
# 
#   # Fit variogram model (try linear)  use show.vgm() to display all possible models
#   v_mod = fit.variogram(v_emp, model = vgm(NA,"Lin",0))
# 
#   # Now do some ordinary krigging to interpolate
#   vario = gstat::krige(formula = Value ~ 1,
#                        locations=test,
#                        model=v_mod,
#                        newdata=crop_square)
# 
# 
#   return(vario)
# }
# 
# 
# 
# 
# for (s in squaresList) {
#   print(s)
#   # 4. Retrieve temperature MERA data
#   # ++++++++++++++++++++++++++++++++++++++++++++++
# 
#   g_mod2 = variogramme(temperature_mera, '20150101', '20150301')
# 
# 
#   tmp2 = subset(tmp, tmp$square == s)
#   pixel_list = unique(tmp2$pixelID)
#   nPixel = length(pixel_list)
# 
#   for (i in pixel_list) {
#     
#     # Retrieve the index of one of the pixels
#     ind = which(tmp$square == s & tmp$pixelID == i)
#     
#     tmp$aspect_slope[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(aspect_slope,
#                                                                                  st_sfc(tmp$geometry[ind[1]],
#                                                                                         crs = crs_modis)))
#     
#     # Create a slope colomn
#     tmp$slope[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(st_as_stars(crop_slope),
#                                                                           st_sfc(tmp$geometry[ind[1]],
#                                                                                  crs = crs_modis)))
#   
# 
#     tmp$temperature[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(g_mod2, tmp$geometry[ind[1]]))
# 
#   }
# 
#   temperature_mean = mean(tmp$temperature,na.rm=T)
#   tmp$anomaly = tmp$temperature - temperature_mean
# }
# 
# 
# tmp$class_aspect[tmp$aspect_slope <= 20] = 'N'
# tmp$class_aspect[tmp$aspect_slope >= 340] = 'N'
# tmp$class_aspect[tmp$aspect_slope > 20 & tmp$aspect_slope < 70] = 'NE'
# tmp$class_aspect[tmp$aspect_slope >= 70 & tmp$aspect_slope <= 110] = 'E'
# tmp$class_aspect[tmp$aspect_slope > 110 & tmp$aspect_slope < 160] = 'SE'
# tmp$class_aspect[tmp$aspect_slope >= 160 & tmp$aspect_slope <= 200] = 'S'
# tmp$class_aspect[tmp$aspect_slope > 200 & tmp$aspect_slope < 250] = 'SW'
# tmp$class_aspect[tmp$aspect_slope >= 250 & tmp$aspect_slope <= 290] = 'W'
# tmp$class_aspect[tmp$aspect_slope > 290 & tmp$aspect_slope < 340] = 'NW'
# 

#save(tmp, file=paste0('anomalyTemp_slope_asp_square1_to_21_in_2015.Rdata'))

load('anomalyTemp_slope_asp_square1_to_21_in_2015.Rdata')

low_slope = subset(tmp, tmp$slope <= 2.5 & tmp$class_aspect == 'N' & tmp$class_aspect == 'S')
high_slope = subset(tmp, tmp$slope > 2.5 & tmp$class_aspect == 'N' & tmp$class_aspect == 'S')

low_slope$class_aspect = relevel(factor(low_slope$class_aspect), ref="S")
high_slope$class_aspect = relevel(factor(high_slope$class_aspect), ref="S")


# Calculate some R-squared between residus and environmentale covariates
m = lm(anomaly~1 + class_aspect, data=low_slope)

summary(m)

predict_lm <- predict(m, low_slope)
predict_lm <- data.frame(anomaly_Pred = predict_lm, anomaly = low_slope$anomaly)


#colour=factor(square)
ggplot(data=low_slope,
       aes(x=anomaly)) +
  #geom_point() +
  #geom_point(aes(colour = CATEGORY)) +
  #geom_boxplot() +
  geom_density() +
  facet_grid(class_aspect ~ .) + # horizontal
  labs(x='Temperature anomaly',
       y='Density of pixels',
       title='Distribution of temperature anomaly January, February 2015 low slope (<= 2.5°)') +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16))


ggsave(filename = 'temp_slope_asp_lowSlope_year2015.png', width=20 , height=20, units='cm')

ggplot(data=low_slope,
       aes(x=anomaly)) +
  #geom_point() +
  #geom_point(aes(colour = CATEGORY)) +
  geom_boxplot() +
  #geom_density() +
  facet_grid(class_aspect ~ .) + # horizontal
  labs(x='Temperature anomaly',
       y='Density of pixels',
       title='Distribution of temperature anomaly January, February 2015 low slope (<= 2.5°)') +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))


ggsave(filename = 'temp_slope_asp_lowSlope_year2015box.png', width=20 , height=20, units='cm')


ggplot(data=predict_lm,
       aes(x=anomaly, 
           y=anomaly_Pred)) +
  geom_point() +
  labs(x='Temperature anomaly',
       y='Temperature anomaly predicted',
       title='Difference between actual and predicted temp anomaly low slope. lm') +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))


ggsave(filename = 'temp_slope_asp_lowSlope_predictions_year2015.png', width=30 , height=20, units='cm')







# Calculate some R-squared between residus and environmentale covariates
m = lm(anomaly~1 + class_aspect, data=high_slope)

summary(m)

predict_lm <- predict(m, high_slope)
predict_lm <- data.frame(anomaly_Pred = predict_lm, anomaly = high_slope$anomaly)


#colour=factor(square)
ggplot(data=high_slope,
       aes(x=anomaly)) +
  #geom_point() +
  #geom_point(aes(colour = CATEGORY)) +
  #geom_boxplot() +
  geom_density() +
  facet_grid(class_aspect ~ .) + # horizontal
  labs(x='Temperature anomaly',
       y='Density of pixels',
       title='Distribution of temperature anomaly January, February 2015 high slope (> 2.5°)') +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))


ggsave(filename = 'temp_slope_asp_highSlope_year2015.png', width=20 , height=20, units='cm')

ggplot(data=high_slope,
       aes(x=anomaly)) +
  #geom_point() +
  #geom_point(aes(colour = CATEGORY)) +
  geom_boxplot() +
  #geom_density() +
  facet_grid(class_aspect ~ .) + # horizontal
  labs(x='Temperature anomaly',
       y='Density of pixels',
       title='Distribution of temperature anomaly January, February 2015 high slope (> 2.5°)') +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))


ggsave(filename = 'temp_slope_asp_highSlope_year2015box.png', width=20 , height=20, units='cm')

ggplot(data=predict_lm,
       aes(x=anomaly,
           y=anomaly_Pred)) +
  geom_point() +
  labs(x='Temperature anomaly',
       y='Temperature anomaly predicted',
       title='Difference between actual and predicted temp anomaly high slope. lm') +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))


ggsave(filename = 'temp_slope_asp_highSlope_predictions_year2015.png', width=30 , height=20, units='cm')



#low_slope$class_aspect = relevel(low_slope$class_aspect, ref='S') Erreur dans relevel.default(low_slope$class_aspect, ref = "S") : 'relevel' only for (unordered) factors