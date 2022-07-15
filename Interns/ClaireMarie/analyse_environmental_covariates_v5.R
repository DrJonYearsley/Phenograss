# Analyse environmental covariates
# 
# cumul precipitation - soil moisture - soil type
#
# Claire-Marie Alla
# 27/07/2021
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
years = c(2017)


# if you have gls_model_2002_2019.Rdata
# this part can be removed
# ++++++++++++++++++++++++++++++++++++++++++++++

# Filename segmented data
for (i in 1:length(squaresList)) {
  
  filename = paste0(input_file_preffix,'_square_',squaresList[i],'_',years[1],'.RData')
  load(file.path(dataDir,filename))
  
  if (i==1) {
    phenology_wide = d_final
  } else {
    phenology_wide = rbind(phenology_wide,
                           d_final)
  }
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++

#load(file.path(dataDir,'gls_model_1_21_2013.Rdata'))
# to big for load in my computer

# phenology_wide$phase1_residual = residuals(gls_phase1)
# phenology_wide$phase1_fit = predict(gls_phase1)

tmp = subset(phenology_wide, phenology_wide$doy >= 0 & phenology_wide$doy <= 59)
#tmp = subset(phenology_wide, !is.na(phase1_fit) & square%in%c(1:21))

# Read in MODIS grid
modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
crs_modis = st_crs(modis)

IR = st_read('~/Stage/Data/Quadrats/country.shp')
squares = st_read('~/Stage/Data/Quadrats/agriclimate_quadrats_Ireland.shp')
squares_modis = st_transform(squares, crs_modis)  # Make sure CRS is identical (required for gstat interpolation)
crs_squares = st_crs(squares)


# Create a geometry column
tmp = st_as_sf(tmp, coords = c("x_MODIS", "y_MODIS"), crs = crs_modis)

# 1. Retrieve elevation data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# elevation = st_mosaic(read_stars(file.path(elevation_dir, 'square_1.tif')),
#                       read_stars(file.path(elevation_dir, 'square_2.tif')),
#                       read_stars(file.path(elevation_dir, 'square_3.tif')),
#                       read_stars(file.path(elevation_dir, 'square_4.tif')),
#                       read_stars(file.path(elevation_dir, 'square_5.tif')),
#                       read_stars(file.path(elevation_dir, 'square_6.tif')),
#                       read_stars(file.path(elevation_dir, 'square_7.tif')),
#                       read_stars(file.path(elevation_dir, 'square_8.tif')),
#                       read_stars(file.path(elevation_dir, 'square_9.tif')),
#                       read_stars(file.path(elevation_dir, 'square_13.tif')),
#                       read_stars(file.path(elevation_dir, 'square_14.tif')),
#                       read_stars(file.path(elevation_dir, 'square_15.tif')),
#                       read_stars(file.path(elevation_dir, 'square_16.tif')),
#                       read_stars(file.path(elevation_dir, 'square_17.tif')),
#                       read_stars(file.path(elevation_dir, 'square_18.tif')),
#                       read_stars(file.path(elevation_dir, 'square_19.tif')),
#                       read_stars(file.path(elevation_dir, 'square_20.tif')),
#                       read_stars(file.path(elevation_dir, 'square_21.tif')))
# #elevation = read_stars(file.path(elevation_dir, paste0('square_', square, '.tif')))
# st_crs(elevation) = crs_modis
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


# # 1.bis 2 : Retrieve slope data
# 
# # Save modis square - not find better to use a raster function (not available in stars)
# fname = paste0("all_squares.tif")
# if (!dir.exists(outputDir)) {
#   dir.create(outputDir)
# }
# write_stars(adrop(elevation), file.path(outputDir,fname))
# 
# # Use the field function
# crop_slope = raster(file.path(outputDir,fname))
# crop_slope = terrain(crop_slope, opt='slope', unit='degrees')


# 2. Retrieve soil data IFS
# ++++++++++++++++++++++++++++++++++++++++++++++
soil_data = st_read(file.path(soilPath,"Soils_IE_WetDry.shp"))
# Convert soil data in MODIS crs
soil_data_modis = st_transform(soil_data, crs = crs_modis)

# Make a join
tmp = st_join(st_as_sf(tmp, crs = crs_modis), soil_data_modis["CATEGORY"])


# Retrieve MERA data
filename = paste0('soilmoisture_2012_2017.RData')
load(file.path(climatePath,filename))

filename2 = paste0('temperature_degrees_2012_2013_2015.RData')
load(file.path(climatePath,filename2))

filename3 = paste0('precipitation_2012_2013_2015_2016_2017.RData')
load(file.path(climatePath,filename3))


# Function : Create variogram
variogramme <- function(covariate_mera, date1, date2) {

  mera_square = subset(covariate_mera, covariate_mera$square == s
                       & covariate_mera$validityDate > date1
                       & covariate_mera$validityDate < date2)


  # Create a geometry column
  mera_square_sf = st_as_sf(mera_square, coords = c("Longitude", "Latitude"),
                            crs = st_crs("EPSG:4326"))

  # Transform into the modis CRS
  mera_square_modis = st_transform(mera_square_sf, crs = crs_modis)

  mera_square_modis$grp = sapply(st_equals(mera_square_modis$geometry), max)
  test = mera_square_modis %>% group_by(grp) %>% summarize(Value = mean(Value,na.rm=T))

  crop_square = st_crop(modis, squares_modis[s,,])

  # Method 2:
  # Interpolate using a model variogram (try a linear variogram)

  # Look at the empirical variogram
  v_emp = variogram(Value ~ 1, data = test)

  # Fit variogram model (try linear)  use show.vgm() to display all possible models
  v_mod = fit.variogram(v_emp, model = vgm(NA,"Lin",0))

  # Now do some ordinary krigging to interpolate
  vario = gstat::krige(formula = Value ~ 1,
                       locations=test,
                       model=v_mod,
                       newdata=crop_square)


  return(vario)
}

variogramme2 <- function(covariate_mera, d) {
  
  mera_square = subset(covariate_mera, covariate_mera$square == s
                       & covariate_mera$date == d)

  # Create a geometry column
  mera_square_sf = st_as_sf(mera_square, coords = c("Longitude", "Latitude"),
                            crs = st_crs("EPSG:4326"))
  
  # Transform into the modis CRS
  mera_square_modis = st_transform(mera_square_sf, crs = crs_modis)
  
  mera_square_modis$grp = sapply(st_equals(mera_square_modis$geometry), max)
  test = mera_square_modis %>% group_by(grp) %>% summarize(Value = mean(Value,na.rm=T))
  
  crop_square = st_crop(modis, squares_modis[s,,])

  # Method 2:
  # Interpolate using a model variogram (try a linear variogram)
  
  # Look at the empirical variogram
  v_emp = variogram(Value ~ 1, data = test)

  if (length(which(v_emp$dist == 0.0)) != 0) {
    return (stars())
  }
  else {

    # Fit variogram model (try linear)  use show.vgm() to display all possible models
    v_mod = fit.variogram(v_emp, model = vgm(NA,"Lin",0))
    
    # Now do some ordinary krigging to interpolate
    vario = gstat::krige(formula = Value ~ 1,
                         locations=test,
                         model=v_mod,
                         newdata=crop_square)
    return(vario)
  }
}


# Convert dates to a date object
precipitation_mera$date = as.Date(as.character(precipitation_mera$validityDate),
                                format="%Y%m%d")

tmp$date = as.Date(as.character(tmp$date), format="%Y-%m-%d")


for (s in squaresList) {
  print(s)

  pixel_list = unique(tmp$pixelID[tmp$square == s])
  nPixel = length(pixel_list)

  # 3. Retrieve soil moisture MERA data
  # ++++++++++++++++++++++++++++++++++++++++++++++

  g_mod = variogramme(soilmoisture_mera, '20170101', '20170301')


  # 4. Retrieve temperature MERA data
  # ++++++++++++++++++++++++++++++++++++++++++++++

  #g_mod2 = variogramme(temperature_mera, '20130501', '20130901')
  
  days = unique(tmp$doy[tmp$square == s])
  days = sort(days)
  
  for (day in days) {
    
    # recuperer la date
    d = which(tmp$square == s & tmp$doy == day)[1]
    
    # 4. Retrieve precipitation MERA data
    # ++++++++++++++++++++++++++++++++++++++++++++++
    
    g_mod3 = variogramme2(precipitation_mera, tmp$date[d])

    if (length(g_mod3) == 0) {
      tmp$precipitation[tmp$square == s & tmp$doy == day] = NA
      
    } else {
      
      #tmp$precipitation[tmp$square == s & tmp$doy == day] = as.numeric(extract(g_mod3, tmp$geometry))
      #tmp$soilmoisture[tmp$square == s & tmp$doy == day] = as.numeric(extract(g_mod, tmp$geometry))
      
      for (i in pixel_list) {

        # Retrieve the index of one of the pixels
        ind = which(tmp$square == s & tmp$pixelID == i & tmp$doy == day)

        if (!is.na(ind[1])) {
          tmp$precipitation[tmp$square == s & tmp$pixelID == i &  tmp$doy == day] = as.numeric(st_extract(g_mod3, tmp$geometry[ind[1]]))
          tmp$soilmoisture[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(g_mod, tmp$geometry[ind[1]]))
        }
      
      }
    }
  }


  # tmp2 = subset(tmp, tmp$square == s)
  # pixel_list = unique(tmp2$pixelID)
  # nPixel = length(pixel_list)
  #
  # for (i in pixel_list) {
  # 
  #   # Retrieve the index of one of the pixels
  #   ind = which(tmp$square == s & tmp$pixelID == i)
  # 
  #   for (ind_elev in ind) {
  # 
  #     # # Create a elevation colomn thanks to geometry
  #     # tmp$elevation[ind_elev] = as.numeric(st_extract(elevation,
  #     #                                                 st_sfc(tmp$geometry[ind[1]],
  #     #                                                        crs = crs_modis)))
  # 
  #     tmp$aspect_slope[ind_elev] = as.numeric(st_extract(aspect_slope,
  #                                                        st_sfc(tmp$geometry[ind[1]],
  #                                                               crs = crs_modis)))
  #     # # Create a slope colomn
  #     # tmp$slope[ind_elev] = as.numeric(st_extract(st_as_stars(crop_slope),
  #     #                                             st_sfc(tmp$geometry[ind[1]],
  #     #                                                    crs = crs_modis)))
  #   }
  # 
  #   tmp$soilmoisture[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(g_mod, tmp$geometry[ind[1]]))
  #   #tmp$temperature[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(g_mod2, tmp$geometry[ind[1]]))
  #   tmp$precipitation[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(g_mod3, tmp$geometry[ind[1]]))
  # 
  # }
}




# Save tmp
#save(tmp, file=paste0(outputDir, '/gls_1_21_59days_2017_soilmoist_precip.Rdata'))
#load(file.path(outputDir,'gls_model_1_21_2012_elev_aspElev_CAT_soilmoist_precip_temp.Rdata'))
#load(file.path(outputDir,'gls_model_1_21_2013_elev_aspElev_CAT_precip_temp.Rdata'))
#load(file.path(outputDir, 'gls_model_2012_2017_precipitation_soilmoist.RData'))


# tmp$class_aspect[tmp$aspect_slope <= 20] = 'N'
# tmp$class_aspect[tmp$aspect_slope >= 340] = 'N'
# tmp$class_aspect[tmp$aspect_slope > 20 & tmp$aspect_slope < 70] = 'NE'
# tmp$class_aspect[tmp$aspect_slope >= 70 & tmp$aspect_slope <= 110] = 'E'
# tmp$class_aspect[tmp$aspect_slope > 110 & tmp$aspect_slope < 160] = 'SE'
# tmp$class_aspect[tmp$aspect_slope >= 160 & tmp$aspect_slope <= 200] = 'S'
# tmp$class_aspect[tmp$aspect_slope > 200 & tmp$aspect_slope < 250] = 'SW'
# tmp$class_aspect[tmp$aspect_slope >= 250 & tmp$aspect_slope <= 290] = 'W'
# tmp$class_aspect[tmp$aspect_slope > 290 & tmp$aspect_slope < 340] = 'NW'


# Average over the days
tmp = tmp %>% group_by(doy, pixelID, square) %>% summarize(precipitation = mean(precipitation,na.rm=T),
                                                           soilmoisture = mean(soilmoisture,na.rm=T),
                                                           CATEGORY = first(CATEGORY))

tmp = tmp %>% group_by(pixelID, square) %>% summarize(precip_cumulate = sum(precipitation,na.rm=T),
                                                      soilmoisture = mean(soilmoisture,na.rm=T),
                                                      CATEGORY = first(CATEGORY))

# 
# n = nrow(tmp2)
# tmp$cum_precip[tmp$square == s] = tmp2$precip_cumulate[n]
# 
# days = unique(tmp$doy)
# days = sort(days)
# 
# c = 0 # counter
# # Accumulate ascending when the precipitation is above 0, 5 and 10 per days
# for (day in days) {
#   
#   if (tmp$precipitation[tmp$doy == day] > 0) {
#     
#     value = tmp$precipitation[tmp$doy == day] - 0
#     
#     tmp$precip_cumulate[tmp$doy == day] = value + c
#     c = tmp$precip_cumulate[tmp$doy == day]
#     
#   } else {
#     tmp$precip_cumulate[tmp$doy == day] = c
#   }
# }


# Calculate some R-squared between residus and environmentale covariates
m = lm(soilmoisture~1 + precip_cumulate, data=tmp)
#m = lm(precipitation~1 + CATEGORY, data=tmp) # 2017

#m = lm(phase1_residual~1 + precipitation, data=tmp) # for residuals of 21 squares (2013)

#m = lm(temperature~1 + class_aspect, data=tmp) # 2012 et 2013
#m = lm(temperature~1 + aspect_slope, data=tmp) # 2012


summary(m)

predict_lm <- predict(m, tmp)
#predict_lm <- data.frame(temperature_Pred = predict_lm, temperature = tmp$temperature)
predict_lm <- data.frame(soilmoisture_Pred = predict_lm, soilmoisture = tmp$soilmoisture, square = tmp$square)


tmp = subset(tmp, tmp$CATEGORY == 'Peat'
             | tmp$CATEGORY == 'Well Drained'
             | tmp$CATEGORY == 'Poorly Drained')

#colour=factor(square)
ggplot(data=tmp,
       aes(x=precip_cumulate,
           y=soilmoisture)) +
  geom_point() +
  geom_point(aes(colour = CATEGORY)) +
  #geom_boxplot() +
  #geom_density() +
  #facet_grid(class_aspect ~ .) +
  #scale_colour_viridis_d(option = "inferno") +
  labs(x='Cumulated precipitations (mm)',
       y='Soil moisture',
       title = 'Distribution of cumulated precipitations and soil moisture in 2 months in 2017') +
  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18))


ggsave(filename = 'gls5_fitted_vs_cumulprecip_over0_soilmoisture_year2017.png', width=30 , height=20, units='cm')

# Display difference between actual and predicted
plot(predict_lm$soilmoisture_Pred - predict_lm$soilmoisture, main = "Difference between actual and predicted soil moisture. lm")
