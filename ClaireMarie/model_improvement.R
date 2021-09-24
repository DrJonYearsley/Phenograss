# 
# Creation  of the dataset of all covariates for the improvement model
#
# Claire-Marie Alla
# 15/07/2021 - 30/07/2021
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
library(rgdal)
library(emmeans)
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
#squaresList = c(1:9, 13:21)
# 2 first months are missing in temperature data
#years = c(2013:2017)

# Test 1
# squaresList = c(20)
# years = c(2013:2017)

# Test 2
squaresList = c(1:9, 13:21)
years = c(2013)



# Filename segmented data
for (i in 1:length(squaresList)) {
  for (y in 1:length(years)) {
    filename = paste0(input_file_preffix,'_square_',squaresList[i],'_',years[y],'.RData')
    load(file.path(dataDir,filename))

    output_smoothed$square = squaresList[i]

    if (y==1 & i==1) {
      phenology = output_smoothed
      d_pheno = d_final
    } else {
      phenology = rbind(phenology,
                        output_smoothed)
      d_pheno = rbind(d_pheno,
                      d_final)
    }
  }
}
phenology = subset(phenology, phenology$phase == 1)

# Create a wide version of phenology
phenology_wide = pivot_wider(data=subset(phenology, warning==FALSE),
                             id_cols = c('pixelID','year', 'x_MODIS','y_MODIS','x_ITM','y_ITM','square'),
                             names_from = 'phase',
                             names_prefix = 'phase',
                             values_from = c(t,slope),
                             values_fn = mean)


# Centre the x and y coordinates
phenology_wide$x_ITM_centre = phenology_wide$x_ITM - mean(phenology_wide$x_ITM, na.rm=TRUE)
phenology_wide$y_ITM_centre = phenology_wide$y_ITM - mean(phenology_wide$y_ITM, na.rm=TRUE)

rm(list='phenology')


# Add environmental covariates
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Read in MODIS grid
modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
crs_modis = st_crs(modis)

IR = st_read('~/Stage/Data/Quadrats/country.shp')
squares = st_read('~/Stage/Data/Quadrats/agriclimate_quadrats_Ireland.shp')
squares_modis = st_transform(squares, crs_modis)  # Make sure CRS is identical (required for gstat interpolation)
crs_squares = st_crs(squares)


# 1. Retrieve slope data

# Save modis square - not find better to use a raster function (not available in stars)
fname = paste0("all_squares_elevation.tif")

# Use the field function
crop_slope = raster(file.path(outputDir,fname))
crop_slope = terrain(crop_slope, opt='slope', unit='degrees')

# phenology_wide$x_MODIS = sapply(phenology_wide$geometry,"[[",1)
# phenology_wide$y_MODIS = sapply(phenology_wide$geometry,"[[",2)

phenology_wide$slope = as.numeric(extract(crop_slope,
                                          SpatialPoints(data.frame(phenology_wide$x_MODIS,
                                                                   phenology_wide$y_MODIS))))



# Create a geometry column for phenology wide
phenology_wide = st_as_sf(phenology_wide, coords = c("x_MODIS", "y_MODIS"),
                          crs = crs_modis)

d_pheno = st_as_sf(d_pheno, coords = c("x_MODIS", "y_MODIS"),
                   crs = crs_modis)
# Convert dates to a date object for temperature and precipitation
d_pheno$date = as.Date(as.character(d_pheno$date), format="%Y-%m-%d")


# 2. Retrieve elevation data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
elevation = read_stars(file.path(outputDir, 'all_squares_elevation.tif'))
st_crs(elevation) = crs_modis

phenology_wide$elevation = as.numeric(st_extract(elevation, st_sfc(phenology_wide$geometry, crs=crs_modis))[[1]])


# 3. Retrieve aspect slope data
aspect_slope = read_stars(file.path(outputDir, 'all_squares_aspect_slope.tif'))
st_crs(aspect_slope) = crs_modis

phenology_wide$aspect_slope = as.numeric(st_extract(aspect_slope, st_sfc(phenology_wide$geometry, crs=crs_modis))[[1]])


phenology_wide$class_aspect[phenology_wide$aspect_slope <= 20] = 'N'
phenology_wide$class_aspect[phenology_wide$aspect_slope >= 340] = 'N'
phenology_wide$class_aspect[phenology_wide$aspect_slope > 20 & phenology_wide$aspect_slope < 70] = 'NE'
phenology_wide$class_aspect[phenology_wide$aspect_slope >= 70 & phenology_wide$aspect_slope <= 110] = 'E'
phenology_wide$class_aspect[phenology_wide$aspect_slope > 110 & phenology_wide$aspect_slope < 160] = 'SE'
phenology_wide$class_aspect[phenology_wide$aspect_slope >= 160 & phenology_wide$aspect_slope <= 200] = 'S'
phenology_wide$class_aspect[phenology_wide$aspect_slope > 200 & phenology_wide$aspect_slope < 250] = 'SW'
phenology_wide$class_aspect[phenology_wide$aspect_slope >= 250 & phenology_wide$aspect_slope <= 290] = 'W'
phenology_wide$class_aspect[phenology_wide$aspect_slope > 290 & phenology_wide$aspect_slope < 340] = 'NW'



# 4. Retrieve soil data IFS
# ++++++++++++++++++++++++++++++++++++++++++++++
soil_data = st_read(file.path(soilPath,"Soils_IE_WetDry.shp"))
# Convert soil data in MODIS crs
soil_data_modis = st_transform(soil_data, crs = crs_modis)

# Make a join
phenology_wide = st_join(st_as_sf(phenology_wide, crs = crs_modis),
                         soil_data_modis["CATEGORY"])





# Retrieve MERA data
filename = paste0('soilmoisture_2012_to_2017.RData')
load(file.path(climatePath,filename))

filename2 = paste0('temperature_degrees_2012_to_2017.RData')
load(file.path(climatePath,filename2))

temperature_mera$date = as.Date(as.character(temperature_mera$validityDate),
                                format="%Y%m%d")

filename3 = paste0('precipitation_2012_to_2017.RData')
load(file.path(climatePath,filename3))

precipitation_mera$date = as.Date(as.character(precipitation_mera$validityDate),
                                  format="%Y%m%d")




# Function : Create variogram
variogramme <- function(mera_square) {

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



# For each square
for (s in squaresList) {
  print(s)

  pixel_list = unique(d_pheno$pixelID[d_pheno$square == s])
  nPixel = length(pixel_list)

  for (y in years) {
    print(y)


    # Subset mera data in the current year
    tmp_temperature_mera = subset(temperature_mera, temperature_mera$square == s
                                  & temperature_mera$date >= paste0(y, '-01-01')
                                  & temperature_mera$date < paste0(y, '-03-05'))

    tmp_precipitation_mera = subset(precipitation_mera, precipitation_mera$square == s
                                    & precipitation_mera$date >= paste0(y, '-01-01')
                                    & precipitation_mera$date < paste0(y, '-03-05'))


    days = unique(d_pheno$doy[d_pheno$year == y & d_pheno$square == s
                              & d_pheno$doy > 0
                              & d_pheno$doy < 60])
    days = sort(days)

    for (day in days) {

      # recuperer la date
      d = which(d_pheno$square == s & d_pheno$year == y
                & d_pheno$doy == day)[1]

      # 6. Retrieve temperature cumulated MERA data
      # ++++++++++++++++++++++++++++++++++++++++++++++

      mera_square_temp = subset(tmp_temperature_mera,
                                tmp_temperature_mera$date == d_pheno$date[d])

      g_mod2 = variogramme(mera_square_temp)

      # 7. Retrieve precipitation cumulated MERA data
      # ++++++++++++++++++++++++++++++++++++++++++++++

      mera_square_precip = subset(tmp_precipitation_mera,
                                  tmp_precipitation_mera$date == d_pheno$date[d])

      g_mod3 = variogramme(mera_square_precip)


      for (i in pixel_list) {

        # Retrieve the index of one of the pixels
        ind = which(d_pheno$square == s & d_pheno$year == y
                    & d_pheno$pixelID == i & d_pheno$doy == day)[1]

        if (length(ind) != 0 & !is.na(ind)) {

          if (length(g_mod2) != 0) {
            d_pheno$temperature[d_pheno$square == s & d_pheno$year == y & d_pheno$pixelID == i
                                & d_pheno$doy == day] = as.numeric(st_extract(g_mod2, d_pheno$geometry[ind])[[1]])
          }

          if (length(g_mod3) != 0) {
            d_pheno$precipitation[d_pheno$square == s & d_pheno$year == y & d_pheno$pixelID == i
                                  & d_pheno$doy == day] = as.numeric(st_extract(g_mod3, d_pheno$geometry[ind])[[1]])
          }

        }

      }
    }

    # 5. Retrieve soil moisture MERA data
    # ++++++++++++++++++++++++++++++++++++++++++++++

    mera_square = subset(soilmoisture_mera, soilmoisture_mera$square == s
                         & (soilmoisture_mera$validityDate > paste0(y, '02-21')
                            & soilmoisture_mera$validityDate < paste0(y, '03-01')))

    g_mod = variogramme(mera_square)


    for (i in pixel_list) {

      # Retrieve the index of one of the pixels
      ind = which(phenology_wide$square == s & phenology_wide$year == y & phenology_wide$pixelID == i)

      if (length(ind) != 0) {

        if (length(g_mod) != 0) {
          # Soil moisture
          phenology_wide$soilmoisture[phenology_wide$square == s
                                      & phenology_wide$year == y
                                      & phenology_wide$pixelID == i] = as.numeric(st_extract(g_mod,phenology_wide$geometry[ind[1]]))
        }

        # Temperature cumulated over 5.5
        phenology_wide$cumul_temp[phenology_wide$square == s
                                  & phenology_wide$year == y
                                  & phenology_wide$pixelID == i] = sum(d_pheno$temperature[d_pheno$square == s
                                                                                           & d_pheno$year == y
                                                                                           & d_pheno$pixel == i
                                                                                           & d_pheno$temperature >= 5.5], na.rm = T)

        # Precipitation cumulated over 0
        phenology_wide$cumul_precip0[phenology_wide$square == s
                                     & phenology_wide$year == y
                                     & phenology_wide$pixelID == i] = sum(d_pheno$precipitation[d_pheno$square == s
                                                                                                & d_pheno$year == y
                                                                                                & d_pheno$pixel == i], na.rm = T)

        # Precipitation cumulated over 1
        phenology_wide$cumul_precip1[phenology_wide$square == s
                                     & phenology_wide$year == y
                                     & phenology_wide$pixelID == i] = sum(d_pheno$precipitation[d_pheno$square == s
                                                                                                & d_pheno$year == y
                                                                                                & d_pheno$pixel == i
                                                                                                & d_pheno$precipitation > 1], na.rm = T)

        # Precipitation cumulated over 5
        phenology_wide$cumul_precip5[phenology_wide$square == s & phenology_wide$year == y
                                     & phenology_wide$pixelID == i] = sum(d_pheno$precipitation[d_pheno$square == s
                                                                                                & d_pheno$year == y
                                                                                                & d_pheno$pixel == i
                                                                                                & d_pheno$precipitation > 5], na.rm = T)
      }

    }

  }

}


# ++++++++++++++++++++++++++++++++++++++++++++++

#save(phenology_wide, file=paste0(dataDir, '/gls_improve_model_square1_to_21_in_2017.Rdata'))
save(phenology_wide, file=paste0(dataDir, '/gls_improve_model_square1_to_21_in_2013.Rdata'))
#save(phenology_wide, file=paste0(dataDir, '/gls_improve_model_square20_in_2013_to_2017.Rdata'))