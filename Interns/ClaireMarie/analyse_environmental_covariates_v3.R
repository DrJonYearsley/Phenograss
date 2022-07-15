# Analyse environmental covariates
#
# Claire-Marie Alla
# 21/06/2021
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


elevation_dir = '~/Stage/Data_created/elevation_SRTM3_square'
modisPath = '~/Stage/Data/MODIS'
soilPath = '~/Stage/Data/soil_data/Soils_IE_WetDry'
climatePath = '~/Stage/Data/Climate'
tempPath = '~/Stage/Data/Climate/Temp_2m_quadrats'

dataDir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created'

input_file_preffix = 'phenology'

# Import data --------
squaresList = c(1:9, 13:21)
years = c(2012,2013)


# if you have gls_model_2002_2019.Rdata
# this part can be removed 
# ++++++++++++++++++++++++++++++++++++++++++++++

# Filename segmented data
for (i in 1:length(squaresList)) {
  for (y in 1:length(years)) {
    filename = paste0(input_file_preffix,'_square_',squaresList[i],'_',years[y],'.RData')
    load(file.path(dataDir,filename))
    
    output_smoothed$square = squaresList[i]
    
    if (y==1 & i==1) {
      phenology = output_smoothed
    } else {
      phenology = rbind(phenology,
                        output_smoothed)
    }
  }
}

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

# Coordinates
locations = subset(phenology_wide, 
                   subset=year==years[1],
                   select=c('pixelID','x_ITM','y_ITM'))

# Calculate the centre coordinates for each square
tmp = aggregate(cbind(x_ITM_centre, y_ITM_centre)~square+year, 
                data=phenology_wide, 
                FUN=mean, 
                na.rm=TRUE)

rm(list='phenology')

# Fit model for phenology dates of phase 1 -----


# Mixed model for phase 1
m_phase1 = lme(t_phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)

summary(m_phase1)


# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time

phenology_wide$dummy = interaction(phenology_wide$year, phenology_wide$square, sep='_')
gls_phase1 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)

summary(gls_phase1)


anova(gls_phase1, type='marginal')

acf(gls_phase1$residuals, lag.max=20, type='partial')  # Makes sense... spatial correlation roughly 1-2 km



library(emmeans)

m_eff = emmeans(gls_phase1, spec='year')
sos = as.data.frame(summary(m_eff))
contrast(m_eff, infer=c(T,T))

as.Date(paste(round(sos$emmean),sos$year),format="%j %Y")

# WOrk out date from day of year
as.Date(paste('52',sos$year),format="%j %Y")

gls_phase1_v2 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                    data=phenology_wide,
                    na.action=na.exclude)

ggplot(data=as.data.frame(summary(m_eff)),
       aes(x=year,
           y=emmean,
           ymin=lower.CL,
           ymax=upper.CL)) +
  geom_pointrange() + 
  labs(x='Year',
       y='Day of Year') +
  theme_bw() + 
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18))


#save(gls_phase1, phenology_wide, file=paste0(dataDir, '/gls_model_2002_2019.Rdata'))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++

#load(file.path(dataDir,'gls_model_2002_2019.Rdata'))
# to big for load in my computer

phenology_wide$phase1_residual = residuals(gls_phase1)
phenology_wide$phase1_fit = predict(gls_phase1)

#square = 20

#tmp = subset(phenology_wide, year==2012 & !is.na(phase1_fit))
tmp = subset(phenology_wide, year==2013 & !is.na(phase1_fit) & square%in%c(1:21))


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
elevation = st_mosaic(read_stars(file.path(elevation_dir, 'square_1.tif')),
                      read_stars(file.path(elevation_dir, 'square_2.tif')),
                      read_stars(file.path(elevation_dir, 'square_3.tif')),
                      read_stars(file.path(elevation_dir, 'square_4.tif')),
                      read_stars(file.path(elevation_dir, 'square_5.tif')),
                      read_stars(file.path(elevation_dir, 'square_6.tif')),
                      read_stars(file.path(elevation_dir, 'square_7.tif')),
                      read_stars(file.path(elevation_dir, 'square_8.tif')),
                      read_stars(file.path(elevation_dir, 'square_9.tif')),
                      read_stars(file.path(elevation_dir, 'square_13.tif')),
                      read_stars(file.path(elevation_dir, 'square_14.tif')),
                      read_stars(file.path(elevation_dir, 'square_15.tif')),
                      read_stars(file.path(elevation_dir, 'square_16.tif')),
                      read_stars(file.path(elevation_dir, 'square_17.tif')),
                      read_stars(file.path(elevation_dir, 'square_18.tif')),
                      read_stars(file.path(elevation_dir, 'square_19.tif')),
                      read_stars(file.path(elevation_dir, 'square_20.tif')),
                      read_stars(file.path(elevation_dir, 'square_21.tif')))
#elevation = read_stars(file.path(elevation_dir, paste0('square_', square, '.tif')))
st_crs(elevation) = crs_modis

filename2 = paste0('temperature_degrees_2012_2013_2015.RData')
load(file.path(climatePath,filename2))

for (s in squaresList) {
  
  # 4. Retrieve temperature MERA data
  # ++++++++++++++++++++++++++++++++++++++++++++++

  # Select the square s in January / February / March
  mera_square = subset(temperature_mera, temperature_mera$square == s
                       & temperature_mera$validityDate > '20130101'
                       & temperature_mera$validityDate < '20140101')

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
  g_mod2 = gstat::krige(formula = Value ~ 1,
                      locations=test,
                      model=v_mod,
                      newdata=crop_square)

  # +++++++++++++++++++++++++++++++++++++++++++++
  
  
  tmp3 = subset(tmp, tmp$square == s)
  pixel_list = unique(tmp3$pixelID)
  nPixel = length(pixel_list)

  for (i in pixel_list) {
    print(i)
  
    # Retrieve the index of one of the pixels
    ind = which(tmp$square == s & tmp$pixelID == i)[1]
  
    # Create a elevation colomn thanks to geometry
    tmp$elevation[tmp$square == s & tmp$pixelID == i] = as.numeric(st_extract(elevation,
                                                                                 st_sfc(tmp$geometry[ind],
                                                                                        crs = crs_modis)))
    
    tmp$temperature[ind] = as.numeric(st_extract(g_mod2, tmp$geometry[ind]))  
  }
}


# Save tmp
save(tmp, file=paste0(outputDir, '/gls_model_2013_square1_20_elevation_temp.Rdata'))
#load(file.path(outputDir,'gls_model_2013_square1_20_elevation_temp.Rdata'))


# Calculate some R-squared between residus and environmentale covariates
#m = lm(phase1_residual~1 + elevation, data=tmp) # for residuals of 21 squares (2017)
#m = lm(phase1_residual~1 + CATEGORY, data=tmp) # for residuals of 21 squares (2017)
#m = lm(soilmoisture~1 + CATEGORY, data=tmp) # 2017

# m = lm(temperature~1 + elevation, data=tmp) # predire l'elevation en fonction de la temperature 2012 ?data_test
# 
# # library(rpart);
# # m = rpart(temperature~1 + elevation, data=tmp)
# 
# predict_lm <- predict(m, tmp2)
# #tmp$temp_residual = residuals(m)
# #tmp$temp_fit = predict(m)
# head(predict_lm)
# 
# 
# # total precipitation with total soil moisture
# #m = lm(soilmoisture~1 + precipitation, data=tmp)
# 
# predict_lm <- data.frame(temperature_Pred = predict_lm, temperature = tmp2$temperature)
# 
# summary(m)
# 
# ggplot(data=tmp2,
#        aes(x=elevation,
#            y=temperature,
#            colour=factor(square))) +
#   geom_point() +
#   geom_abline(slope=1, 
#               intercept=0, 
#               colour='darkred', 
#               linetype='dashed',
#               size=1.5) +
#   scale_colour_viridis_d('square') +
#   labs(x='Elevation (m)',
#        y='Temperature (°C)') +
#   theme_bw() + 
#   theme(axis.title = element_text(size=18),
#         axis.text = element_text(size=18))
# 
# 
# ggsave(filename = 'gls_fitted_vs_elevation_temp_year2012.png', width=30 , height=20, units='cm')
# 
# # Display difference between actual and predicted
# plot(predict_lm$temperature_Pred - predict_lm$temperature, main = "Difference between actual and predicted. lm")
# 
# # ggplot(data = predict_lm, ) + geom_point() + facet_grid(~)
# # ou geom_boxplot() + 
# # labs(y = "Taux de croissance du produit intérieur brut\n",
# # x = "\nRatio dette publique / produit intérieur brut (%)",
# # title = "Données Reinhart et Rogoff corrigées, 1946-2009\n") +
