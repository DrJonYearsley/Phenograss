# 
# Make the accumulation of temperature of 2 months (Jan, Feb.)
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

modisPath = '~/Stage/Data/MODIS'
climatePath = '~/Stage/Data/Climate'

dataDir = '~/Stage/Data/MODIS/Phenophase_estimates'
outputDir = '~/Stage/Data_created'

input_file_preffix = 'phenology'

# Import data --------
squaresList = c(1:9, 13:21)
years = c(2015)
# Test sur 2 jours
#days = c(0:90)

# Filename segmented data
for (i in 1:length(squaresList)) {
  
    filename = paste0(input_file_preffix,'_square_',squaresList[i],'_',years[1],'.RData')
    load(file.path(dataDir,filename))
    
    if (i==1) {
      phenology = d_final
    } else {
      phenology = rbind(phenology,
                        d_final)
    }
}



# Read in MODIS grid
modis = read_stars(file.path(modisPath, 'modis_grid_ireland.tif'))
crs_modis = st_crs(modis)

squares = st_read('~/Stage/Data/Quadrats/agriclimate_quadrats_Ireland.shp')
squares_modis = st_transform(squares, crs_modis)  # Make sure CRS is identical (required for gstat interpolation)
crs_squares = st_crs(squares)


# Select the 2 months January and February (0 to 90 days)
tmp = subset(phenology, phenology$doy >= 0 & phenology$doy <= 59) # 2 5 7 8 10...

# Create a geometry column
tmp = st_as_sf(tmp, coords = c("x_MODIS", "y_MODIS"), crs = crs_modis)


filename2 = paste0('temperature_degrees_2012_2013_2015.RData')
load(file.path(climatePath,filename2))

temperature_mera = subset(temperature_mera, temperature_mera$validityDate >= '20150101'
                          & temperature_mera$validityDate < '20160101')


# Convert dates to a date object
temperature_mera$date = as.Date(as.character(temperature_mera$validityDate),
                                format="%Y%m%d")

tmp$date = as.Date(as.character(tmp$date), format="%Y-%m-%d")


# Function : Create variogram
variogramme <- function(covariate_mera, d) {

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




for (s in squaresList) {
  print(s)
  
  days = unique(tmp$doy[tmp$square == s])
  days = sort(days)
  
  for (day in days) {
    
    # recuperer la date
    d = which(tmp$square == s & tmp$doy == day)[1]
    
    # 4. Retrieve temperature MERA data
    # ++++++++++++++++++++++++++++++++++++++++++++++
    
    g_mod2 = variogramme(temperature_mera, tmp$date[d])
    
    if (length(g_mod2) == 0) {
      tmp$temperature[tmp$square == s & tmp$doy == day] = NA
      
    } else {
      
      tmp2 = subset(tmp, tmp$square == s)
      pixel_list = unique(tmp2$pixelID)
      nPixel = length(pixel_list)
      
      for (i in pixel_list) {
        
        # Retrieve the index of one of the pixels
        ind = which(tmp$square == s & tmp$pixelID == i & tmp$doy == day)[1]
        if (is.na(ind)) {
          tmp$temperature[tmp$square == s & tmp$pixelID == i &  tmp$doy == day] = NA
          #tmp$precipitation[tmp$square == s & tmp$pixelID == i &  tmp$doy == day] = NA
        } else {
          tmp$temperature[tmp$square == s & tmp$pixelID == i &  tmp$doy == day] = as.numeric(st_extract(g_mod2, tmp$geometry[ind]))
          #tmp$precipitation[tmp$square == s & tmp$pixelID == i &  tmp$doy == day] = as.numeric(st_extract(g_mod3, tmp$geometry[ind]))
        }
        
      } 
      
    }
  }
  
  # Average over the days
  tmp2 = subset(tmp, tmp$square == s & !is.na(tmp$temperature))
  tmp2 = tmp2 %>% group_by(doy) %>% summarize(temperature = mean(temperature,na.rm=T))
  
  
  # Optimized in improvement_model.R
  c = 0 # counter
  # Accumulate ascending when the temperature is above 5.5 per days
  for (day in days) {
    
    if (tmp2$temperature[tmp2$doy == day] > 5.5) {
      
      value = tmp2$temperature[tmp2$doy == day] - 5.5
      
      tmp2$temp_cumulate[tmp2$doy == day] = value + c
      c = tmp2$temp_cumulate[tmp2$doy == day]
      
    } else {
      tmp2$temp_cumulate[tmp2$doy == day] = c
    }
  }
  
  
  ggplot(data=tmp2,
         aes(x=doy,
             y=temp_cumulate)) +
    geom_line() +
    #geom_boxplot() +
    #geom_density() +
    labs(x='Day of the year',
         y='Cumulate temperature (°C)') +
    theme_bw() +
    theme(axis.title = element_text(size=12),
          axis.text = element_text(size=12))
  
  
  #ggsave(filename = paste0('temp_cumulate_square', s, '_year2015.png'))
  
  # print(mean(tmp$temperature,na.rm=T))
  # print(min(tmp$temperature,na.rm=T))
  # print(max(tmp$temperature,na.rm=T))
  
  n = nrow(tmp2)
  tmp$cum_temp[tmp$square == s] = tmp2$temp_cumulate[n]
}


# Display distribution of temperature cumulated
ggplot(data=tmp,
       aes(x=square,
           y=cum_temp)) +
  geom_point() +
  #geom_boxplot() +
  labs(x='squares',
       y='Cumulate temperature in 2015 (°C)') +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12))

ggsave(filename = paste0('temp_cumulate_with_square_year2015.png'))

