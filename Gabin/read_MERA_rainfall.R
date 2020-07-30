# Read precipitation data
#
# Convert preciptation to mm
# Calculate daily precipitation in each 10 km quadrate
#
# Jon Yearsley
# July 2020
#
# ***************************************


setwd("~/WorkFiles/MEGA/Projects/GrasslandPhenology/Data/")
#setwd("~/MEGAsync/Projects/GrasslandPhenology/Data/")

rm(list=ls())
library(rgdal)
library(raster)
library(ggplot2)
library(data.table)

output_file = 'TPrecip_2017_squares.RData'

files = list.files(path = '~/WorkFiles/Data/Climate/MERA/TPrecip',
                   pattern='TotalPrecip_2017',
                   full.names = TRUE)

# files = list.files(path = '~/Data/MERA',
#                    pattern='TotalPrecip_2013',
#                    full.names = TRUE)

# Import data
# rain = read.table('TotalPrecip_2017_01.txt', header=TRUE, sep='', 
#                   col.names=c('Latitude', 'Longitude', 'Rain',
#                 'dataDate', 'dataTime', 'validityDate','ValidityTime'))

# Import the quadrat and Ireland data
quadrats = readOGR(dsn='Quadrats/', layer='agriclimate_quadrats_Ireland')
ireland = readOGR(dsn='Quadrats/', layer='country')

# Calculate bounding box in degrees (WGS84)
ireland_wgs84 = spTransform(ireland, CRSobj=CRS("+init=epsg:4326"))
quadrats_wgs84 = spTransform(quadrats, CRSobj=CRS("+init=epsg:4326"))

box = bbox(ireland_wgs84)

# Using fread is faster
rain = fread(files[[1]])
names(rain) = c('Latitude', 'Longitude', 'Value', 'dataDate',
                'dataTime', 'validityDate', 'validityTime')


# Find values of Latitude within the file
ind = which(grepl('Latitude', rain$Latitude))

# Create a subset of rain with first day's results
rain_sub = rain[1:(ind[1]-1),]


# Convert variables and remove the headers for each day
rain_sub$Latitude = as.numeric(rain_sub$Latitude)
rain_sub$Longitude = as.numeric(rain_sub$Longitude)
ind_long = rain_sub$Longitude>180
rain_sub$Longitude[ind_long] =  rain_sub$Longitude[ind_long] - 360

coordinates(rain_sub) = ~ Longitude + Latitude
crs(rain_sub) = CRS("+init=epsg:4326")

quadrats_wgs84$ID = c(1:nrow(quadrats_wgs84))
tmp = over(rain_sub, quadrats_wgs84)
ind_quadrats = which(!is.na(tmp[,1]))

sq_info = rbind(tmp[ind_quadrats,c(1,2,5)])

for (f in 1:length(files)) {
  if (f>1) {
    rain = fread(files[[f]])
    names(rain) = c('Latitude', 'Longitude', 'Value', 'dataDate',
                    'dataTime', 'validityDate', 'validityTime')
    # Find values of Latitude within the file
    ind = which(grepl('Latitude', rain$Latitude))
  }
  
  for (i in 0:length(ind)) {
    if (i==0) {
      offset = 0
    } else {
      offset = ind[i]
    }
    
    timeStr = rain$validityTime[1 + offset]
    if (nchar(timeStr)<4) {
      timeStr = paste0('0',timeStr)
    }
    
    tmp_data = data.frame(rain[ind_quadrats + offset,],
                          date = as.Date(rain$validityDate[1 + offset], "%Y%m%d"),
                          sq_info,
                          time = strptime(paste0(rain$validityDate[1 + offset],' ',
                                   strtrim(timeStr,2),':00:00'),
                             format="%Y%m%d %H:%M:%S"))
    
    
    if (i==0 & f==1) {
      rain_final = tmp_data
    } else {
      rain_final = rbind(rain_final, tmp_data)
    }
  }
}

rain_final$Latitude = as.numeric(rain_final$Latitude)
rain_final$Longitude = as.numeric(rain_final$Longitude)
rain_final$Value = as.numeric(rain_final$Value)



# Convert daily rainfall from kg/m-3 to mm
rain_final$Rain_mm = rain_final$Value * 1000 / 997

rain_final$doy = as.numeric(format(rain_final$date,'%j'))

# Save data as text file
#write.csv(rain_final, file = 'TPrecip_squares.csv')
save(rain_final, file = output_file)

