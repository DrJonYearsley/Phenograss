# Read MERA temperature data
#
# Convert Kelvin to degrees C
# Crop to quadrats for phenology
#
# Jon Yearsley
# July 2020
#
# ***************************************

setwd("~/WorkFiles/MEGA/Projects/GrasslandPhenology/Data/")

rm(list=ls())
library(rgdal)
library(raster)
library(ggplot2)
library(data.table)

output_file = 'Temp2m_2017_squares.RData'

files = list.files(path = '~/WorkFiles/Data/Climate/MERA/T2m',
                   pattern='Temp2m_2017',
                   full.names = TRUE)



# Import data

# Import the quadrat and Ireland data
quadrats = readOGR(dsn='Quadrats/', layer='agriclimate_quadrats_Ireland')
ireland = readOGR(dsn='Quadrats/', layer='country')

# Calculate bounding box in degrees (WGS84)
ireland_wgs84 = spTransform(ireland, CRSobj=CRS("+init=epsg:4326"))
quadrats_wgs84 = spTransform(quadrats, CRSobj=CRS("+init=epsg:4326"))

box = bbox(ireland_wgs84)

# Using fread is faster
temp = fread(files[[1]], select=c(1,2,3,6,7))
names(temp) = c('Latitude', 'Longitude', 'Value', 
                'validityDate', 'validityTime')


# Find values of Latitude within the file
ind = which(grepl('Latitude', temp$Latitude))

# Create a subset of temp with first day's results
temp_sub = temp[1:(ind[1]-1),]


# Convert variables and remove the headers for each day
temp_sub$Latitude = as.numeric(temp_sub$Latitude)
temp_sub$Longitude = as.numeric(temp_sub$Longitude)
ind_long = temp_sub$Longitude>180
temp_sub$Longitude[ind_long] =  temp_sub$Longitude[ind_long] - 360

coordinates(temp_sub) = ~ Longitude + Latitude
crs(temp_sub) = CRS("+init=epsg:4326")

quadrats_wgs84$ID = c(1:nrow(quadrats_wgs84))
tmp = over(temp_sub, quadrats_wgs84)
ind_quadrats = which(!is.na(tmp[,1]))

sq_info = rbind(tmp[ind_quadrats,c(1,2,5)])

for (f in 1:length(files)) {
  if (f>1) {
    temp = fread(files[[f]], select=c(1,2,3,6,7))
    names(temp) = c('Latitude', 'Longitude', 'Value', 
                    'validityDate', 'validityTime')
    # Find values of Latitude within the file
    ind = which(grepl('Latitude', temp$Latitude))
  }
  
  for (i in 0:length(ind)) {
    if (i==0) {
      offset = 0
    } else {
      offset = ind[i]
    }
    
    timeStr = temp$validityTime[1 + offset]
    if (nchar(timeStr)<4) {
      timeStr = paste0('0',timeStr)
    }
    tmp_data = data.frame(temp[ind_quadrats + offset,],
                          date = as.Date(temp$validityDate[1 + offset], "%Y%m%d"),
                          sq_info,
                          time = strptime(paste0(temp$validityDate[1 + offset],' ',
                                                 strtrim(timeStr,2),':00:00'),
                                          format="%Y%m%d %H:%M:%S"))
    
    if (i==0 & f==1) {
      temp_final = tmp_data
    } else {
      temp_final = rbind(temp_final, tmp_data)
    }
  }
}

temp_final$Latitude = as.numeric(temp_final$Latitude)
temp_final$Longitude = as.numeric(temp_final$Longitude)
temp_final$Value = as.numeric(temp_final$Value)
temp_final$Temp = temp_final$Value - 273.15
temp_final$doy = as.numeric(format(temp_final$date,'%j'))



# Save data as text file
#write.csv(temp_final, file = 'TPrecip_squares.csv')
save(temp_final, file = output_file)

