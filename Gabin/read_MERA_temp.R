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


# Import first column of the file
# Using fread is faster
temp = fread(files[[1]], select=c(1))
names(temp) = c('Latitude')


# Find values of Latitude within the file
ind = which(grepl('Latitude', temp$Latitude))

# Import just one set of results
temp_sub = fread(files[[1]], 
                 skip=0, 
                 nrows=nrow_read,
                 select=c(1,2),
                 showProgress=FALSE)
names(temp_sub) = c('Latitude', 'Longitude')


# Convert variables and calculate subset overlapping quadrats
ind_long = temp_sub$Longitude>180
temp_sub$Longitude[ind_long] =  temp_sub$Longitude[ind_long] - 360

coordinates(temp_sub) = ~ Longitude + Latitude
crs(temp_sub) = CRS("+init=epsg:4326")

quadrats_wgs84$ID = c(1:nrow(quadrats_wgs84))
tmp = over(temp_sub, quadrats_wgs84)
ind_quadrats = which(!is.na(tmp[,1]))

sq_info = rbind(tmp[ind_quadrats,c(1,2,5)])


# Import data by reading sections of the files
for (f in 1:length(files)) {
  print(paste('Reading file ',files[[f]]))
  
  if (f>1) {
    temp = fread(files[[f]], select=1)
    names(temp) = c('Latitude')
    # Find values of Latitude within the file
    ind = which(grepl('Latitude', tmp$Latitude))
  }
  for (i in 0:length(ind)) {
    if (i==0) {
      offset = 0
    } else {
      offset = ind[i]
    }
    
    temp_sub = fread(files[[f]], 
                     skip=offset, 
                     nrows=nrow_read,
                     select=c(1,2,3,6,7),
                     colClasses=c(rep('numeric',3),rep('character',4)))
    names(temp_sub) = c('Latitude', 'Longitude', 'Value', 
                        'validityDate', 'validityTime')
    
    
    timeStr = temp_sub$validityTime[1]
    if (nchar(timeStr)<4) {
      timeStr = paste0('0',timeStr)
    }
    tmp_data = data.frame(temp_sub[ind_quadrats,],
                          date = as.Date(temp_sub$validityDate[1], 
                                         "%Y%m%d"),
                          sq_info,
                          time = strptime(paste0(temp_sub$validityDate[1],' ',
                                                 strtrim(timeStr,2),':00:00'),
                                          format="%Y%m%d %H:%M:%S"))
    
    if (i==0 & f==1) {
      temp_final = tmp_data
    } else {
      temp_final = rbind(temp_final, tmp_data)
    }
  }
}

temp_final$Temp = temp_final$Value - 273.15
temp_final$doy = as.numeric(format(temp_final$date,'%j'))



# Save data as text file
#write.csv(temp_final, file = 'TPrecip_squares.csv')
save(temp_final, file = output_file)

