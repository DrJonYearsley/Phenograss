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

# Import the quadrat and Ireland data
quadrats = readOGR(dsn='Quadrats/', layer='agriclimate_quadrats_Ireland')
ireland = readOGR(dsn='Quadrats/', layer='country')

# Calculate bounding box in degrees (WGS84)
ireland_wgs84 = spTransform(ireland, CRSobj=CRS("+init=epsg:4326"))
quadrats_wgs84 = spTransform(quadrats, CRSobj=CRS("+init=epsg:4326"))

# Using fread is faster
rain = fread(files[[1]], select=1, showProgress=FALSE)
names(rain) = c('Latitude')


# Find values of Latitude within the file
ind = which(grepl('Latitude', rain$Latitude))

# Calculate number of grid squares (rows) to read in
nrow_read = ind[1] - 1


# Clear rain data
rm(list='rain')

# Import just one set of results
rain_sub = fread(files[[1]], 
                 skip=0, 
                 nrows=nrow_read,
                 select=c(1,2),
                 showProgress=FALSE)
names(rain_sub) = c('Latitude', 'Longitude')



# Convert variables 
ind_long = rain_sub$Longitude>180
rain_sub$Longitude[ind_long] =  rain_sub$Longitude[ind_long] - 360

coordinates(rain_sub) = ~ Longitude + Latitude
crs(rain_sub) = CRS("+init=epsg:4326")

quadrats_wgs84$ID = c(1:nrow(quadrats_wgs84))
tmp = over(rain_sub, quadrats_wgs84)
ind_quadrats = which(!is.na(tmp[,1]))

sq_info = rbind(tmp[ind_quadrats,c(1,2,5)])

for (f in 1:length(files)) {
  print(paste('Reading file ',files[[f]]))
        
  if (f>1) {
    rain = fread(files[[f]], select=1)
    names(rain) = c('Latitude')
    # Find values of Latitude within the file
    ind = which(grepl('Latitude', rain$Latitude))
  }
        
  for (i in 0:length(ind)) {
  
    if (i==0) {
      offset = 0
    } else {
      offset = ind[i]
    }
    
    
    rain_sub = fread(files[[f]], 
                     skip=offset, 
                     nrows=nrow_read,
                     select=c(1,2,3,6,7),
                     colClasses=c(rep('numeric',3),rep('character',4)))
    names(rain_sub) = c('Latitude', 'Longitude','Value',
                        'validityDate', 'validityTime')
    
    
    timeStr = rain_sub$validityTime[1]
    if (nchar(timeStr)<4) {
      timeStr = paste0('0',timeStr)
    }
    
    tmp_data = data.frame(rain_sub[ind_quadrats,],
                          date = as.Date(rain_sub$validityDate[1], 
                                         "%Y%m%d"),
                          sq_info,
                          time = strptime(paste0(rain_sub$validityDate[1],' ',
                                   strtrim(timeStr,2),':00:00'),
                             format="%Y%m%d %H:%M:%S"))
    
    
    if (i==0 & f==1) {
      rain_final = tmp_data
    } else {
      rain_final = rbind(rain_final, tmp_data)
    }
  }
}


# Convert daily rainfall from kg/m-3 to mm
rain_final$Rain_mm = rain_final$Value * 1000 / 997

rain_final$doy = as.numeric(format(rain_final$date,'%j'))

# Save data as text file
#write.csv(rain_final, file = 'TPrecip_squares.csv')
save(rain_final, file = output_file)

