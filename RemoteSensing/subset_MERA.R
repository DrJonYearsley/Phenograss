# Read in MERA data and subset for the 20 10Km squares
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# June 2021
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd('~/git_repos/Phenograss/RemoteSensing/')

# library(rgdal)
# library(sp)

library(data.table)
library(sf)
#library(stars)
#library(gstat)

base = "/Volumes"
meraDir = file.path(base,"MERA_Data")
squareDir = file.path(base,"MODIS_data/Quadrats")
modisDir = file.path(base, "MODIS_data/MODIS")

meraFilePrefix = "T2m/Temp2m"
outputDir = file.path(base,"MODIS_data/MERA")
output_prefix = "Temp_subset"

# Define padding (in degrees) around square
pad = 0.05


months = c("01","02","03","04","05","06","07","08","09","10","11","12")
#months = c("01","02","03","04","05","06","07","08","09","10")
year=2015

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Import the square quadrats
squares = st_read(file.path(squareDir, 'agriclimate_quadrats_Ireland.shp'))
nSquares = nrow(squares)
crs_modis = st_crs(squares)

# Convert squares to wgs84
squares_wgs84 = st_transform(squares, crs=4326)


# Import first 30k lines
tmp = fread(file.path(meraDir, paste0(meraFilePrefix,"_",year,'_01.txt')),
            select=1,
            nrows=300000,
            header=TRUE)
names(tmp) = "x1"
tmp$x1 = as.numeric(tmp$x1)
ind = which(is.na(tmp$x1))

chunkSize = ind[1]-1

# Import the data for each month
for (m in months) {
  # Import MERA data 
  meraFilename = paste0(meraFilePrefix,"_",year,'_',m,'.txt')
  
  print(paste0("Importing file ",meraFilename))
  
  fileIO = file(file.path(meraDir,meraFilename), 'rt')

  # Read header
  file_head = read.table(fileIO, nrows=1, header=FALSE)
  colNames = file_head[1,]
  
  # Import data in chunks
  import=TRUE
  counter = 1
  while (import) {
    now = Sys.time()
    file_chunk = read.table(fileIO, 
                            nrows = chunkSize, 
                            header=FALSE,
                            col.names = colNames,
                            colClasses = "numeric")
    
    counter = counter+1
    
    # Throw away the line with text. 
    # If end of file there will be an error and discard will be NULL
    discard = tryCatch(read.table(fileIO, nrows=1, header=FALSE),
                       error = function(e) {return(NULL)})
    
    # Stop import if the "discard" line doesn't exist
    if(length(discard)==0) {
      import=FALSE
    } else if (any(file_head!=discard)) {
      warning(paste0("Chunk ",counter,"didn''t end as expected"))
    }
    
    
    # Convert Longitude values
    file_chunk$Longitude[file_chunk$Longitude>180] =  file_chunk$Longitude[file_chunk$Longitude>180]-360
    
    for (s in 1:nSquares) {
      bbox = st_bbox(squares_wgs84[s,])
      d_sub = subset(file_chunk, Longitude>bbox[1]-pad & Longitude<bbox[3]+pad & 
                       Latitude>bbox[2]-pad & Latitude<bbox[4]+pad)
      
      # Record square ID
      d_sub$square = s
      
      
      if (!exists('mera_hourly')) {
        mera_hourly = d_sub
      } else {
        mera_hourly = data.table::rbindlist(list(mera_hourly, d_sub))
      }
      rm(list=c('d_sub'))
    }
    rm(list=c('file_chunk'))

    print(paste0('Read chunk ',counter,' in ',signif(difftime(Sys.time(), now),3),' secs.'))
  }
  # Close the connection to the file
  close(fileIO)
  
  # Save the spatially sampled data set
  save(mera_hourly, file = file.path(outputDir,paste0(output_prefix,'_hourly_',year,'_',m,'.RData')))
  
  
  rm(list=c('mera_hourly'))
}


# Create daily average data -----
d = vector('list',length=length(months))
f = paste0(names(mera_hourly)[3],'~',names(mera_hourly)[1],'+',names(mera_hourly)[2],'+',names(mera_hourly)[6],'+ square')
for (m in 1:length(months)) {
  load(file.path(outputDir,paste0(output_prefix,'_hourly_',year,'_',months[m],'.RData')))

  # Average over times in each day
  d[[m]] = aggregate(formula(f), 
                  data=mera_hourly, 
                  FUN=mean,
                  na.rm=TRUE)  
  
}

mera_daily = data.table::rbindlist(d)

write.csv(mera_daily, 
          file=file.path(outputDir,paste0(output_prefix,'_',year,'.csv')),
          row.names=FALSE, quote=FALSE)
save(mera_daily, file = file.path(outputDir,paste0(output_prefix,'_',year,'.RData')))



# library(automap)
# v_mod_ok = autofitVariogram(Value. ~ 1, as(d_modis, "Spatial"))
# g = gstat(formula = Value. ~ 1, model = v_mod_ok, data = d_modis)
# z = predict(crop_modis, model=g)
# 
# ggplot() +  
#   geom_sf(data=squares[s,]) +
#   geom_sf(data=d_modis) 
# 
