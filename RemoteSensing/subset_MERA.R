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

library(sf)
#library(stars)
#library(gstat)

meraDir = "/Volumes/MODIS_data/MERA"
squareDir = "/Volumes/MODIS_data/Quadrats"
modisDir = "/Volumes/MODIS_data/MODIS"

meraFilePrefix = "TPrecip/TotalPrecip"
outputDir = "/Volumes/MODIS_data/MERA"
output_prefix = "TotalPrecip_subset"

# Define padding (in degrees) around square
pad = 0.05


months = c("01","02","03","04","05","06","07","08","09","10","11","12")
year =2012

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Import the square quadrats
squares = st_read(file.path(squareDir, 'agriclimate_quadrats_Ireland.shp'))
nSquares = nrow(squares)
crs_modis = st_crs(squares)

# Convert squares to wgs84
squares_wgs84 = st_transform(squares, crs=4326)

for (m in months) {
  # Import MERA data
  meraFilename = paste0(meraFilePrefix,"_",year,'_',m,'.txt')
  
  print(paste0("Importing file ",meraFilename))
  d = read.table(file.path(meraDir,meraFilename),
                 sep='',
                 header=TRUE)
  
  # Remove text data from the data frame
  d$Latitude. = as.numeric(d$Latitude.) # Text will be converted to NA
  d$Longitude. = as.numeric(d$Longitude.) # Text will be converted to NA
  d$Value. = as.numeric(d$Value.)
  d$validityDate. = as.numeric(d$validityDate.)
  d$validityTime = as.numeric(d$validityTime)
  d = na.omit(d)
  
  # Convert Longitude values
  d$Longitude.[d$Longitude.>180] =  d$Longitude.[d$Longitude.>180]-360
  
  
  # # Read in MODIS grid
  # modis = read_stars(file.path(modisDir,'modis_grid_ireland.tif'))
  # st_crs(modis) =  crs_squares  # Make sure modis and squares have same CRS
  
  
  for (s in 1:nSquares) {
    bbox = st_bbox(squares_wgs84[s,])
    d_sub = subset(d, Longitude.>bbox[1]-pad & Longitude.<bbox[3]+pad & 
                     Latitude.>bbox[2]-pad & Latitude.<bbox[4]+pad)
    
    d_sub_agg = aggregate(Value.~Longitude.+Latitude.+validityDate., 
                          data=d_sub, 
                          FUN=mean,
                          na.rm=TRUE)
    d_sub_agg$square = s
    
    if (!exists('mera_subset')) {
      mera_subset = d_sub_agg
    } else {
      mera_subset = rbind(mera_subset, d_sub_agg)
    }
    # # Convert this subset into a spatial object
    # d_wgs84 = st_as_sf(d_sub_agg, 
    #                coords=c("Longitude.", "Latitude."),
    #                crs=4326)
    # # Transform onto modis grid
    # d_modis = st_transform(d_wgs84, crs_squares)
    # 
    # v = variogram(Value. ~ 1, data = d_modis)
    # v.fit = fit.variogram(v, vgm(0, "Lin", NA))
    # g = gstat(formula = Value. ~ 1, data = d_modis)
    # 
    # # Now interpolate d_modis to put it on the modis grid
    # crop_modis = st_crop(modis, squares[s,])
    # d_crop = predict(crop_modis, model=g)
  }
  
}
names(mera_subset) = c("Longitude",'Latitude','ValidityDate','DailyMean','SquareID')


write.csv(mera_subset, 
          file=file.path(outputDir,paste0(output_prefix,'_',year,'.csv')),
          row.names=FALSE, quote=FALSE)
save(mera_subset, file = file.path(outputDir,paste0(output_prefix,'_',year,'.RData')))



# library(automap)
# v_mod_ok = autofitVariogram(Value. ~ 1, as(d_modis, "Spatial"))
# g = gstat(formula = Value. ~ 1, model = v_mod_ok, data = d_modis)
# z = predict(crop_modis, model=g)
# 
# ggplot() +  
#   geom_sf(data=squares[s,]) +
#   geom_sf(data=d_modis) 
# 
