# Script to process MODIS data for Ireland suitable for empirical dynamical modelling
#
#
# This script performs these tasks
#  1. read in MODIS data (ndvi, evi, quality and day of acquisition)
#  2. crop to a predefined area, 
#  3. rescale NDVI and EVI 
#  4. remove poor quality pixels
#  5. read in CORINE (vector product)
#  6. reproject CORINE onto MODIS CRS
#  7. subset MODIS data to be only pixels for different land covers (e.g. pasture = CORINE==231)
#  8. saves the processed raster data for NDVI, EVI and the quality data
#  9. convert MODIS data into a data frame and save data frame
#  10. Make data frame into a spatial object
#  11. Reproject spatial data frame into ITM (espg:2157)

# The final data that is saved to file will have a spatial resolution of approximately 250m 
#        (i.e. the raw MODIS resolution, albeit transformed onto Irish Grid TM75)
#
# Jon Yearsley Nov 2019 
# ****************************************************


library(rgdal)
library(raster)
library(gdalUtils)

rm(list=ls())
#setwd("F:/")
setwd('~/WorkFiles/PeopleStuff/GrasslandPhenology/')

# Directories containing the input and output MODIS data 
inputDir = c('./Data/MODIS/MYD13Q1.006','./Data/MODIS/MOD13Q1.006')
outputDir = './EDM'
outputFileSuffix = 'modis_multi_landcover_data_'

save_rasters = FALSE  # If true save rasters for each MODIS file
yearStr = 'A201[67]' # Some text (or reg experession) that specifies the year of the data (e.g. 'A20[0-9][0-9]' specifies years 2000-2019) ###modified to 2018###
minQuality = 1 # Minimum quality to use: 0 = use only best quality pixels, 1=use reasonable pixels
scalingFactor = 0.0001 # Scale factor to apply to NDVI and EVI data from MODIS
allIreland=TRUE     # If TRUE then calculate for all of Ireland otherwise use 10km sq

# define a 10 km square within Ireland
sq_10km = extent(-504000,-493000, 5893000, 5904400)


# Land cover type to contrast
corineInclude = list(pasture=c(231), 
                     forest=c(311,312,313), 
                     urban=c(111, 112, 121,122,123,124))

corinePath = './Data/CORINE2018_Ireland/CLC18_IE_ITM'

# Read in CORINE data. make sure this is vector data (shapefile)
corine = readOGR(dsn=file.path(corinePath,'CLC18_IE_ITM.shp'), layer='CLC18_IE_ITM')

# Load MODIS data
regexp = paste(yearStr,'[[:graph:]]+.hdf$',sep='')
hdf.files = list.files(path=inputDir,pattern=regexp,recursive=T, full.names=TRUE)
nFiles = length(hdf.files) # Calculate number of files to import


# Extract date of the file
file.dates = sapply(hdf.files, 
               FUN=function(x){
                 tmp=regexpr(pattern = '/[0-9]{4}.[0-9]{2}.[0-9]{2}/',x)
                 dat=substring(x, first=tmp[1]+1, last=tmp[1]+attr(tmp,'match.length')-2 )
                 return(dat)}, 
               simplify=TRUE, USE.NAMES=FALSE) 
r.file.date = strptime(file.dates, format = "%Y.%m.%d", tz = "")

# Define extent of Irleand (roughly) in MODIS CRS 
ir = extent(-7.5E5, -3.3E5,5.7E6, 6.17E6)

 
# Find MODIS CRS and reproject CORINE onto MODIS
sds <- get_subdatasets(hdf.files[1])
evi = crop(raster(sds[grep("250m 16 days EVI", sds)], as.is=T), ir)*scalingFactor^2
modis_crs = crs(evi)


# Reproject CORINE
corine_modis = spTransform(corine, CRS=modis_crs)


# Extract land covers using CORINE
lc_modis = vector('list', length=length(corineInclude))
names(lc_modis) = names(corineInclude)
for (i in 1:length(corineInclude)) {
  lc_modis[[i]] = subset(corine_modis, CODE_18%in%corineInclude[[i]])
}

# Create list to hold data from different land covers
d = vector('list', length=length(corineInclude))
names(d) = names(corineInclude)
coord_info = vector('list', length=length(corineInclude))
names(coord_info) = names(corineInclude)

for (f in 1:length(hdf.files)) {
  print(paste('File',f,'from',length(hdf.files),':',hdf.files[f]))
  
  # Read in the MODIS data and crop to Ireland 
  sds <- get_subdatasets(hdf.files[f])
  
  # These lines read in the individual wavelength bands. Used to check the scaling factor for ndvi and evi
  if (allIreland) {
    evi = crop(raster(sds[grep("250m 16 days EVI", sds)], as.is=T),  ir)*scalingFactor^2
    QC = crop(raster(sds[grep("16 days pixel reliability",sds)]),  ir)
    doy = crop(raster(sds[grep("16 days composite day of the year",sds)]), ir)
    # Reliability, 0=good, 1=OK but use with care, 2=snow/icd, 3=cloudy,-1=no data 
  } else {
    evi = crop(raster(sds[grep("250m 16 days EVI", sds)], as.is=T),  sq_10km)*scalingFactor^2
    QC = crop(raster(sds[grep("16 days pixel reliability",sds)]),  sq_10km)
    doy = crop(raster(sds[grep("16 days composite day of the year",sds)]), sq_10km)
    # Reliability, 0=good, 1=OK but use with care, 2=snow/icd, 3=cloudy,-1=no data 
  }
  # Keep only good quality data (reliability=0 or 1) and reproject onto Irish grid
  evi[QC<0 | QC>minQuality] <- NA
  QC[QC<0 | QC>minQuality] <- NA
  
  # Set negative evi & ndvi to zero
  evi[evi<0] <- 0
  
  
  for (i in 1:length(corineInclude)) {
    # Extract cells corresponding to  CORINE land cover
    evi_lc = mask(evi, lc_modis[[i]])
    
    # Create a data frame
    if (f==1) {
      # Create spatial data to convert coords to ITM
      xy_modis <- SpatialPoints(coords = coordinates(evi_lc), 
                                proj4string = modis_crs)
      xy_ITM = spTransform(xy_modis, CRSobj = CRS("+init=epsg:2157"))
      
      coord_info[[i]] = as.data.frame(cbind(coordinates(xy_ITM), coordinates(xy_modis)))
      names(coord_info[[i]]) = c('x_ITM','y_ITM','x_MODIS','y_MODIS')
      
      
      d[[i]] = cbind(coord_info[[i]],
                data.frame(year=format(r.file.date[f],"%Y"),
                           doy= getValues(doy), 
                           evi=getValues(evi_lc), 
                           QC=getValues(QC)))
    } else {
      d[[i]] = rbind(d[[i]], 
                cbind(coord_info[[i]],
                      data.frame(year=format(r.file.date[f],"%Y"),
                                 doy= getValues(doy), 
                                 evi=getValues(evi_lc), 
                                 QC=getValues(QC))  ))
    }
  }
}  
  
# Remove missing data
for (i in 1:length(corineInclude)) {
  d[[i]] = subset(d[[i]], is.finite(evi))
}

# Save data to a file
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

if (allIreland) {
  regionStr='allIreland'
} else {
  regionStr = '10kmsq'
}


fname.df = file.path(outputDir,paste0(outputFileSuffix,'_',regionStr,'_',yearStr,'.RData') )
save(d, r.file.date, hdf.files,corineInclude, file = fname.df)

