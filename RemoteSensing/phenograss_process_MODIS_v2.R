# Script to process MODIS data for Ireland. 
# This version imports preprocessed CORINE data that 
# gives the fraction of a MODIS pixels that is pasture
#
# This script performs these tasks
#  1. read in MODIS data (ndvi, evi, quality and day of acquisition)
#  2. read in squares (roughty 10 km square) for subsets within Ireland
#  3. read in preprocessed CORINE data
#  For each square: 
#  4. Crop MODIS data to the extent of the square
#  5. rescale NDVI and EVI 
#  6. remove poor quality pixels
#  7. remove pixels with less than pastureThreshold fraction of pasture
#  8. Set NDVI and EVI of less than 0 to 0
#  9. Remove pixels with missing EVI data
#  10. convert MODIS data into a data frame and save data frame
#  11. Make data frame into a spatial object
#  12. Include the ITM (espg:2157) coords of each pixel in the data frame
#  13. Save the data frame for each square in a separate file

# The final data that is saved to file will have a spatial resolution of approximately 250m 
#        (i.e. the raw MODIS resolution, albeit transformed onto Irish Grid TM75)
#
# Jon Yearsley Jon.Yearsey@ucd.ie
# Jan 2020 
# ****************************************************

library(sp)
library(rgdal)
library(raster)
library(gdalUtils)

rm(list=ls())

# setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology')
setwd('~/Research/Phenograss')

# Directories containing the input and output MODIS data 
# inputDir = c('./Data/MODIS/MYD13Q1.006','./Data/MODIS/MOD13Q1.006')
inputDir = c('/Volumes/Untitled/MODIS/MYD13Q1.006','/Volumes/Untitled/Data/MODIS/MOD13Q1.006')
quadratPath = './Data/Quadrats'
outputDir = './Data/MODIS_squares'
outputSuffix = 'pasture'
save_rasters = FALSE  # If true save rasters for each MODIS file
yearStr = 'A2005' # Some text (or reg experession) that specifies the year of the data (e.g. 'A20[0-9][0-9]' specifies years 2000-2019) 
minQuality = 1 # Minimum quality to use: 0 = use only best quality pixels, 1=use reasonable pixels
scalingFactor = 0.0001 # Scale factor to apply to NDVI and EVI data from MODIS
corinePath = './Data/CORINE_Ireland'
corineFilename = 'corine2018_pasturecover.gri'
pastureThreshold = 0.7 # The minimum fraction of a pixel that is pasture
squareList = c(1:21)  # List of squares to analyse


# Import squares from shapefile
squares = readOGR(dsn=quadratPath,
                  layer='agriclimate_quadrats_Ireland')
nSquares = length(squares)

# Load MODIS grid data
modis = raster('./Data/MODIS/modis_grid_ireland.grd')

# Save MODIS CRS
modis_crs = crs(modis)

# Read in CORINE preprocessed data.
corine_modis = raster(file.path(corinePath,corineFilename))


# Define extent of Irleand (roughly) in MODIS CRS 
#ir = extent(-7.5E5, -3.3E5,5.7E6, 6.17E6)



extractbit = function(x,n1,n2) {
  # Extract the nth bit from the number x
  return(bitwAnd(bitwShiftR(x,n1),1))
}


# Search for relevant MODIS data files
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

# Create a list to contain data from each square
d = vector('list', length=length(squareList))

for (f in 1:nFiles) {
  print(hdf.files[f])
  
  # Read in the MODIS data and crop to Ireland 
  sds <- get_subdatasets(hdf.files[f])
  
  if (grepl('MOD13',hdf.files[f])) {
    satellite = 'Terra'
  } else if (grepl('MYD13',hdf.files[f])){
    satellite = 'Aqua'
  } else {
    satellite = NA
  }
  
  # These lines read in the individual wavelength bands. Used to check the scaling factor for ndvi and evi
  # ndvi is (nir-red)/(nir+red)
  #  red = crop(raster(sds[grep("250m 16 days red reflectance", sds)], as.is=T), ir)
  #  nir = crop(raster(sds[grep("250m 16 days NIR reflectance", sds)], as.is=T), ir)
  
  
  ndvi = raster(sds[grep("250m 16 days NDVI", sds)], as.is=T)
  evi = raster(sds[grep("250m 16 days EVI", sds)], as.is=T)
  QC = raster(sds[grep("16 days pixel reliability",sds)])
  qualbit = raster(sds[grep("16 days VI Quality",sds)])
  doy = raster(sds[grep("16 days composite day of the year",sds)])
  # Reliability, 0=good, 1=OK but use with care, 2=snow/icd, 3=cloudy,-1=no data 
  
  # # Extract a few quality controls 
  # shadow = calc(qualbit, fun=function(x) {floor(x/2^15) %% 2})
  
  for (s in squareList) {
    # Crop to one of the squares
    evi_crop = crop(evi, squares[s,]) * scalingFactor^2
    ndvi_crop = crop(ndvi, squares[s,]) * scalingFactor^2
    QC_crop = crop(QC, squares[s,])
    doy_crop = crop(doy, squares[s,])
    qualbit_crop = crop(qualbit, squares[s,])
    
    corine_crop = crop(corine_modis, squares[s,])
    
    
    # Keep only good quality data (reliability=0 or 1) and reproject onto Irish grid
    # and only areas with more than pastureThreshold
    maskRaster = QC_crop>=0 & QC_crop<=minQuality & corine_crop>pastureThreshold
    evi_pasture = mask(evi_crop, maskRaster, maskvalue=FALSE, updatevalue=NA)
    ndvi_pasture = mask(ndvi_crop, maskRaster, maskvalue=FALSE, updatevalue=NA)
    doy_pasture = mask(doy_crop, maskRaster, maskvalue=FALSE, updatevalue=NA)
    QC_pasture = mask(QC_crop, maskRaster, maskvalue=FALSE, updatevalue=NA)
    
    
    # Set negative evi & ndvi to zero
    ndvi_pasture[ndvi_pasture<0] <- 0
    evi_pasture[evi_pasture<0] <- 0
    
    # Create spatial data to convert coords to ITM
    xy_modis <- SpatialPoints(coords = coordinates(ndvi_pasture), 
                              proj4string = modis_crs)
    xy_ITM = spTransform(xy_modis, CRSobj = CRS("+init=epsg:2157"))
    
    coord_info = as.data.frame(cbind(coordinates(xy_ITM), coordinates(xy_modis)))
    names(coord_info) = c('x_ITM','y_ITM','x_MODIS','y_MODIS')
    
    # Create a data frame
    if (f==1) {
      d[[s]] = cbind(coord_info,
                data.frame(satellite = satellite,
                           square = s,
                           year=format(r.file.date[f],"%Y"),
                           doy= getValues(doy_pasture), 
                           evi=getValues(evi_pasture), 
                           ndvi=getValues(ndvi_pasture), 
                           QC=getValues(QC_pasture)))
      d[[s]] = subset(d[[s]], is.finite(evi))
    } else {
      tmp = cbind(coord_info,
                  data.frame(satellite = satellite,
                             square = s,
                             year=format(r.file.date[f],"%Y"),
                             doy= getValues(doy_pasture), 
                             evi=getValues(evi_pasture), 
                             ndvi=getValues(ndvi_pasture), 
                             QC=getValues(QC_pasture)))
      d[[s]] = rbind(d[[s]], subset(tmp, is.finite(evi)))
    }
  }
  if (save_rasters) {
    # Write the rasters to a new file
    fname.ndvi = paste0(outputDir,'/NDVI_',outputSuffix,'_',format(r.file.date[f],"%Y_%m_%d")) 
    fname.evi = paste0(outputDir,'/EVI_',outputSuffix,'_',format(r.file.date[f],"%Y_%m_%d")) 
    fname.doy = paste0(outputDir,'/doy_',outputSuffix,'_',format(r.file.date[f],"%Y_%m_%d")) 
    fname.qc = paste0(outputDir,'/QC_pasture_',format(r.file.date[f],"%Y_%m_%d")) 
    
    # Save data as a geo-Tiff
    writeRaster(ndvi_pasture,file=fname.ndvi,format='GTiff',overwrite=TRUE)
    writeRaster(evi_pasture,file=fname.evi,format='GTiff',overwrite=TRUE)
    writeRaster(doy,file=fname.doy,format='GTiff',overwrite=TRUE)
    writeRaster(QC,file=fname.qc,format='GTiff',overwrite=TRUE)  
  }
}


# Save the data frame for each square in a separate file
for (t in squareList) {
  if (nrow(d[[t]]>0)) {
    # Create a date for each observation
    d[[t]]$date = strptime(paste0(d[[t]]$year,'-',d[[t]]$doy), format="%Y-%j")
    
    d_sq = d[[t]]
    # Save data to a file
    fname.df = file.path(outputDir,paste0('modis_pasture_',yearStr,'_square',t,'.RData') )
    save(d_sq, pastureThreshold, r.file.date, hdf.files, file = fname.df)
  }
}



