# Script to process MODIS data for Ireland. 
#
# This script performs these tasks
#  1. read in MODIS data (ndvi, evi, quality and day of acquisition)
#  2. crop to a predefined area (square define in Irish Grid TM75)
#  3. rescale NDVI and EVI 
#  4. remove poor quality pixels
#  5. read in CORINE (vector product)
#  6. reproject CORINE onto MODIS CRS
#  7. subset MODIS data to be only pixels over pasture (CORINE==231)
#  8. saves the processed raster data for NDVI, EVI and the quality data
#  9. convert MODIS data into a data frame and save data frame
#  10. Make data frame into a spatial object
#  11. Reproject spatial data frame into ITM (espg:2157)

# The final data that is saved to file will have a spatial resolution of approximately 250m 
#        (i.e. the raw MODIS resolution, albeit transformed onto Irish Grid TM75)
#
# Jon Yearsley Nov 2019 
# ****************************************************

library(sp)
library(rgdal)
library(raster)
library(gdalUtils)

rm(list=ls())
setwd('F:/') ##### The directory is the only thing we have different #####

# Directories containing the input and output MODIS data 
inputDir = c('../Data/MODIS/MYD13Q1.006','../Data/MODIS/MOD13Q1.006')
outputDir = '../Data/MODIS'
outputSuffix = 'pasture'
save_rasters = FALSE  # If true save rasters for each MODIS file
yearStr = 'A2017' # Some text (or reg experession) that specifies the year of the data (e.g. 'A20[0-9][0-9]' specifies years 2000-2019) ###modified to 2018###
corineInclude = c(231)  # Specify corine codes to include (pasture = 231)
minQuality = 1 # Minimum quality to use: 0 = use only best quality pixels, 1=use reasonable pixels
scalingFactor = 0.0001 # Scale factor to apply to NDVI and EVI data from MODIS

####
# Define the square quadrat for DUBLIN and the size of the square
square_centre_TM75 = data.frame(east=317500, north=257100) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for WEXFORD and the size of the square
square_centre_TM75 = data.frame(east=288825, north=122375) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for WICKLOW and the size of the square
square_centre_TM75 = data.frame(east=320175, north=184075) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for OFFALY and the size of the square
square_centre_TM75 = data.frame(east=213500, north=200100) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for ROSCOMMON and the size of the square
square_centre_TM75 = data.frame(east=192050, north=243650) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for TIPPERARY and the size of the square
square_centre_TM75 = data.frame(east=224525, north=138825) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for CAVAN and the size of the square
square_centre_TM75 = data.frame(east=246200, north=297425) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for LONGFORD and the size of the square
square_centre_TM75 = data.frame(east=229650, north=278475) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for WESTMEATH and the size of the square
square_centre_TM75 = data.frame(east=244975, north=246650) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for ANTRIM and the size of the square
square_centre_TM75 = data.frame(east=316675, north=398150) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for ARMAGH and the size of the square
square_centre_TM75 = data.frame(east=292400, north=346025) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for DERRY and the size of the square
square_centre_TM75 = data.frame(east=267250, north=414850) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for DONEGAL and the size of the square
square_centre_TM75 = data.frame(east=219550, north=418950) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for FERMANAGH and the size of the square
square_centre_TM75 = data.frame(east=193150, north=356425) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for TYRONE and the size of the square
square_centre_TM75 = data.frame(east=222225, north=382600) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for CORK and the size of the square
square_centre_TM75 = data.frame(east=134700, north=058975) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for KERRY and the size of the square
square_centre_TM75 = data.frame(east=085300, north=100375) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for LIMERICK and the size of the square
square_centre_TM75 = data.frame(east=156100, north=146100) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for Galway and the size of the square
square_centre_TM75 = data.frame(east=136150, north=244250) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for MAYO and the size of the square
square_centre_TM75 = data.frame(east=103700, north=283775) # Define in Irish Grid TM75 (espg:29903)

# Define the square quadrat for SLIGO and the size of the square
square_centre_TM75 = data.frame(east=157150, north=320900) # Define in Irish Grid TM75 (espg:29903)

square_size_m = 10000  # units = m. Each pixel is roughly 250 m
####

# Define location of CORINE data
corinePath = '../Data/CORINE_Ireland/CLC18_IE_ITM'

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
ndvi = crop(raster(sds[grep("250m 16 days NDVI", sds)], as.is=T), ir)*scalingFactor^2
modis_crs = crs(ndvi)

sq_df = data.frame(east = square_size_m * c(-0.5, -0.5, 0.5, 0.5) + square_centre_TM75$east, 
                   north = square_size_m * c(-0.5, 0.5, 0.5, -0.5)+ square_centre_TM75$north)

sq = SpatialPoints(sq_df, proj4string=CRS("+init=epsg:29903"))

# Convert square to MODIS CRS and extract extent
sq_modis = extent(spTransform(sq, CRS=modis_crs))


# Reproject CORINE
corine_modis = spTransform(corine, CRS=modis_crs)
# Extract pasture (code 231 in the vector data)
pasture_modis = subset(corine_modis, CODE_18%in%corineInclude)

extractbit = function(x,n1,n2) {
  # Extract the nth bit from the number x
  return(bitwAnd(bitwShiftR(x,n1),1))
}

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
  
  
  ndvi = crop(raster(sds[grep("250m 16 days NDVI", sds)], as.is=T), sq_modis)*scalingFactor^2
  evi = crop(raster(sds[grep("250m 16 days EVI", sds)], as.is=T),  sq_modis)*scalingFactor^2
  QC = crop(raster(sds[grep("16 days pixel reliability",sds)]),  sq_modis)
  qualbit = crop(raster(sds[grep("16 days VI Quality",sds)]),  sq_modis)
  doy = crop(raster(sds[grep("16 days composite day of the year",sds)]), sq_modis)
  # Reliability, 0=good, 1=OK but use with care, 2=snow/icd, 3=cloudy,-1=no data 

  # Extract a few quality controls 
  shadow = calc(qualbit, fun=function(x) {floor(x/2^15) %% 2})
  
  # Keep only good quality data (reliability=0 or 1) and reproject onto Irish grid
  maskRaster = QC>=0 & QC<=minQuality
  evi = mask(evi, maskRaster, maskvalue=FALSE, updatevalue=NA)
  ndvi = mask(ndvi, maskRaster, maskvalue=FALSE, updatevalue=NA)
  QC = mask(QC, maskRaster, maskvalue=FALSE, updatevalue=NA)


  # Set negative evi & ndvi to zero
  ndvi[ndvi<0] <- 0
  evi[evi<0] <- 0
  


  # Extract cells corresponding to pasture in CORINE
  ndvi_pasture = mask(ndvi, pasture_modis)
  evi_pasture = mask(evi, pasture_modis)
  
  # Create a data frame
  if (f==1) {
    # Create spatial data to convert coords to ITM
    xy_modis <- SpatialPoints(coords = coordinates(ndvi_pasture), 
                              proj4string = modis_crs)
    xy_ITM = spTransform(xy_modis, CRSobj = CRS("+init=epsg:2157"))
    
    coord_info = as.data.frame(cbind(coordinates(xy_ITM), coordinates(xy_modis)))
    names(coord_info) = c('x_ITM','y_ITM','x_MODIS','y_MODIS')
    
    
    d = cbind(coord_info,
              data.frame(satellite = satellite,
                         year=format(r.file.date[f],"%Y"),
                         doy= getValues(doy), 
                         evi=getValues(evi_pasture), 
                         ndvi=getValues(ndvi_pasture), 
                         QC=getValues(QC)))
    d = subset(d, is.finite(evi))
  } else {
    tmp = cbind(coord_info,
                data.frame(satellite = satellite,
                           year=format(r.file.date[f],"%Y"),
                           doy= getValues(doy), 
                           evi=getValues(evi_pasture), 
                           ndvi=getValues(ndvi_pasture), 
                           QC=getValues(QC))  )
    d = rbind(d, subset(tmp, is.finite(evi)))
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

# Remove missing data
d = subset(d, is.finite(evi))


# Create a date for each observation
d$date = strptime(paste0(d$year,'-',d$doy), format="%Y-%j")



# Save data to a file
fname.df = file.path(outputDir,paste0('modis_pasture_data_',yearStr,'.RData') )
save(d, r.file.date, hdf.files, file = fname.df)

