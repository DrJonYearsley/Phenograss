# Script to process MODIS data for Ireland. 
#
# This script performs these tasks
#  1. read in MODIS data (ndvi, evi, quality and day of acquisition)
#  2. crop to a predefined area, 
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


library(rgdal)
library(raster)
library(gdalUtils)

rm(list=ls())
#setwd("F:/")
setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology/RScript')

# Directories containing the input and output MODIS data 
inputDir = c('../Data/MODIS/MYD13Q1.006','../Data/MODIS/MOD13Q1.006')
outputDir = './MODIS_Testing'
outputSuffix = 'pasture'
save_rasters = FALSE  # If true save rasters for each MODIS file
yearStr = 'A20' # Some text (or reg experession) that specifies the year of the data (e.g. 'A20[0-9][0-9]' specifies years 2000-2019) ###modified to 2018###
corineInclude = c(231)  # Specify corine codes to include (pasture = 231)
minQuality = 1 # Minimum quality to use: 0 = use only best quality pixels, 1=use reasonable pixels
scalingFactor = 0.0001 # Scale factor to apply to NDVI and EVI data from MODIS

corinePath = '../Data/CORINE2018_Ireland/CLC18_IE_ITM'

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

# define a 10 km square within Ireland
sq_10km = extent(-504000,-493000, 5893000, 5904400)



# Find MODIS CRS and reproject CORINE onto MODIS
sds <- get_subdatasets(hdf.files[1])
ndvi = crop(raster(sds[grep("250m 16 days NDVI", sds)], as.is=T), ir)*scalingFactor^2
modis_crs = crs(ndvi)


# Reproject CORINE
corine_modis = spTransform(corine, CRS=modis_crs)
# Extract pasture (code 231 in the vector data)
pasture_modis = subset(corine_modis, CODE_18%in%corineInclude)

for (f in 1:length(hdf.files)) {
  # Read in the MODIS data and crop to Ireland 
  sds <- get_subdatasets(hdf.files[f])
  
 # These lines read in the individual wavelength bands. Used to check the scaling factor for ndvi and evi
  # ndvi is (nir-red)/(nir+red)
  #  red = crop(raster(sds[grep("250m 16 days red reflectance", sds)], as.is=T), ir)
  #  nir = crop(raster(sds[grep("250m 16 days NIR reflectance", sds)], as.is=T), ir)
  
  
  ndvi = crop(raster(sds[grep("250m 16 days NDVI", sds)], as.is=T), sq_10km)*scalingFactor^2
  evi = crop(raster(sds[grep("250m 16 days EVI", sds)], as.is=T),  sq_10km)*scalingFactor^2
  QC = crop(raster(sds[grep("16 days pixel reliability",sds)]),  sq_10km)
  doy = crop(raster(sds[grep("16 days composite day of the year",sds)]), sq_10km)
  # Reliability, 0=good, 1=OK but use with care, 2=snow/icd, 3=cloudy,-1=no data 

  # Keep only good quality data (reliability=0 or 1) and reproject onto Irish grid
  ndvi[QC<0 | QC>minQuality] <- NA 
  evi[QC<0 | QC>minQuality] <- NA
  QC[QC<0 | QC>minQuality] <- NA
  
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
              data.frame(year=format(r.file.date[f],"%Y"),
                         doy= getValues(doy), 
                         evi=getValues(evi_pasture), 
                         ndvi=getValues(ndvi_pasture), 
                         QC=getValues(QC)))
  } else {
    d = rbind(d, 
              cbind(coord_info,
                    data.frame(year=format(r.file.date[f],"%Y"),
                               doy= getValues(doy), 
                               evi=getValues(evi_pasture), 
                               ndvi=getValues(ndvi_pasture), 
                               QC=getValues(QC))  ))
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
d = subset(d, is.finite(QC))

# Save data to a file
fname.df = file.path(outputDir,paste0('modis_pasture_data_',yearStr,'.RData') )
save(d, r.file.date, hdf.files, file = fname.df)

