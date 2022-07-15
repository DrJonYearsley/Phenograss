# Script to process MODIS data for Ireland. 
# This version imports preprocessed CORINE data that 
# gives the fraction of a MODIS pixels that is pasture
# This version processes all of Ireland, not just the 21 10 km squares
#
# This script performs these tasks
#  1. read in MODIS data (ndvi, evi, quality and day of acquisition)
#  2. read in squares (roughty 10 km square) for subsets within Ireland
#  3. read in preprocessed CORINE data
#  For all of Ireland 
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

library(terra)

rm(list=ls())

# setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology')
#setwd('~/Research/Phenograss')
# setwd('~/Research/Phenograss')
setwd('/media/jon/MODIS_data')

# Directories containing the input and output MODIS data 
# inputDir = c('./Data/MODIS/MYD13Q1.006','./Data/MODIS/MOD13Q1.006')
inputDir = c('/media/jon/MODIS_data/MODIS/MYD13Q1.061','/media/jon/MODIS_data/MODIS/MOD13Q1.061')
outputDir = './MODIS_AllIreland'
outputSuffix = 'pasture'
save_rasters = FALSE  # If true save rasters for each MODIS file
yearStr = paste0('A',c(2022)) # Some text (or reg experession) that specifies the year of the data (e.g. 'A20[0-9][0-9]' specifies years 2000-2019) 
minQuality = 1 # Minimum quality to use: 0 = use only best quality pixels, 1=use reasonable pixels
scalingFactor = 0.0001 # Scale factor to apply to NDVI and EVI data from MODIS
corinePath = './Corine'
corineFilename = 'corine2018_pasturecover_All_Ireland.grd'
pastureThreshold = 0.7 # The minimum fraction of a pixel that is pasture



# Load MODIS grid data
modis = rast('./MODIS/modis_grid_ireland.grd')

# Save MODIS CRS
modis_crs = crs(modis, proj=TRUE)

# Read in CORINE preprocessed data.
corine_modis = rast(file.path(corinePath,corineFilename))


# Define extent of Irleand (roughly) in MODIS CRS 
#ir = extent(-7.5E5, -3.3E5,5.7E6, 6.17E6)



extractbit = function(x,n1,n2) {
  # Extract the nth bit from the number x
  return(bitwAnd(bitwShiftR(x,n1),1))
}

for (y in yearStr) {
  # Search for relevant MODIS data files
  regexp = paste(y,'[[:graph:]]+.hdf$',sep='')
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
  
  
  
  for (f in 1:nFiles) {
    print(hdf.files[f])
    
    # Read in the MODIS data and crop to Ireland 
    sds <- rast(hdf.files[f])
    
    if (grepl('MOD13',hdf.files[f])) {
      satellite = 'Terra'
    } else if (grepl('MYD13',hdf.files[f])){
      satellite = 'Aqua'
    } else {
      satellite = NA
    }
    
    # Crop MODIS data to be same as corine
    sds_cropped = crop(sds, corine_modis)
    
    
    
    
    # Keep only good quality data (reliability=0 or 1) and reproject onto Irish grid
    # and only areas with more than pastureThreshold
    maskRaster = sds_cropped$`"250m 16 days pixel reliability"`>=0 & 
      sds_cropped$`"250m 16 days pixel reliability"`<=minQuality & 
      corine_modis>pastureThreshold
    
    evi_pasture = mask(sds_cropped$`"250m 16 days EVI"`, maskRaster, maskvalues=FALSE, updatevalue=NA) * scalingFactor^2
    doy_pasture = mask(sds_cropped$`"250m 16 days composite day of the year"`, maskRaster, maskvalue=FALSE, updatevalue=NA)
    QC_pasture = mask(sds_cropped$`"250m 16 days VI Quality"`, maskRaster, maskvalue=FALSE, updatevalue=NA) 
    
    
    # Set negative evi  to zero
    evi_pasture[evi_pasture<0] <- 0
    
    # Create spatial data to convert coords to ITM
    xy_modis <- vect(crds(evi_pasture, na.rm=FALSE), 
                     crs = modis_crs)
    xy_ITM = project(xy_modis, "epsg:2157")
    
    coord_info = cbind(crds(xy_ITM, df=TRUE), crds(xy_modis, df=TRUE))
    names(coord_info) = c('x_ITM','y_ITM','x_MODIS','y_MODIS')
    
    # Create a data frame
    if (f==1) {
      d = cbind(coord_info,
                data.frame(satellite = satellite,
                           year=format(r.file.date[f],"%Y"),
                           doy= values(doy_pasture, mat=FALSE), 
                           evi=values(evi_pasture, mat=FALSE), 
                           # ndvi=getValues(ndvi_pasture), 
                           QC=values(QC_pasture,mat=FALSE)))
      d = subset(d, is.finite(evi))
    } else {
      tmp = cbind(coord_info,
                  data.frame(satellite = satellite,
                             year=format(r.file.date[f],"%Y"),
                             doy= values(doy_pasture,mat=FALSE), 
                             evi=values(evi_pasture,mat=FALSE), 
                             # ndvi=getValues(ndvi_pasture), 
                             QC=values(QC_pasture,mat=FALSE)))
      d = rbind(d, subset(tmp, is.finite(evi)))
    }
    
    
  }
  
  
  if (save_rasters) {
    # Write the rasters to a new file
    # fname.ndvi = paste0(outputDir,'/NDVI_',outputSuffix,'_',format(r.file.date[f],"%Y_%m_%d")) 
    fname.evi = paste0(outputDir,'/EVI_',outputSuffix,'_',format(r.file.date[f],"%Y_%m_%d"),'.tif') 
    fname.doy = paste0(outputDir,'/doy_',outputSuffix,'_',format(r.file.date[f],"%Y_%m_%d"),'.tif') 
    fname.qc = paste0(outputDir,'/QC_pasture_',format(r.file.date[f],"%Y_%m_%d"),'.tif') 
    
    # Save data as a geo-Tiff
    # writeRaster(ndvi_pasture,file=fname.ndvi,format='GTiff',overwrite=TRUE)
    writeRaster(evi_pasture,file=fname.evi,filetype='GTiff',overwrite=TRUE)
    writeRaster(doy,file=fname.doy,filetype='GTiff',overwrite=TRUE)
    writeRaster(QC,file=fname.qc,filetype='GTiff',overwrite=TRUE)  
    
    
    rm(list=c('doy_pasture','evi_pasture','QC_pasture','sds','sds_cropped','tmp','maskRaster'))
    gc(verbose=FALSE)
  }
  
  
  
  # Save the data frame for modis data
  
  # Create a date for each observation
  d$date = strptime(paste0(d$year,'-',d$doy), format="%Y-%j")
  
  # Save data to a file
  fname.df = file.path(outputDir,paste0('modis_pasture_',y,'.RData') )
  save(d, pastureThreshold, r.file.date, hdf.files, file = fname.df)
  
  rm(list=c('d'))
  gc(verbose=FALSE)
}

