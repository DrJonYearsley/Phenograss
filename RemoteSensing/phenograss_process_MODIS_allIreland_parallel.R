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


rm(list=ls())


library(terra)
library(foreach)
library(doParallel)


# Register cluster with 2 nodes
cl<-makeCluster(5)
registerDoParallel(cl)


# setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology')
#setwd('~/Research/Phenograss')
# setwd('~/Research/Phenograss')
setwd('/media/jon/MODIS_data')

# Directories containing the input and output MODIS data 
# inputDir = c('./Data/MODIS/MYD13Q1.006','./Data/MODIS/MOD13Q1.006')
inputDir = c('/media/jon/MODIS_data/MODIS/MYD13Q1.006','/media/jon/MODIS_data/MODIS/MOD13Q1.006')
outputDir = './MODIS_AllIreland'
outputSuffix = 'pasture'
yearStr = 'A2008' # Some text (or reg experession) that specifies the year of the data (e.g. 'A20[0-9][0-9]' specifies years 2000-2019) 
minQuality = 1 # Minimum quality to use: 0 = use only best quality pixels, 1=use reasonable pixels
scalingFactor = 0.0001 # Scale factor to apply to NDVI and EVI data from MODIS
corinePath = './Corine'
corineFilename = 'corine2018_pasturecover_All_Ireland.grd'
pastureThreshold = 0.7 # The minimum fraction of a pixel that is pasture

# Read in CORINE preprocessed data.
corine_modis = rast(file.path(corinePath,corineFilename))
corine_mask = corine_modis>pastureThreshold

# Load MODIS grid data
modis = crop(rast('./MODIS/modis_grid_ireland.grd'), corine_modis)

# Save MODIS CRS
modis_crs = crs(modis, proj=TRUE)


# Create spatial data to convert coords to ITM
xy_modis <- vect(crds(modis, na.rm=FALSE), 
                 crs = modis_crs)
xy_ITM = project(xy_modis, "epsg:2157")

coord_info = cbind(crds(xy_ITM, df=TRUE), crds(xy_modis, df=TRUE))
names(coord_info) = c('x_ITM','y_ITM','x_MODIS','y_MODIS')


rm(list=c('modis','xy_modis','xy_ITM','corine_modis'))

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



import_modis = function(filename, modis_crs, corine_mask, scalingFactor, coord_info) {

  print(hdf.files[f])
  flush.console()
  
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
  sds_cropped = crop(sds, corine_mask)
  rm(list='sds')
  
  # Keep only good quality data (reliability=0 or 1) and reproject onto Irish grid
  # and only areas with more than pastureThreshold
  maskRaster = corine_mask & sds_cropped$`"250m 16 days pixel reliability"`>=0 & 
    sds_cropped$`"250m 16 days pixel reliability"`<=minQuality 
  
  evi_pasture = mask(sds_cropped$`"250m 16 days EVI"`, maskRaster, maskvalues=FALSE, updatevalue=NA) * scalingFactor^2
  doy_pasture = mask(sds_cropped$`"250m 16 days composite day of the year"`, maskRaster, maskvalue=FALSE, updatevalue=NA)
  QC_pasture = mask(sds_cropped$`"250m 16 days VI Quality"`, maskRaster, maskvalue=FALSE, updatevalue=NA) 
  rm(list=c('sds_cropped'))
  
  # Set negative evi  to zero
  evi_pasture[evi_pasture<0] <- 0
  
  
  
  # Create a data frame
    dout = cbind(coord_info,
              data.frame(satellite = satellite,
                         year=format(r.file.date[f],"%Y"),
                         doy= values(doy_pasture, mat=FALSE), 
                         evi=values(evi_pasture, mat=FALSE), 
                         QC=values(QC_pasture,mat=FALSE)))
    dout = subset(dout, is.finite(evi))

    print(paste("Finished with ",hdf.files[f]))
    flush.console()
    
    return(dout)
}





# Parallelize the segmentation of the EVI data for each pixel in a square  ------
d <- foreach (f = icount(nFiles), .packages=c('terra'), .inorder=FALSE, .combine='rbind') %dopar% {
  import_modis(hdf.files[f], modis_crs, corine_mask, scalingFactor, coord_info)
}




# Save the data frame for modis data

# Create a date for each observation
d$date = strptime(paste0(d$year,'-',d$doy), format="%Y-%j")

# Save data to a file
fname.df = file.path(outputDir,paste0('modis_pasture_',yearStr,'.RData') )
save(d, pastureThreshold, r.file.date, hdf.files, file = fname.df)


# Disconnect the compute nodes
stopCluster(cl) 
