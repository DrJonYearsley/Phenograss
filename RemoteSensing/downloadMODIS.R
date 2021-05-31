# Use MODIS pacakge to download data on BDRF and test NBAR calculations
#
# This MODIS package will only download files that are required. 
# Previously downloaded files will not be re-downloaded

# Jon Yearsley (jon.yearsley@ucd.ie)
# Aug 2019
# *****************************************


#setwd('~/WorkFiles/PeopleStuff/GrasslandPhenology/Data')
setwd('~/MEGAsync/Projects/GrasslandPhenology/Data')
library(MODIS)
library(rgdal)
library(gdalUtils)

period = as.Date(c("2020/11/01", "2021/03/01"))

tile = c(17,03)  # The horizontal and vertical tile coords for Ireland (h17v03)

## Change options for MODIS package
# # Permanently set MODIS options
MODISoptions(localArcPath = "~/Research/Phenograss/Data/MODIS/",
             MODISserverOrder = c( "LPDAAC"))
# 
# # Move the location of the data archive
# orgStruc(to = '~/WorkFiles/PeopleStuff/GrasslandPhenology/Data', move=TRUE)
# orgStruc(to = '~/Research/Phenograss/Data/', move=TRUE)

# # Download MODIS data for Ireland

# Albedo parameters
# a = getHdf(product='MCD43A1', begin=period[1], end=period[2],tileH=tile[1], tileV=tile[2], quiet=FALSE)

# 250m 8 day surface reflectances
# b = getHdf(product='M.D09Q1', begin=period[1], end=period[2],tileH=tile[1], tileV=tile[2], quiet=FALSE)

# 500m NADIR corrected
# c = getHdf(product='MCD43A4', begin=period[1], end=period[2],tileH=tile[1], tileV=tile[2], quiet=FALSE)

# 1km satellite viewing angles
#d = getHdf(product='M.D09GA', begin=period[1], end=period[2],tileH=tile[1], tileV=tile[2], quiet=FALSE)

# 250m 16 day composite VI data
e = getHdf(product='M.D13Q1', begin=period[1], end=period[2],tileH=tile[1], tileV=tile[2], quiet=FALSE)

# # Load raster MODIS HDF file
# 
# # Loaf albedo parameters
# sds1 <- get_subdatasets('MODIS/MCD43A1.006/2018.05.09/MCD43A1.A2018129.h17v03.006.2018157004043.hdf')
# params = readGDAL(sds1[11])
# 
# 
# # Load MOD09 surface reflectances
# sds2 <- get_subdatasets('MODIS/MOD09Q1.006/2018.05.09/MOD09Q1.A2018129.h17v03.006.2018151145036.hdf')
# band1 = readGDAL(sds2[1])
# 
# # Load MCD nbar surface reflectances
# sds3 <- get_subdatasets('MODIS/MCD43A4.006/2018.05.09/MCD43A4.A2018129.h17v03.006.2018157004043.hdf')
# nbar_band1 = readGDAL(sds3[8])
# 
# k_iso = d['band1']
# k_vol = d['band2']
# k_geo = d['band3']

