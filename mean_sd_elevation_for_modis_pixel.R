#Calculate mean elevation and standard deviation for every modis pixel
#17.02.2020 
#by Dana Looschelders
#contains test loop for elevation raster splitted into a 100 files
#and code for loop that splits the raster into modis pixel size
#***************************************************************************
#Load libraries
library(raster)
library(rgdal)
library(sp)
library(sf)
library(ggplot2)
library(gstat)
library(devtools)
library(MODIS)
library(MODISTools)
library(gdalUtils)
#library(SpaDES)
#*****************************************************************************
setwd("C:/00 Dana/Uni/Internship/Work/remote sensing/MODIS_Data")

#download shapefile for Ireland
Ireland=getData('GADM',country="IRL",level=0, download=T)

#get MODIS CRS
Ir_ex=raster("modis_grid_ireland.grd")
modis_crs=crs(Ir_ex)

#convert shapefile to modis crs
Ireland_modis=spTransform(Ireland, modis_crs)

#set new working directory for .hgt files
setwd("F:/Digital_Elevation-20200209T175429Z-001/Digital_Elevation") #for all files

#read in one modis file for reference regarding resolution etc
sds=getSds(HdfName = "MOD13Q1.A2016001.h17v03.006.2016029070523.hdf")
modis_layer=raster(readGDAL(sds$SDS4gdal[1], as.is=TRUE))

#read in SRTM.hgt files
# file N33W177.hgt --> lat: 33 to 34 N and lon: 177 to 178 West
hgt.files = list.files(pattern=".hgt") #read in files
dat=sapply(X = hgt.files, FUN = raster) #list files
list2env(x = dat, envir = .GlobalEnv) #unlist files as dataframes in global environment

#get elevation crs
crs_height=crs(dat[[1]])
#crop the modis file to size of Ireland
modis_crop=crop(modis_layer, Ireland_modis)
#mosaic the files as RasterLayer 
name=names(dat) #not an elegant solution but mosaic won't take lists
mosaic_elevation_IE=mosaic(N51W005.hgt, N51W006.hgt, N51W008.hgt, N51W009.hgt, N51W010.hgt, 
                  N51W011.hgt, N52W005.hgt, N52W006.hgt, N52W007.hgt, N52W008.hgt, 
                  N52W009.hgt, N52W010.hgt, N52W011.hgt, N53W005.hgt, N53W006.hgt,
                  N53W007.hgt, N53W008.hgt, N53W009.hgt, N53W010.hgt, N53W011.hgt,
                  N54W005.hgt, N54W006.hgt, N54W007.hgt, N54W008.hgt, N54W009.hgt,
                  N54W010.hgt, N54W011.hgt, N55W005.hgt, N55W006.hgt, N55W007.hgt,
                  N55W008.hgt, N55W009.hgt, fun=mean)

#crop to shapefile of Ireland to reduce size (crop to same size as modis so pixel correspond)
crop_elevation_IE=crop(mosaic_elevation_IE, Ireland)

#get number of rows, columns and cells for cropped modis grid file  
col_modis=ncol(modis_crop)
row_modis=nrow(modis_crop)
cells_modis=ncell(modis_crop)
ncell(modis_crop)==col_modis*row_modis #check that cell numbers are correct

load("sourcecode_splitRaster.r") 
#only needed when SpaDES won't work
#run the code to implement the function splitRaster from the SpaDES package 
#(package won't run on my laptop for whatever reason)

#**************************************************************************
#TEST
#split the subset rasters of a modis layer and the mosaic alevation file 
#in 100 pieces and calculate stats for them 
files2=splitRaster(r=crop_elevation_IE, nx=10, ny=10)
files3=splitRaster(r=modis_crop, nx=10, ny=10)

sampleStratified(x = files2[[80]], size=10)
#prepare dataframe for for loop
nfiles=length(files2)
data=data.frame(mean=rep(NA, 100), sd=rep(NA, 100), min_lon=rep(NA,100), 
                min_lat=rep(NA, 100), x_modis=rep(NA,100), y_modis=rep(NA,100),
                sum_nas=rep(NA, 100))
#use for loop to calculate and assign values for mean, sd, no of NAs, and coordiantes
for (f in 1:nfiles){
   #dummy_SS_mean=sampleStratified(dummy, size=10, na.rm=T, xy=F)
   data[f,1]=mean(values(files2[[f]]), na.rm=T) #mean height value
   data[f,2]=sd(values(files2[[f]]), na.rm=T) #standard deviation of height
   coords=bbox(files2[[f]]) #get bbox coords in lat/lon
   data[f,3]=coords[1,1] 
   data[f,4]=coords[2,1]
   coords_modis=bbox(files3[[f]]) #get bbox coords in modis crs
   data[f,5]=coords_modis[1,1]
   data[f,6]=coords_modis[2,1]
   data[f,7]=sum(is.na(values(files2[[f]]))) #calculate NAs in subset
}
save(data, file="data.rdata") #save dataframe as rdata file

#*****************************************************************************
#split into modis pixel size
files_height=splitRaster(r=crop_elevation_IE, nx=row_modis, ny=col_modis)
files_modis=splitRaster(r=modis_crop, nx=row_modis, ny=col_modis)

#write
nfiles=length(files_height)
data=data.frame(mean=rep(NA, cells_modis), sd=rep(NA, cells_modis), min_lon=rep(NA,cells_modis), 
                min_lat=rep(NA, cells_modis), x_modis=rep(NA,cells_modis), y_modis=rep(NA,cells_modis),
                sum_nas=rep(NA, cells_modis))
for (f in 1:nfiles){
   #dummy_SS_mean=sampleStratified(dummy, size=10, na.rm=T, xy=F)
   data[f,1]=mean(values(files_height[[f]]), na.rm=T) #mean height value
   data[f,2]=sd(values(files_height[[f]]), na.rm=T) #standard deviation of height
   coords=bbox(files_height[[f]]) #get bbox coords in lat/lon
   data[f,3]=coords[1,1] 
   data[f,4]=coords[2,1]
   coords_modis=bbox(files_modis[[f]]) #get bbox coords in modis crs
   data[f,5]=coords_modis[1,1]
   data[f,6]=coords_modis[2,1]
   data[f,7]=sum(is.na(values(files_height[[f]]))) #calculate NAs in subset
}
save(data, file="data.rdata") #save dataframe as rdata file