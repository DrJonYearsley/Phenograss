library(rgdal)
library(raster)
library(sp)
library(tidyverse)
library(ggplot2)
library(plyr)
library(maptools)
library(rgeos)
setwd("F:/Soils_IE_WetDry")
load("soil.new.rdata")

setwd("C:/00 Dana/Uni/Internship/Work/remote sensing/MODIS_Data")
Ir_grid=raster("modis_grid_ireland.grd") #read in modis grid

res.modis=res(Ir_grid)
setwd("F:/Soils_IE_WetDry")
moist=readOGR("Soils_IE_WetDry.shp")

crs.moist=crs(moist)
ex=extent(moist)

#subset dataframe to Carlow for simplicity
sub_carlow <- subset(moist, COUNTY == "CARLOW")
ex.carlow=extent(sub_carlow)

#transformation to visualize with ggplot
sub_carlow@data$id = rownames(sub_carlow@data)
carlow.points = fortify(sub_carlow, region="id")

carlow.df = join(carlow.points, sub_carlow@data, by="id")

ggplot(carlow.df) + 
  aes(long,lat,group=group,fill=CATEGORY) + 
  geom_polygon() +
  coord_equal() +
  scale_fill_brewer("Soil moisture",palette ="Set1" )

#*****************************************************
#still in progress:
#****************************************************
#get soil type for every modis pixel
#create a raster with modis resolution (use Ir_grid) and extent of soil
r=raster(ext=ex.carlow, res=res.modis) 
rows=nrow(r)
cols=ncol(r)

#ideas:
#1)
#use aggregate to create raster with mode of soil category
agg_carlow=aggregate(sub_carlow["CATEGORY"], by=r, FUN=mode) #use aggregate with mode to get most abundent soil type
agg_carlow=aggregate(sub_carlow, by=sub_carlow@data$CATEGORY,fact=c(rows, cols), FUN=mode)
agg_carlow=aggregate(sub_carlow, by="CATEGOY", sums=list(list(mode, "CATEGORY")))

plot(agg_carlow)
#2)
#aggregate dataframe: create a dataframe and add column to indicate modis pixel (with number)

#3)
#use rasterize to create raster object from soil data and use same method as with elevation data (extremely slow)
?rasterize
sub_carlow@data=sub_carlow@data$CATEGORY
data=sub_carlow@data
sub_carlow@data=NULL
  as.integer((sub_carlow@data$CATEGORY))
r=raster(nrow=10, ncol=10)
raster=rasterize(sub_carlow, r, update=TRUE, 
                 field=sub_carlow@data[,7], fun=mean, updateValue="NA")
rasterr=rasterize(sub_carlow, r)
class(sub_carlow)

#4.
#use dataframe for xyz to raster -> doesn't work as is grid is not regular

#for the whole of Ireland
#visualization
moist@data$id = rownames(moist@data)
moist.points = fortify(moist, region="id")

moist.df = join(moist.points, moist@data, by="id")

ggplot(moist.df) + 
  aes(long,lat,group=group,fill=CATEGORY) + 
  geom_polygon() +
  coord_equal() +
  scale_fill_brewer("Soil moisture",palette ="Set1" )

save.image(file="soil.new.rdata")