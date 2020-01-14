# create_quadrat_polygons.R
#
# Reads in the eatsings and northings (TM75) of centre of 
# each quadrat, converts coordinates to MODIS CRS and then 
# defines a 10 km (roughly) box around the centre. 
# The corners of this box are saved as a shapefile in the MODIS CRS.
# 
# Jon Yearsley  jon.yearsley@ucd.ie
# Jan 2020
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())


library(sp)
library(rgdal)
library(raster)
library(readxl)

rm(list=ls())

setwd('/home/jon/WorkFiles/PeopleStuff/GrasslandPhenology')
outputDir = './Data/Quadrats'

# Quadrat definitions file
quadratFile = './Data/Quadrats/QuadratsAgroclimateRegionsIreland.xlsx'

# Define a square quadrat and the size of the square
square_size_m = 10000  # units = m. Each pixel is roughly 250 m

# Read in data  on quadrats
quadrats = read_xlsx(path=quadratFile, 
                     skip=14)

quadrats$Easting = as.numeric(quadrats$Easting)
quadrats$Northing = as.numeric(quadrats$Northing)
nQuadrat = nrow(quadrats)



# Load MODIS grid data
modis = raster('./Data/MODIS/modis_grid_ireland.grd')

# Save MODIS CRS
modis_crs = crs(modis)


# Define spatial data frame for quadrat centers
sq_df = SpatialPointsDataFrame(coords = quadrats[,c(3,4)], 
                               data=quadrats, 
                               proj4string=CRS("+init=epsg:29903"))
sq_df_modis = spTransform(sq_df, modis_crs)


for (q in 1:nQuadrat) {
  sq_width = SpatialPoints(data.frame(east=quadrats$Easting[q]+c(0.5,-0.5)*square_size_m,
                                      north=quadrats$Northing[q])
                           , proj4string=CRS("+init=epsg:29903"))
  sq_height = SpatialPoints(data.frame(east=quadrats$Easting[q],
                                       north=quadrats$Northing[q]+c(0.5,-0.5)*square_size_m)
                            , proj4string=CRS("+init=epsg:29903"))
  
  # Convert square and side length into  MODIS
  sq_width_modis = coordinates(spTransform(sq_width, modis_crs))[,1]
  sq_height_modis = coordinates(spTransform(sq_height, modis_crs))[,2]
  
  # Create quadrat in the MODIS CRS
  quadrat_points_modis = data.frame(x=c(sq_width_modis, rev(sq_width_modis)),
                                    y=rep(sq_height_modis, each=2))
  
  sq_name = paste0('quadrat',q,'_',quadrats$County[q],'_Clust',quadrats$Cluster[q])
  
  # Define final points
  # sq_modis = SpatialPointsDataFrame(quadrat_points_modis, 
  #                                   data=data.frame(coordinates(quadrat_points_modis)), 
  #                                   proj4string=modis_crs)
  
  # Define spatial polygon for the quadrat
  poly = SpatialPolygons(list(Polygons(list(Polygon(coords=quadrat_points_modis)),
                                       ID=sq_name)),
                         proj4string=modis_crs)
  
  sq_modis2 = SpatialPolygonsDataFrame(poly, data=quadrats[1,], match.ID=FALSE)
  
  
  writeOGR(sq_modis2, dsn=file.path(outputDir,paste0(sq_name,'.shp')), 
           layer=sq_name, 
           driver='ESRI Shapefile',
           overwrite_layer=TRUE) 
}


