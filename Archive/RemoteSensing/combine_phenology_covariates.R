# A script to combine the phenology estimates with the environmental covariates
#
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# March 2022
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
# 
# library(segmented)
library(ggplot2)
library(tidyr)
library(viridisLite)
library(terra)

dataDir = '/media/jon/MODIS_data/PhenologyOutput_toAnalyse/'
envFile = '/media/jon/MODIS_data/Data_created/all_envData.RData'
pastureFile = "/media/jon/MODIS_data/Corine/corine2018_pasturecover_All_Ireland.grd"

input_file_preffix = 'phenology'

# Import data --------
squares = c(2:12,16:18)  # Squares 1, 14  have too little data to use squares 13, 15, 20, 21 in north west removed


# Import all the data
for (i in 1:length(squares)) {
  filename = list.files(path=dataDir, 
                        pattern=paste0(input_file_preffix,"_square_",squares[i],".RData"),
                        full.names = TRUE)
  load(filename)  
  
  if (i==1) {
    phenology_long = phenology
  } else {
    phenology_long = rbind(phenology_long, phenology)
  }
  
  rm(list="phenology")
}
  
# Create a wide version of phenology
phenology_wide = pivot_wider(data=subset(phenology_long, warning==FALSE),
                             id_cols = c('pixelID','year', 'x_MODIS','y_MODIS','x_ITM','y_ITM','square'),
                             names_from = 'phase',
                             names_prefix = 'phase',
                             values_from = c(t),
                             values_fn = mean)

# Centre the x and y coordinates on the whole of Ireland
phenology_wide$x_ITM_centre = phenology_wide$x_ITM - mean(phenology_wide$x_ITM, na.rm=TRUE)
phenology_wide$y_ITM_centre = phenology_wide$y_ITM - mean(phenology_wide$y_ITM, na.rm=TRUE)

# Coordinates
years=unique(phenology_wide$year)
locations = subset(phenology_wide, 
                   subset=year==years[1],
                   select=c('pixelID','x_ITM','y_ITM'))

# Calculate the centre coordinates for each square
tmp = aggregate(cbind(x_ITM_centre, y_ITM_centre)~square, 
                data=phenology_wide, 
                FUN=mean, 
                na.rm=TRUE)

phenology_wide$x_square = NA
phenology_wide$y_square = NA
ind = match(phenology_wide$square, tmp$square)
phenology_wide$x_square = tmp$x_ITM_centre[ind]
phenology_wide$y_square = tmp$x_ITM_centre[ind]
head(phenology_wide)

rm(list='phenology_long')







# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load up environmental data

# Load data on pasture landcover
pasture = terra::rast(pastureFile)
# Add proportion of pasture for each pixel
y = vect(phenology_wide[,c(3,4)], 
         geom=c('x_MODIS','y_MODIS'), 
         crs(pasture))

tmp = terra::extract(pasture, y, xy=FALSE)
phenology_wide$p_pasture = tmp$layer

# Import some environmental data
load(envFile)

# Assign a pixel ID to env data
all_squares_df$pixelID = NA
for (s in unique(phenology_wide$square)) {
  pheno_sub = phenology_wide[phenology_wide$square==s,c(1:7)]
  pixelID_List = unique(pheno_sub$pixelID)
  inds = match(pixelID_List, pheno_sub$pixelID)
  for (p in 1:length(pixelID_List)) {
    ind = abs(all_squares_df$x_MODIS - pheno_sub$x_MODIS[inds[p]])<50 & 
      abs(all_squares_df$y_MODIS - pheno_sub$y_MODIS[inds[p]])<50 
    all_squares_df$pixelID[ind] = pixelID_List[p]
  }
}

nam = names(all_squares_df)
nam[7] = "square"
names(all_squares_df) = nam

d = merge(phenology_wide, all_squares_df)

# Add a categorical aspect variable
d$aspect_cat = NA
delta = 10
d$aspect_cat[d$ASPECT>360-delta | d$ASPECT<delta ] = "N"
d$aspect_cat[d$ASPECT>45 & d$ASPECT<90+45 ] = "E"
d$aspect_cat[d$ASPECT>180-delta & d$ASPECT<180+delta ] = "S"
d$aspect_cat[d$ASPECT>270-delta & d$ASPECT<270+delta ] = "W"
d$aspect_cat = as.factor(d$aspect_cat)

d$year = as.factor(d$year)

# Calculate the average SOS, POS and EOS for each square for each year and then work out an anomaly
# This could be useful to look at the effect of environment (e.g. aspect and slope)
agg_median = aggregate(cbind(x_ITM,y_ITM, phase1,phase2,phase3, ELEVATION, p_pasture)~year+square, 
                data=d, 
                FUN=median, 
                na.rm=TRUE)


agg_mean = aggregate(cbind(x_ITM,y_ITM, phase1,phase2,phase3, ELEVATION, p_pasture)~year+square, 
                       data=d, 
                       FUN=mean, 
                       na.rm=TRUE)


for (s in unique(d$square)) {
  d$anom1[d$square==s] = d$phase1[d$square==s] - agg_median$phase1[match(d$year[d$square==s], agg_median$year[agg_median$square==s])]
  d$anom2[d$square==s] = d$phase2[d$square==s] - agg_median$phase2[match(d$year[d$square==s], agg_median$year[agg_median$square==s])]
  # d$anom3[d$square==s] = d$phase3[d$square==s] - agg$phase3[match(d$year[d$square==s], agg$year[agg$square==s])]
}



save(phenology_wide, d, agg_median, agg_mean,file=file.path(dataDir,'combine_phenology_data_final_report.RData'))




