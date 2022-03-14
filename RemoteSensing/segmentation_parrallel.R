# Segmentation modelling of phenology
#
# Apply segmentation method developed by Gabin
# Run it in parrallel across the pixels from a square
#
# Jon Yearsley Jon.Yearsley@ucd.ie
# Aug 2020
# ****************************************

setwd('~/git_repos/Phenograss/RemoteSensing/')

rm(list=ls())

# dataDir = '~/Research/Phenograss/Data/MODIS_squares'
# outputDir = '~/Research/Phenograss/Data/PhenologyOutput_test/'
dataDir = '/media/jon/MODIS_data/MODIS_squares'
outputDir = '/media/jon/MODIS_data/PhenologyOutput_test/'


library(mgcv)
library(segmented)
library(ggplot2)
library(foreach)
library(doParallel)


# Register cluster with 2 nodes
cl<-makeCluster(10)
registerDoParallel(cl)


# Set basic parameters
square = c(13:21)
square = c(1:21)
# square = 20


year = 2015
knots = -1
min_obs = 15      # Minimum number of rows for trying to segment data
#starting_breaks = c(50, 100, 200, 300)  # Initial guess for break points in doy
nSegBreaks = 5    # Number of breakpoints to use for segmentation
print_pixel=FALSE


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++ Start of function definitions ++++++++++++++++++++++++++

# Import some functions defines in segmentation_functions.R
source("segmentation_functions.R")


# Create data frame to contain parameters
params = data.frame(knots = knots, min_obs=min_obs, nSegBreaks=nSegBreaks)

# Parallel execution function ----------
# Function to package segmentation input and output so 
#    it's suitable for the foreach loop
segment_within_pixel = function(input_data, pixel_info, output_labels, params) {
  
  # Perform segmentation on the input data
  segments = segmentEVI(input_data, 
                        nBreaks = params$nSegBreaks,
                        knots=params$knots, 
                        use_raw = TRUE,
                        sd_filter_threshold = 6)
  
  # Now put the segmented data into a single line
  modelStr = c('raw', 'smooth')
  
  for (model in c(1,2)) {
    # Model 1 corresponds to using the raw data, model 2 corresponds to using the smoothed data
    
    # Record estimates if they have been made 
    if (!is.null(segments[[model]])) {
      # Calculate estimated segment break dates and 95% CI
      CI = confint(segments[[model]])
      
      # calculate slopes
      slope_est = slope(segments[[model]])$doy[,1]
      
      
      # Record estimates (break point and following slope)
      tmp = data.frame(pixelID = pixel_info$pixelID,
                            x_ITM = pixel_info$x_ITM,
                            y_ITM = pixel_info$y_ITM,
                            x_MODIS = pixel_info$x_MODIS,
                            y_MODIS = pixel_info$y_MODIS,
                            year = pixel_info$year,
                            square = pixel_info$square,
                            model=modelStr[model],
                            phase = c(1:params$nSegBreaks), 
                            t=CI[,1],
                            lowerCI=CI[,2],
                            upperCI=CI[,3],
                            slope=slope_est[-1],
                            slopeprev=slope_est[1:params$nSegBreaks])

    } else {
      tmp = data.frame(pixelID = pixel_info$pixelID,
                       x_ITM = pixel_info$x_ITM,
                       y_ITM = pixel_info$y_ITM,
                       x_MODIS = pixel_info$x_MODIS,
                       y_MODIS = pixel_info$y_MODIS,
                       year = pixel_info$year,
                       square = pixel_info$square,
                       model=modelStr[model],
                       phase = c(1:params$nSegBreaks), 
                       t=NA,
                       lowerCI=NA,
                       upperCI=NA,
                       slope=NA,
                       slopeprev=NA)
    }
    
    # Attach tmp of the two models together
    if (model==1) {
      output = tmp
    } else {
      output = rbind(output, tmp)
    }
  }
  
  return(output)    
}
# End of parallel execution loop


# +++++++++++++ End of function definitions ++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++









# Loop around 10 km squares and perform segmentation for each pixel time series  --------

for (s in square) {
  print(paste0('Analysing square ',s,' at ',Sys.time()))
  
  # Import data...
  # Load the focal year and years either side
  file_prev = list.files(pattern=paste0(year-1,'_square',s,'.RData'),
                         path = dataDir,
                         full.names = TRUE)
  file = list.files(pattern=paste0(year,'_square',s,'.RData'),
                    path = dataDir,
                    full.names = TRUE)
  file_next = list.files(pattern=paste0(year+1,'_square',s,'.RData'),
                         path = dataDir,
                         full.names = TRUE)
  if (length(file)!=1 | length(file_prev)!=1 | length(file_next)!=1) {
    error('File names not found correctly')
  }
  
  # Load data from previous year after day 265
  load(file_prev)
  d_prev = d_sq[d_sq$doy>265,]
  d_prev$doy = 265 - d_prev$doy
  
  # Load data from previous year before day 100
  load(file_next)
  d_next = d_sq[d_sq$doy<100,]
  d_next$doy = 365 + d_next$doy
  
  # Load data from the focal year
  load(file)
  
  # Combine 3 year's data into one data frame
  d_final = rbind(d_prev, d_sq, d_next)
  
  # Add unique pixel ID
  x_values = unique(d_final$x_MODIS)
  y_values = unique(d_final$y_MODIS)
  
  d_final$pixelID = paste0('x',match(d_final$x_MODIS,x_values),
                           'y',match(d_final$y_MODIS,y_values))
  
  
  
  # Remove scaling factor from evi if evi is too small
  scalingFactor = 0.0001
  if (mean(d_final$evi, na.rm=T)<1e-6) {
    d_final$evi = d_final$evi/scalingFactor^2
    d_final$ndvi = d_final$ndvi/scalingFactor^2
  }
  
  
  
  
  
  
  # *****************************************************
  # Perform segmentation on all pixels in the square ----
  
  pixel_list = unique(d_final$pixelID)
  nPixel = length(pixel_list)
  
  # Create data frame with pixel info
  ind = match(pixel_list, d_final$pixelID)
  pixel_data = data.frame(pixelID =  d_final$pixelID[ind],
                          x_ITM = d_final$x_ITM[ind],
                          y_ITM = d_final$y_ITM[ind],
                          x_MODIS = d_final$x_MODIS[ind],
                          y_MODIS = d_final$y_MODIS[ind])
  
  
  
  # Define output labels and pixel data frame
  # Labels for the final data frame variables
  # t_ is the time of the segmentation
  # slope_ is the regression slope after t_
  # slopeprev_ is the regression slope before t_
  label = c('t_', 'lowerCI_', 'upperCI_','slope_','slopeprev_')
  modelStr2 = c('_raw','_smooth')
  
  output_labels = paste0(rep(paste0(rep(label,times=nSegBreaks),
                                         rep(c(1:nSegBreaks),each=length(label))), times=2),
                              rep(modelStr2,each=nSegBreaks*length(label)))
  
  
  # Create data frame for the output
  pixel_data$year = year
  pixel_data$square = s
  
  
  
  
  # Parallelize the segmentation of the EVI data for each pixel in a square  ------
  segment_output <- foreach (p = icount(nPixel), .packages=c('segmented','mgcv'), .inorder=FALSE, .combine='rbind') %dopar% {
   segment_within_pixel(d_final[d_final$pixelID==pixel_list[p],], 
                        pixel_data[pixel_data$pixelID==pixel_list[p],],
                        label, params)
  }

  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Finished estimating segmentation for all pixels
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  
  
  
  
  
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Loop around all pixels and move the phenophases ----
  # so that: 
  #    phase 1 
  #         starts in the year 
  #         at the first positive slope
  #    phase 2
  #         is the first maximum (i.e. slope is positive before and negative after)
  #    phase 3
  #         is within the year
  #         has a negative slope before
  #         a slope after that is less steep (or maybe positive) compared to before
  
  
  
  segment_output$phase = NA            # Remove draft phase information
  segment_output$warning = NA
  
  nModels = length(unique(segment_output$model))
  
  for (m in c(1:nModels)) {  # Loop around the two models (raw and smoothed)
    for (p in 1:nPixel) {
      # Estimate phenophases for smoothed data
      ind = which(segment_output$pixelID==pixel_list[p] & segment_output$model==unique(segment_output$model)[m])
      if (all(is.finite(segment_output$t[ind]))) {
        
        # Call function to assign phenophases
        phenophase = assignPhenophase(segment_output[ind,])
        
        # Put the phenophases into the final data frame
        segment_output$phase[ind] = phenophase$phase
        segment_output$warning[ind] = phenophase$warning
      }
    }
  }
  
  # Flag phenophases with very wide 95% CI
  ind = (segment_output$upperCI - segment_output$lowerCI)>30
  segment_output$wideCI = ind
  
  # Save data -------------------
  filename = paste0('phenology_square_',s,'_',year,'.RData')
  save(knots,year,s,nSegBreaks,segment_output,d_final,
       file=file.path(outputDir,filename))
} # Finish looping around all the squares


# Disconnect the compute nodes
stopCluster(cl) 
  
  

