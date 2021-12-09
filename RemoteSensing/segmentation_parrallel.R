# Segmentation modelling of phenology
#
# Apply segmentation method developed by Gabin
# Run it in parrallel across the pixels from a square
#
# Jon Yearsley Jon.Yearsley@ucd.ie
# Aug 2020
# ****************************************

rm(list=ls())

dataDir = '~/Research/Phenograss/Data/MODIS_squares'
outputDir = '~/Research/Phenograss/Data/PhenologyOutput_test/'


library(mgcv)
library(segmented)
library(ggplot2)
library(foreach)
library(doParallel)


# Register cluster with 2 nodes
cl<-makeCluster(2)
registerDoParallel(cl)


# Set basic parameters
square = c(1:9,13:21)
square = c(10:12)
square = c(1:21)
square = 5

year = 2014
knots = -1
min_obs = 15      # Minimum number of rows for trying to segment data
#starting_breaks = c(50, 100, 200, 300)  # Initial guess for break points in doy
nSegBreaks = 4    # Number of breakpoints to use for segmentation
print_pixel=FALSE


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++ Start of function definitions ++++++++++++++++++++++++++


# *************************************
# Function to perform segmentation on the whole time series ----
segmentEVI = function(d_sub, 
                      nBreaks = 4, 
                      useHighQuality=FALSE, 
                      knots=-1,
                      use_raw = FALSE,
                      sd_filter_threshold=4) {
  # This functions fits a segmented linear model to raw data and smoothed data
  # the smoothed data gives fewer NA's and more consistent results. 
  
  require(segmented)
  require(mgcv)
  
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++
  # Perform segmentation on data smoothed using a GAM
  
  
  
  
  
  m_gam = tryCatch(gam(evi~s(doy, bs='cr',k=knots), 
                       data=d_sub,
                       gamma=1),  # Gamma controls over/under fitting
                   error = function(e) {
                     NULL
                   },
                   warning = function(e) {
                     NULL
                   })
  
  
  if (!is.null(m_gam)) {
    # Remove EVI data points that are more than 6 SE below from the prediction from m_gam
    tmp = predict(m_gam, se.fit = TRUE, newdata = d_sub)
    evi_deviation = ((d_sub$evi-tmp$fit)/tmp$se.fit)
    filter_ind = evi_deviation>-sd_filter_threshold
    # # Option to visualise the filtered data
    # plot(d_sub$doy, d_sub$evi)
    # points(d_sub$doy[filter_ind], d_sub$evi[filter_ind], pch=20)
    
    
    # Smooth data after removing data more than sd_filter_threshold se below prediction  
    m_gam2 = gam(evi~s(doy, bs='cr',k=knots), 
                 data=d_sub[filter_ind,],
                 gamma=1)  # Gamma controls over/under fitting
    
    # Add smoothed predictions to the data frame
    d_sub$evi_smooth = NA
    d_sub$evi_smooth[filter_ind] = predict(m_gam2)
    
    
    # Segmenting the smoothed predictions
    m_smooth = lm(evi_smooth~doy, data=d_sub[filter_ind,])
    
    m_seg_smooth = tryCatch(segmented(m_smooth, 
                                      seg.Z = ~doy,
                                      npsi = nBreaks,
                                      control=seg.control(it.max=50, 
                                                          fix.npsi=TRUE, 
                                                          n.boot=15, 
                                                          display=FALSE)),
                            error = function(e) {
                              NULL
                            },
                            warning = function(e) {
                              NULL
                            })
  } else {
    m_seg_smooth = NULL
  }
  
  
  
  if (use_raw) {
    # +++++++++++++++++++++++++++++++++++++++++++++++++
    # Perform segmentation on raw evi data 
    
    m_raw = lm(evi~doy, data=d_sub[filter_ind,])
    
    m_seg_raw = tryCatch(segmented(m_raw, 
                                   seg.Z = ~doy,
                                   npsi=nBreaks,
                                   control=seg.control(it.max=50, 
                                                       fix.npsi=TRUE, 
                                                       n.boot=15, 
                                                       display=FALSE)),
                         error = function(e) {
                           NULL
                         },
                         warning = function(e) {
                           NULL
                         })
  } else {
    m_seg_raw = NULL
  }       
  
  
  return(list(m_seg_raw, m_seg_smooth, filtered=filter_ind, d_sub=d_sub, m_gam2=m_gam2)) 
}



# Function to assign phenophases by applying conditions for SOS, POS and EOS -------
assignPhenophase = function(input_data) {
  # A function to assign phenophases based upon criteria for Start of season (SOS), 
  # peak of season (POS) and end of season (EOS)
  
  # input data is a data frame with variable:
  #     t            time of break point estimate  (doy)
  #     slope        slope of segment following the breakpoint
  #     slopeprev    slope of segment preceeding the breakpoint
  
  # output_data is a data frame with variables:
  #     phase        estimated phenophases (1=SOS, 2=POS, 3=EOS)
  #     warning      TRUE if SOS is after 1st June or EOS after end of year
  
  
  
  output_data = data.frame(phase = rep(NA, times=nrow(input_data)),
                           warning = rep(FALSE, times=nrow(input_data)))
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Phenophase 1 
  #         starts in the year 
  #         at the first positive slope
  #         a warning if it is after 1st June
  
  test_cond_phase1 = input_data$t>0 & input_data$slope>0
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Phenophase 2 
  #       is the first maximum (i.e. slope is positive before and negative after)
  
  test_cond_phase2 = input_data$t>0 & 
    input_data$t<366 & 
    input_data$slopeprev>0 &
    input_data$slope<0      
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Phenophase 3 
  #         
  #         has a negative slope before
  #         a slope after that is less steep (or maybe positive) compared to before
  #         A warning if it is not within the year
  
  test_cond_phase3 = input_data$t>0 & 
    input_data$slopeprev<0 & 
    input_data$slopeprev< input_data$slope
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Find possible break points for each phenophase
  phase1_ind = NA
  phase2_ind = NA
  phase3_ind = NA
  
  # Phenophase 1 (SOS) ++++++++++++++++++++++++++++++++++++++++++
  if (any(test_cond_phase1)) {
    phase1_ind = which(test_cond_phase1)[1]  # Take first valid phase 1 (SOS) break point
    
    
    # Adjust any preceding breakpoints
    output_data$phase[1:phase1_ind] = phase1_ind-rev(c(1:phase1_ind))
  } 
  
  
  
  
  
  # Phenophase 2 (POS) ++++++++++++++++++++++++++++++++++++++++++
  # Pick phenophase 2 if it is after phenophase 1
  if (any(test_cond_phase2) & is.na(phase1_ind)) {
    # No phenophase 1 defined
    phase2_ind = which(test_cond_phase2)[1]
  } else if (!is.na(phase1_ind)) {
    if (any(test_cond_phase2[-c(1:phase1_ind)]) ) {
      # Phenophase 1 defined
      phase2_ind = which(test_cond_phase2[-c(1:phase1_ind)])[1] + phase1_ind
    }
  }
  
  # Adjust any breakpoints between phenophases 1 and 2 to be a value of 1.5
  if (!is.na(phase2_ind)) {
    output_data$phase[phase2_ind] = 2     # Set first valid phase 2 breakpoint (POS)
    
    if (!is.na(phase1_ind) & phase2_ind>(phase1_ind+1)) {
      # Adjust breakpoints between phase 1 and 2 to be 1.5
      output_data$phase[(phase1_ind+1):(phase2_ind-1)] = 1.5
    }
    if (is.na(phase1_ind)) {
      # If no phase 1 set all phases before 2 to NA
      output_data$phase[1:(phase2_ind-1)] = NA
    }
  }
  
  
  
  
  
  
  # Phenophase 3 (EOS) ++++++++++++++++++++++++++++++++++++++++++
  # Pick phenophase 3 if it is after phenophase 1 and phenophase 2
  if (any(test_cond_phase3) & is.na(phase2_ind) & is.na(phase1_ind)) {
    # No phenophases 1 & 2 defined
    phase3_ind = which(test_cond_phase3)[1]
  } 
  if (!is.na(phase2_ind)) {
    # Phenophase 2 defined
    if (any(test_cond_phase3[-c(1:phase2_ind)])) {
      phase3_ind = which(test_cond_phase3[-c(1:phase2_ind)])[1] + phase2_ind
    }
  } else if (!is.na(phase1_ind)) {
    # No phenophase 2 but phenophase 1 defined
    if (any(test_cond_phase3[-c(1:phase1_ind)])) {
      phase3_ind = which(test_cond_phase3[-c(1:phase1_ind)])[1] + phase1_ind
    }
  }
  
  if (!is.na(phase3_ind)) {
    output_data$phase[phase3_ind] = 3     # Set first valid phase 3 breakpoint (EOS)
    
    if (!is.na(phase2_ind) & phase3_ind>(phase2_ind+1)) {
      # Adjust any breakpoints between phenophases 2 and 3  to be a value of 2.5
      output_data$phase[(phase2_ind+1):(phase3_ind-1)] = 2.5
    }
    if (is.na(phase2_ind) & !is.na(phase1_ind) & phase3_ind - phase1_ind>1) {
      # If no phase 2,  set all phases between phase 1 and phase 3 to NA
      output_data$phase[(phase1_ind+1):(phase3_ind-1)] = NA
    }
    if (is.na(phase2_ind) & is.na(phase1_ind)) {
      # If no phase 1 or phase 2,  set all phases before 3 to NA
      output_data$phase[1:(phase3_ind-1)] = NA
    }
    if (phase3_ind<nrow(output_data)) {
      # Set all breaks after phase 3 to be phase 4
      output_data$phase[-c(1:phase3_ind)] = 4
    }
  }
  
  
  
  # Flag warnings +++++++++++++++++++++++  
  
  # If phase one is later than day 152 (1st June on non leap years) then raise a warning
  if (is.finite(phase1_ind)) {
    if (input_data$t[phase1_ind]>=152) {
    output_data$warning[phase1_ind] = TRUE
    }
  }
  
  # If phase three is later than day 365 (31st December on non leap years) then raise a warning
  if (is.finite(phase3_ind)) {
    if (input_data$t[phase3_ind]>365) {
    output_data$warning[phase3_ind] = TRUE
    }
  }
  
  
  return(output_data)
}




# Parallel execution function ----------
segment_within_pixel = function(input_data, pixel_info, output_labels) {
  
  # Perform segmentation on the input data
  segments = segmentEVI(input_data, 
                        nBreaks = nSegBreaks,
                        knots=knots, 
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
                            phase = c(1:nSegBreaks), 
                            t=CI[,1],
                            lowerCI=CI[,2],
                            upperCI=CI[,3],
                            slope=slope_est[-1],
                            slopeprev=slope_est[1:nSegBreaks])

    } else {
      tmp = data.frame(pixelID = pixel_info$pixelID,
                       x_ITM = pixel_info$x_ITM,
                       y_ITM = pixel_info$y_ITM,
                       x_MODIS = pixel_info$x_MODIS,
                       y_MODIS = pixel_info$y_MODIS,
                       year = pixel_info$year,
                       square = pixel_info$square,
                       model=modelStr[model],
                       phase = c(1:nSegBreaks), 
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





# Import data --------

for (s in square) {
  print(paste0('Analysing square ',s,' at ',Sys.time()))
  
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
  modelStr = c('_raw','_smooth')
  
  output_labels = paste0(rep(paste0(rep(label,times=nSegBreaks),
                                         rep(c(1:nSegBreaks),each=length(label))), times=2),
                              rep(modelStr,each=nSegBreaks*length(label)))
  
  
  # Create data frame for the output
  pixel_data$year = year
  pixel_data$square = s
  
  
  # Segment the EVI data for each pixel in a square  
  segment_output <- foreach (p = icount(nPixel), .packages=c('segmented','mgcv'), .inorder=FALSE, .combine='rbind') %dopar% {
   segment_within_pixel(d_final[d_final$pixelID==pixel_list[p],], pixel_data[pixel_data$pixelID==pixel_list[p],],label)
  }

  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Finished estimating segmentation for all pixels
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Reformat dataframe to put
  # all phenophase dates in one column
  # library(tidyr)
  # output_tmp = pivot_longer(output, 
  #                           cols = -c(1:6),
  #                           names_to = c('t','phase','smooth'),
  #                           names_sep = '_',
  #                           values_to = 'doy')
  # 
  # # This pulls the 95% confidence interval out into their own columns
  # output_long = pivot_wider(output_tmp, 
  #                           id_cols=-doy, 
  #                           names_from = t,
  #                           values_from = doy)
  # 
  
  
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
  
  for (m in c(1:2)) {  # Loop around the two models (raw and smoothed)
    for (p in 1:nPixel) {
      # Estimate phenophases for smoothed data
      ind = which(segment_output$pixelID==pixel_list[p] & segment_output$model==modelStr[m])
      if (all(is.finite(segment_output$t[ind]))) {
        
        # Call function to assign phenophases
        phenophase = assignPhenophase(segment_output[ind,])
        
        # Put the phenophases into the final data frame
        segment_output$phase[ind] = phenophase$phase
        segment_output$warning[ind] = phenophase$warning
      }
    }
  }
  
  
  # Save data -------------------
  filename = paste0('phenology_square_',s,'_',year,'.RData')
  save(knots,year,s,nSegBreaks,output_long,d_final,
       file=file.path(outputDir,filename))
} # Finish looping around all the squares


# Disconnect the compute nodes
stopCluster(cl) 
  
  
  