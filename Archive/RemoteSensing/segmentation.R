# Segmentation modelling of phenology
#
# Apply segmentation method developed by Gabin
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


# Set basic parameters
square = c(1:9,13:21)
square = c(10:12)
square = c(1:21)

year = 2014
knots = -1
min_obs = 15      # Minimum number of rows for trying to segment data
#starting_breaks = c(50, 100, 200, 300)  # Initial guess for break points in doy
nSegBreaks = 4    # Number of breakpoints to use for segmentation
print_pixel=FALSE


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++ Start of function definitions ++++++++++++++++++++++++++

# *************************************
# Function to perform segmentation on spring, autumn and summer separately ----
segmentEVI2 = function(d_sub, 
                       nBreaks = 2, 
                       knots=-1,
                       use_raw = FALSE,
                       sd_filter_threshold=4) {
  # This functions fits a segmented linear model to raw data and smoothed data
  # the smoothed data gives fewer NA's and more consistent results. 
  
  # This function differs from segmentEVI by splitting the year into 
  # two parts: sprint and autumn (e.g. spring is  doy<180,  autumn is doy>150) 
  # and performing a segmentation on each part separately
  
  require(segmented)
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++
  # Smooth time series in order to identify outlier data
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
  }
  
  # Fit a gam and find the time of the maximum  evi
  m_gam2 = tryCatch(gam(evi~s(doy, bs='cr',k=knots), 
                        data=d_sub[filter_ind, ],
                        gamma=1),  # Gamma controls over/under fitting
                    error = function(e) {
                      NULL
                    },
                    warning = function(e) {
                      NULL
                    }) 
  if (!is.null(m_gam2)) {
    pred_df = data.frame(doy=c(1:365))
    pred_df$fit = predict(m_gam2, newdata=pred_df)
    peak_doy = pred_df$doy[pred_df$fit==max(pred_df$fit)]
  } else {
    peak_doy = NA
  }
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++
  # Perform segmentation on raw evi data 
  if (is.na(peak_doy)) {
    doy_divider = 180
  } else {
    doy_divider = peak_doy
  }
  
  
  # Estimate start of season
  m_raw_spring = lm(evi~doy, data=d_sub[filter_ind & d_sub$doy<doy_divider,])   # DOY=150 is 30th May
  spring_break=nBreaks
  m_seg_spring = NULL
  while(spring_break>0 & is.null(m_seg_spring)) {
    m_seg_spring = tryCatch(segmented(m_raw_spring, 
                                      seg.Z = ~doy,
                                      npsi=spring_break,
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
    spring_break = spring_break-1
  }
  
  
  # Estimate end of season
  m_raw_autumn = lm(evi~doy, data=d_sub[filter_ind & d_sub$doy>doy_divider,])   # DOY=150 is 30th May
  autumn_break=nBreaks
  m_seg_autumn = NULL
  while(autumn_break>0 & is.null(m_seg_autumn)) {
    m_seg_autumn = tryCatch(segmented(m_raw_autumn, 
                                      seg.Z = ~doy,
                                      npsi=autumn_break,
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
    autumn_break = autumn_break - 1
  }
  
  
  return(list(spring=m_seg_spring, autumn=m_seg_autumn, 
              peak_doy=peak_doy, filtered=filter_ind, d_sub=d_sub, 
              m_gam2=m_gam2, pred_df=pred_df)) 
}



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
  
  
  
  # Labels for the final data frame variables
  # t_ is the time of the segmentation
  # slope_ is the regression slope after t_
  # slopeprev_ is the regression slope before t_
  label = c('t_', 'lowerCI_', 'upperCI_','slope_','slopeprev_')
  modelStr = c('_raw','_smooth')
  
  # Create data frame for the output
  output = pixel_data
  output$year = year
  output$square = s
  
  for (i in 1:nSegBreaks) {
    # Add in data holder for raw data segmentation
    d_add = data.frame(t = as.numeric(NA), 
                       lowerCI=as.numeric(NA), 
                       upperCI=as.numeric(NA), 
                       slope=as.numeric(NA), 
                       slopeprev = as.numeric(NA))
    names(d_add) = paste0(label,i,modelStr[1])    
    output = cbind(output, d_add)
    
    # Add in data holder for smoothed data segmentation
    names(d_add) = paste0(label,i,modelStr[2])    
    output = cbind(output, d_add)
  }
  
  
  
  
  # Loop around all the pixels
  for (p in 1:nPixel) {
    if (print_pixel) {
      print(paste0('Estimating pixel ',pixel_list[p]))
    }

    
    segments = segmentEVI(d_final[d_final$pixelID==pixel_list[p],], 
                          nBreaks = nSegBreaks,
                          knots=knots, 
                          use_raw = TRUE,
                          sd_filter_threshold = 6)
    
    # # Save the smoothed data
    # d_pixel = d_final[d_final$pixelID==pixel_list[p],]
    # d_pixel$evi_smooth = segments$d_sub$evi_smooth
    # 
    # # Plot the points used in the analysis and filtered points
    # plot(d_pixel$doy, d_pixel$evi)
    # points(d_pixel$doy[segments$filtered], d_pixel$evi[segments$filtered], pch=20, col="red")

    
    
    for (model in c(1,2)) {
      # Model 1 corresponds to using the raw data, model 2 corresponds to using th smoothed data
      
      # Record estimates if they have been made 
      if (!is.null(segments[[model]])) {
        # Calculate estimated segment break dates and 95% CI
        CI = confint(segments[[model]])
        
        # calculate slopes
        slope_est = slope(segments[[model]])$doy[,1]
        
        # Locate pixel in data frame
        ind = which(output$pixelID==pixel_list[p])
        
        # Record estimates (break point and following slope)
        for (i in 1:nSegBreaks) {
          varNames = paste0(label,i,modelStr[model])
          output[ind,varNames] = c(CI[i,], slope_est[i+1], slope_est[i]) 
        }
      }
    }
  }
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Finished estimating segmentation for all pixels
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  
  
  
  
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Reformat dataframe to put
  # all phenophase dates in one column
  library(tidyr)
  output_tmp = pivot_longer(output, 
                            cols = -c(1:6),
                            names_to = c('t','phase','smooth'),
                            names_sep = '_',
                            values_to = 'doy')
  
  # This pulls the 95% confidence interval out into their own columns
  output_long = pivot_wider(output_tmp, 
                            id_cols=-doy, 
                            names_from = t,
                            values_from = doy)
  
  
  
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
  
  
  
  output_long$phase = NA            # Remove draft phase information
  output_long$warning = NA
  
  for (p in 1:nPixel) {
    # Estimate phenophases for smoothed data
    ind = which(output_long$pixelID==pixel_list[p] & output_long$smooth=="smooth")
    if (all(is.finite(output_long$t[ind]))) {
      
      phenophase = assignPhenophase(output_long[ind,])
      
      # Put the phenophases into the final data frame
      output_long$phase[ind] = phenophase$phase
      output_long$warning[ind] = phenophase$warning
    }
    
    
    # Estimate phenophases for raw data
    ind = which(output_long$pixelID==pixel_list[p] & output_long$smooth=="raw")
    if (all(is.finite(output_long$t[ind]))) {
      
      phenophase = assignPhenophase(output_long[ind,])
      
      # Put the phenophases into the final data frame
      output_long$phase[ind] = phenophase$phase
      output_long$warning[ind] = phenophase$warning
    }
  }
  
  
  # Save data -------------------
  filename = paste0('phenology_square_',s,'_',year,'.RData')
  save(knots,year,s,nSegBreaks,output_long,d_final,
       file=file.path(outputDir,filename))
} # Finish looping around all the squares


 
  
  
  