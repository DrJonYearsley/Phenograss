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
  
  
  return(list(spring=m_seg_spring, autumn=m_seg_autumn, peak_doy=peak_doy, filtered=filter_ind, d_sub=d_sub, m_gam2=m_gam2, pred_df=pred_df)) 
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
    
    # If a lone prediction is 
    
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


# +++++++++++++ End of function definitions ++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Import data --------
square = c(1:9,13:21)
square = c(1)


year = 2012
knots = -1
min_obs = 15      # Minimum number of rows for trying to segment data
starting_breaks = c(50, 100, 200, 300)

for (s in square) {
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
  
  for (i in 1:length(starting_breaks)) {
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
    print(paste0('Estimating pixel ',pixel_list[p]))
    
    # Subset data from this pixel
    p = sample.int(nPixel, size=1)
    d_pixel = d_final[d_final$pixelID==pixel_list[p],]
    
    segments = segmentEVI2(d_pixel, 
                           nBreaks = 3,
                           knots=knots, 
                           use_raw = TRUE,
                           sd_filter_threshold = 6)
    
    #d_pixel$evi_smooth = segments$d_sub$evi_smooth
    
    plot(d_pixel$doy, d_pixel$evi)
    plot(segments$spring,  add=T, col="blue")
    plot(segments$autumn,  add=T, col="red")
    points(d_pixel$doy[segments$filtered], d_pixel$evi[segments$filtered], pch=20, col="red")
    points(segments$pred_df$doy,segments$pred_df$fit, pch=20)
    
    
    
    for (model in c(1,2)) {
      # Record estimates if they have been made 
      if (!is.null(segments[[model]])) {
        # Calculate estimated segment break dates and 95% CI
        CI = confint(segments[[model]])
        
        # calculate slopes
        slope_est = slope(segments[[model]])$doy[,1]
        
        # Locate pixel in data frame
        ind = which(output$pixelID==pixel_list[p])
        
        # Record estimates (break point and following slope)
        for (i in 1:length(starting_breaks)) {
          varNames = paste0(label,i,modelStr[model])
          output[ind,varNames] = c(CI[i,], slope_est[i+1], slope_est[i]) 
        }
      }
    }
  }
  
  
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
  
  
  # Just extract smoothed data
  output_smoothed = subset(output_long, smooth=='smooth')
  
  # Loop around all pixels and move the phenophases so that: 
  #    phase 1 
  #         starts in the year 
  #         at the first positive slope
  #    phase 2
  #         is the first maximum (i.e. slope is positive before and negative after)
  #    phase 3
  #         is within the year
  #         has a negative slope before
  #         a slope after that is less steep (or maybe positive) compared to before
  
  
  output_smoothed$phase = as.numeric(output_smoothed$phase)
  output_smoothed$warning = FALSE
  for (p in 1:nPixel) {
    ind = which(output_smoothed$pixelID==pixel_list[p])
    if (all(is.finite(output_smoothed$t[ind]))) {
      
      # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Phenophase 1 
      #         starts in the year 
      #         at the first positive slope
      
      test_cond_phase1 = output_smoothed$t[ind]>0 & output_smoothed$slope[ind]>0
      
      
      # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Phenophase 2 
      #       is the first maximum (i.e. slope is positive before and negative after)
      
      test_cond_phase2 = output_smoothed$t[ind]>0 & 
        output_smoothed$t[ind]<366 & 
        output_smoothed$slopeprev[ind]>0 &
        output_smoothed$slope[ind]<0      
      
      # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Phenophase 3 
      #         is within the year
      #         has a negative slope before
      #         a slope after that is less steep (or maybe positive) compared to before
      
      test_cond_phase3 = output_smoothed$t[ind]>0 & 
        output_smoothed$t[ind]<366 & 
        output_smoothed$slopeprev[ind]<0 & 
        output_smoothed$slopeprev[ind] < output_smoothed$slope[ind]
      
      
      # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Find possible break points for each phenophase
      phase1_ind = NA
      phase2_ind = NA
      phase3_ind = NA
      if (any(test_cond_phase1)) {
        phase1_ind = which(test_cond_phase1)[1]  # Take first valid phase 1 (SOS) break point
        
        
        # Adjust any preceding breakpoints
        output_smoothed$phase[ind[1:phase1_ind]] = phase1_ind-rev(c(1:phase1_ind)) + 1
      } 
      
      # Pick phenophase 2 if it is after phenophase 1
      if (any(test_cond_phase2) & is.na(phase1_ind)) {
        # No phenophase 1 defined
        phase2_ind = which(test_cond_phase2)[1]
      } else if (any(test_cond_phase2[-c(1:phase1_ind)]) ) {
        # Phenophase 1 defined
        phase2_ind = which(test_cond_phase2[-c(1:phase1_ind)])[1]
      }
      
      # Adjust any breakpoints between phenophases 1 and 2 to be a value of 1.5
      if (!is.na(phase2_ind)) {
        if (!is.na(phase1_ind) & phase2_ind>(phase1_ind+1)) {
          output_smoothed$phase[ind[(phase1_ind+1):(phase2_ind-1)]] = 1.5
        }
        if (is.na(phase1_ind)) {
          output_smoothed$phase[ind[1:(phase2_ind-1)]] = NA
        }
      }
      
      
      # Pick phenophase 3 if it is after phenophase 1 and phenophase 2
      if (any(test_cond_phase3) & is.na(phase2_ind) & is.na(phase1_ind)) {
        # No phenophases 1 & 2 defined
        phase3_ind = which(test_cond_phase3)[1]
      } 
      if (!is.na(phase2_ind)) {
        # Phenophase 2 defined
        if (any(test_cond_phase3[-c(1:phase2_ind)])) {
          phase3_ind = which(test_cond_phase3[-c(1:phase2_ind)])[1]
        }
      } else if (!is.na(phase1_ind)) {
        # No phenophase 2 but phenophase 1 defined
        if (any(test_cond_phase3[-c(1:phase1_ind)])) {
          phase3_ind = which(test_cond_phase3[-c(1:phase1_ind)])[1]
        }
      }
      # Adjust any breakpoints between phenophases 2 and 3  to be a value of 2.5
      
      
      # Correct phase
      output_smoothed$phase[ind] = output_smoothed$phase[ind] - output_smoothed$phase[ind[phase1_ind]] + 1
      
      # If previous phase has higher slope then remove phases
      if (output_smoothed$slopeprev[ind[phase1_ind]]>output_smoothed$slope[ind[phase1_ind]]) {
        output_smoothed$warning[ind] = TRUE
        output_smoothed$phase[ind] = NA
      }
      
      # If phase one is later than day 152 (1st June on non leap years) then raise a warning
      if (output_smoothed$t[ind[phase1_ind]]>=152) {
        output_smoothed$warning[ind] = TRUE
      }
      
    } else {  # No breakpoint meets condition
      output_smoothed$phase[ind] = NA
      output_smoothed$warning[ind] = TRUE
    }
    
    
  }
}


# Save data -------------------
filename = paste0('phenology_square_',s,'_',year,'.RData')
save(knots,year,s,starting_breaks,output_smoothed,d_final,
     file=file.path(outputDir,filename))


}
