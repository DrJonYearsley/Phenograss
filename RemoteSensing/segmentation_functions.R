# Define two functions for segmenting an EVI time series
#
# The functions are:
#       segmentEVI:        Take an EVI time series, clean it and 
#                          segment it into linear components
#       assignPhenophase:  Take the output of segmentEVI and identify three 
#                          phenophases that obey qualitative conditions for each phase. 
#                          The three phases are: Start of season (SOS), 
#                         peak of season (POS) and end of season (EOS)
#
#  Jon yearsley (Jon.Yearsley@ucd.ie)
#  Dec 2021
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





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
  
  
  # Fit initial GAM
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
    # Remove EVI data points that are more than sd_filter_threshold SE below from the prediction from m_gam
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
    m_smooth = lm(evi_smooth~doy, data=d_sub[filter_ind,])  # Create base lm model
    
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
  
  # Outputs are segmented model using raw data, segmented model using smoothed data, 
  # indices for data kept in analysis, the original data, the gam used for smoothing
  return(list(m_seg_raw, m_seg_smooth, filtered=filter_ind, d_sub=d_sub, m_gam2=m_gam2)) 
}









# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    output_data$phase[1:phase1_ind] = 2-rev(c(1:phase1_ind))
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



