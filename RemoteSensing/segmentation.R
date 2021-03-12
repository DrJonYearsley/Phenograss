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


#library(quantreg)
library(mgcv)
library(segmented)
library(ggplot2)

# Import data --------
square = c(1:9,13:21)
#square = c(13:21)


year = 2012
knots = -1
min_obs = 15      # Minimum number of rows for trying to segment data
starting_breaks = c(50,100, 200, 300)

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
  
  
  # file_prev = paste0('modis_pasture_A',year-1,'_square',s,'.RData')
  # file = paste0('modis_pasture_A',year,'_square',s,'.RData')
  # file_next = paste0('modis_pasture_A',year+1,'_square',s,'.RData')
  
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
  
  
  # *************************************
  # Function to perform segmentation ----
  segmentEVI = function(d_sub, 
                        starting_breaks = c(100, 200,300), 
                        useHighQuality=FALSE, 
                        knots=-1,
                        use_raw = FALSE) {
    # This functions fits a segmented linear model to raw data and smoothed data
    # the smoothed data gives fewer NA's and more consistent results. 
    
    require(segmented)
    
    
    
    if (use_raw) {  
      # +++++++++++++++++++++++++++++++++++++++++++++++++
      # Perform segmentation on raw evi data
      m_raw = lm(evi~doy, data=d_sub)
      
      m_seg_raw = tryCatch(segmented(m_raw, 
                                     seg.Z = ~doy,
                                     psi = list(doy=starting_breaks),
                                     control=seg.control(display=FALSE)),
                           error = function(e) {
                             NA
                           },
                           warning = function(e) {
                             NA
                           })
    } else {
      m_seg_raw = NA
    }       
    
    # +++++++++++++++++++++++++++++++++++++++++++++++++
    # Perform segmentation on data smoothed using a GAM
    
    m_gam = tryCatch(gam(evi~s(doy, bs='cr',k=knots), 
                         data=d_sub),
                     error = function(e) {
                       NA
                     },
                     warning = function(e) {
                       NA
                     })
    
    
    if (!is.na(m_gam)) {
      # Remove EVI data points that are more than 6 SE below from the prediction
      tmp = predict(m_gam, se.fit = TRUE)
      evi_deviation = (d_sub$evi-tmp$fit)/tmp$se.fit
      
      # Smooth data after removing data more than 6 se below prediction  
      m_gam = gam(evi~s(doy, bs='cr',k=knots), data=d_sub[evi_deviation>-6,])
      
      # Add smoothed predictions to the data frame
      d_sub$evi_smooth = NA
      d_sub$evi_smooth[evi_deviation>-6] = predict(m_gam)
      
      # If a lone prediction is 
      
      # Segmenting the smoothed predictions
      m_smooth = lm(evi_smooth~doy, data=d_sub[evi_deviation>-6,])
      
      m_seg_smooth = tryCatch(segmented(m_smooth, 
                                        seg.Z = ~doy,
                                        psi = list(doy=starting_breaks),
                                        control=seg.control(display=FALSE)),
                              error = function(e) {
                                NA
                              },
                              warning = function(e) {
                                NA
                              })
    } else {
      m_seg = NA
      m_seg_smooth = NA
    }
    
    
    return(list(m_seg_raw,m_seg_smooth)) 
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
    d_add = data.frame(t = NA, lowerCI=NA, upperCI=NA, slope=NA, slopeprev = NA)
    names(d_add) = paste0(label,i,modelStr[1])    
    output = cbind(output, d_add)
    
    # Add in data holder for smoothed data segmentation
    d_add = data.frame(t = NA, lowerCI=NA, upperCI=NA, slope=NA, slopeprev = NA)
    names(d_add) = paste0(label,i,modelStr[2])    
    output = cbind(output, d_add)
  }
  
  
  
  
  # Loop around all the pixels
  for (p in 1:nPixel) {
    print(paste0('Estimating pixel ',pixel_list[p]))
    
    # Subset data from this pixel
    d_sub = d_final[d_final$pixelID==pixel_list[p],]
    
    segments = segmentEVI(d_sub, 
                          starting_breaks = starting_breaks,
                          useHighQuality=FALSE,
                          knots=knots)
    
    for (model in c(1,2)) {
      # Record estimates if they have been made 
      if (!is.na(segments[[model]][1])) {
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
  
  # Loop around all pixels and move the phenophases so that phase 1 starts in the year
  # at the first positive slope
  output_smoothed$phase = as.integer(output_smoothed$phase)
  output_smoothed$warning = FALSE
  for (p in 1:nPixel) {
    ind = which(output_smoothed$pixelID==pixel_list[p])
    if (all(is.finite(output_smoothed$t[ind]))) {
      
      # Phenophase 1 needs a growing segment in the focal year that 
      # shows increased growth from the last phase
      test_cond = output_smoothed$t[ind]>0 & output_smoothed$slope[ind]>0
      
      if (any(test_cond)) {
        phase1_ind = which(test_cond)[1]
        
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

# # ************************************************
# # Define phenophase 1 to be the date in the year when the 
# # slope following the date is positive.
# 
# # If phenophase 1 is <0 move all phenophases by 1
# output_long$phase = as.integer(output_long$phase)
# ind = output_long$phase==1 & is.finite(output_long$t) & output_long$t<0 & output_long$smooth=='smooth'
# ind_shift = output_long$pixelID%in%output_long$pixelID[ind] & output_long$smooth=='smooth'
# output_long$phase[ind_shift] = output_long$phase[ind_shift] - 1
# 
# ind = output_long$phase==1 & is.finite(output_long$t) & output_long$t<0 & output_long$smooth=='raw'
# ind_shift = output_long$pixelID%in%output_long$pixelID[ind] & output_long$smooth=='raw'
# output_long$phase[ind_shift] = output_long$phase[ind_shift] - 1
# 
# output_long$phase = as.factor(output_long$phase)





# 
# 
# # Plot estimate of phenology dates ---
# ggplot(data=output_smoothed, 
#        aes(x=t,
#            fill=factor(phase))) +
#   geom_histogram(position='dodge',
#                  bins = 50) +
#   scale_fill_brewer('Phase',
#                     palette='Dark2') +
#   labs(x='Day of Year',
#        y = 'Number of Pixels',
#        title=paste('Square',square,'Year=',year)) +
#   theme_bw() + 
#   theme(axis.text = element_text(size=18),
#         axis.title = element_text(size=20),
#         legend.text = element_text(size=14),
#         legend.title = element_text(size=14))
# 
# ggsave(width=11, height=6,
#        filename = paste0('phenophases_square',square,'_',year,'.png'))
# 
# 
# ggplot(data=subset(output_smoothed,phase==0), 
#        aes(x=x_MODIS,
#            y=y_MODIS,
#            fill = t)) +
#   geom_tile() + 
#   scale_fill_viridis_c('Day of Year',option='magma') +
#   labs(title=paste('Phase 0: Square',square,'Year=',year)) +
#   theme_bw()
# 
# ggplot(data=subset(output_smoothed,phase==1), 
#        aes(x=x_MODIS,
#            y=y_MODIS,
#            fill = t)) +
#   geom_tile() + 
#   scale_fill_viridis_c('Day of\nYear',option='magma') +
#   labs(title=paste('Phase 1: Square',square,'Year=',year)) +
#   theme_bw()
# 
# 
# # *****************************************
# # Visualise a single time series  --------
# 
# p = sample(pixel_list, size=1)
# 
# # Select unusual pixels
# ind = which(output_smoothed$phase==1 & output_smoothed$t>100 & output_smoothed$t<130)
# p = sample(output_smoothed$pixelID[ind],size=1)
# 
# d_sub = d_final[d_final$pixelID==p,]
# 
# segments = segmentEVI(d_sub, 
#                       starting_breaks = starting_breaks,
#                       useHighQuality=FALSE,
#                       knots = knots)
# 
# # Make predictions across the year
# d_pred = data.frame(doy=c(-100:465), 
#                     evi=NA, 
#                     evi_smooth=NA, 
#                     evi_smooth_lwr=NA, 
#                     evi_smooth_upr=NA)
# 
# m_gam = gam(evi~s(doy, bs='tp',k=knots), data=d_sub)
# d_pred$evi = predict(m_gam, newdata=d_pred)
# 
# 
# 
# 
# 
# if (is.na(segments[[2]][1])) {
#   d_pred$evi_smooth = NA
#   breaks_smooth = NA
# } else {
#   tmp = predict(segments[[2]], newdata=d_pred, interval='confidence')
#   d_pred$evi_smooth = tmp[,1]
#   d_pred$evi_smooth_lwr = tmp[,2]
#   d_pred$evi_smooth_upr = tmp[,3]
#   breaks_smooth = confint(segments[[2]])
# }
# 
# breaks_df = as.data.frame(rbind(breaks_smooth,
#                                 breaks_smooth))
# breaks_df$phase = rep(c(1:length(starting_breaks)), times=2)
# breaks_df$y = rep(range(d_sub$evi), each=length(starting_breaks))
# names(breaks_df) = c('x','xmin','xmax','phase','y')
# 
# ggplot() +
#   geom_ribbon(data=d_pred, aes(x=doy,
#                                y=evi_smooth,
#                                ymin=evi_smooth_lwr, 
#                                ymax=evi_smooth_upr),
#               alpha=0.3, fill='blue') +
#   geom_line(data=d_pred, aes(x=doy, y=evi_smooth),colour='blue') +
#   geom_ribbon(data=breaks_df,
#               aes(xmin=xmin,
#                   xmax=xmax,
#                   y=y,
#                   group=factor(phase)), alpha=0.3) +
#   geom_vline(xintercept=breaks_smooth[,1],colour='blue') +
#   geom_point(data=d_sub,aes(x=doy, y=evi,fill=factor(QC)), shape=21, size=2) +
#   geom_path(data=d_pred, aes(x=doy, y=evi)) +
#   scale_colour_brewer(palette = 'Dark2') +
#   scale_fill_discrete('MODIS\nQuality\nFlag',
#                       type = c('black','white'), 
#                       label=c('V. Good','Good')) +
#   labs(x='Day of Year',
#        y='EVI',
#        title=paste('Pixel',p)) +
#   theme_bw() +
#   theme(axis.text = element_text(size=18),
#         axis.title = element_text(size=20),
#         legend.text = element_text(size=14),
#         legend.title = element_text(size=14))
# 
# ggsave(width=11, height=6,filename = paste0('pheno_pix',p,'_square',square,'_',year,'.png'))
