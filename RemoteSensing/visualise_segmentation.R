# Visualise phenology dates from segmentation
#
# Jon Yearsley
# Nov 2020
# 
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(segmented)
library(ggplot2)
library(mgcv)

setwd('~/Research/Phenograss/Data/PhenologyOutput/')

dataDir = '.'

input_file_prefix = 'phenology'

# Import data --------
square = 20
year = 2019

# starting_breaks = c(50, 100, 200, 300)


# Filename segmented data
filename = paste0(input_file_prefix,'_square_',square,'_',year,'.RData')
load(file.path(dataDir,filename))

d_final$evi = d_final$evi / 0.0001^2

pixel_list = unique(d_final$pixelID)
nPixel = length(pixel_list)


# Plot estimate of phenology dates ---
ggplot(data=subset(output_smoothed, warning==FALSE), 
       aes(x=t,
           fill=factor(phase))) +
  geom_histogram(position='dodge',
                 bins = 50) +
  scale_fill_brewer('Phase',
                    palette='Dark2') +
  labs(x='Day of Year',
       y = 'Number of Pixels',
       title=paste('Square',square,'Year=',year)) +
  theme_bw() + 
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))

ggsave(width=11, height=6,
       filename = paste0('phenophases_square',square,'_',year,'.png'))



# Plot map of phenophase dates ----
ggplot(data=subset(output_smoothed,phase==0), 
       aes(x=x_MODIS,
           y=y_MODIS,
           fill = t)) +
  geom_tile() + 
  scale_fill_viridis_c('Day of Year',option='magma') +
  labs(title=paste('Phase 0: Square',square,'Year=',year)) +
  theme_bw()

ggplot(data=subset(output_smoothed,phase==1 & warning==FALSE), 
       aes(x=x_MODIS,
           y=y_MODIS,
           fill = t)) +
  geom_tile() + 
  coord_equal() +
  scale_fill_viridis_c('Day of\nYear',option='magma', limits=c(0,150)) +
  labs(x = 'X coord (MODIS CRS)',
       y = 'Y coord (MODIS CRS)',
    title=paste('Phase 1: Square',square,'Year=',year)) +
  theme_bw()

ggsave(width=11, height=6,
       filename = paste0('phenophase1_map_square',square,'_',year,'.png'))

ggplot(data=subset(output_smoothed,phase==2 & warning==FALSE), 
       aes(x=x_MODIS,
           y=y_MODIS,
           fill = t)) +
  geom_tile() + 
  scale_fill_viridis_c('Day of\nYear',option='magma') +
  labs(title=paste('Phase 2: Square',square,'Year=',year)) +
  theme_bw()

# *****************************************
# Visualise a single time series  --------

p = sample(pixel_list, size=1)

# # Select unusual pixels
# ind = which(output_smoothed$phase==1 & output_smoothed$t>200 & output_smoothed$t<1300)
# p = sample(output_smoothed$pixelID[ind],size=1)

d_sub = d_final[d_final$pixelID==p,]


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
  
  m_gam = gam(evi~s(doy, bs='cr',k=knots), data=d_sub)
  
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
  
  
  return(list(m_seg_raw,m_seg_smooth, d_sub[evi_deviation>-6,], m_gam)) 
}




segments = segmentEVI(d_sub, 
                      starting_breaks = starting_breaks,
                      useHighQuality=FALSE,
                      knots = knots)

# Make predictions across the year
d_pred = data.frame(doy=c(-100:465), 
                    evi=NA, 
                    evi_smooth=NA, 
                    evi_smooth_lwr=NA, 
                    evi_smooth_upr=NA)

m_gam = gam(evi~s(doy, bs='cr',k=knots), data=d_sub)
d_pred$evi = predict(m_gam, newdata=d_pred)

d_pred$evi = predict(segments[[4]], newdata=d_pred)



if (is.na(segments[[2]][1])) {
  d_pred$evi_smooth = NA
  breaks_smooth = NA
} else {
  tmp = predict(segments[[2]], newdata=d_pred, interval='confidence')
  d_pred$evi_smooth = tmp[,1]
  d_pred$evi_smooth_lwr = tmp[,2]
  d_pred$evi_smooth_upr = tmp[,3]
  breaks_smooth = confint(segments[[2]])
}

breaks_df = as.data.frame(rbind(breaks_smooth,
                                breaks_smooth))
breaks_df$phase = rep(c(1:length(starting_breaks)), times=2)
breaks_df$y = rep(range(d_sub$evi), each=length(starting_breaks))
names(breaks_df) = c('x','xmin','xmax','phase','y')

ggplot() +
  geom_ribbon(data=d_pred, aes(x=doy,
                               y=evi_smooth,
                               ymin=evi_smooth_lwr, 
                               ymax=evi_smooth_upr),
              alpha=0.3, fill='blue') +
  geom_line(data=d_pred, aes(x=doy, y=evi_smooth),colour='blue') +
  geom_ribbon(data=breaks_df,
              aes(xmin=xmin,
                  xmax=xmax,
                  y=y,
                  group=factor(phase)), alpha=0.3) +
  geom_vline(xintercept=breaks_smooth[,1],colour='blue') +
  geom_point(data=segments[[3]],aes(x=doy, y=evi,fill=factor(QC)), shape=21, size=2) +
  geom_path(data=d_pred, aes(x=doy, y=evi)) +
  scale_colour_brewer(palette = 'Dark2') +
  scale_fill_discrete('MODIS\nQuality\nFlag',
                      type = c('black','white'), 
                      label=c('V. Good','Good')) +
  labs(x='Day of Year',
       y='EVI',
       title=paste('Pixel',p)) +
  theme_bw() +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))

ggsave(width=11, height=6,filename = paste0('pheno_pix',p,'_square',square,'_',year,'.png'))

