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

#setwd('~/Research/Phenograss/Data/PhenologyOutput_test/')
setwd('/media/jon/MODIS_data/PhenologyOutput_test/')
# Load functions to perform segmentation
source('~/git_repos/Phenograss/RemoteSensing/segmentation_functions.R')

dataDir = '.'

input_file_prefix = 'phenology'

# Import data --------
square = 14
year = 2012

# starting_breaks = c(50, 100, 200, 300)


# Filename segmented data
filename = paste0(input_file_prefix,'_square_',square,'_',year,'.RData')
#filename = paste0(input_file_prefix,'_square_',square,'.RData')
load(file.path(dataDir,filename))

#d_final$evi = d_final$evi / 0.0001^2

pixel_list = unique(d_final$pixelID)
nPixel = length(pixel_list)


# Plot estimate of phenology dates ---
ggplot(data=subset(segment_output, model=="smooth" & warning==FALSE & phase%in%c(1,2,3) & !wideCI), 
       aes(x=t,
           fill=factor(phase))) +
  geom_histogram(position='identity',
                 binwidth = 10,
                 alpha = 0.8,
                 colour='darkgrey') +
  scale_fill_brewer('Phenohase',
                    palette='Dark2',
                    labels=c('SOS','POS','EOS')) +
  # Label 1st April, 1st July and 1st Oct
  scale_x_continuous(breaks=c(1,91,182,274,365),
    labels=c(1,91,182,274,365),
                     limits = c(0,370)) +
  labs(x='Day of Year',
       y = 'Number of Pixels') +
  theme_bw() + 
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.position = 'top',
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(width=11, height=6,
       filename = paste0('phenophases_square',square,'_',year,'.png'))



# # Plot map of phenophase dates ----
# ggplot(data=subset(segment_output,phase==1 & model=="smooth"), 
#        aes(x=x_MODIS,
#            y=y_MODIS,
#            fill = t)) +
#   geom_tile() + 
#   scale_fill_viridis_c('Day of Year',option='magma',limits=c(0,100)) +
#   labs(title=paste('Phase 1: Square',square,'Year=',year)) +
#   theme_bw()
# 



# tmp = subset(output_smoothed,phase==1 & warning==FALSE)
# range(tmp$t, na.rm=TRUE)
# nColour = 9
# vals = round(seq(10,
#                  80,
#                  length.out=nColour-1), digits=2)
# col_lab = paste('< ',vals[1])
# tmp$sos = col_lab[1]
# for (i in 2:(nColour-1)) {
#   ind = tmp$t>=vals[i-1] & tmp$t<vals[i]
#   col_lab = c(col_lab, paste0(vals[i-1],' - ',vals[i]))
#   tmp$sos[ind] = col_lab[i]
# }
# ind = tmp$t>=vals[nColour-1]
# col_lab = c(col_lab, paste0('> ',vals[i]))
# tmp$sos[ind] = col_lab[i+1]
# tmp$sos[is.na(tmp$t)] = NA
# tmp$sos = factor(tmp$sos, ordered=T, levels=rev(col_lab))
# 
# ggplot(data=tmp,
#        aes(x=x_MODIS,
#            y=y_MODIS,
#            fill=sos)) +
#   geom_tile(colour='darkblue') +
#   coord_equal() +
#   scale_fill_brewer('Start of Season\n(day of year)',
#                     palette='Greens', 
#                     direction=1,
#                     na.value='darkgray') +
#   labs(x='X Coord (MODIS CRS)',
#        y='Y Coord (MODIS CRS)') +
#   theme_bw()
# 
# ggsave(width=11, height=6,
#        filename = paste0('phenophase1_map_square',square,'_',year,'.png'))





# ggplot(data=subset(output_smoothed,phase==2 & warning==FALSE), 
#        aes(x=x_MODIS,
#            y=y_MODIS,
#            fill = t)) +
#   geom_tile() + 
#   scale_fill_viridis_c('Day of\nYear',option='magma') +
#   labs(title=paste('Phase 2: Square',square,'Year=',year)) +
#   theme_bw()

# *****************************************
# Visualise a single time series  --------

p = sample(pixel_list, size=1)

p="x40y7"
# which(p==pixel_list)

type = c("raw","smooth")   # Either "raw" or "smooth"


# # Select unusual pixels
# ind = which(output_smoothed$phase==1 & output_smoothed$t>200 & output_smoothed$t<1300)
# p = sample(output_smoothed$pixelID[ind],size=1)

d_sub = d_final[d_final$pixelID==p,]
seg_sub = segment_output[segment_output$pixelID==p,]
pl = vector('list',length=2)

for (t in type) {
  
  
  if (t=="raw") {
    segment_ind = 1
  } else {
    segment_ind = 2
  }
  
  # Look at smoothed output
  segments = segmentEVI(d_sub, 
                        nBreaks=5,
                        useHighQuality=FALSE,
                        use_raw = TRUE,
                        knots = knots,
                        sd_filter_threshold=6)
  
  # Make predictions across the year
  d_pred = data.frame(doy=c(-100:465), 
                      evi_smoothed=NA, 
                      evi_segmented=NA, 
                      evi_segmented_lwr=NA, 
                      evi_segmented_upr=NA)
  
  
  
  if (is.null(segments[[segment_ind]])) {
    # No fitted model
    d_pred$evi_smooth = NA
    breaks_smooth = array(NA, dim=c(nSegBreaks,3))
    
    pl[[which(t==type)]] = ggplot() +
      # Plot original data points
      geom_point(data=segments[[4]],aes(x=doy, y=evi,fill=factor(QC)), shape=21, size=2) +
      # Highlight data points that were removed
      geom_point(data=segments[[4]][!segments[[3]],],aes(x=doy, y=evi), shape=21, size=2, fill="red")
    
  } else {
    # Create the smoothed data
    d_pred$evi_smoothed = predict(segments[[5]], newdata=d_pred)
    
    # Create the segmented data
    tmp = predict(segments[[segment_ind]], newdata=d_pred, interval='confidence')
    d_pred$evi_segmented = tmp[,1]
    d_pred$evi_segmented_lwr = tmp[,2]
    d_pred$evi_segmented_upr = tmp[,3]
    
    if (t=="smooth") {
      breaks = seg_sub[seg_sub$model=="smooth" & seg_sub$phase%in%c(1,2,3),]
      breaks_df = rbind(breaks, breaks)
      breaks_df$y = rep(c(0,1), each=3)

      #       breaks_smooth = confint(segments[[segment_ind]])
      # 
      # breaks_df = as.data.frame(rbind(breaks_smooth,
      #                                 breaks_smooth))
      # breaks_df$phase = rep(c(1:nSegBreaks), times=2)
      # breaks_df$y = rep(range(d_sub$evi), each=nSegBreaks)
      # breaks_df$y = rep(c(0,1), each=nSegBreaks)
      # names(breaks_df) = c('x','xmin','xmax','phase','y')
    }
    pl[[which(t==type)]] = ggplot() +
      geom_ribbon(data=d_pred, aes(x=doy,
                                   y=evi_segmented,
                                   ymin=evi_segmented_lwr, 
                                   ymax=evi_segmented_upr),
                  alpha=1, fill='grey') +
      geom_line(data=d_pred, aes(x=doy, y=evi_segmented),colour='black') +
      

    
      # Plot original data points
      geom_point(data=segments[[4]],aes(x=doy, y=evi), shape=21, size=2,fill='darkgrey') +
      # Highlight data points that were removed
      geom_point(data=segments[[4]][!segments[[3]],],aes(x=doy, y=evi), shape=22, size=3, fill="red")
  }
  
  if (t=="smooth") {
    # Add in lines for the phenophases
    pl[[which(t==type)]] = pl[[which(t==type)]] +
      geom_ribbon(data=breaks_df,
                aes(xmin=lowerCI,
                    xmax=upperCI,
                    y=y,
                    group=factor(phase),
                    fill=factor(phase)), 
                alpha=0.3) + geom_segment(data=breaks_df, aes(x=t, y=0, xend=t, yend=1, colour=factor(phase)))
  }
  
  # Add styling to the plot
  pl[[which(t==type)]] = pl[[which(t==type)]] + 
    # scale_fill_discrete('MODIS\nQuality\nFlag',
    #                     type = c('black','white'), 
    #                     label=c('V. Good','Good')) +
    lims(y=c(0,1)) +
    labs(x='Day of Year',
         y='EVI') +
    scale_fill_brewer('Phenophase',
                      palette='Dark2',
                      labels=c('SOS','POS','EOS')) +
    scale_colour_brewer('Phenophase',
                      palette='Dark2',
                      labels=c('SOS','POS','EOS')) +
    scale_x_continuous(breaks=c(-100,1,91,182,274,365,450),
                       labels=c(-100,1,91,182,274,365,450),
                       limits = c(-110,470)) +
    theme_bw() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=24),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          legend.position = 'top',
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  # if (t=="smooth") {
  #   # Add smoothed line to plot if appropriate
  #   pl[[which(t==type)]] = pl[[which(t==type)]] + geom_path(data=d_pred, aes(x=doy, y=evi_smoothed))   # Add in smoothed result if relevant
  # }
  # 
}

# pl[[1]]  # Raw
pl[[2]]  # Smooth



ggsave(width=11, height=6,filename = paste0('pheno_pix',p,'_square',square,'_',year,'.png'))

