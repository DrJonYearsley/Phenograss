# Segmentation modelling of phenology
#
# Apply segementation method developed by Gabin
#
# Jon Yearsley Jon.Yearsley@ucd.ie
# Aug 2020
# ****************************************

rm(list=ls())

dataDir = '~/MEGAsync/Projects/GrasslandPhenology/Data/'

library(quantreg)
library(mgcv)
library(segmented)
library(ggplot2)

# Import data --------

square = 9
year = 2017

# Load the focal year and years either side
file_prev = paste0('modis_pasture_A',year-1,'_square',square,'.RData')
file = paste0('modis_pasture_A',year,'_square',square,'.RData')
file_next = paste0('modis_pasture_A',year+1,'_square',square,'.RData')

load(file.path(dataDir,file_prev))
d_prev = d_sq[d_sq$doy>265,]
d_prev$doy = 265 - d_prev$doy

load(file.path(dataDir,file_next))
d_next = d_sq[d_sq$doy<100,]
d_next$doy = 365 + d_next$doy

load(file.path(dataDir,file))

d_final = rbind(d_prev, d_sq, d_next)

# Add unique pixel ID
x_values = unique(d_final$x_MODIS)
y_values = unique(d_final$y_MODIS)

d_final$pixelID = paste0('x',match(d_final$x_MODIS,x_values),
                         'y',match(d_final$y_MODIS,y_values))



# Test routine on a pixel
pixel = sample(unique(d_final$pixelID),1)
d_sub = d_final[d_final$pixelID==pixel,]

ind = d_sub$QC==0

# Try smoothing the data using a GAM
m_gam1 = gam(evi~s(doy, bs='cr'), data=d_sub)
m_gam2 = gam(evi~s(doy, bs='cr'), data=d_sub[ind,])

# Add smoothed values to the data frame
d_sub$evi_smooth1 = predict(m_gam1)
d_sub$evi_smooth2 = NA
d_sub$evi_smooth2[ind] = predict(m_gam2)


# Plot time series
p =ggplot(data=d_sub,
       aes(x=doy, 
           y=evi)) +
  geom_point(aes(colour=factor(QC))) + 
  geom_line(aes(y=evi_smooth1), colour='blue') +
  geom_line(aes(y=evi_smooth2), colour='orange') +
  labs(title=pixel) +
  theme_bw()

p



# Try segmenting the curve
m1 = lm(evi~doy, data=d_sub)
m2 = lm(evi_smooth1~doy, data=d_sub)
m3 = lm(evi_smooth2~doy, data=d_sub[ind,])

starting_breaks = c(100, 200,300)

m_seg1 = segmented(m1, 
                  seg.Z = ~doy,
                  psi = list(doy=starting_breaks),
                  control=seg.control(display=TRUE))

m_seg2 = segmented(m2, 
                   seg.Z = ~doy,
                   psi = list(doy=starting_breaks),
                   control=seg.control(display=TRUE))

m_seg3 = segmented(m3, 
                   seg.Z = ~doy,
                   psi = list(doy=starting_breaks),
                   control=seg.control(display=TRUE))

# Make predictions across the year
d_pred = data.frame(doy=c(-100:465), evi=NA, evi_smooth=NA)
d_pred$evi = predict(m_seg1, newdata=d_pred)
d_pred$evi_smooth1 = predict(m_seg2, newdata=d_pred)
d_pred$evi_smooth2 = predict(m_seg3, newdata=d_pred)

breaks1 = confint(m_seg1)
breaks2 = confint(m_seg2)
breaks3 = confint(m_seg3)

p + 
  geom_line(data=d_pred, colour='red') +
  geom_line(data=d_pred, aes(y=evi_smooth1), colour='blue') +
  geom_line(data=d_pred, aes(y=evi_smooth2), colour='orange') +
  geom_vline(xintercept=breaks1[,1], colour='red') +
  geom_vline(xintercept=breaks2[,1],colour='blue') +
  geom_vline(xintercept=breaks3[,1],colour='orange')


summary(m_seg2)
