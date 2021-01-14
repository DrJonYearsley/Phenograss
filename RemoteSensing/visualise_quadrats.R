# Visualise phenology from quadrats
#
#
#
# **************************************************

rm(list=ls())
setwd('~/MEGAsync/Projects/GrasslandPhenology/RemoteSensing/')

library(ggplot2)

# Draw map of the squares
library(rgdal)

IR = readOGR(dsn='../Data/Quadrats/country.shp')
squares = readOGR('../Data/Quadrats/agriclimate_quadrats_Ireland.shp')

IR_modis = spTransform(IR, CRSobj = proj4string(squares))
plot(IR_modis)
plot(squares, add=T, col='black')

squares$ID = c(1:21)
plot(subset(squares, ID==8), add=T, col='red')




# *****************************************

load ("~/MEGAsync/Projects/GrasslandPhenology/Data/modis_pasture_A2018_square8.RData")

x_values = unique(d_sq$x_MODIS)
y_values = unique(d_sq$y_MODIS)

d_sq$pixelID = paste0('x',match(d_sq$x_MODIS,x_values),
                      'y',match(d_sq$y_MODIS,y_values))


table(d_sq$pixelID)

# Aggregate across all pixels in a square
d_ag = aggregate(cbind(evi,ndvi) ~ doy+satellite, data=d_sq, FUN=median, subset = (QC==1))





ggplot(data=d_sq, 
       aes(x=doy, 
           y=evi,
           colour=factor(QC))) +
  geom_point() +
  geom_smooth()


ggplot(data=d_ag, 
       aes(x=doy, 
           y=evi)) +
  geom_point() +
  geom_smooth()




library(gstat)
library(sp)
library(spacetime)

# Estimate spatial variogram of the data
d_sq = d_sq[order(d_sq$doy),]
d_spatial = d_sq
coordinates(d_spatial) = ~x_MODIS + y_MODIS
d_spatial$month = format(d_spatial$date, "%m")
months = unique(d_spatial$month)
for (d in 1:length(months)) {
  vg = variogram(evi~1, 
                 data=subset(d_spatial, month%in%months[d]),
                 cutoff=5000)
 
  if (d==1) {
    vg_df = as.data.frame(vg)
    vg_df$month=months[d]
  } else {
    tmp = as.data.frame(vg)
    tmp$month=months[d]
    vg_df = rbind(vg_df, tmp)
  }
  
# # Fit Variogram model
# vmod = fit.variogram(vg, vgm("Exp"))
}
  
#   # Plot the results
# plot(vg, vmod, main=paste0('Days ',start,' to ',start+40),add=T)


ggplot(data=vg_df,
       aes(x=dist, y=gamma, colour=month)) +
  geom_path(size=2) +
  scale_color_viridis_d('Month', option = "plasma") +
  labs(x='Lag (m)', 
       y='Semi-Variance', 
       title='Semivariograms for square 8 in 2018') +
  theme_bw() + theme(axis.title = element_text(size=16),
                     axis.text = element_text(size=16))

ggsave('variograms_square8_2018.png')
