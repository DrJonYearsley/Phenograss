# Fourier decomposition of the MODIS data
#
# Jon Yearsley
# July 2020
# ********************************************

rm(list=ls())

library(ggplot2)
library(rgdal)
library(timeSeries)
library(spectral)
library(sp)

dataDir = '~/MEGAsync/Projects/GrasslandPhenology/Data'

squares = 9
years = c(2017)

# Load HNV map
hnv = readOGR(dsn=file.path(dataDir,'HNV'),layer='HNVf-ED')


# Load MODIS data
flag = TRUE
for (f in c(1:length(years))) {
  for (s in 1:length(squares)) {
    filename = paste0('MODIS_pasture_A',years[f],'_square',squares[s],'.RData')
    load(file.path(dataDir,filename))
    
    if (flag) {
      d = d_sq
      flag=FALSE
    } else {
      d = rbind(d, d_sq)
    }
  }
}

# Calculate unique pixel ID's
x_values = unique(d$x_MODIS)
y_values = unique(d$y_MODIS)

d$pixelID = paste0('x',match(d$x_MODIS,x_values),
                      'y',match(d$y_MODIS,y_values))

pixel_list = unique(d$pixelID)

# Add in julian day
d$julian = as.numeric(julian(d$date))

# Visualise the spatial distribution of pixels
tmp = aggregate(evi~x_MODIS+y_MODIS, data=d, FUN=function(x){return(quantile(x, na.rm=TRUE, probs = c(0.01, 0.99)))})
tmp$evi_min = tmp$evi[,1]
tmp$evi_max = tmp$evi[,2]
tmp$evi_range = tmp$evi[,2] - tmp$evi[,1]

ggplot(data=tmp,
       aes(x=x_MODIS, y=y_MODIS, colour=evi_min)) +
  geom_point() +  
  scale_color_viridis_c() +
  labs(title='Min EVI') +
  theme_bw()


ggplot(data=tmp,
       aes(x=x_MODIS, y=y_MODIS, colour=evi_max)) +
  geom_point() + 
  scale_color_viridis_c() +
  labs(title='Max EVI') +
  theme_bw()

ggplot(data=tmp,
       aes(x=x_MODIS, y=y_MODIS, colour=evi_range)) +
  geom_point() + 
  scale_color_viridis_c(option='magma') +
  labs(title='Max-Min EVI') +
  theme_bw()


px=sample(pixel_list, size=5)

ggplot(data=d[d$pixelID %in% px,],
       aes(x=julian, y=evi, colour=as.factor(pixelID))) +
  geom_point() + 
  geom_smooth(method='loess', formula=y~x, se=FALSE, span=0.4) +
  scale_colour_brewer(palette='Set1') +
  labs(title=paste0('Pixel', paste0(px,collapse=', '))) +
  theme_bw()




# Perform FFT on a random pixel
x = d$evi[d$pixelID%in%px[1]]
t = d$julian[d$pixelID%in%px[1]]
# Reorder x
x = x[order(t)]
ff = fft(analyticFunction(x))

f = data.frame(power=Mod(ff)[c(2:floor(length(x)/2))])
f$smooth = filter(f$power,rep(1,3)/3, method='convolution')
f$x = c(1:nrow(f))

ggplot(data=f, aes(x=x, y=power)) +
  geom_point() +
  geom_smooth(method='loess', formula=y~x, span=0.1) +
  geom_line(aes(y=smooth)) +
  theme_bw()

ll = lsp(x, times=t, type='freq', alpha=0.05)


     