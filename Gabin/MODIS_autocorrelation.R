# Look at Morane's I
#
#
# ************************

rm(list=ls())
dataDir = '~/MEGAsync/Projects/GrasslandPhenology/Data'

load(file.path(dataDir,'modis_pasture_A2017_square1.RData'))

# Calculate unique pixel ID's
x_values = unique(d_sq$x_MODIS)
y_values = unique(d_sq$y_MODIS)

d_sq$pixelID = paste0('x',match(d_sq$x_MODIS,x_values),
                      'y',match(d_sq$y_MODIS,y_values))

pixel_list = unique(d_sq$pixelID)

# Calculate distances between all pixels 
# (use ITM because it is a rectangular CRS)
ind = match(pixel_list,d_sq$pixelID)

pixels = cbind(d_sq$x_ITM[ind], d_sq$y_ITM[ind])
rownames(pixels) <- pixel_list


# Distance is in km
#d = as.matrix(dist(pixels, method="euclidean"))
d = as.matrix(dist(cbind(d_sq$x_ITM, d_sq$y_ITM), method="euclidean"))

ind2 = match(d_sq$pixelID, pixel_list)

# Calculate Moran's I for different distance bins

moran = data.frame(bins = seq(from=0.1,to=10,length=10)*1000,
                   I = NA,
                   I_lowerCI = NA,
                   I_upperCI=NA)

for (i in 1:nrow(moran)) {
  if (i==nrow(moran)) {
    lower = moran$bins[i]
    upper = max(d+10)
  } else {
    lower = moran$bins[i]
    upper = moran$bins[i+1]
  }
  w = which(d>=lower & d<upper, arr.ind=TRUE)
  I = cor.test(x=d_sq$evi[w[,1]], 
               y=d_sq$evi[w[,2]],
               method='pearson',
               conf.level=0.95)
  
  moran$I[i] = I$estimate
  moran$I_lowerCI[i] = I$conf.int[1]
  moran$I_upperCI[i] = I$conf.int[2]
  moran$N[i] = nrow(w)
  }

# Visualise Moran's I
library(ggplot2)

ggplot(data=moran,
       aes(x=bins,
           y=I,
           ymin=I_lowerCI,
           ymax=I_upperCI)) +
  geom_abline(intercept=0, slope=0, colour='blue') +
  geom_path(aes(y=-1/(N-1)), colour='grey') +
  geom_ribbon(alpha=0.3) +
  geom_point(size=2) +
  geom_path() +
  theme_bw() +
  labs(x='Minimum of Distance Class (km)', y='Moran\'s I')
