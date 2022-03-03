# Lyons Farm validation
#
# Compare remote sensing data to field data from Lyon's farm
#
# Jon Yearsley (Jon.Yearsley@ucd.ie)
# Dec 2021
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

# Lyons farm data
lyonsDataFile = "~/Research/Phenograss/Data/LyonsFarm/DM_production_PRG.csv"

# Modis data dir
modisDir = "~/Research/Phenograss/Data/MODIS/MODIS_v6_gri_format/"

# Long and lat of Lyon's farm grass plots (taken from google maps!)
lyonsFarm = data.frame(lon=-6.532056, lat=53.298288)

library(ggplot2)
library(terra)
library(mgcv)

# Import Lyons Farm data ------
lyons = read.csv(lyonsDataFile)


# Convert date to POSIX data
lyons$Date = as.Date(lyons$Date, format="%d/%m/%Y")

lyons$year = as.factor(as.numeric(format(lyons$Date,"%Y")))


# Add doy of year
lyons$doy = as.numeric(format(lyons$Date, "%j"))

# Define polygon to extract modis data from Lyon's farm
crdref <- "+proj=longlat +datum=WGS84"
pts <- vect(lyonsFarm, crs="epsg:4326")

# Create a buffer around this point
region = buffer(x = pts, width=300)


# Read MODIS data for Ireland for years in Lyon's farm data and crop the data

# Use a regulat expression to find terra and aqua files from relevant years
modisFiles = list.files(path=modisDir, 
                        pattern=paste0("MODIS_MOL[AT]_(", 
                                       paste0(unique(lyons$Year), collapse="|"),
                                       ")[[:graph:]]+.grd"), 
                        full.names=TRUE)

# Read one Modis file and define extent to crop raster
mod = rast(modisFiles[1])
# Layer names
names(mod)

# Project Lyons buffer region onto modis CRS
region_modis = project(region, y=mod)
pts_modis = project(pts, mod)

# Extract EVI layer and crop a little to be near the flight site
tmp = crop(mod, ext(region_modis))

writeRaster(tmp, filename = "tmp.tiff", overwrite=TRUE)


# Set up data frame and import MODIS data
evi = data.frame(x_MODIS=NA,
                 y_MODIS=NA, 
                 lon=lyonsFarm$lon, 
                 lat = lyonsFarm$lat, 
                 evi=NA, 
                 ndvi=NA, 
                 year=NA, 
                 doy=rep(NA, times=length(modisFiles)), 
                 QA = NA)


for (m in 1:length(modisFiles)) {
  mod = rast(modisFiles[m])
  evi_point = extract(mod, pts_modis, xy=TRUE)
  
  evi$year[m] = as.numeric(regmatches(modisFiles[m], 
                               regexpr(text=modisFiles[m],
                                       pattern=paste0("(", 
                                                      paste0(unique(lyons$Year), collapse="|"),
                                                      ")"))))
  
  evi$x_MODIS[m] = evi_point$x
  evi$y_MODIS[m] = evi_point$y
  evi$doy[m] = evi_point$DOY
  evi$evi[m] = evi_point$EVI
  evi$ndvi[m] = evi_point$NDVI
  evi$QA[m] = evi_point$PIXEL_RELIABILITY

  
  evi_region = crop(mod,ext(region_modis))
  xy = crds(evi_region)

  
  tmp = data.frame(pixelID = letters[1:nrow(xy)],
                   x_MODIS = xy[,1],
                   y_MODIS = xy[,2],
                   year=evi$year[m],
                   doy= values(evi_region, dataframe=TRUE)$DOY,
                   evi = values(evi_region, dataframe=TRUE)$EVI,
                   ndvi = values(evi_region, dataframe=TRUE)$NDVI,
                   QA = values(evi_region, dataframe=TRUE)$PIXEL_RELIABILITY)
  
  if (m==1) {
    evi2 = tmp
  } else {
    evi2 = rbind(evi2, tmp)
  }
}

evi2$year = as.factor(evi2$year)
evi$year = as.factor(evi$year)

# ++++++++++++++++++++++++++++++++++++++++++
# Visualise the data

ggplot(data=lyons,
       aes(x=doy,
           y=Pregraze.cover,
           colour=as.factor(Year),
           group=as.factor(Year))) + 
  geom_line() + 
  geom_point(size=2) +
  # geom_smooth(method="loess", 
  #             formula="y~x",
  #             span=0.5) +
  scale_colour_brewer("Year", palette="Dark2") +
  labs(x="Day of Year",
       y = "Pre-Graze Cover (kg DM/ha)") +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=24))

ggsave('~/Desktop/lyons_timeseries.png')


ggplot(data=subset(evi2, !outlier),
       aes(x=doy,
           y=evi,
           colour=year)) + 
  geom_point(size=2) + 
  geom_smooth(method="loess", 
              formula="y~x",
              span=0.4) +
  scale_colour_brewer("Year", palette="Dark2") +
  labs(x="Day of Year",
       y="Vegetation Index (EVI)") +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=24))

ggsave("~/Desktop/evi_timeseries.png")



# Add EVI data to Lyons Farm by taking an average +/- window around doy
window = 5
lyons$evi = NA
for (i in 1:nrow(lyons)) {
  sub = subset(evi2, year==lyons$year[i] & abs(doy-lyons$doy[i])<window & !outlier)
  if (nrow(sub)>0) {
    lyons$evi[i] = quantile(sub$evi, 0.75)
  }
}


# Visualise
ggplot(data=lyons,
       aes(x=Pregraze.cover,
           y=evi,
           colour=as.factor(year))) +
  geom_point(size=2) + 
  labs(x="Pre-graze cover (kg DM/ha)",
       y="Vegetation Index (EVI)") +
  scale_colour_brewer("Year", palette="Dark2") +
  theme_bw() + 
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=24))

ggsave("~/Desktop/evi_lyons.png")


summary(lm(evi~Pregraze.cover, data=lyons))

# ++++++++++++++++++++++++++++++++++++++++++++++++
# Analyse the MODIS and Lyons farm data ----------


# Smooth evi time series and remove outlier points

evi_gam = gam(evi~ s(doy, bs="tp", by=year), 
              data=evi2, 
              gamma=0.7, 
              na.action=na.exclude)
summary(evi_gam)
plot(evi_gam, residuals=TRUE, pch=20, col="black",  se=TRUE, seWithMean=TRUE)

tmp = predict(evi_gam, newdata=evi2, type='response', se.fit=TRUE)
evi2$pred = tmp$fit
evi2$predse = tmp$se.fit


# Calculate standadised residuals
evi2$resid_std = residuals(evi_gam) / evi2$predse

evi2$outlier = evi2$resid_std < -6


# Refit gam exlcuding outliers
evi_gam2 = gam(evi~ s(doy, bs="tp", by=year), 
               data=evi2, 
               subset=!evi2$outlier,
               gamma=0.7, 
               na.action=na.exclude)
summary(evi_gam2)
plot(evi_gam2, residuals=TRUE, pch=20, se=TRUE, seWithMean=TRUE)

tmp = predict(evi_gam2, newdata=evi2, type='response', se.fit=TRUE)
evi2$pred2 = tmp$fit
evi2$predse2 = tmp$se.fit



# Add gam predcition into the lyons farm
tmp = predict(evi_gam2, newdata=lyons, se.fit=TRUE)
lyons$gam_pred=tmp$fit
lyons$gam_predse=tmp$se.fit


lyons_gam = gam(Pregraze.cover ~ s(doy,bs="tp",by=year), data=lyons)
plot(lyons_gam)


test = data.frame(year=rep(c("2015","2016"), each=200), doy=rep(c(1:200), times=2))
test$lyons = predict(lyons_gam, newdata=test)
test$modis = predict(evi_gam2, newdata=test)


# Visualise -------
ggplot(data=test,
       aes(x=lyons,
           y=modis,
           colour=year)) + 
  geom_point() + 
  theme_bw()
 



ggplot(data=test,
       aes(x=doy,
           y=lyons,
           colour=year)) +
  geom_point() +
  geom_point(aes(y=modis))
  