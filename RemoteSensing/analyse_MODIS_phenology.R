# Analyse MODIS phenology
#
# Analyse MODIS phenology and add in environmental correlates
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Sept 2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
# 
# library(segmented)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)
library(raster)
library(sp)
library(rgdal)
library(gstat)
library(DHARMa)

dataDir = '/media/jon/MODIS_data/PhenologyOutput_toAnalyse/'
envFile = '/media/jon/MODIS_data/Data_created/all_envData.RData'
pastureFile = "/media/jon/MODIS_data/Corine/corine2018_pasturecover_All_Ireland.gri"

input_file_preffix = 'phenology'

# Import data --------
squares = c(2:14,16:21)  # Squares 1 and 15 have too little data to use


# Import all the data
for (i in 1:length(squares)) {
  filename = list.files(path=dataDir, 
                        pattern=paste0(input_file_preffix,"_square_",squares[i],".RData"),
                        full.names = TRUE)
  load(filename)  
  
  if (i==1) {
    phenology_long = phenology
  } else {
    phenology_long = rbind(phenology_long, phenology)
  }
  
  rm(list="phenology")
}
  
# Create a wide version of phenology
phenology_wide = pivot_wider(data=subset(phenology_long, warning==FALSE),
                             id_cols = c('pixelID','year', 'x_MODIS','y_MODIS','x_ITM','y_ITM','square'),
                             names_from = 'phase',
                             names_prefix = 'phase',
                             values_from = c(t),
                             values_fn = mean)

# Centre the x and y coordinates on the whole of Ireland
phenology_wide$x_ITM_centre = phenology_wide$x_ITM - mean(phenology_wide$x_ITM, na.rm=TRUE)
phenology_wide$y_ITM_centre = phenology_wide$y_ITM - mean(phenology_wide$y_ITM, na.rm=TRUE)

# Coordinates
years=unique(phenology_wide$year)
locations = subset(phenology_wide, 
                   subset=year==years[1],
                   select=c('pixelID','x_ITM','y_ITM'))

# Calculate the centre coordinates for each square
tmp = aggregate(cbind(x_ITM_centre, y_ITM_centre)~square, 
                data=phenology_wide, 
                FUN=mean, 
                na.rm=TRUE)

phenology_wide$x_square = NA
phenology_wide$y_square = NA
ind = match(phenology_wide$square, tmp$square)
phenology_wide$x_square = tmp$x_ITM_centre[ind]
phenology_wide$y_square = tmp$x_ITM_centre[ind]
head(phenology_wide)

rm(list='phenology_long')



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate spatial autocorrelation in a typical square

for (s in c(2:14,16:21)) {
  tmp = subset(phenology_wide, year==2014 & square%in% s)
  
  ind = is.finite(tmp$phase1)
  
  vg.obs = variogram(phase1~1, data=tmp[ind,], locations=~x_ITM_centre+y_ITM_centre)
  vg.fit = fit.variogram(vg.obs, model=vgm(psill=1,
                                           model="Sph",
                                           range=NA,
                                           kappa=NA, 
                                           nugget=NA))
  
  print(paste(s,sum(ind),round(vg.fit$range[2]), round(vg.fit$psill[1])))
}


# Spherical model seems best and a range of roughly 1000 is broadly appropriate and a nugget of 100
plot(vg.obs, vg.fit,cutoff=2000)


# Look across squares in a singe year
tmp2 = subset(phenology_wide, year==2008 & square%in%c(1))
ind = is.finite(tmp2$phase1)
vg.obs = variogram(phase1~1+square, 
                   data=tmp2[ind,], 
                   locations=~x_ITM_centre+y_ITM_centre)
vg.fit = fit.variogram(vg.obs, model=vgm(psill=NA,
                                         model="Exp",
                                         range=500,
                                         kappa=NA, 
                                         nugget=NA))

vg.fit
plot(vg.obs, vg.fit,cutoff=2000)


# ++++++++++++++++++++++++++++++++++++++++++++++++
# Explore data with graphs ----------


# Plot square averages
ggplot(data=agg, 
      aes(x=y_ITM,
          y=phase2,
          colour=year)) +
  geom_smooth(method='lm', formula='y~x') +
  geom_point() + 
  scale_color_viridis_d()

ggplot(data=agg, 
       aes(x=year,
           y=phase2,
           colour=y_ITM)) +
  geom_smooth(method='loess', formula='y~x', se=FALSE) +
  geom_point() + 
  scale_color_viridis_c()

# Create subset
tmp=subset(d, year%in%c(2002:2009) & 
             square==2 &
             SLOPE_PERCENT>0 & 
             !is.na(aspect_cat))

ggplot(data=tmp,
       aes(y=y_MODIS,
           x=x_MODIS,
           colour=ASPECT)) +
  geom_point() +
theme_bw()

ggplot(data=tmp,
       aes(y=anom1,
           x=y_MODIS,
           colour=aspect_cat)) +
  geom_point() +
  geom_smooth(method='lm',
              se=T)


ggplot(data=tmp,
       aes(y=anom1,
           x=aspect_cat,
           colour=aspect_cat)) +
  geom_boxplot() 

ggplot(data=tmp,
       aes(y=anom3,
           x=aspect_cat,
           colour=aspect_cat)) +
  geom_boxplot() 


m = lm(t_phase1~factor(year)+aspect_cat*SLOPE_PERCENT + ELEVATION + SOIL_TYPE  + square, 
       data=subset(d,p_pasture>0.99 & year%in%c(2002:2009)))

d$dummy = interaction(d$year, d$square, sep='_')
gls_phase1 = gls(t_phase1~1+factor(year) + aspect_cat*SLOPE_PERCENT + ELEVATION + SOIL_TYPE + (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=subset(d, p_pasture>0.99 & year%in%c(2002:2009)),
                 na.action=na.exclude)


summary(gls_phase1)

aggregate(p_pasture~SOIL_TYPE, data=d, FUN=median)

# ++++++++++++++++++++++++++++++++++++++++++++++++
# Fit model for phenology dates of phase 1 -----


# Look at one square for a couple of year and fit the spatial model

test = subset(phenology_wide, year%in%c(2005:2009) & square==13)


# Fit a spatial model using a spherical semi-variogram

# Generate a blocking variable
test$dummy = interaction(test$year, test$square, sep='_')
test$dummy = as.factor(test$dummy)
test$year = as.factor(test$year)


gls_phase1 = gls(phase1~1 + year + (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = c(500, 0.5),   # Range and then nugget
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=TRUE,
                                      fixed=F),
                 data=test,
                 na.action=na.exclude,
                 method="REML")

summary(gls_phase1)  # range is 500 - 1000 and nugget is 0.5-0.7


# # Fit a similar model using glmmTMB but with an exponential spatial structure
# # but doesn't work that well
# library(glmmTMB)
# test$pos = numFactor(test$x_ITM_centre, test$y_ITM_centre)
# tmb_phase1 = glmmTMB(phase1~1 + (x_ITM_centre + y_ITM_centre) + exp(pos+0|dummy),
#                  data=test,
#                  na.action=na.exclude)
# 
# summary(tmb_phase1)


# +++++++++++++++++++++++++++++++++++++


# ++++++++++++++++++++++++++++++++++++++
# Calculate average across each square and then 
# fit a model to the square averages

phenology_wide$year = as.factor(phenology_wide$year)
pheno_ave = aggregate(cbind(x_ITM_centre, y_ITM_centre, phase1, phase2, phase3)~year+square, 
                      data=phenology_wide, 
                      FUN=mean, 
                      na.rm=TRUE)


lm_phases = lm(cbind(phase1,phase2,phase3)~1+year+x_ITM_centre+y_ITM_centre,
               data=pheno_ave)
lm_phase1 = lm(phase1~1+year+x_ITM_centre+y_ITM_centre,
               data=pheno_ave)
lm_phase2 = lm(phase2~1+year+x_ITM_centre+y_ITM_centre,
               data=pheno_ave)
lm_phase3 = glm(phase3~1+year+x_ITM_centre+y_ITM_centre,
               data=pheno_ave,
               family=Gamma)

plot(lm_phase3)

summary(lm_phase1)
anova(lm_phase1)





m_phase1 = lme(phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=test,
               na.action=na.exclude)
m_phase2 = lme(phase2~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=test,
               na.action=na.exclude)
m_phase3 = lme(phase3~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
                                                   random= ~1 | pixelID,
                                                   data=test,
                                                   na.action=na.exclude)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Mixed model for phase 1, 2 and 3
m_phase1 = lme(phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)


m_phase2 = lme(phase2~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)
m_phase3 = lme(phase3~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)


summary(m_phase1)
plot(ACF(m_phase1, lag.max=20, type='partial')[-1,], alpha=0.05)  
plot(ACF(m_phase2, lag.max=20, type='partial')[-1,], alpha=0.05)  
plot(ACF(m_phase3, lag.max=20, type='partial')[-1,], alpha=0.05)  
# Makes sense... spatial correlation roughly 1-2 km
# Roughly the same for all 3 phenophases



# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time
phenology_wide$dummy = interaction(phenology_wide$year, phenology_wide$square, sep='_')
phenology_wide$dummy = as.factor(phenology_wide$dummy)

gls_phase1 = gls(phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = 1000,
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=subset(phenology_wide, square%in%c(3:14) & p_pasture>0.99),
                 na.action=na.exclude)

summary(gls_phase1)
ACF(gls_phase1, lag.max=20, resType='response',form=)

anova(gls_phase1, type='marginal')

acf(gls_phase1$residuals, lag.max=20, type='partial')  # Makes sense... spatial correlation roughly 1-2 km


res = residuals(gls_phase1)

library(gstat)


library(emmeans)

m_eff = emmeans(gls_phase1, spec='year')
sos = as.data.frame(summary(m_eff))
contrast(m_eff, infer=c(T,T))

as.Date(paste(round(sos$emmean),sos$year),format="%j %Y")

# WOrk out date from day of year
as.Date(paste('52',sos$year),format="%j %Y")

gls_phase1_v2 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                 data=phenology_wide,
                 na.action=na.exclude)

ggplot(data=as.data.frame(summary(m_eff)),
       aes(x=year,
           y=emmean,
           ymin=lower.CL,
           ymax=upper.CL)) +
  geom_pointrange() + 
  labs(x='Year',
       y='Day of Year') +
  theme_bw() + 
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18))


save(gls_phase1, phenology_wide, file='gls_model_2002_2019.Rdata')





# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load up environmental data

# Load data on pasture landcover
pasture = raster::raster(pastureFile)
# Add proportion of pasture for each pixel
y = SpatialPoints(phenology_wide[,c(3,4)], proj4string = crs(pasture))

tmp = extract(pasture, y)
phenology_wide$p_pasture = tmp

# Import some environmental data
load(envFile)

# Assign a pixel ID to env data
all_squares_df$pixelID = NA
for (s in unique(phenology_wide$square)) {
  pheno_sub = phenology_wide[phenology_wide$square==s,c(1:7)]
  pixelID_List = unique(pheno_sub$pixelID)
  inds = match(pixelID_List, pheno_sub$pixelID)
  for (p in 1:length(pixelID_List)) {
    ind = abs(all_squares_df$x_MODIS - pheno_sub$x_MODIS[inds[p]])<50 & 
      abs(all_squares_df$y_MODIS - pheno_sub$y_MODIS[inds[p]])<50 
    all_squares_df$pixelID[ind] = pixelID_List[p]
  }
}

nam = names(all_squares_df)
nam[7] = "square"
names(all_squares_df) = nam

d = merge(phenology_wide, all_squares_df)

# Add a categorical aspect variable
d$aspect_cat = NA
delta = 10
d$aspect_cat[d$ASPECT>360-delta | d$ASPECT<delta ] = "N"
d$aspect_cat[d$ASPECT>45 & d$ASPECT<90+45 ] = "E"
d$aspect_cat[d$ASPECT>180-delta & d$ASPECT<180+delta ] = "S"
d$aspect_cat[d$ASPECT>270-delta & d$ASPECT<270+delta ] = "W"
d$aspect_cat = as.factor(d$aspect_cat)

d$year = as.factor(d$year)

# Calculate the average SOS, POS and EOS for each square for each year and then work out an anomaly
# This could be useful to look at the effect of environment (e.g. aspect and slope)
agg = aggregate(cbind(x_ITM,y_ITM, phase1,phase2,phase3)~year+square, data=d, FUN=median, na.rm=TRUE)

for (s in unique(d$square)) {
  d$anom1[d$square==s] = d$phase1[d$square==s] - agg$phase1[match(d$year[d$square==s], agg$year[agg$square==s])]
  d$anom2[d$square==s] = d$phase2[d$square==s] - agg$phase2[match(d$year[d$square==s], agg$year[agg$square==s])]
  d$anom3[d$square==s] = d$phase3[d$square==s] - agg$phase3[match(d$year[d$square==s], agg$year[agg$square==s])]
}







# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load pre-save model fit data ------

load(file.path(dataDir,'gls_model_2002_2019.Rdata'))
phenology_wide$phase1_residual = residuals(gls_phase1)
phenology_wide$phase1_fit = predict(gls_phase1)

library(sf)
library(rgdal)
# Produce a map of fitted values ---
tmp = subset(phenology_wide, year==2019 & !is.na(phase1_fit) & square%in%c(1:21))
IR = readOGR(dsn='~/MEGAsync/Projects/GrasslandPhenology/Data/Quadrats/country.shp')
squares = readOGR('~/MEGAsync/Projects/GrasslandPhenology/Data/Quadrats/agriclimate_quadrats_Ireland.shp')

tmp = aggregate(phase1_fit~square, data=subset(phenology_wide, year==2019), FUN=mean, na.rm=TRUE)
nColour = 7
vals = round(seq(49,
                 55,
                 length.out=nColour), 
             digits=2)
tmp$fit = NA
col_lab = c()
for (i in 2:nColour) {
  ind = tmp$phase1_fit>=vals[i-1] & tmp$phase1_fit<vals[i]
  col_lab = c(col_lab, paste0(vals[i-1],' - ',vals[i]))
  tmp$fit[ind] = col_lab[i-1]
}

squares@data$fit = NA
squares@data$fit[tmp$square] = tmp$fit

squares_sf = st_as_sf(squares)
IR_sf = st_as_sf(IR)

ggplot() +
  geom_sf(data=IR_sf, fill='white') +
  geom_sf(data=squares_sf, aes(fill=factor(fit, ordered=T, levels=rev(col_lab)))) +
  scale_fill_brewer('Start of Season\n(day of year)',
                      palette='Greens', 
                      direction=1,na.value='darkgray') +
  theme_void()

ggsave(filename = 'gls_fitted_map_year2019.png', height=10, units='cm')






# Calculate some R-squared between predictions and data -----
m = lm(t_phase1~factor(square) + phase1_fit, 
       data=subset(phenology_wide, square=20))

summary(m)

ggplot(data=subset(phenology_wide, square==20),
       aes(x=phase1_fit,
           y=t_phase1,
           colour=factor(year))) +
  geom_point() + 
  geom_abline(slope=1, 
              intercept=0, 
              colour='darkred', 
              linetype='dashed',
              size=1.5) +
  scale_colour_viridis_d('Year') +
  labs(x='Fitted Value (day of year)',
       y='Pixel Start of Season (day of year)') +
  theme_bw() + 
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18))


ggsave(filename = 'gls_fitted_vs_data_square20_year2019.png', height=20, units='cm')



# Plot model fit ----
tmp = subset(phenology_wide, square==20 & year==2019)
nColour = 5
vals = round(seq(52.5,
           52.6,
           length.out=nColour-1), digits=2)
col_lab = paste('< ',vals[1])
tmp$fit = col_lab[1]
for (i in 2:(nColour-1)) {
  ind = tmp$phase1_fit>=vals[i-1] & tmp$phase1_fit<vals[i]
  col_lab = c(col_lab, paste0(vals[i-1],' - ',vals[i]))
  tmp$fit[ind] = col_lab[i]
}


ind = tmp$phase1_fit>=vals[i]
col_lab = c(col_lab, paste0('> ',vals[i]))
tmp$fit[ind] = col_lab[i+1]
tmp$fit[is.na(tmp$phase1_fit)] = NA
tmp$fit = factor(tmp$fit, ordered=T, levels=rev(col_lab))

ggplot(data=tmp,
       aes(x=x_MODIS,
           y=y_MODIS,
           fill=fit)) +
  geom_tile(colour='darkblue') +
  coord_equal() +
  scale_fill_brewer('Start of Season\n(day of year)',
                    palette='Greens', 
                    direction=1,na.value='darkgray') +
  labs(x='X Coord (MODIS CRS)',
       y='Y Coord (MODIS CRS)') +
  theme_bw()

ggsave(filename = 'gls_fit_square20_year2019.png',height=10,units='cm')



# Plot residuals ----
nColour = 7
vals = seq(-50,50,by=100/(nColour-2))
col_lab = '< -50'
tmp$residual = col_lab[1]
for (i in 2:(nColour-1)) {
  ind = tmp$phase1_residual>=vals[i-1] & tmp$phase1_residual<vals[i]
  col_lab = c(col_lab, paste0(vals[i-1],' - ',vals[i]))
  tmp$residual[ind] = col_lab[i]
}
ind = tmp$phase1_residual>=vals[i]
col_lab = c(col_lab, paste0('> ',vals[i]))
tmp$residual[ind] = col_lab[i+1]
tmp$residual[is.na(tmp$phase1_residual)] = NA
tmp$residual = factor(tmp$residual, ordered=T, levels=rev(col_lab))

ggplot(data=tmp,
       aes(x=x_MODIS,
           y=y_MODIS,
           fill=residual)) +
  geom_tile(colour='darkblue') +
  coord_equal() +
  scale_fill_brewer('Residuals (days)',
                    palette='RdBu', 
                    direction=1,
                    na.value='darkgray') +
  labs(x='X Coord (MODIS CRS)',
       y='Y Coord (MODIS CRS)') +
  theme_bw()

ggsave(filename = 'gls_residual_square20_year2019.png',height=10,units='cm')

# # INLA version of spatial model
# Mesh <- inla.mesh.2d(locations, 
#                      max.edge = c(10, 10))



# Permutation test ---------
# Test whether SOS is more variable after 2011

library(emmeans)
# Calculate fitted effect across years
m_eff = emmeans(gls_phase1, spec='year')
sos = as.data.frame(summary(m_eff))


var_stat = function(d, y1, y2) {
  # Function to evaluate variation before and after a set date
  return(sd(d$emmean[d$year%in%y2]) / sd(d$emmean[d$year%in%y1]))
}

mean_stat = function(d, y) {
  # Function to evaluate variation before and after a set date
  return(mean(d$emmean[d$year>y]) / mean(d$emmean[d$year<=y]))
}

y1 = c(2002:2010)
y2=c(2011:2019)
nPerm = 10000
stat_array = array(NA, dim=c(nPerm,1))
ind = sos$year%in%c(y1,y2)
for (p in 1:nPerm) {
  d_perm = data.frame(year=sample(sos$year[ind]), emmean=sos$emmean[ind])
  stat_array[p] = var_stat(d_perm, y1,y2)
}

stat_obs = var_stat(sos, y1, y2)

# Empirical p-value
(sum(stat_array>=stat_obs)+1)/(nPerm+1)

hist(stat_array, 100)
