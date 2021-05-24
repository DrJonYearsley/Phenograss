# Analyse MODIS phenology
#
#
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Jan 2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(segmented)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)

dataDir = '~/Research/Phenograss/Data/PhenologyOutput'

input_file_preffix = 'phenology'

# Import data --------
squares = c(1:9,13:21)
years = c(2002:2019)


# Filename segmented data
for (i in 1:length(squares)) {
  for (y in 1:length(years)) {
    filename = paste0(input_file_preffix,'_square_',squares[i],'_',years[y],'.RData')
    load(file.path(dataDir,filename))
    
    output_smoothed$square = squares[i]
    
    if (y==1 & i==1) {
      phenology = output_smoothed
    } else {
      phenology = rbind(phenology,
                        output_smoothed)
    }
  }
}

# Create a wide version of phenology
phenology_wide = pivot_wider(data=subset(phenology, warning==FALSE),
                             id_cols = c('pixelID','year', 'x_MODIS','y_MODIS','x_ITM','y_ITM','square'),
                             names_from = 'phase',
                             names_prefix = 'phase',
                             values_from = c(t,slope),
                             values_fn = mean)

# Centre the x and y coordinates
phenology_wide$x_ITM_centre = phenology_wide$x_ITM - mean(phenology_wide$x_ITM, na.rm=TRUE)
phenology_wide$y_ITM_centre = phenology_wide$y_ITM - mean(phenology_wide$y_ITM, na.rm=TRUE)

# Coordinates
locations = subset(phenology_wide, 
                   subset=year==years[1],
                   select=c('pixelID','x_ITM','y_ITM'))

# Calculate the centre coordinates for each square
tmp = aggregate(cbind(x_ITM_centre, y_ITM_centre)~square+year, 
                data=phenology_wide, 
                FUN=mean, 
                na.rm=TRUE)

rm(list='phenology')

# # Try an INLA model
# library(inlabru)
# library(INLA)


# Fit model for phenology dates of phase 1 -----


# Mixed model for phase 1
m_phase1 = lme(t_phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)

summary(m_phase1)

# # Fit this same model using INLA
# m_inla = inla(t_phase1~1+as.factor(year)+f(pixelID, model='iid'),
#               data=phenology_wide,
#               family='gaussian',
#               control.compute = list(dic = TRUE))


# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time

phenology_wide$dummy = interaction(phenology_wide$year, phenology_wide$square, sep='_')
gls_phase1 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)

summary(gls_phase1)


anova(gls_phase1, type='marginal')

acf(gls_phase1$residuals, lag.max=20, type='partial')  # Makes sense... spatial correlation roughly 1-2 km




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
