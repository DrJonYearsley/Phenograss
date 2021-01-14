# Script to compare phenology across years
#
#
# Jon Yearsley
#
# ++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(segmented)
library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)

dataDir = '~/MEGAsync/Projects/GrasslandPhenology/Data/PhenologyOutput/'

input_file_preffix = 'phenology'

# Import data --------
squares = c(1:9)
years = c(2014:2017)

#starting_breaks = c(50, 100, 200, 300)


# Filename segmented data
for (s in 1:length(squares)) {
  for (y in 1:length(years)) {
    filename = paste0(input_file_preffix,'_square_',squares[s],'_',years[y],'.RData')
    load(file.path(dataDir,filename))
    
    output_smoothed$square = squares[s]
    
    if (y==1 & s==1) {
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

# Calculate the centre coordinates for each square

tmp = aggregate(cbind(x_ITM_centre, y_ITM_centre)~square+year, 
                data=phenology_wide, 
                FUN=mean, 
                na.rm=TRUE)



# Fit a variogram to the phenology estimates for each year ---
library(gstat)
library(sp)

# Look at phase 1
d_spatial = phenology_wide
coordinates(d_spatial) = ~x_ITM + y_ITM

range_estimates = array(NA, dim=length(years))

for (y in 1:length(years)) {
  vg = variogram(t_phase1~1, 
                 data=subset(d_spatial, is.finite(t_phase1) & year%in%years[y]),
                 cutoff=5000)
  
  if (y==1) {
    vg_df = as.data.frame(vg)
    vg_df$year=years[y]
  } else {
    tmp = as.data.frame(vg)
    tmp$year=years[y]
    vg_df = rbind(vg_df, tmp)
  }
  
  # # Fit Variogram model
  vmod = fit.variogram(vg, vgm("Exp"))
  
  range_estimates[y] = vmod$range[2]
}

range_estimates

ggplot(data=vg_df,
       aes(x=dist,
           y=gamma,
           colour=factor(year))) +
  geom_point()


# Fit model for phenology dates of phase 1 -----


# Mixed model for phase 1
m_phase1 = lme(t_phase1~as.factor(year),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)

# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time
gls_phase1 = gls(t_phase1~1+factor(year) + (x_ITM_centre + y_ITM_centre),
               correlation = corExp(value = 500,
                                    form=~x_ITM_centre + y_ITM_centre | year, 
                                    nugget=FALSE,
                                    fixed=T),
               data=phenology_wide,
               na.action=na.exclude)


gls_phase1 = gls(t_phase1~ (x_ITM_centre + y_ITM_centre),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)


summary(gls_phase1)


# Fit model for phenology dates of phase 2 - phase 1
gls_phase2 = gls(I(t_phase2 - t_phase1)~1+factor(year),
                 correlation = corExp(value = 500,
                                      form=~x_ITM + y_ITM | year, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)

summary(gls_phase2)


# Fit model for phenology dates of phase 3 - phase 1
gls_phase3 = gls(I(t_phase3 - t_phase1)~1+factor(year),
                 correlation = corExp(value = 500,
                                      form=~x_ITM + y_ITM | year, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phenology_wide,
                 na.action=na.exclude)

summary(gls_phase3)


# +++++++++++++++++++++++++++++++++++++++++
# Plot predictions from phase 1 model -----

phenology_wide$t1_pred = predict(gls_phase1)

ggplot(data=phenology_wide,
       aes(x=y_MODIS,
           y=t_phase1)) +
  geom_point()
  

ggplot(data=phenology_wide,
       aes(x=x_ITM,
           y=y_ITM,
           colour=t1_pred)) +
  geom_point(size=0.05) +
  scale_colour_viridis_c('Day of Year',option='magma') +
  theme_bw()

ggplot(data=phenology_wide,
       aes(x=x_ITM,
           y=y_ITM,
           colour=t_phase1)) +
  geom_point(size=0.05) +
  scale_colour_viridis_c('Day of Year',option='magma') +
  theme_bw()
