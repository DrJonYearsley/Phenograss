# Analyse MODIS phenology
#
#  Script modified to run on sonic (e.g. remove visualisations)
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Jan 2021
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

library(nlme)
library(tidyr)

dataDir = '~/data/Phenophases'

input_file_preffix = 'phenology'
output_file = 'gls_model_squares1_21_years_2002_2019.Rdata'

# Import data --------
squares = c(1:21)
years = c(2002:2019)

all_data_file = file.path(dataDir, 
                          paste0(input_file_preffix,
                                 '_years_',
                                 years[1],'_',years[2],
                                 '_squares_',
                                 min(squares),'_',max(squares),
                                 '.RData'))

if (exists(all_data_file)) {
  load(all_data_file)
} else {
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
  
  save(tmp, phenology_wide, file=all_data_file)
}



# +++++++++++++++++++++++++++++++++++++++++++++++
# Fit model for phenology dates of phase 1 -----


# Mixed model for phase 1
m_phase1 = lme(t_phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
               random= ~1 | pixelID,
               data=phenology_wide,
               na.action=na.exclude)



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



save(gls_phase1, file=output_file)
