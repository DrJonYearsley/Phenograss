# Analyse MODIS phenology
#
# Analyse MODIS phenology and add in environmental correlates
# for the EPA final report
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# May 2022
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())


library(ggplot2)
library(nlme)
library(tidyr)
library(viridisLite)
library(terra)
library(gstat)
library(DHARMa)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set parameters -------------
dataDir = '/media/jon/MODIS_data/PhenologyOutput_toAnalyse/'         # Folder with phenophase estimate data
envFile = '/media/jon/MODIS_data/Data_created/all_envData.RData'     # Folder with environmental data
pastureFile = "/media/jon/MODIS_data/Corine/corine2018_pasturecover_All_Ireland.grd"   # Raster of % pasture cover for Ireland

input_file_prefix = 'phenology' # Prefix for input files
calculate_variograms = FALSE    # If TRUE calculate spatial variograms for the phenopahses
squares = c(2:13,16:18)         # The squares to include in the analysis
yearList = c(2003:2019)         # The years to include in the analysis

phenophase = 'SOS'              # Phenophase to analyse: SOS, POS, EOS, LOS

# Note: The number of squares is reduced in order to keep calculations 
# within the workstation's memory limit (32GB)








# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import data --------

if (file.exists(file.path(dataDir, 'combine_phenology_data_final_report.RData'))) {
  # If data have already been pre-processed and placed into one file
  
  load(file.path(dataDir, 'combine_phenology_data_final_report.RData'))
} else {
  # If data haven't been pre-processed then import all individual files and combine them
    for (i in 1:length(squares)) {
    filename = list.files(path=dataDir, 
                          pattern=paste0(input_file_prefix,"_square_",squares[i],".RData"),
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
  phenology_wide = pivot_wider(data=subset(phenology_long, warning==FALSE & upperCI-lowerCI<50),
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
  
  
  # ++++++++++++++++++++++++++++++++++++++++++  
  # Add in pasture cover data ---------------
  
  
  # Load data on pasture landcover
  pasture = terra::rast(pastureFile)
  # Add proportion of pasture for each pixel
  y = vect(phenology_wide[,c(3,4)], geom=c('x_MODIS','y_MODIS'), crs(pasture))
  
  tmp = extract(pasture, y, xy=FALSE)
  phenology_wide$p_pasture = tmp$layer  
}




# Create a variable with the phenophase data
if (phenophase=="SOS") {
  phenology_wide$phenophase = phenophase_wide$phase1
} else if (phenophase=="POS") {
  phenology_wide$phenophase = phenophase_wide$phase2
} else if (phenophase=="EOS") {
  phenology_wide$phenophase = phenophase_wide$phase3
} else if (phenophase=="LOS") {
  phenology_wide$phenophase = phenophase_wide$phase3 - phenophase_wide$phase1
}  




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Spatial variogram calculation ------------
if (calculate_variograms) { 
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Calculate spatial autocorrelation in typical squares
  
  # Look across squares in a singe year
  squareList = c(16, 6, 4, 10, 8, 11)
  yearList= 2017
  for (s in 1:length(squareList)) {
    tmp2 = subset(phenology_wide, year==yearList & square%in%squareList[s])
    ind = is.finite(tmp2$phenophase)
    vg.obs = variogram(phenophase~1, 
                       data=tmp2[ind,], 
                       locations=~x_ITM_centre+y_ITM_centre)
    vg.fit = fit.variogram(vg.obs, model=vgm(psill=NA,
                                             model="Exp",
                                             range=500,
                                             kappa=NA, 
                                             nugget=NA))
    
    
    # Place the empirical and fitted variograms in a data frame
    if (s==1) {
      vg_df = data.frame(type='Empirical',square=squareList[s], year=yearList, dist=vg.obs$dist, value = vg.obs$gamma)
      
      tmp = variogramLine(vg.fit, maxdist = max(vg.obs$dist))
      
      vg_df = rbind(vg_df, data.frame(type='Fitted',square=squareList[s], year=yearList, dist=tmp$dist, value = tmp$gamma))
    } else {
      vg_df = rbind(vg_df, data.frame(type='Empirical',square=squareList[s], year=yearList, dist=vg.obs$dist, value = vg.obs$gamma))
      
      tmp = variogramLine(vg.fit, maxdist = max(vg.obs$dist))
      vg_df = rbind(vg_df, data.frame(type='Fitted',square=squareList[s], year=yearList, dist=tmp$dist, value = tmp$gamma))
    }
  }
  
  
  ggplot(data=subset(vg_df, type=='Empirical'),
         aes(x=dist,
             y=value,
             colour=factor(square))) +
    geom_point(size=4) + 
    geom_path(data=subset(vg_df, type=='Fitted'), size=1.5) +
    lims(x=c(0,3000)) +
    labs(x='Distance between pixels (m)',
         y=paste0('Semi-variance for ',phenophase)) +
    scale_colour_brewer('10 km Square',palette='Set2') +
    theme_bw() + 
    theme(axis.text = element_text(size=18),
          axis.title = element_text(size=20),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14),
          legend.position = 'right',
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  ggsave(file=paste0('figure_2_4_variogram_',phenophase,'_',yearList,'.png'), dpi=300, width=6, height=5, units='in')
  ggsave(file=paste0('figure_2_4_variogram_',phenophase,'_',yearList,'.tiff'), dpi=300, width=6, height=5, units='in')
  
  
  
}









# ++++++++++++++++++++++++++++++++++++++++++++++++
# Fit model for phenology dates of phase -----


# Initial data wrangling
phenology_wide$year = as.factor(phenology_wide$year)
phenology_wide$dummy = interaction(phenology_wide$year, phenology_wide$square, sep='_')
phenology_wide$dummy = as.factor(phenology_wide$dummy)

# Define a time factor for pre 1st Jan 2011 and post 1st Jan 2011
phenology_wide$timeFac = factor(phenology_wide$year%in%c(2003:2010))





# 
# # ++++++++++++++++++++++++++++++++++++++
# # Fit models to average phenology for each square --------
# # Calculate average across each square and then 
# # fit a model to the square averages
# pheno_ave = aggregate(cbind(x_ITM_centre, y_ITM_centre, phase1, phase2, phase3)~year+square,
#                       data=phenology_wide,
#                       FUN=mean,
#                       na.rm=TRUE)
# 
# 
# lm_phases = lm(cbind(phase1,phase2,phase3)~1+year+x_ITM_centre+y_ITM_centre + I(x_ITM_centre^2) + I(y_ITM_centre^2) + x_ITM_centre:y_ITM_centre,
#                data=pheno_ave)
# summary(lm_phases)

# lm_phase1 = lm(phase1~1+year+x_ITM_centre+y_ITM_centre,
#                data=pheno_ave)
# lm_phase2 = lm(phase2~1+year+x_ITM_centre+y_ITM_centre,
#                data=pheno_ave)
# lm_phase3 = lm(phase3^7~1+year+x_ITM_centre+y_ITM_centre,
#                data=pheno_ave)
# 
# # Do some validation
# sim1=simulateResiduals(lm_phase1)
# plot(sim1)
# sim2=simulateResiduals(lm_phase2)
# plot(sim2)
# sim3=simulateResiduals(lm_phase3)
# plot(sim3)
# # Looks OK
# 
# 
# summary(lm_phases)
# 
# 
# summary(lm_phase3)
# anova(lm_phase3)
# 


# 
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Mixed model for phase 1, 2 and 3
# lme_phase1 = lme(phase1~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
#                random= ~1 | pixelID,
#                data=phenology_wide,
#                na.action=na.exclude)
# 
# 
# lme_phase2 = lme(phase2~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
#                random= ~1 | pixelID,
#                data=phenology_wide,
#                na.action=na.exclude)
# 
# lme_phase3 = lme(phase3~as.factor(year)+ (x_ITM_centre + y_ITM_centre),
#                random= ~1 | pixelID,
#                data=phenology_wide,
#                na.action=na.exclude)
# 
# 
# 
# summary(lme_phase1)
# plot(ACF(lme_phase1, lag.max=20, type='partial')[-1,], alpha=0.05)  
# plot(ACF(lme_phase2, lag.max=20, type='partial')[-1,], alpha=0.05)  
# plot(ACF(lme_phase3, lag.max=20, type='partial')[-1,], alpha=0.05)  
# # Makes sense... spatial correlation roughly 1-2 km (one unit is 250m)
# # Roughly the same for all 3 phenophases
# 


# Fit a generalised least squares with spatial autocorrelation
# spatial structure is fixed to reduce computation time


phase_sub = subset(phenology_wide, 
                    p_pasture>0.9 & year%in% yearList, 
                    select=c('x_ITM_centre','y_ITM_centre','dummy','year','phenophase'))

gls_phase = gls(phenophase~1 + year + (x_ITM_centre+y_ITM_centre + I(x_ITM_centre^2) + I(y_ITM_centre^2) + x_ITM_centre:y_ITM_centre),
                 correlation = corExp(value = 500,
                                      form=~x_ITM_centre + y_ITM_centre| dummy, 
                                      nugget=FALSE,
                                      fixed=T),
                 data=phase_sub,
                 na.action=na.exclude)




save(gls_phase, phase_sub, file=paste0('gls_',phenophase,'_model_2003_2019_quadratic_final_report.Rdata'))









# Visualise phenophase variation across years ----------

library(emmeans)
# Calculate fitted effect across years
m_eff = emmeans(gls_phase, spec='year')
phase_emm = as.data.frame(summary(m_eff))



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

ggsave(filename=paste0(phenophase,'_time_series_2003_2019.png'), width=15, height=4)





# +++++++++++++++++++++++++++++++++++++++++++++
# Permutation test ---------
# Test whether phenophase is more variable after 2011

var_stat = function(d, y1, y2) {
  # Function to evaluate variation before and after a set date
  return(sd(d$emmean[d$year%in%y2]) / sd(d$emmean[d$year%in%y1]))
}

mean_stat = function(d, y) {
  # Function to evaluate variation before and after a set date
  return(mean(d$emmean[d$year>y]) / mean(d$emmean[d$year<=y]))
}

y1 = c(2003:2010)
y2=c(2011:2021)
nPerm = 10000
stat_array = array(NA, dim=c(nPerm,1))
ind = phase_emm$year%in%c(y1,y2)
for (p in 1:nPerm) {
  d_perm = data.frame(year=sample(phase_emm$year[ind]), emmean=phase_emm$emmean[ind])
  stat_array[p] = var_stat(d_perm, y1,y2)
}

stat_obs = var_stat(phase_emm, y1, y2)

# Empirical p-value
(sum(stat_array>=stat_obs)+1)/(nPerm+1)

hist(stat_array, 100)
