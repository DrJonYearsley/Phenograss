# Script to import germination data and preproces it 
# for use in a survival analysis
# 
# Processing requires a dataframe with each row
# being an object (seed). Treatment, Variety and SeedID
# are recorded and then the day of germination, and whether 
# germination was ever recorded. 
# If a seed never germinated then germination day is set 
# to the maximum day. These seeds will be treated as 
# right-censored in the survival analysis.
#
# Survial analysis is performed by script germination_analysis.R
#
# Jon Yearsley (24th March 2020)
# Jon.Yearsley@ucd.ie
# **************************************************



setwd('~/MEGA/MainFolder/Projects/GrasslandPhenology/Data/')
rm(list=ls())

library(readxl)
library(tidyr)

# Read raw germination data (skip first line)
d = read_excel('Raw Data Germination Trial.xlsx', 
               skip=1, 
               col_names=TRUE)

# Relabel Treatment
# Remove waterlogging treatment because it is not relevant
# relabel treatments
d$Treatment[d$Treatment=='CON'] = 'Ambient'
d$Treatment[d$Treatment=='eCO2'] = '+CO2+Temp'
d$Treatment[d$Treatment=='WAT'] = 'Ambient'
d$Treatment[d$Treatment=='eCO2W'] = '+CO2+Temp'

# Pivot raw data into long format
d_long = pivot_longer(data=d, 
                   names_to='Date', 
                   values_to='NumSeeds', 
                   cols=c(4:ncol(d)))

# Look at expt design
table(d$Variety, d$Treatment)

# Extract Day as a numerical variable
d_long$Day = as.numeric(substring(d_long$Date, first = 5))

# Combine data in d_long from different plantsID's
d_long_ag = aggregate(NumSeeds~Variety+Treatment+Day, 
                      data=d_long, 
                      FUN=sum)

# Count total numbers germinating per variety per treatment 
# Do this by comparing proportions to total seed germinated

# Check: Moy and all the semi-natural and wild genotypes were sown 
# 5 seeds per pot and all the other cultivars were sown 10 
# seeds per pot. 
d_agg = aggregate(NumSeeds~Variety+Treatment, data=d_long, FUN=sum)


props = as.data.frame(read_excel('Percentage Cumulative Germination.xlsx'))
props = props[,c(1,2,ncol(props))]
props$Treatment[props$Treatment=='eCO2'] = '+CO2+Temp'
names(props) = c('Treatment', 'Variety', 'Frac')

# Calculate the total number of seeds as variable SeedTot
d_tot = merge(d_agg, props)
d_tot$SeedTot = round(d_tot$NumSeeds / d_tot$Frac * 100)

# Set Treatment aand Variety as factors
d_tot$Treatment = as.factor(d_tot$Treatment)
d_tot$Variety = as.factor(d_tot$Variety)

# Make the final data frame for the survival analysis
varList = levels(d_tot$Variety)
treatList = levels(d_tot$Treatment)
nTreatment = nlevels(d_tot$Treatment)
nVariety = nlevels(d_tot$Variety)


# Loop over all treatments and varieties
seedCount = 0
for (t in 1:nTreatment) {
  for (v in 1:nVariety) {
    tmp1 = subset(d_tot, Variety==varList[v] & Treatment==treatList[t])
    tmp2 = subset(d_long_ag, Variety==varList[v] & 
                    Treatment==treatList[t] & 
                    NumSeeds>0)
    
    
    d_tmp = data.frame(SeedID=seedCount+c(1:tmp1$SeedTot),
                       Treatment=treatList[t],
                       Variety=varList[v],
                       Day = NA)
    
    for (r in 1:nrow(tmp2)) {
      if (r==1) {
        tmp3 = rep(tmp2$Day[r], each=tmp2$NumSeeds[r])
      } else {
        tmp3 = c(tmp3, rep(tmp2$Day[r], each=tmp2$NumSeeds[r]))
      } 
    }
    
    # Fill tmp3 data into d_tmp
    d_tmp$Day[1:length(tmp3)] = tmp3
    
    if (t==1 & v==1) {
      germination = d_tmp
    } else {
      germination = rbind(germination, d_tmp)
    }
    seedCount = seedCount + tmp1$SeedTot
  }
}

# Add in the Date as a POSIXct 
germination$Date = as.Date(germination$Day-1, origin=as.Date('16/12/2019', format='%d/%m/%Y'))

# Add in booblean for whether a seed germinated
germination$germinated = !is.na(germination$Day)

# Code non-germinating seeds to have Day=maximum day
germination$Day[!germination$germinated] = max(germination$Day, na.rm=T)

# Add in Start of observation period (i.e. 1 Day earlier)
germination$Day0 = germination$Day-1


# Add in Variety groupings
wildList = c('Wild6', 'Wild7','Wild4')
seminaturalList = c('Semi-natural6','Semi-natural7','Semi-natural11')

germination$Status = 'Cultivar'
germination$Status[germination$Variety%in%seminaturalList]  = 'Seminatural'
germination$Status[germination$Variety%in%wildList]  = 'Wild'
germination$Status = as.factor(germination$Status)

# Ploidy
tetraploidList = c('Carraig','Dunluce','Abergain','Aspect')
germination$ploidy = 'Diploid'
germination$ploidy[germination$Variety%in%tetraploidList] = 'Tetraploid'
germination$ploidy = as.factor(germination$ploidy)

# Heading
earlyList = c('Moy','Lilora')
intermediateList = c('Carraig','Dunluce','Solomon')
lateList = c('Aberchoice','Abergain','Aspect')

germination$heading = NA
germination$heading[germination$Variety%in%earlyList] = 'Early'
germination$heading[germination$Variety%in%intermediateList] = 'Intermediate'
germination$heading[germination$Variety%in%lateList] = 'Late'
germination$heading = as.factor(germination$heading)


save(germination, file='germination_data.RData')


