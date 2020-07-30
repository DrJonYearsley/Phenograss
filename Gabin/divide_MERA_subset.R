# Read subsetted MERA data and output a csv file for each square
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# July 2020
# *************************************************************

rm(list=ls())

setwd('~/WorkFiles/MEGA/Projects/GrasslandPhenology/Data/')
load('Temp2m_2017_squares.RData')


squares = unique(temp_final$ID)

for (s in squares) {
  temp_data = subset(temp_final, ID==s)
  
  write.csv(temp_data, 
            file=paste0('Temp2m_2017_square_',s,'.csv'), 
            row.names = FALSE, 
            quote=FALSE)
}

