# Compare precipitation to soil moisture
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# July2021
# ++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd('/Volumes/MODIS_data/MERA/DailyData_subsetted/')

year = 2017


files = list.files(pattern=paste0(year,'.RData'))

load(files[grepl("Precip",files)])

mera = mera_daily

load(files[grepl("Moist",files)])
