# Double check the MERA  data
# Against daily data for Mace Head
#
# Station Name: MACE HEAD		
# Station Height: 21 M 		
# Latitude:53.326  	Longitude: -9.901	
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# July 2021
# +++++++++++++++++++++++++++++++++++++++++++

setwd('/Volumes/MERA_Data')
library(data.table)

mace_latlong = c(53.326, 360-9.901)

# Read in data from MAce Head Observatory
mace = read.csv('Mace_head_daily_data.csv', skip=24)

mace$date = as.Date(mace$date, format="%d-%b-%Y")
mace$year = as.numeric(format(mace$date,"%Y"))
mace$month = as.numeric(format(mace$date,"%m"))

mace_sub = subset(mace, year==2015 & month==11)



# Read in MERA data 
mera = fread('TPrecip/TotalPrecip_2015_11.txt')

# Remove header rows within the data
ind = mera$Latitude=="Latitude"
mera = mera[!ind,]

mera = as.data.frame(lapply(mera,as.numeric))
nam = names(mera)
nam[3] = 'Value'
names(mera) = nam

# Subset Mace Head
delta = abs(mera$Latitude - mace_latlong[1]) +  abs(mera$Longitude - mace_latlong[2])
mace_mera = mera[delta==min(delta),]

# Correlate rain from the two data sets
cor(mace_mera$Value, mace_sub$rain)            # Corresponds to dataDate
cor(mace_mera$Value[-nrow(mace_mera)], mace_sub$rain[-1])   # Corresponds to validityDate


plot(mace_mera$Value, mace_sub$rain)                        # Corresponds to dataDate
plot(mace_mera$Value[-nrow(mace_mera)], mace_sub$rain[-1])  # Corresponds to validityDate

