#read new data in R
#load necessary librarys
library(readxl) #library to read .xlsx files (otherwise .csv files are needed)
library(ggplot2) #library for graphing data
library(tidyverse)#library to organise data

#packages for mixed models
#library(lme4)
#library(nlme)
#library(arm)
#library(car)
#library(MASS)

#read in data directly from github folder
Biomass=read_excel("data_rosemount/Biomass Data.xlsx")
str(Biomass) #check if data was read in correctly

D50=read_excel("data_rosemount/D50 Data.xlsx")
str(D50) #check if data was read in correctly

Tiller_first=read_excel("data_rosemount/Date of 1st Tiller Appearance.xlsx")
str(Tiller_first) #check if data was read in correctly

Senescence=read_excel("data_rosemount/Date of Leaf Senescence.xlsx")
str(Senescence) #check if data was read in correctly

Flower=read_excel("data_rosemount/Flowering Plants and Tillers.xlsx")
str(Flower) #check if data was read in correctly

Pheno_Aug=read_excel("data_rosemount/Leaf Phenology August.xlsx")
str(Pheno_Aug) #check if data was read in correctly

Pheno_Jul=read_excel("data_rosemount/Leaf Phenology July.xlsx")
str(Pheno_Jul) #check if data was read in correctly

Pheno_Jun=read_excel("data_rosemount/Leaf Phenology June.xlsx")
str(Pheno_Jun) #check if data was read in correctly

Pheno_May=read_excel("data_rosemount/Leaf Phenology May.xlsx")
str(Pheno_May) #check if data was read in correctly

Pheno_Sep=read_excel("data_rosemount/Leaf Phenology September.xlsx")
str(Pheno_Sep) #check if data was read in correctly

Phyllo_12=read_excel("data_rosemount/Leaf Phyllochron 1-2.xlsx")
str(Phyllo_12) #check if data was read in correctly

Phyllo_13=read_excel("data_rosemount/Leaf Phyllochron 1-3.xlsx")
str(Phyllo_13) #check if data was read in correctly

Phyllo_14=read_excel("data_rosemount/Leaf Phyllochron 1-4.xlsx")
str(Phyllo_14) #check if data was read in correctly

Phyllo_23=read_excel("data_rosemount/Leaf Phyllochron 2-3.xlsx")
str(Phyllo_23) #check if data was read in correctly

Phyllo_24=read_excel("data_rosemount/Leaf Phyllochron 2-4.xlsx")
str(Phyllo_24) #check if data was read in correctly

Phyllo_34=read_excel("data_rosemount/Leaf Phyllochron 3-4.xlsx")
str(Phyllo_34) #check if data was read in correctly

Phyllo_34=read_excel("data_rosemount/Leaf Phyllochron 3-4.xlsx")
str(Phyllo_34) #check if data was read in correctly

Phyllo_Germ_Cum=read_excel("data_rosemount/Percentage Cumulative Germination.xlsx")
str(Phyllo_Germ_Cum) #check if data was read in correctly
#column name is a date and therefore not read correctly 
#create column names
dates=seq.Date(from=as.Date("2019-12-16"), to=as.Date("2020-01-14"), by = "day") #create date sequence
columnnames=c("Treatment", "Variety", as.character(dates)) #create vector with column names
Phyllo_Germ_Cum=read_excel("data_rosemount/Percentage Cumulative Germination.xlsx", 
                           col_names = as.character(columnnames)) #read in data and assign creates vector with column names
str(Phyllo_Germ_Cum) #check again

Height=read_excel("data_rosemount/Plant Height Data.xlsx")
str(Height) #check if data was read in correctly

Growth_rate=read_excel("data_rosemount/Growth Rate Data.xlsx", 
                       col_types = c("numeric", "text","text", 
                                     "numeric","numeric","numeric",
                                     "numeric","numeric", "numeric"))
#this time coltypes need tobe specified, otherwise reads the growth rate as string
str(Growth_rate) #check if data was read in correctly
#there is extra data in the spreadsheet -> remove those rows
Growth_rate=Growth_rate[1:336,]

Max_Growth=read_excel("data_rosemount/Maximum Plant Height Data.xlsx")
str(Max_Growth) #check if data was read in correctly

Germ_raw=read_excel("data_rosemount/Raw Data Germination Trial.xlsx")
str(Germ_raw) #check if data was read in correctly

Tiller=read_excel("data_rosemount/Tiller Number Data.xlsx")
str(Tiller) #check if data was read in correctly

Germ_Tiller=read_excel("data_rosemount/Time to 1st Tiller from Germination.xlsx")
str(Germ_Tiller) #check if data was read in correctly

Germ_Sene=read_excel("data_rosemount/Time to Leaf Senescence from Germination.xlsx")
str(Germ_Sene) #check if data was read in correctly

Germ_Cum=read_excel("data_rosemount/Total Cumulative Germination.xlsx")
str(Germ_Cum) #check if data was read in correctly

Germ_total_num=read_excel("data_rosemount/Total Number of Germinated Seeds Data.xlsx")
str(Germ_total_num) #check if data was read in correctly

Transplanted=read_excel("data_rosemount/Transplanted Pots.xlsx")
str(Transplanted) #check if data was read in correctly

Tussock=read_excel("data_rosemount/Tussock Diameter Data.xlsx")
str(Tussock) #check if data was read in correctly

#combine some parameters for mixed model 
dat=cbind.data.frame(Biomass[,c(1:4,10)], Germ_Sene$Senescence, Tiller[,10], 
                     Tussock$`Tussock Diameter 18/05/2020`,
                     Growth_rate[,5:length(Growth_rate)])
str(dat)
colnames(dat)[5:13]=c("Biomass", "Senescence", "Tiller", "Tussock", "Growth_May", "Growth_Jun", "Growth_Jul", "Growth_Aug", "Growth_Sep")

#find fitting distribution
qqp(dat$Biomass, "lnorm") #lognormal
qqp(dat$Biomass, "norm") #normal
shapiro.test(dat$Biomass) #not normaly distributed

qqp(dat$Tussock, "lnorm") #lognormal
qqp(dat$Tussock, "norm") #normal
shapiro.test(dat$Tussock) #p-value 0.013 -> prob normally distributed

qqp(dat$Senescence, "lnorm") #lognormal
qqp(dat$Senescence, "norm") #normal
shapiro.test(dat$Senescence) #not normally distributed

qqp(dat$Tiller, "lnorm") #lognormal
qqp(dat$Tiller, "norm") #normal
shapiro.test(dat$Tiller) #not normally distributed

qqp(dat$Growth_May, "lnorm") #lognormal
qqp(dat$Growth_May, "norm") #normal
shapiro.test(dat$Growth_May) #p value 0.03

qqp(dat$Growth_Jun, "lnorm") #lognormal
qqp(dat$Growth_Jun, "norm") #normal
shapiro.test(dat$Growth_Jun) #p value 0.058

qqp(dat$Growth_Jul, "lnorm") #lognormal
qqp(dat$Growth_Jul, "norm") #normal
shapiro.test(dat$Growth_Jul) #p value 0.003

qqp(dat$Growth_Aug, "lnorm") #lognormal
qqp(dat$Growth_Aug, "norm") #normal
shapiro.test(dat$Growth_Aug) #not normally distributed

qqp(dat$Growth_Sep, "lnorm") #lognormal
qqp(dat$Growth_Sep, "norm") #normal
shapiro.test(dat$Growth_Sep) #not normally distributed

#most of the data is not normally distributed
