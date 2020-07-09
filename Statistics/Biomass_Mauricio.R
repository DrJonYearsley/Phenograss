### Mauricio Biomass 09 July 2020

#Remove older files and clear WD
rm(list=ls())

#Set Working Directory
setwd('E:')

#load necessary librarys to read new data in R
library(readxl) #library to read .xlsx files (otherwise .csv files are needed)
library(ggplot2) #library for graphing data
library(tidyverse)#library to organise data

#other packages that are useful
library(lme4)
library(nlme)
library(arm)
library(car)
library(MASS)

# download data and replace spaces in names with '.'
Biomass=read_excel("data_rosemount/Biomass Data.xlsx", .name_repair = 'universal')
str(Biomass) #check if data was read in correctly

#Put data in long format
Biomass_long = pivot_longer(Biomass,
                             cols = -c(1:5),
                             names_to = 'Month',
                             values_to = 'Grams',
                             values_drop_na = TRUE)

#Check whether data was imported correctly
str(Biomass_long)

# Visualise the data
format = theme(axis.title = element_text(size=14),
            axis.text = element_text(size=12))


#Box plot with Treatments ###OBS: Need to put X-axis in order by month
ggplot(data=Biomass_long, 
       aes(x=Month, y=Grams, fill=Treatment))+ #specifies the data used for x and y
  geom_boxplot(width=0.2)+ #adds tiny boxplot in violin 
  scale_fill_brewer(name='Treatment', palette='Set2')+
  xlab("Month")+ #adds title for x axis
  ylab("Dried Biomass (g)")+ #adds title for y axis
  ggtitle("Biomass")+ #adds plot title
  theme_bw()+
  format


#Violin plot with Genetic Status
ggplot(data=Biomass_long, 
       aes(x=Genetic, y=Grams, fill=Treatment))+ #specifies the data used for x and y
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ #adds violin shape without trimmed edges
  scale_fill_brewer(name='Treatment', palette='Set2')+
  xlab("Genetic Status")+ #adds title for x axis
  ylab("Dried Biomass (g)")+ #adds title for y axis
  ggtitle("Biomass")+ #adds plot title
  theme_bw()+
  format


#Violin plot by Variety with boxplot ###OBS: Boxplot are not inside the violin
ggplot(data=Biomass_long,
       aes(x=Variety, y=Grams, fill=Treatment))+ #specifies the data used for x and y
  geom_violin(trim=F)+ #adds violin shape without trimmed edges
  geom_boxplot(width=0.1)+ #adds boxplot to violin
  scale_fill_brewer(name='Treatment', palette='Set2')+ # legend and colour of violins/boxplots
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ #rotates the axix labels on x axis
  xlab("Variety")+ #adds title for x axis
  ylab("Dried Biomass (g)")+ #adds title for y axis
  ggtitle("Dried Biomass (g) per Variety") #adds plot title


#Violin plot by Variety with boxplot displaying Genetic Status in different colours
ggplot(data=Biomass_long,
       aes(x=Variety, y=Grams, fill=Genetic))+ #specifies the data used for x and y
  geom_violin(trim=F)+ #adds violin shape without trimmed edges
  geom_boxplot(width=0.1)+ #adds boxplot to violin
  scale_fill_brewer(name='Genetic Status', palette='Set2')+ # legend and colour of violins/boxplots
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ #rotates the axix labels on x axis
  xlab("Variety")+ #adds title for x axis
  ylab("Dried Biomass (g)")+ #adds title for y axis
  ggtitle("Dried Biomass (g) per Variety") #adds plot title


#Histogram plot with Frequency
hist(Biomass_long$Grams, #with frequency
     main = "Histogram of Biomass (Frequency)") #plot title

#Histogram plot with Density
hist(Biomass_long$Grams, prob=TRUE, 
     main = 'Histogram of Biomass (Density)', col='green') #histogram 
lines(density(Biomass_long$Grams), col="black") #density as red line

#Summary to get all relevant statistical parameters
summary(Biomass_long)

#Do it for the other parameters and perhaps subset for an especific time/treament ###===###