#Tiller number data analysis
#Dana Looschelders
#Changelog:
  #07.02.20 add new data
#Comments marked with #********* contain instructions for neccessary changes to script

library(ggplot2) #load packages
library(tidyr)
library(dplyr)
library(tidyverse)

#*********************************************************************
#choose your working directory (should contain a csv.file with data)
#**********************************************************************

setwd("C:/00 Dana/Uni/Internship/Work") #set working directory (place of files)

#import data
Tiller=read.table("Tiller Number Data Collection.csv", sep=";", dec=".", header=T)
str(Tiller) #check data
#delete empty cloumns off data.frame
Tiller=Filter(function(x)!all(is.na(x)), Tiller)
str(Tiller)#check again

#calculate mean Tiller number for every variety in every chamber
#function for statistics for every week (just add data and run function)

#*************************************************************************************************
#names: names for columns in output dataframe (calles stats_Tiller)
#no_col: number of columns to process (DON'T count first column as dataframe is created with it)
#***************************************************************************************************
str(Tiller)
str(Height)
stats.Tiller=function(no_col=5, names=c("week 2", "week 3", "week 4", "week 5", "week 6")){
  cols=seq(from=6, by=1, to=6+no_col-1)
  stats_Tiller=aggregate(Tiller[,5], by=list(Varietys=Tiller$Variety, Chamber=Tiller$Chamber), mean)
  names(stats_Tiller)[3]="week 1"
  for (i in cols) {
    dummy_stats=aggregate(Tiller[,i], by=list(Varietys=Tiller$Variety, Chamber=Tiller$Chamber), mean)
    col=i-2
    stats_Tiller[,col]=dummy_stats$x
    names(stats_Tiller)[col]=names[i-5]
  }
  print(stats_Tiller)
}

stats_Tiller=stats.Tiller()
write.csv(x = stats_Tiller, file = "stats_Tiller.csv")

#create time series
#***********************************************************************************************************
#to modify for new data: add Chamber_x$'weex x' to function rbind and and the number for the week to cbind
#************************************************************************************************************

get.tiller.time.series=function(Chambernr=1, name=One){
Chamber_x=stats_Tiller[stats_Tiller$Chamber==Chambernr,]
Chamber_try=rbind(Chamber_x$`week 1`, Chamber_x$`week 2`, Chamber_x$`week 3`, Chamber_x$`week 4`, Chamber_x$`week 5`, Chamber_x$`week 6`)
colnames(Chamber_try)=Chamber_x$Varietys
Chamber_try=cbind(Chamber_try, week=c(1,2,3,4,5,6))
Chamber_try=as.data.frame(Chamber_try)
names(Chamber_try)[8]="semi_natural11"
names(Chamber_try)[9]="semi_natural16"
names(Chamber_try)[10]="semi_natural17"
#collapse the measures for the different varieties into one column in order to plot the time series
name <- Chamber_try %>%
  select(Aberchoice, Abergain, Aspect, Carraig, Dunluce, Lilora, Moy, semi_natural11, semi_natural16, semi_natural17,Solomon, Wild4, Wild6, Wild7, week) %>%
  gather(key = "variable", value = "value", -week)
print(name)
ggplot(name, aes(x = week, y = value)) + #basic plot
  geom_line(aes(color = variable), size = 1) + #adds tiller number as lines
  geom_point(aes(shape=variable), size=3)+ #adds points
  theme_minimal()+ #defines design of plot
  scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))+ #changes shape of points
  ylab(label="Tiller number ")+ #labels y axis
  xlab(label="Week")+ #labels x axis
  labs(title="Chamber 1 Lolium perenne Tiller number 6th Feb 2020", #title of plot
       subtitle = "bar: biomass cut")+ #subtitle of plot
  geom_vline(xintercept=2)+ #adds vertical line to indicate biomass cut 
  geom_vline(xintercept=5)
}

get.tiller.time.series(Chambernr = 8, name=Eight)

Chambers=c(1,2,7,8)

#function to save plots for all four chambers automatically in working directory
save.tiller.plot=function(){
  for (i in Chambers){
    title=paste("Tiller number",i,".pdf") 
    get.tiller.time.series(Chambernr=i)
    ggsave(file=title, width=12, height=8)
  }
}
save.tiller.plot() #call function to execute
