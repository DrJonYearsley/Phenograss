library(ggplot2) #load packages
library(tidyr)
library(tidyverse)

#*********************************************************************
#choose your working directory (should contain a csv.file with data)
#**********************************************************************
setwd("C:/00 Dana/Uni/Internship/Work") #set working directory (place of files)
#import data
Biomass=read.table("Dried Biomass Data Collection.csv", sep=";", dec=".", header=T)
str(Biomass) #check data
#delete empty cloumns off data.frame
Biomass=Filter(function(x)!all(is.na(x)), Biomass)
str(Biomass)#check again

#calculate mean Biomass for every variety in every Treatment
#function for statistics for every week (just add data and run function)
#names: names for columns in output dataframe (calles stats_Biomass)
#no_col: number of columns to process (DON'T count first column as dataframe is created with it)

#*************************************************************************************************
#names: names for columns in output dataframe (calles stats_Tiller)
#no_col: number of columns to process (DON'T count first column as dataframe is created with it)
#***************************************************************************************************

stats.Biomass=function(no_col=1, names=c("cut 2")){
  cols=seq(from=6, by=1, to=6+no_col-1)
  stats_Biomass=aggregate(Biomass[,5], by=list(Varietys=Biomass$Variety, Treatment=Biomass$Treatment), mean)
  names(stats_Biomass)[3]="cut 1"
  for (i in cols) {
    dummy_stats=aggregate(Biomass[,i], by=list(Varietys=Biomass$Variety, Treatment=Biomass$Treatment), mean)
    col=i-2
    stats_Biomass[,col]=dummy_stats$x
    names(stats_Biomass)[col]=names[i-5]
  }
  print(stats_Biomass)
}

stats_Biomass=stats.Biomass()
str(stats_Biomass)
write.csv(x = stats_Biomass, file="stats_Biomass.csv")
#create time series
#***********************************************************************************************************
#to modify for new data: add Treatment_x$'weex x' to function rbind and and the number for the week to cbind
#************************************************************************************************************

get.Biomass.time.series=function(Treatmentnr=1){
  Treatment=stats_Biomass[stats_Biomass$Treatment==Treatmentnr,]
  Treatment_dummy=rbind(Treatment$`cut 1`, Treatment$`cut 2`)
  colnames(Treatment_dummy)=Treatment$Varietys
  Treatment_dummy=cbind(Treatment_dummy,week=c(1,2))
  Treatment_dummy=as.data.frame(Treatment_dummy)
  #the semi-natural varieties need to be renamed because the - messes with R (an underscore is fine)
  names(Treatment_dummy)[8]="semi_natural11" 
  names(Treatment_dummy)[9]="semi_natural16"
  names(Treatment_dummy)[10]="semi_natural17"
  #collapse the measures for the different varieties into one column in order to plot the time series
  df <- Treatment_dummy %>% #creates new dataframe
    select(Aberchoice, Abergain, Aspect, Carraig, Dunluce, Lilora, Moy, semi_natural11, semi_natural16, semi_natural17,Solomon, Wild4, Wild6, Wild7, week) %>%
    gather(key = "variable", value = "value", -week) 
  #plot as time series
  ggplot(df, aes(x = week, y = value)) + #basic plot 
    geom_line(aes(color = variable), size = 1) + #adds lines
    geom_point(aes(shape=variable), size=3)+ #adds points
    theme_minimal()+ #defines design of plot
    scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))+ #changes shape of points
    ylab(label="Plant Biomass [cm]")+ #labels y axes
    xlab(label="Cut")+ #labels x axes
    labs(title=paste("Treatment",Treatmentnr,"Lolium perenne Biomass 27. Feb 2020")) #title of plot #subtitle of plot
}
get.Biomass.time.series(Treatmentnr=1)

Treatments=c(1,2,7,8)

#function to save plots for all four Treatments
save.Biomass.plot=function(){
  for (i in Treatments){
    title=paste("Biomass",i,".pdf") 
    get.Biomass.time.series(Treatmentnr=i)
    ggsave(title, width=12, height=8)
  }
}
save.Biomass.plot() #call function to execute

#for a single date plot with ggplot2
#example with 30th Jan
H_plot3=ggplot(stats_Biomass_30, aes(y=avg_Biomass_30 ,x=Varietys)) #create plot
H_bars3=geom_bar(stat="identity") #create bar plot layer
#add mean per Treatment
mean_Biomass_Treatment3=aggregate(stats_Biomass_30$avg_Biomass_30, 
                                by = list(Treatment=stats_Biomass_30$Treatment), mean) #calculate mean per Treatment
str(mean_Biomass_Treatment) #check data
H_lines3=geom_hline(data=mean_Biomass_Treatment3, aes(yintercept=x), color="red") #create line layer for mean
H_plot3+H_bars3+facet_grid(.~Treatment)+
  labs(title="Lolium perenne Biomass 30. Jan 2020", 
       subtitle = "Treatments: 1,2,7,8 \n red bar: mean per Treatment") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))+
  H_lines3+ylab(label="Average Biomass [cm]")+
  xlab(label="Ryegrass Varieties")


#plot with growth
#re-arrange dataframe in order to display growth as stacked bars
stats_Biomass_dummy=data.frame(Dummy=rep(1,112), "Biomass"=stats_Biomass$avg_Biomass_17, "Varieties"=stats_Biomass$Varietys, "Treatment"=stats_Biomass$Treatment)
stats_Biomass_dummy[57:112,2]=stats_Biomass$Growth
stats_Biomass_dummy[57:112,3]=stats_Biomass$Varietys
stats_Biomass_dummy[57:112,4]=stats_Biomass$Treatment
stats_Biomass_dummy[57:112,1]=rep(2,56)
str(stats_Biomass_dummy)
#plot with stacked bars
H_plot3=ggplot(stats_Biomass_dummy, aes(fill=Dummy,y=Biomass ,x=Varieties)) #create plot
H_bars3=geom_bar(stat="identity", position="stack") #create bar plot layer
H_lines=geom_hline(data=mean_Biomass_Treatment, aes(yintercept=x), color="grey") #create line layer for mean
H_lines2=geom_hline(data=mean_Biomass_Treatment2, aes(yintercept=x), color="black") #create line layer for mean
H_plot3+H_bars3+facet_grid(.~Treatment)+
  labs(title="Lolium perenne growth 23rd Jan 2020",
       subtitle = "Treatments: 1,2,7,8\n darkblue: growth until 17th \n lightblue: growth until 23d \n grey bar: mean per Treatment on 17th \n black bar: mean per Treatment on 23rd") +
  theme_minimal()+theme(axis.text.x = element_text(angle = 90))+
  H_lines+H_lines2+ylab(label="Average Biomass [cm]")+
  xlab(label="Ryegrass Varieties")+guides(fill=FALSE)

#Export stats
write.csv(stats_Biomass, "Biomass_statistics.csv")

